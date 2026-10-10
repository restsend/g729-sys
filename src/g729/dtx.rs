// G.729 Annex B discontinuous transmission (DTX), ported from bcg729/src/dtx.c.

use crate::g729::basic_operations::*;
use crate::g729::cng::compute_comfort_noise_excitation_vector;
use crate::g729::codebooks::SID_GAIN_CODEBOOK;
use crate::g729::compute_lp::auto_correlation_2_lp;
use crate::g729::fixed_point_math::g729_log2_q0q16;
use crate::g729::interpolate_q_lsp::interpolate_q_lsp;
use crate::g729::ld8k::*;
use crate::g729::lp2lsp_conversion::lp2lsp_conversion;
use crate::g729::lsp_quantization::noise_lsp_quantization;
use crate::g729::q_lsp_2_lp::q_lsp_2_lp;

const SID_FRAME: u8 = 2;
const UNTRANSMITTED_FRAME: u8 = 0;

/// DTX state; kept inline (no heap) for `no_std`.
pub struct DtxChannelContext {
    autocorrelation_coefficients: [[i32; NB_LSP_COEFF + 1]; 7],
    autocorrelation_coefficients_scale: [i8; 7],
    previous_vad_flag: u8,
    pseudo_random_seed: u16,
    previous_residual_energy: i32,
    previous_residual_energy_scale: i8,
    sid_lp_coefficient_autocorrelation: [i32; NB_LSP_COEFF + 1],
    current_sid_gain: i16,
    previous_decoded_log_energy: i8,
    decoded_log_energy: i8,
    smoothed_sid_gain: i16,
    reflection_coefficients: [i32; NB_LSP_COEFF],
    count_fr: u8,
    q_lsp_coefficients: [i16; NB_LSP_COEFF],
}

impl Default for DtxChannelContext {
    fn default() -> Self {
        Self::new()
    }
}

impl DtxChannelContext {
    pub fn new() -> Self {
        let mut ctx = DtxChannelContext {
            autocorrelation_coefficients: [[0; NB_LSP_COEFF + 1]; 7],
            autocorrelation_coefficients_scale: [0; 7],
            previous_vad_flag: 1,
            pseudo_random_seed: CNG_DTX_RANDOM_SEED_INIT,
            previous_residual_energy: 0,
            previous_residual_energy_scale: 0,
            sid_lp_coefficient_autocorrelation: [0; NB_LSP_COEFF + 1],
            current_sid_gain: 0,
            previous_decoded_log_energy: 0,
            decoded_log_energy: 0,
            smoothed_sid_gain: 0,
            reflection_coefficients: [0; NB_LSP_COEFF],
            count_fr: 0,
            q_lsp_coefficients: [0; NB_LSP_COEFF],
        };
        // Set the past autocorrelation[0] to 1 to avoid arithmetic problems.
        for i in 0..7 {
            ctx.autocorrelation_coefficients[i][0] = ONE_IN_Q30;
            ctx.autocorrelation_coefficients_scale[i] = 30;
        }
        ctx
    }

    /// Reflection coefficients (Q31) of the filter used for the last transmitted
    /// SID frame, used to build the RFC3389 payload.
    pub fn reflection_coefficients(&self) -> &[i32; NB_LSP_COEFF] {
        &self.reflection_coefficients
    }

    /// Decoded frame mean energy in dB, used to build the RFC3389 payload.
    pub fn decoded_log_energy(&self) -> i8 {
        self.decoded_log_energy
    }
}

/// Rescale and sum several autocorrelation vectors.
fn sum_autocorrelation_coefficients(
    autocorrelation_coefficients: &[[i32; NB_LSP_COEFF + 1]],
    autocorrelation_coefficients_scale: &[i8],
    nb_elements: usize,
    auto_correlation_coefficients_result: &mut [i32; NB_LSP_COEFF + 1],
    autocorrelation_coefficients_scale_results: &mut i8,
) {
    let mut auto_correlation_sum_buffer = [0i64; NB_LSP_COEFF + 1];
    let mut max: i64 = 0;
    let mut rescaled_autocorrelation_coefficients = [[0i32; NB_LSP_COEFF + 1]; 7];
    let mut right_shift_to_normalise: i8 = 0;

    let mut min_scale = autocorrelation_coefficients_scale[0];
    for &scale in &autocorrelation_coefficients_scale[1..nb_elements] {
        if scale < min_scale {
            min_scale = scale;
        }
    }

    for j in 0..nb_elements {
        let rescaling = autocorrelation_coefficients_scale[j] - min_scale;
        for i in 0..NB_LSP_COEFF + 1 {
            rescaled_autocorrelation_coefficients[j][i] =
                shr32(autocorrelation_coefficients[j][i], rescaling as u32);
        }
    }

    for (i, sum) in auto_correlation_sum_buffer.iter_mut().enumerate() {
        *sum = rescaled_autocorrelation_coefficients[0][i] as i64;
        for row in rescaled_autocorrelation_coefficients
            .iter()
            .take(nb_elements)
            .skip(1)
        {
            *sum = add64(*sum, row[i] as i64);
        }
        if sum.abs() > max {
            max = sum.abs();
        }
    }

    if max > MAXINT32 as i64 {
        loop {
            max = shr64(max, 1);
            right_shift_to_normalise += 1;
            if max <= MAXINT32 as i64 {
                break;
            }
        }
        for i in 0..NB_LSP_COEFF + 1 {
            auto_correlation_coefficients_result[i] = shr64(
                auto_correlation_sum_buffer[i],
                right_shift_to_normalise as u32,
            ) as i32;
        }
    } else {
        for i in 0..NB_LSP_COEFF + 1 {
            auto_correlation_coefficients_result[i] = auto_correlation_sum_buffer[i] as i32;
        }
    }

    *autocorrelation_coefficients_scale_results = min_scale - right_shift_to_normalise;
}

/// Residual energy quantization (B4.2.1).
fn residual_energy_quantization(
    residual_energy: i32,
    residual_energy_scale: i8,
    decoded_log_energy: &mut i8,
) -> u8 {
    // -479849 is log2(aw / (NCur*80)).
    let mut acc = sub32(
        g729_log2_q0q16(residual_energy),
        add32(479849, (residual_energy_scale as i32) << 16),
    );
    acc = shr32(acc, 1);
    acc = mult16_32_q15(INV_LOG2_10_Q15, acc);

    if acc < -26214 {
        *decoded_log_energy = -12;
        0
    } else if acc < 45875 {
        let mut acc = acc + 19661;
        if acc < 0 {
            acc = 0;
        } else {
            acc = mult16_32_q13(20480, acc);
        }
        let steps = shr32(acc, 15) as u8;
        *decoded_log_energy = (-2 + 4 * steps as i32) as i8;
        1 + steps
    } else if acc < 216268 {
        let mut acc = acc - 49152;
        if acc < 0 {
            acc = 0;
        } else {
            acc = mult16_32_q12(20480, acc);
        }
        let steps = shr32(acc, 15) as u8;
        *decoded_log_energy = (16 + 2 * steps as i32) as i8;
        6 + steps
    } else {
        *decoded_log_energy = 66;
        31
    }
}

/// LP coefficient autocorrelation (eq B.13).
fn compute_lpc_coefficient_autocorrelation(
    lp_coefficients: &[i16; NB_LSP_COEFF],
    lp_autocorrelation: &mut [i32; NB_LSP_COEFF + 1],
) {
    lp_autocorrelation[0] = (4096 * 4096) >> 4;
    for &coefficient in lp_coefficients.iter() {
        lp_autocorrelation[0] = mac16_16_q4(lp_autocorrelation[0], coefficient, coefficient);
    }

    for j in 1..NB_LSP_COEFF + 1 {
        lp_autocorrelation[j] = shl(lp_coefficients[j - 1] as i32, 9);
        for k in 0..NB_LSP_COEFF - j {
            lp_autocorrelation[j] = mac16_16_q3(
                lp_autocorrelation[j],
                lp_coefficients[k],
                lp_coefficients[k + j],
            );
        }
    }
}

/// Compare two LPC filters (eq B.12); returns 1 when they differ significantly.
fn compare_lpc_filters(
    lp_coefficients_autocorrelation: &[i32; NB_LSP_COEFF + 1],
    autocorrelation_coefficients: &[i32; NB_LSP_COEFF + 1],
    residual_energy: i32,
    threshold: i32,
) -> u8 {
    let mut acc: i64 = 0;
    for i in 0..NB_LSP_COEFF + 1 {
        acc = mac64(
            acc,
            lp_coefficients_autocorrelation[i],
            autocorrelation_coefficients[i],
        );
    }

    if acc >= mult32_32(residual_energy, threshold) {
        1
    } else {
        0
    }
}

/// Save the current autocorrelation vector in the DTX context (B4.1.1).
pub fn update_dtx_context(
    ctx: &mut DtxChannelContext,
    autocorrelation_coefficients: &[i32],
    autocorrelation_coefficients_scale: i8,
) {
    for i in (1..7).rev() {
        ctx.autocorrelation_coefficients[i] = ctx.autocorrelation_coefficients[i - 1];
        ctx.autocorrelation_coefficients_scale[i] = ctx.autocorrelation_coefficients_scale[i - 1];
    }
    ctx.autocorrelation_coefficients[0][..NB_LSP_COEFF + 1]
        .copy_from_slice(&autocorrelation_coefficients[..NB_LSP_COEFF + 1]);
    ctx.autocorrelation_coefficients_scale[0] = autocorrelation_coefficients_scale;
}

/// Called on every frame; updates the VAD-flag history and, for NOISE frames,
/// produces the SID parameters and the comfort-noise excitation.
#[allow(clippy::too_many_arguments)]
pub fn encode_sid_frame(
    ctx: &mut DtxChannelContext,
    previous_lsp_coefficients: &mut [i16; NB_LSP_COEFF],
    previous_q_lsp_coefficients: &mut [i16; NB_LSP_COEFF],
    vad_flag: u8,
    previous_q_lsf: &mut [[i16; NB_LSP_COEFF]; MA_MAX_K],
    excitation_vector: &mut [i16],
    q_lp_coefficients: &mut [i16; 2 * NB_LSP_COEFF],
    bit_stream: &mut [u8],
    bit_stream_length: &mut u8,
) {
    let mut summed_autocorrelation_coefficients = [0i32; NB_LSP_COEFF + 1];
    let mut summed_autocorrelation_coefficients_scale: i8 = 0;
    let mut lp_coefficients = [0i16; NB_LSP_COEFF];
    let mut lsp_coefficients = [0i16; NB_LSP_COEFF];
    let mut reflection_coefficients = [0i32; NB_LSP_COEFF];
    let mut residual_energy: i32 = 0;
    let frame_type: u8;
    let quantized_residual_energy: u8;
    let mut decoded_log_energy: i8 = 0;
    let mut parameters = [0u8; 3];
    let mut interpolated_q_lsp = [0i16; NB_LSP_COEFF];

    if vad_flag == 1 {
        ctx.pseudo_random_seed = CNG_DTX_RANDOM_SEED_INIT;
        ctx.previous_vad_flag = 1;
        return;
    }

    // A NOISE frame: sum the autocorrelation of the current and previous frames.
    sum_autocorrelation_coefficients(
        &ctx.autocorrelation_coefficients,
        &ctx.autocorrelation_coefficients_scale,
        2,
        &mut summed_autocorrelation_coefficients,
        &mut summed_autocorrelation_coefficients_scale,
    );

    auto_correlation_2_lp(
        &summed_autocorrelation_coefficients,
        &mut lp_coefficients,
        &mut reflection_coefficients,
        &mut residual_energy,
    );

    if ctx.previous_vad_flag == 1 {
        // First noise frame after speech: always send a SID frame (B.10).
        frame_type = SID_FRAME;
        quantized_residual_energy = residual_energy_quantization(
            residual_energy,
            summed_autocorrelation_coefficients_scale,
            &mut decoded_log_energy,
        );
    } else {
        let mut flag_chang = 0;

        let mean_energy: i32;
        let mean_energy_scale: i8;
        if summed_autocorrelation_coefficients_scale < ctx.previous_residual_energy_scale {
            mean_energy_scale = summed_autocorrelation_coefficients_scale;
            mean_energy = add32(
                shr(residual_energy, 1),
                svshr32(
                    ctx.previous_residual_energy,
                    (ctx.previous_residual_energy_scale - summed_autocorrelation_coefficients_scale)
                        as i32
                        + 1,
                ),
            );
        } else {
            mean_energy_scale = ctx.previous_residual_energy_scale;
            mean_energy = add32(
                svshr32(
                    residual_energy,
                    (summed_autocorrelation_coefficients_scale - ctx.previous_residual_energy_scale)
                        as i32
                        + 1,
                ),
                shr(ctx.previous_residual_energy, 1),
            );
        }
        quantized_residual_energy =
            residual_energy_quantization(mean_energy, mean_energy_scale, &mut decoded_log_energy);

        if compare_lpc_filters(
            &ctx.sid_lp_coefficient_autocorrelation,
            &summed_autocorrelation_coefficients,
            residual_energy,
            THRESHOLD1_IN_Q20,
        ) != 0
        {
            flag_chang = 1;
        }

        if abs(ctx.previous_decoded_log_energy as i32 - decoded_log_energy as i32) > 2 {
            flag_chang = 1;
        }

        // Decide whether a new SID frame is transmitted (eq B.11).
        ctx.count_fr += 1;
        if ctx.count_fr < 3 {
            frame_type = UNTRANSMITTED_FRAME;
        } else {
            if flag_chang == 1 {
                frame_type = SID_FRAME;
            } else {
                frame_type = UNTRANSMITTED_FRAME;
            }
            ctx.count_fr = 3;
        }
    }

    if frame_type == SID_FRAME {
        let mut sid_lp_autocorrelation_coefficients = [0i32; NB_LSP_COEFF + 1];
        let mut sid_lp_autocorrelation_coefficients_scale: i8 = 0;
        let mut past_average_lp_coefficients = [0i16; NB_LSP_COEFF];
        let mut past_average_reflection_coefficients = [0i32; NB_LSP_COEFF];
        let mut past_average_residual_energy: i32 = 0;

        ctx.count_fr = 0;

        // Past average filter over the last 6 frames (B4.2.2).
        sum_autocorrelation_coefficients(
            &ctx.autocorrelation_coefficients[1..],
            &ctx.autocorrelation_coefficients_scale[1..],
            6,
            &mut sid_lp_autocorrelation_coefficients,
            &mut sid_lp_autocorrelation_coefficients_scale,
        );

        auto_correlation_2_lp(
            &sid_lp_autocorrelation_coefficients,
            &mut past_average_lp_coefficients,
            &mut past_average_reflection_coefficients,
            &mut past_average_residual_energy,
        );

        compute_lpc_coefficient_autocorrelation(
            &past_average_lp_coefficients,
            &mut ctx.sid_lp_coefficient_autocorrelation,
        );

        ctx.decoded_log_energy = decoded_log_energy;

        // Select the filter to encode (eq B.17).
        if compare_lpc_filters(
            &ctx.sid_lp_coefficient_autocorrelation,
            &summed_autocorrelation_coefficients,
            residual_energy,
            THRESHOLD3_IN_Q20,
        ) == 0
        {
            if !lp2lsp_conversion(&past_average_lp_coefficients, &mut lsp_coefficients) {
                lsp_coefficients.copy_from_slice(previous_q_lsp_coefficients);
            }
            ctx.reflection_coefficients
                .copy_from_slice(&past_average_reflection_coefficients);
        } else {
            compute_lpc_coefficient_autocorrelation(
                &lp_coefficients,
                &mut ctx.sid_lp_coefficient_autocorrelation,
            );
            if !lp2lsp_conversion(&lp_coefficients, &mut lsp_coefficients) {
                lsp_coefficients.copy_from_slice(previous_q_lsp_coefficients);
            }
            ctx.reflection_coefficients
                .copy_from_slice(&reflection_coefficients);
        }

        previous_lsp_coefficients.copy_from_slice(&lsp_coefficients);

        noise_lsp_quantization(
            previous_q_lsf,
            &lsp_coefficients,
            &mut ctx.q_lsp_coefficients,
            &mut parameters,
        );

        ctx.previous_decoded_log_energy = decoded_log_energy;
        ctx.current_sid_gain = SID_GAIN_CODEBOOK[quantized_residual_energy as usize];
    }

    ctx.previous_residual_energy = residual_energy;
    ctx.previous_residual_energy_scale = summed_autocorrelation_coefficients_scale;

    // Target gain smoothing (eq B.19).
    if ctx.previous_vad_flag == 1 {
        ctx.smoothed_sid_gain = ctx.current_sid_gain;
    } else {
        ctx.smoothed_sid_gain = sub16(ctx.smoothed_sid_gain, shr16(ctx.smoothed_sid_gain, 3));
        ctx.smoothed_sid_gain = add16(ctx.smoothed_sid_gain, shr16(ctx.current_sid_gain, 3));
    }

    compute_comfort_noise_excitation_vector(
        ctx.smoothed_sid_gain,
        &mut ctx.pseudo_random_seed,
        excitation_vector,
    );

    interpolate_q_lsp(
        previous_q_lsp_coefficients,
        &ctx.q_lsp_coefficients,
        &mut interpolated_q_lsp,
    );
    previous_q_lsp_coefficients.copy_from_slice(&ctx.q_lsp_coefficients);

    q_lsp_2_lp(&interpolated_q_lsp, &mut q_lp_coefficients[0..NB_LSP_COEFF]);
    q_lsp_2_lp(
        &ctx.q_lsp_coefficients,
        &mut q_lp_coefficients[NB_LSP_COEFF..],
    );

    if frame_type == SID_FRAME {
        *bit_stream_length = 2;
        // The reference packs the 2 MSB of L2 with "%0x03" (maps 3 to 0);
        // kept verbatim for bit-exactness.
        bit_stream[0] = ((parameters[0] & 0x01) << 7)
            | ((parameters[1] & 0x1F) << 2)
            | ((parameters[2] >> 2) % 0x03);
        bit_stream[1] = ((parameters[2] & 0x03) << 6) | ((quantized_residual_energy & 0x1F) << 1);
    } else {
        *bit_stream_length = 0;
    }

    ctx.previous_vad_flag = 0;
}
