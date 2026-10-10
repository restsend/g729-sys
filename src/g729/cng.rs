// G.729 Annex B comfort noise generation (CNG), ported from bcg729/src/cng.c.

use crate::g729::basic_operations::*;
use crate::g729::codebooks::{
    L1, L1_SUBSET_INDEX, L2L3, L2_SUBSET_INDEX, L3_SUBSET_INDEX, NOISE_MA_PREDICTOR,
    NOISE_MA_PREDICTOR_SUM, SID_GAIN_CODEBOOK,
};
use crate::g729::decode_adaptative_code_vector::compute_adaptative_codebook_vector;
use crate::g729::decode_lsp::compute_q_lsf;
use crate::g729::fixed_point_math::{
    g729_cos_q13q15, g729_exp2_q11q16, g729_inv_sqrt_q0q31, g729_sqrt_q0q7,
};
use crate::g729::interpolate_q_lsp::interpolate_q_lsp;
use crate::g729::ld8k::*;
use crate::g729::lp2lsp_conversion::lp2lsp_conversion;
use crate::g729::q_lsp_2_lp::q_lsp_2_lp;
use crate::g729::utils::pseudo_random;

pub struct CngChannelContext {
    pub received_sid_gain: i16,
    pub smoothed_sid_gain: i16,
    pub q_lsp: [i16; NB_LSP_COEFF],
}

static SID_Q_LSP_INITIAL_VALUES: [i16; NB_LSP_COEFF] = [
    31441, 27566, 21458, 13612, 4663, -4663, -13612, -21458, -27566, -31441,
];

pub fn init_bcg729_cng_channel() -> CngChannelContext {
    CngChannelContext {
        received_sid_gain: 0,
        smoothed_sid_gain: 0,
        q_lsp: SID_Q_LSP_INITIAL_VALUES,
    }
}

/// Comfort noise excitation generation (B4.4/B4.5).
///
/// `excitation_vector` is the full buffer; the current frame starts at
/// `L_PAST_EXCITATION`.
pub fn compute_comfort_noise_excitation_vector(
    target_gain: i16,
    random_generator_seed: &mut u16,
    excitation_vector: &mut [i16],
) {
    for subframe_index in (0..L_FRAME).step_by(L_SUBFRAME) {
        let base = L_PAST_EXCITATION + subframe_index;
        let mut gaussian_random_excitation = [0i16; L_SUBFRAME];
        let mut eg: i32 = 0;
        let mut gg: i32;

        let mut ea: i32 = 0;
        let mut ei: i32;

        let mut gf: i32;

        let mut sign = [0i16; 4];
        let mut position = [0i16; 4];
        let mut delta_scale_factor: u8 = 0;
        let mut delta: i64;

        let mut random_number_buffer = pseudo_random(random_generator_seed);
        let mut frac_pitch_delay = (random_number_buffer & 0x0003) as i16 - 1;
        if frac_pitch_delay == 2 {
            frac_pitch_delay = 0;
        }
        random_number_buffer >>= 2;
        let int_pitch_delay = (random_number_buffer & 0x003F) as i16 + 40;
        random_number_buffer >>= 6;

        position[0] = ((random_number_buffer & 0x0007) * 5) as i16;
        random_number_buffer >>= 3;
        sign[0] = (random_number_buffer & 0x0001) as i16;
        random_number_buffer >>= 1;
        position[1] = ((random_number_buffer & 0x0007) * 5 + 1) as i16;
        random_number_buffer >>= 3;
        sign[1] = (random_number_buffer & 0x0001) as i16;
        random_number_buffer = pseudo_random(random_generator_seed);
        position[2] = ((random_number_buffer & 0x0007) * 5 + 2) as i16;
        random_number_buffer >>= 3;
        sign[2] = (random_number_buffer & 0x0001) as i16;
        random_number_buffer >>= 1;
        position[3] = ((random_number_buffer & 0x0001) + 3) as i16;
        random_number_buffer >>= 1;
        position[3] += ((random_number_buffer & 0x0007) * 5) as i16;
        random_number_buffer >>= 3;
        sign[3] = (random_number_buffer & 0x0001) as i16;
        // Adaptive gain Ga (eq B.22, max 0.5).
        let ga: i16 = ((pseudo_random(random_generator_seed) & 0x1fff) << 1) as i16;

        for sample in gaussian_random_excitation.iter_mut() {
            let mut tmp_buffer: i32 = 0;
            for _ in 0..12 {
                tmp_buffer = add32(
                    tmp_buffer,
                    pseudo_random(random_generator_seed) as i16 as i32,
                );
            }
            *sample = shr32(tmp_buffer, 7) as i16;
            eg = mac16_16(eg, *sample, *sample);
        }

        gg = mult16_32_q15(GAUSSIAN_EXCITATION_COEFF_FACTOR, g729_inv_sqrt_q0q31(eg));
        gg = mult16_32_q15(target_gain, gg);

        for sample in gaussian_random_excitation.iter_mut() {
            if *sample < 0 {
                *sample = (-saturate(pshr(mult16_32_q15(-*sample, gg), 2), MAXINT16 as i32)) as i16;
            } else {
                *sample = pshr(mult16_32_q15(*sample, gg), 2) as i16;
            }
        }

        compute_adaptative_codebook_vector(
            excitation_vector,
            frac_pitch_delay,
            int_pitch_delay,
            base,
        );

        for i in 0..L_SUBFRAME {
            excitation_vector[base + i] = saturate(
                mult16_16_p15(excitation_vector[base + i], ga),
                MAXINT16 as i32,
            ) as i16;
        }

        for i in 0..L_SUBFRAME {
            excitation_vector[base + i] = saturate(
                add32(
                    excitation_vector[base + i] as i32,
                    gaussian_random_excitation[i] as i32,
                ),
                MAXINT16 as i32,
            ) as i16;
        }

        for i in 0..L_SUBFRAME {
            ea = mac16_16(ea, excitation_vector[base + i], excitation_vector[base + i]);
        }

        ei = 0;
        for i in 0..4 {
            if sign[i] == 0 {
                ei = sub32(ei, excitation_vector[base + position[i] as usize] as i32);
            } else {
                ei = add32(ei, excitation_vector[base + position[i] as usize] as i32);
            }
        }

        let k: i32 = mult16_32(
            target_gain,
            shr32(mult16_16(L_SUBFRAME as i16, target_gain), 3),
        );

        delta = (ei as i64).wrapping_mul(ei as i64)
            + (((k as i64).wrapping_sub((ea as i64) << 3)) >> 1);

        if delta < 0 {
            excitation_vector[base..base + L_SUBFRAME]
                .copy_from_slice(&gaussian_random_excitation[..L_SUBFRAME]);

            ei = 0;
            for i in 0..4 {
                if sign[i] == 0 {
                    ei = sub32(ei, excitation_vector[base + position[i] as usize] as i32);
                } else {
                    ei = add32(ei, excitation_vector[base + position[i] as usize] as i32);
                }
            }
            delta = (ei as i64).wrapping_mul(ei as i64) + mult16_32_p15(COEFF_K, k) as i64;
        }

        while delta >= 0x0000000080000000i64 {
            delta >>= 1;
            delta_scale_factor += 1;
        }
        if delta_scale_factor % 2 == 1 {
            delta >>= 1;
            delta_scale_factor += 1;
        }

        let delta = g729_sqrt_q0q7(delta as u32);

        ei = svshr32(ei, delta_scale_factor as i32 / 2 - 7);

        gf = sub32(delta, ei);
        let x2: i32 = -add32(delta, ei);
        if abs(x2) < abs(gf) {
            gf = x2;
        }
        gf = svshr32(gf, 2 + 7 - delta_scale_factor as i32 / 2);

        for i in 0..4 {
            let idx = base + position[i] as usize;
            if sign[i] == 0 {
                excitation_vector[idx] = sub32(excitation_vector[idx] as i32, gf) as i16;
            } else {
                excitation_vector[idx] = add32(excitation_vector[idx] as i32, gf) as i16;
            }
        }
    }
}

/// Decode a SID or missing frame and update the comfort-noise state (B4.4/B4.5).
#[allow(clippy::too_many_arguments)]
pub fn decode_sid_frame(
    cng_channel_context: &mut CngChannelContext,
    previous_frame_is_active_flag: u8,
    bit_stream: Option<&[u8]>,
    bit_stream_length: u8,
    excitation_vector: &mut [i16],
    previous_q_lsp: &mut [i16; NB_LSP_COEFF],
    lp: &mut [i16],
    pseudo_random_seed: &mut u16,
    previous_l_code_word: &mut [[i16; NB_LSP_COEFF]; MA_MAX_K],
    rfc3389_payload_flag: u8,
) {
    let mut interpolated_q_lsp = [0i16; NB_LSP_COEFF];

    if let Some(bit_stream) = bit_stream {
        if rfc3389_payload_flag != 0 {
            let mut lp_coefficients = [0i32; NB_LSP_COEFF + 1];
            let mut lp_coefficients_q12 = [0i16; NB_LSP_COEFF];
            let mut previous_iteration_lp_coefficients = [0i32; NB_LSP_COEFF + 1];
            let mut cn_filter_order = bit_stream_length.saturating_sub(1) as usize;
            let mut k = [0i16; NB_LSP_COEFF];

            if cn_filter_order > NB_LSP_COEFF {
                cn_filter_order = NB_LSP_COEFF;
            }

            let byte0 = bit_stream[0] as i32;
            let mut received_sid_gain_log = -byte0 + 90;
            if received_sid_gain_log > 66 {
                received_sid_gain_log = 66;
            }
            received_sid_gain_log = mult16_16(received_sid_gain_log as i16, 680);
            let received_sid_gain_energy = g729_exp2_q11q16(received_sid_gain_log as i16);
            if received_sid_gain_energy > 0 {
                cng_channel_context.received_sid_gain =
                    shr32(g729_sqrt_q0q7(received_sid_gain_energy as u32), 12) as i16;
                if cng_channel_context.received_sid_gain < SID_GAIN_CODEBOOK[0] {
                    cng_channel_context.received_sid_gain = SID_GAIN_CODEBOOK[0];
                }
            } else {
                cng_channel_context.received_sid_gain = SID_GAIN_CODEBOOK[0];
            }

            for (i, k_i) in k.iter_mut().take(cn_filter_order).enumerate() {
                let b = bit_stream.get(i + 1).copied().unwrap_or(0);
                *k_i = mult16_16(add16(b as i16, 127), 258) as i16;
            }
            for k_i in k
                .iter_mut()
                .skip(cn_filter_order)
                .take(NB_LSP_COEFF - cn_filter_order)
            {
                *k_i = 0;
            }

            // Rebuild the LP coefficients (G.711 Appendix II 5.2.1.3).
            lp_coefficients[0] = ONE_IN_Q27;
            lp_coefficients[1] = -shl(k[0] as i32, 12);
            for i in 2..NB_LSP_COEFF + 1 {
                previous_iteration_lp_coefficients[1..i].copy_from_slice(&lp_coefficients[1..i]);
                lp_coefficients[i] = -shl(k[i - 1] as i32, 16);
                for j in 1..i {
                    lp_coefficients[j] = mac32_32_q31(
                        lp_coefficients[j],
                        lp_coefficients[i],
                        previous_iteration_lp_coefficients[i - j],
                    );
                }
                lp_coefficients[i] = shr(lp_coefficients[i], 4);
            }

            for i in 0..NB_LSP_COEFF {
                lp_coefficients_q12[i] =
                    saturate(pshr(lp_coefficients[i + 1], 15), MAXINT16 as i32) as i16;
            }

            if !lp2lsp_conversion(&lp_coefficients_q12, &mut cng_channel_context.q_lsp) {
                cng_channel_context.q_lsp.copy_from_slice(previous_q_lsp);
            }
        } else {
            // Regular G.729 SID payload on 2 bytes.
            let mut current_q_lsf = [0i16; NB_LSP_COEFF];
            let byte0 = bit_stream[0];
            let byte1 = bit_stream.get(1).copied().unwrap_or(0);
            let l0 = ((byte0 >> 7) & 0x01) as usize;
            let l1_index = ((byte0 >> 2) & 0x1F) as usize;
            let l2_index = (((byte0 & 0x03) << 2) | ((byte1 >> 6) & 0x03)) as usize;

            cng_channel_context.received_sid_gain =
                SID_GAIN_CODEBOOK[((byte1 >> 1) & 0x1F) as usize];

            for i in 0..NB_LSP_COEFF / 2 {
                current_q_lsf[i] = add16(
                    L1[L1_SUBSET_INDEX[l1_index]][i],
                    L2L3[L2_SUBSET_INDEX[l2_index]][i],
                );
            }
            for i in NB_LSP_COEFF / 2..NB_LSP_COEFF {
                current_q_lsf[i] = add16(
                    L1[L1_SUBSET_INDEX[l1_index]][i],
                    L2L3[L3_SUBSET_INDEX[l2_index]][i],
                );
            }
            compute_q_lsf(
                &mut current_q_lsf,
                previous_l_code_word,
                l0,
                &NOISE_MA_PREDICTOR,
                &NOISE_MA_PREDICTOR_SUM,
            );

            for (q_lsp, &q_lsf) in cng_channel_context
                .q_lsp
                .iter_mut()
                .zip(current_q_lsf.iter())
            {
                *q_lsp = g729_cos_q13q15(q_lsf);
            }
        }
    }

    interpolate_q_lsp(
        previous_q_lsp,
        &cng_channel_context.q_lsp,
        &mut interpolated_q_lsp,
    );
    previous_q_lsp.copy_from_slice(&cng_channel_context.q_lsp);

    q_lsp_2_lp(&interpolated_q_lsp, &mut lp[0..NB_LSP_COEFF]);
    q_lsp_2_lp(&cng_channel_context.q_lsp, &mut lp[NB_LSP_COEFF..]);

    // Target gain smoothing (eq B.19).
    if previous_frame_is_active_flag != 0 {
        cng_channel_context.smoothed_sid_gain = cng_channel_context.received_sid_gain;
    } else {
        cng_channel_context.smoothed_sid_gain = sub16(
            cng_channel_context.smoothed_sid_gain,
            shr16(cng_channel_context.smoothed_sid_gain, 3),
        );
        cng_channel_context.smoothed_sid_gain = add16(
            cng_channel_context.smoothed_sid_gain,
            shr16(cng_channel_context.received_sid_gain, 3),
        );
    }

    compute_comfort_noise_excitation_vector(
        cng_channel_context.smoothed_sid_gain,
        pseudo_random_seed,
        excitation_vector,
    );
}
