// G.729 Annex B voice activity detection (VAD), ported from bcg729/src/vad.c.

use crate::g729::basic_operations::*;
use crate::g729::codebooks::LOW_BAND_FILTER;
use crate::g729::fixed_point_math::g729_log2_q0q16;
use crate::g729::ld8k::*;
use crate::g729::utils::get_min_in_array;

const VOICE: u8 = 1;
const NOISE: u8 = 0;

/// VAD state; kept inline (no heap) so it works in `no_std`.
pub struct VadChannelContext {
    ef_buffer: [Word16; N0],
    frame_count: usize,
    update_count: u32,
    nb_valid_init_frame: u32,
    init_ef_sum: Word32,
    init_zc_sum: Word32,
    init_lsf_sum: [Word32; NB_LSP_COEFF],
    mean_ef: Word16,
    mean_el: Word16,
    mean_zc: Word16,
    mean_lsf: [Word16; NB_LSP_COEFF],
    svd_m1: u8,
    svd_m2: u8,
    count_inert: u8,
    second_stage_vad_smoothing_flag: u8,
    smoothing_counter: u8,
    previous_frame_ef: Word16,
    noise_continuity_counter: u8,
}

impl Default for VadChannelContext {
    fn default() -> Self {
        Self::new()
    }
}

impl VadChannelContext {
    pub fn new() -> Self {
        VadChannelContext {
            ef_buffer: [MAX_16; N0],
            frame_count: 0,
            update_count: 0,
            nb_valid_init_frame: 0,
            init_ef_sum: 0,
            init_zc_sum: 0,
            init_lsf_sum: [0; NB_LSP_COEFF],
            mean_ef: 0,
            mean_el: 0,
            mean_zc: 0,
            mean_lsf: [0; NB_LSP_COEFF],
            svd_m1: VOICE,
            svd_m2: VOICE,
            count_inert: 0,
            second_stage_vad_smoothing_flag: 1,
            smoothing_counter: 0,
            previous_frame_ef: 0,
            noise_continuity_counter: 0,
        }
    }
}

/// Annex B multi-boundary initial decision (B3.5).
fn multi_boundary_initial_voice_activity_decision(
    delta_s: Word32,
    delta_ef: Word16,
    delta_el: Word16,
    delta_zc: Word16,
) -> u8 {
    let delta_ef32 = mult16_16(10, delta_ef); // Q11
    let delta_el32 = mult16_16(10, delta_el); // Q11
    let delta_s = mult16_32(830, delta_s); // Q28 (830 = 1/(4*pi^2) in Q15)

    if delta_s > add32(mult16_32_q12(delta_zc, 58720), 228170) {
        return VOICE;
    }
    if delta_s > add32(mult16_32_q12(delta_zc, -152520), 311141) {
        return VOICE;
    }
    if delta_ef32 < add32(mult16_32_q15(delta_zc, -51200), -10240) {
        return VOICE;
    }
    if delta_ef32 < add32(mult16_32_q15(delta_zc, 40960), -12288) {
        return VOICE;
    }
    if delta_ef32 < -9626 {
        return VOICE;
    }
    if delta_ef32 < add32(mult16_32_q12(275, delta_s), -24986) {
        return VOICE;
    }
    if delta_s > 241592 {
        return VOICE;
    }
    if delta_ef32 < add32(mult16_32_q15(delta_zc, 51200), -14336) {
        return VOICE;
    }
    if delta_ef32 < add32(mult16_32_q15(delta_zc, -59578), -9868) {
        return VOICE;
    }
    if delta_ef32 < -10854 {
        return VOICE;
    }
    if delta_el32 < add32(mult16_32_q13(875, delta_s), -31744) {
        return VOICE;
    }
    // The reference uses deltaEl32 on both sides here; kept verbatim.
    if delta_el32 > add32(mult16_32_q15(30427, delta_el32), 2341) {
        return VOICE;
    }
    if delta_el32 < add32(mult16_32_q14(-24576, delta_el32), -18432) {
        return VOICE;
    }
    if delta_el32 < add32(mult16_32_q15(23406, delta_el32), -4389) {
        return VOICE;
    }

    NOISE
}

/// Annex B VAD decision: returns 1 for an active voice frame, 0 for silence.
///
/// `signal_current_frame` must start one element before the current frame
/// (the decision accesses indices `[-1, L_FRAME[`).
pub fn bcg729_vad(
    ctx: &mut VadChannelContext,
    reflection_coefficient: Word32,
    lsf_coefficients: &[Word16; NB_LSP_COEFF],
    auto_correlation_coefficients: &[Word32],
    auto_correlation_coefficients_scale: i8,
    signal_current_frame: &[Word16],
) -> u8 {
    let ef: Word16;
    let emin: Word16;
    let el: Word16;
    let zc: Word16;
    let delta_s: Word32;
    let mut acc: Word32;
    let mut ivd: u8;

    // Full-band energy Ef/10 (B3.1), Q11.
    acc = sub32(
        g729_log2_q0q16(auto_correlation_coefficients[0].wrapping_add(1)),
        (auto_correlation_coefficients_scale as i32) << 16,
    );
    acc = shr32(sub32(acc, LOG2_240_Q16), 1);
    acc = mult16_32_q15(INV_LOG2_10_Q15, acc);
    ef = pshr(acc, 4) as Word16;

    ctx.ef_buffer[ctx.frame_count % N0] = ef;

    // Low-band energy El/10 (B3.1), Q11.
    acc = mult16_32_q15(LOW_BAND_FILTER[0], auto_correlation_coefficients[0]);
    for i in 1..NB_LSP_COEFF + 3 {
        acc = mac16_32_q14(acc, LOW_BAND_FILTER[i], auto_correlation_coefficients[i]);
    }
    if acc <= 0 {
        acc = 1;
    }
    acc = sub32(
        g729_log2_q0q16(acc),
        (auto_correlation_coefficients_scale as i32) << 16,
    );
    acc = shr32(sub32(acc, LOG2_240_Q16), 1);
    acc = mult16_32_q15(INV_LOG2_10_Q15, acc);
    el = pshr(acc, 4) as Word16;

    // Zero-crossing rate, Q15 (1/80 per crossing).
    let mut zc_acc: Word16 = 0;
    for i in 0..L_FRAME {
        if mult16_16(signal_current_frame[i], signal_current_frame[i + 1]) < 0 {
            zc_acc = add16(zc_acc, 410);
        }
    }
    zc = zc_acc;

    // B3.2: initialisation of the background-noise running averages.
    if ctx.frame_count == NI {
        if ctx.nb_valid_init_frame > 0 {
            let mean_en = div32(ctx.init_ef_sum, ctx.nb_valid_init_frame as i32) as Word16;
            ctx.mean_zc = div32(ctx.init_zc_sum, ctx.nb_valid_init_frame as i32) as Word16;
            for i in 0..NB_LSP_COEFF {
                ctx.mean_lsf[i] =
                    div32(ctx.init_lsf_sum[i], ctx.nb_valid_init_frame as i32) as Word16;
            }
            ctx.mean_ef = sub16(mean_en, 2048);
            ctx.mean_el = sub16(mean_en, 2458);
        } else {
            ctx.frame_count = 0;
        }
    }

    if ctx.frame_count < NI {
        if ef < 3072 {
            ivd = NOISE;
        } else {
            ivd = VOICE;
            ctx.nb_valid_init_frame += 1;
            ctx.init_ef_sum = add32(ctx.init_ef_sum, ef as Word32);
            ctx.init_zc_sum = add32(ctx.init_zc_sum, zc as Word32);
            for i in 0..NB_LSP_COEFF {
                ctx.init_lsf_sum[i] = add32(ctx.init_lsf_sum[i], lsf_coefficients[i] as Word32);
            }
        }

        ctx.frame_count += 1;
        ctx.previous_frame_ef = ef;
        ctx.svd_m2 = ctx.svd_m1;
        ctx.svd_m1 = ivd;

        return ivd;
    }

    emin = get_min_in_array(&ctx.ef_buffer, N0); // B3.3

    // B3.4: spectral distortion and energy/zero-crossing deviations.
    let mut delta_s_acc: Word32 = 0;
    for i in 0..NB_LSP_COEFF {
        let acc16 = sub16(lsf_coefficients[i], ctx.mean_lsf[i]);
        delta_s_acc = mac16_16_q13(delta_s_acc, acc16, acc16);
    }
    delta_s = delta_s_acc;

    let delta_ef = sub16(ctx.mean_ef, ef);
    let delta_el = sub16(ctx.mean_el, el);
    let delta_zc = sub16(ctx.mean_zc, zc);

    if ef < 3072 {
        ivd = NOISE;
    } else {
        ivd = multi_boundary_initial_voice_activity_decision(delta_s, delta_ef, delta_el, delta_zc);
    }

    // B3.6: voice activity decision smoothing.
    if ivd == VOICE {
        ctx.count_inert = 0;
    }

    if ivd == NOISE && ctx.count_inert < 6 {
        ctx.count_inert += 1;
        ivd = VOICE;
    }

    if ivd == NOISE && ctx.svd_m1 != 0 && delta_ef > 410 && ef > 3072 {
        ivd = VOICE;
    }

    if ctx.second_stage_vad_smoothing_flag == 1
        && ivd == NOISE
        && ctx.svd_m1 == VOICE
        && ctx.svd_m2 == VOICE
        && abs16(sub16(ef, ctx.previous_frame_ef)) <= 614
    {
        ivd = VOICE;
        ctx.smoothing_counter += 1;
        if ctx.smoothing_counter <= 4 {
            ctx.second_stage_vad_smoothing_flag = 1;
        } else {
            ctx.second_stage_vad_smoothing_flag = 0;
            ctx.smoothing_counter = 0;
        }
    } else {
        ctx.second_stage_vad_smoothing_flag = 1;
    }

    if ivd == NOISE {
        ctx.noise_continuity_counter += 1;
    }

    if ivd == VOICE && ctx.noise_continuity_counter > 10 && sub16(ef, ctx.previous_frame_ef) <= 614
    {
        ivd = NOISE;
        ctx.noise_continuity_counter = 0;
        ctx.count_inert = 6;
    }

    if ivd == VOICE {
        ctx.noise_continuity_counter = 0;
    }

    // B3.7: update of the running averages.
    if sub16(ef, 614) < ctx.mean_ef && reflection_coefficient < 1610612736 {
        let (beta_e, beta_e_complement, beta_zc, beta_zc_complement, beta_lsf, beta_lsf_complement);
        ctx.update_count += 1;
        let update_count = ctx.update_count;
        if update_count < 20 {
            beta_e = 24576;
            beta_e_complement = 8192;
            beta_zc = 26214;
            beta_zc_complement = 6554;
            beta_lsf = 19661;
            beta_lsf_complement = 13107;
        } else if update_count < 30 {
            beta_e = 31130;
            beta_e_complement = 1638;
            beta_zc = 30147;
            beta_zc_complement = 2621;
            beta_lsf = 21299;
            beta_lsf_complement = 11469;
        } else if update_count < 40 {
            beta_e = 31785;
            beta_e_complement = 983;
            beta_zc = 30802;
            beta_zc_complement = 1966;
            beta_lsf = 22938;
            beta_lsf_complement = 9830;
        } else if update_count < 50 {
            beta_e = 32440;
            beta_e_complement = 328;
            beta_zc = 31457;
            beta_zc_complement = 1311;
            beta_lsf = 24756;
            beta_lsf_complement = 8192;
        } else if update_count < 60 {
            beta_e = 32604;
            beta_e_complement = 164;
            beta_zc = 32440;
            beta_zc_complement = 328;
            beta_lsf = 24576;
            beta_lsf_complement = 8192;
        } else {
            beta_e = 32702;
            beta_e_complement = 66;
            beta_zc = 32604;
            beta_zc_complement = 164;
            beta_lsf = 24576;
            beta_lsf_complement = 8192;
        }

        ctx.mean_ef = add16(
            mult16_16_q15(ctx.mean_ef, beta_e) as Word16,
            mult16_16_q15(ef, beta_e_complement) as Word16,
        );
        ctx.mean_el = add16(
            mult16_16_q15(ctx.mean_el, beta_e) as Word16,
            mult16_16_q15(el, beta_e_complement) as Word16,
        );
        ctx.mean_zc = add16(
            mult16_16_q15(ctx.mean_zc, beta_zc) as Word16,
            mult16_16_q15(zc, beta_zc_complement) as Word16,
        );
        for i in 0..NB_LSP_COEFF {
            ctx.mean_lsf[i] = add16(
                mult16_16_q15(ctx.mean_lsf[i], beta_lsf) as Word16,
                mult16_16_q15(lsf_coefficients[i], beta_lsf_complement) as Word16,
            );
        }
    }

    if ctx.frame_count > N0
        && ((ctx.mean_ef < emin && delta_s < 819) || (ctx.mean_ef > add16(emin, 2048)))
    {
        ctx.mean_ef = emin;
        ctx.update_count = 0;
    }

    ctx.frame_count += 1;
    ctx.previous_frame_ef = ef;
    ctx.svd_m2 = ctx.svd_m1;
    ctx.svd_m1 = ivd;

    ivd
}
