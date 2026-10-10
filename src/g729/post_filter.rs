use crate::g729::basic_operations::*;
use crate::g729::fixed_point_math::g729_sqrt_q0q7;
use crate::g729::ld8k::*;
use crate::g729::lp_synthesis_filter::lp_synthesis_filter;
use crate::g729::utils::count_leading_zeros;

pub fn init_post_filter() -> (
    [i16; MAXIMUM_INT_PITCH_DELAY + L_FRAME],
    [i16; MAXIMUM_INT_PITCH_DELAY + L_FRAME],
    [i16; 1 + L_SUBFRAME],
    [i16; NB_LSP_COEFF + L_SUBFRAME],
    i16,
) {
    (
        [0; MAXIMUM_INT_PITCH_DELAY + L_FRAME],
        [0; MAXIMUM_INT_PITCH_DELAY + L_FRAME],
        [0; 1 + L_SUBFRAME],
        [0; NB_LSP_COEFF + L_SUBFRAME],
        /* previous gain for adaptative gain control: 1 in Q12 */
        4096,
    )
}

/*****************************************************************************/
/* postFilter: filter the reconstructed speech according to spec A.4.2       */

/*      -(i) LPCoefficients: 10 LP coeff for current subframe in Q12         */
/*      -(i) reconstructedSpeech: output of LP Synthesis, 50 values in Q0    */

/*      -(i) intPitchDelay: the integer part of Pitch Delay in Q0            */

/*      -(o) postFilteredSignal: 40 values in Q0                             */

/*****************************************************************************/
#[allow(clippy::too_many_arguments)]
pub fn post_filter(
    residual_signal_buffer: &mut [i16],
    scaled_residual_signal_buffer: &mut [i16],
    long_term_filtered_residual_signal_buffer: &mut [i16],
    short_term_filtered_residual_signal_buffer: &mut [i16],
    previous_adaptative_gain: &mut i16,
    lp_coefficients: &[i16],
    reconstructed_speech: &[i16],
    mut int_pitch_delay: i16,
    subframe_index: usize,
    post_filtered_signal: &mut [i16],
) {
    let mut lp_gamma_n_coefficients = [0_i16; NB_LSP_COEFF]; /* in Q12 */
    let mut correlation_max: i32 = MININT32;
    let mut best_int_pitch_delay: i16 = 0;
    let mut residual_signal_energy: i32 = 0; /* in Q-4 */
    let mut delayed_residual_signal_energy: i32 = 0; /* in Q-4 */
    let mut maximum_three: i32;
    let mut correlation_max_word16: i16 = 0;
    let mut residual_signal_energy_word16: i16 = 0;
    let mut delayed_residual_signal_energy_word16: i16 = 0;
    let mut lp_gamma_d_coefficients = [0_i16; NB_LSP_COEFF]; /* in Q12 */
    let mut hf = [0_i16; 22]; /* the truncated impulse response to short term filter Hf in Q12 */
    let mut rh1: i32;
    let mut tilt_compensated_signal = [0_i16; L_SUBFRAME]; /* in Q0 */
    let mut gain_scaling_factor: i16; /* in Q12 */
    let mut short_term_filtered_residual_signal_square_sum: u32 = 0;

    /********************************************************************/

    /********************************************************************/
    /*** Compute LPGammaN and LPGammaD coefficients : LPGamma[0] = LP[0]*Gamma^(i+1) (i=0..9) ***/
    /* GAMMA_XX constants are in Q15 */
    lp_gamma_n_coefficients[0] = mult16_16_p15(lp_coefficients[0], GAMMA_N1) as i16;
    lp_gamma_n_coefficients[1] = mult16_16_p15(lp_coefficients[1], GAMMA_N2) as i16;
    lp_gamma_n_coefficients[2] = mult16_16_p15(lp_coefficients[2], GAMMA_N3) as i16;
    lp_gamma_n_coefficients[3] = mult16_16_p15(lp_coefficients[3], GAMMA_N4) as i16;
    lp_gamma_n_coefficients[4] = mult16_16_p15(lp_coefficients[4], GAMMA_N5) as i16;
    lp_gamma_n_coefficients[5] = mult16_16_p15(lp_coefficients[5], GAMMA_N6) as i16;
    lp_gamma_n_coefficients[6] = mult16_16_p15(lp_coefficients[6], GAMMA_N7) as i16;
    lp_gamma_n_coefficients[7] = mult16_16_p15(lp_coefficients[7], GAMMA_N8) as i16;
    lp_gamma_n_coefficients[8] = mult16_16_p15(lp_coefficients[8], GAMMA_N9) as i16;
    lp_gamma_n_coefficients[9] = mult16_16_p15(lp_coefficients[9], GAMMA_N10) as i16;

    /*** Compute the residual signal as described in spec 4.2.1 eq79 ***/

    let residual_signal_offset = MAXIMUM_INT_PITCH_DELAY + subframe_index;

    for i in 0..L_SUBFRAME {
        let mut acc = sshl(reconstructed_speech[NB_LSP_COEFF + i] as i32, 12); /* reconstructedSpeech in Q0 shifted to set acc in Q12 */
        for j in 0..NB_LSP_COEFF {
            acc = mac16_16(
                acc,
                lp_gamma_n_coefficients[j],
                reconstructed_speech[NB_LSP_COEFF + i - j - 1],
            );
        }
        residual_signal_buffer[residual_signal_offset + i] =
            saturate(pshr(acc, 12), MAX_INT16 as i32) as i16; /* shift back acc to Q0 and saturate it to avoid overflow when going back to 16 bits */
        scaled_residual_signal_buffer[residual_signal_offset + i] =
            pshr(residual_signal_buffer[residual_signal_offset + i] as i32, 2) as i16;
        /* shift acc to Q-2 and saturate it to get the scaled version of the signal */
    }

    /*** Compute the maximum correlation on scaledResidualSignal delayed by intPitchDelay +/- 3 to get the best delay. Spec 4.2.1 eq80 ***/
    /* using a scaled(Q-2) signals gives correlation in Q-4. */
    if int_pitch_delay > (MAXIMUM_INT_PITCH_DELAY - 3) as i16 {
        int_pitch_delay = (MAXIMUM_INT_PITCH_DELAY - 3) as i16;
    }

    for i in (int_pitch_delay - 3)..=(int_pitch_delay + 3) {
        let mut correlation: i32 = 0;

        let delayed_residual_signal_offset = residual_signal_offset as isize - i as isize;

        for j in 0..L_SUBFRAME {
            correlation = mac16_16(
                correlation,
                scaled_residual_signal_buffer
                    [(delayed_residual_signal_offset + j as isize) as usize],
                scaled_residual_signal_buffer[residual_signal_offset + j],
            );
        }

        if correlation > correlation_max {
            correlation_max = correlation;
            best_int_pitch_delay = i;
        }
    }

    if correlation_max < 0 {
        correlation_max = 0;
    }

    /*** Compute the signal energy ∑r(n)*r(n) and delayed signal energy ∑rk(n)*rk(n) which shall be used to compute gl spec 4.2.1 eq81, eq 82 and eq83 ***/

    let delayed_residual_signal_offset =
        residual_signal_offset as isize - best_int_pitch_delay as isize;

    for i in 0..L_SUBFRAME {
        residual_signal_energy = mac16_16(
            residual_signal_energy,
            scaled_residual_signal_buffer[residual_signal_offset + i],
            scaled_residual_signal_buffer[residual_signal_offset + i],
        );
        delayed_residual_signal_energy = mac16_16(
            delayed_residual_signal_energy,
            scaled_residual_signal_buffer[(delayed_residual_signal_offset + i as isize) as usize],
            scaled_residual_signal_buffer[(delayed_residual_signal_offset + i as isize) as usize],
        );
    }

    /*** Scale correlationMax, residualSignalEnergy and delayedResidualSignalEnergy to the best fit on 16 bits ***/

    maximum_three = correlation_max;
    if maximum_three < residual_signal_energy {
        maximum_three = residual_signal_energy;
    }
    if maximum_three < delayed_residual_signal_energy {
        maximum_three = delayed_residual_signal_energy;
    }

    if maximum_three > 0 {
        let leading_zeros = count_leading_zeros(maximum_three) as i16;
        if leading_zeros < 16 {
            correlation_max_word16 = shr32(correlation_max, (16 - leading_zeros) as u32) as i16;
            residual_signal_energy_word16 =
                shr32(residual_signal_energy, (16 - leading_zeros) as u32) as i16;
            delayed_residual_signal_energy_word16 =
                shr32(delayed_residual_signal_energy, (16 - leading_zeros) as u32) as i16;
        } else {
            correlation_max_word16 = correlation_max as i16;
            residual_signal_energy_word16 = residual_signal_energy as i16;
            delayed_residual_signal_energy_word16 = delayed_residual_signal_energy as i16;
        }
    }

    /* g = gl/2 (as γp=0.5)= (eq83) correlationMax/(2*delayedResidualSignalEnergy) */

    /*** eq82 -> (correlationMax^2)/(residualSignalEnergy*delayedResidualSignalEnergy)<0.5 ***/

    if (mult16_16(correlation_max_word16, correlation_max_word16) < shr(mult16_16(residual_signal_energy_word16, delayed_residual_signal_energy_word16), 1)) /* eq82 */
        || ((correlation_max_word16 == 0) && (delayed_residual_signal_energy_word16 == 0))
    {
        long_term_filtered_residual_signal_buffer[1..1 + L_SUBFRAME].copy_from_slice(
            &residual_signal_buffer[residual_signal_offset..residual_signal_offset + L_SUBFRAME],
        );
    } else {
        /* eq82 gives long term filter enabled, */
        let g0: i16;
        let g1: i16;
        /* eq83: gl = correlationMax/delayedResidualSignalEnergy bounded in ]0,1] */

        if correlation_max > delayed_residual_signal_energy {
            g0 = 21845; /* 2/3 in Q15 */
            g1 = 10923; /* 1/3 in Q15 */
        } else {
            g1 = div32(
                shl32(correlation_max_word16 as i32, 15),
                add32(
                    shl32(delayed_residual_signal_energy_word16 as i32, 1),
                    correlation_max_word16 as i32,
                ),
            ) as i16; /* g1 in Q15 */
            g0 = sub16(32767, g1); /* g0 = 1 - g1 in Q15 */
        }

        let delayed_residual_signal_offset =
            residual_signal_offset as isize - best_int_pitch_delay as isize;
        for i in 0..L_SUBFRAME {
            long_term_filtered_residual_signal_buffer[1 + i] = saturate(
                pshr(
                    add32(
                        mult16_16(g0, residual_signal_buffer[residual_signal_offset + i]),
                        mult16_16(
                            g1,
                            residual_signal_buffer
                                [(delayed_residual_signal_offset + i as isize) as usize],
                        ),
                    ),
                    15,
                ),
                MAX_INT16 as i32,
            ) as i16;
        }
    }

    /********************************************************************/

    /********************************************************************/

    /* compute hf the truncated (to 22 coefficients) impulse response of the filter A(z/γn)/A(z/γd) described in spec 4.2.2 eq84 */

    /* GAMMA_XX constants are in Q15 */
    lp_gamma_d_coefficients[0] = mult16_16_p15(lp_coefficients[0], GAMMA_D1) as i16;
    lp_gamma_d_coefficients[1] = mult16_16_p15(lp_coefficients[1], GAMMA_D2) as i16;
    lp_gamma_d_coefficients[2] = mult16_16_p15(lp_coefficients[2], GAMMA_D3) as i16;
    lp_gamma_d_coefficients[3] = mult16_16_p15(lp_coefficients[3], GAMMA_D4) as i16;
    lp_gamma_d_coefficients[4] = mult16_16_p15(lp_coefficients[4], GAMMA_D5) as i16;
    lp_gamma_d_coefficients[5] = mult16_16_p15(lp_coefficients[5], GAMMA_D6) as i16;
    lp_gamma_d_coefficients[6] = mult16_16_p15(lp_coefficients[6], GAMMA_D7) as i16;
    lp_gamma_d_coefficients[7] = mult16_16_p15(lp_coefficients[7], GAMMA_D8) as i16;
    lp_gamma_d_coefficients[8] = mult16_16_p15(lp_coefficients[8], GAMMA_D9) as i16;
    lp_gamma_d_coefficients[9] = mult16_16_p15(lp_coefficients[9], GAMMA_D10) as i16;

    hf[0] = 4096; /* 1 in Q12 as LPGammaNCoefficients and LPGammaDCoefficient doesn't contain the first element which is 1 and past values of hf are 0 */
    for i in 1..11 {
        let mut acc = sshl(lp_gamma_n_coefficients[i - 1] as i32, 12);
        for j in 0..NB_LSP_COEFF {
            if j < i {
                acc = msu16_16(acc, lp_gamma_d_coefficients[j], hf[i - j - 1]);
            }
        }
        hf[i] = saturate(pshr(acc, 12), MAX_INT16 as i32) as i16; /* get result back in Q12 and saturate on 16 bits */
    }
    for i in 11..22 {
        let mut acc: i32 = 0;
        for j in 0..NB_LSP_COEFF {
            acc = msu16_16(acc, lp_gamma_d_coefficients[j], hf[i - j - 1]);
        }
        hf[i] = saturate(pshr(acc, 12), MAX_INT16 as i32) as i16; /* get result back in Q12 and saturate on 16 bits */
    }

    /* hf is then used to compute k'1 spec 4.2.3 eq87: k'1 = -rh1/rh0 */

    rh1 = mult16_16(hf[0], hf[1]);
    for i in 1..21 {
        rh1 = mac16_16(rh1, hf[i], hf[i + 1]); /* rh1 in Q24 */
    }

    if rh1 < 0 {
        tilt_compensated_signal[..L_SUBFRAME]
            .copy_from_slice(&long_term_filtered_residual_signal_buffer[1..1 + L_SUBFRAME]);
    } else {
        let mut rh0 = mult16_16(hf[0], hf[0]);
        for &h in &hf[1..22] {
            rh0 = mac16_16(rh0, h, h); /* rh0 in Q24 */
        }
        rh1 = mult16_32_q15(GAMMA_T, rh1); /* GAMMA_T in Q15, rh1 in Q24*/
        let tilt_compensation_gain: i16 =
            saturate(div32(rh1, pshr(rh0, 12)), MAX_INT16 as i32) as i16;

        /* compute filter Ht (spec A.4.2.3 eqA14) = 1 + gain*z(-1) */
        for i in 0..L_SUBFRAME {
            tilt_compensated_signal[i] = msu16_16_q12(
                long_term_filtered_residual_signal_buffer[1 + i] as i32,
                tilt_compensation_gain,
                long_term_filtered_residual_signal_buffer[1 + i - 1],
            ) as i16;
        }
    }

    long_term_filtered_residual_signal_buffer[0] =
        long_term_filtered_residual_signal_buffer[1 + L_SUBFRAME - 1];

    /********************************************************************/

    /********************************************************************/

    lp_synthesis_filter(
        &tilt_compensated_signal,
        &lp_gamma_d_coefficients,
        short_term_filtered_residual_signal_buffer,
    );

    for i in 0..NB_LSP_COEFF {
        short_term_filtered_residual_signal_buffer[i] =
            short_term_filtered_residual_signal_buffer[L_SUBFRAME + i];
    }

    /********************************************************************/
    /* Adaptive Gain Control spec A.4.2.4                               */

    /********************************************************************/

    /*** compute G(gain scaling factor) according to eqA15 : G = Sqrt((∑s(n)^2)/∑sf(n)^2 ) ***/

    for i in 0..L_SUBFRAME {
        let val = short_term_filtered_residual_signal_buffer[NB_LSP_COEFF + i];
        let prod = mult16_16(val, val);
        let prod_q4 = shr(prod, 4);
        short_term_filtered_residual_signal_square_sum = uadd32(
            short_term_filtered_residual_signal_square_sum,
            prod_q4 as u32,
        );
    }

    /* the reset of previousAdaptativeGain is not mentionned in the spec but in ITU code only */
    if short_term_filtered_residual_signal_square_sum == 0 {
        *previous_adaptative_gain = 0;
        post_filtered_signal[..L_SUBFRAME].copy_from_slice(
            &short_term_filtered_residual_signal_buffer[NB_LSP_COEFF..NB_LSP_COEFF + L_SUBFRAME],
        );
    } else {
        let mut current_adaptative_gain: i16;

        let mut reconstructed_speech_square_sum: u32 = 0;
        for i in 0..L_SUBFRAME {
            let val = reconstructed_speech[NB_LSP_COEFF + i];
            let prod = mult16_16(val, val);
            let prod_q4 = shr(prod, 4);
            reconstructed_speech_square_sum =
                uadd32(reconstructed_speech_square_sum, prod_q4 as u32);
        }

        if reconstructed_speech_square_sum == 0 {
            gain_scaling_factor = 0;
        } else {
            let mut fraction_result: u32; /* stores  ∑s(n)^2)/∑sf(n)^2 in Q10 on a 32 bit unsigned */

            /* Compute ∑s(n)^2)/∑sf(n)^2  result shall be in Q10 */

            let numerator_shift = reconstructed_speech_square_sum.leading_zeros() as i16;
            reconstructed_speech_square_sum =
                ushl(reconstructed_speech_square_sum, numerator_shift as u32);

            /* normalise denominator to get the result directly in Q10 if possible */

            let scaled_short_term_filtered_residual_signal_square_sum: u32 =
                if 10 - numerator_shift >= 0 {
                    short_term_filtered_residual_signal_square_sum >> (10 - numerator_shift)
                } else {
                    short_term_filtered_residual_signal_square_sum << (numerator_shift - 10)
                };

            if scaled_short_term_filtered_residual_signal_square_sum == 0 {
                fraction_result = udiv32(
                    reconstructed_speech_square_sum,
                    short_term_filtered_residual_signal_square_sum,
                );

                if numerator_shift - 10 >= 0 {
                    fraction_result >>= numerator_shift - 10;
                } else {
                    fraction_result <<= 10 - numerator_shift;
                }
            } else {
                fraction_result = udiv32(
                    reconstructed_speech_square_sum,
                    scaled_short_term_filtered_residual_signal_square_sum,
                ); /* result in Q10 */
            }

            gain_scaling_factor =
                saturate(g729_sqrt_q0q7(fraction_result), MAX_INT16 as i32) as i16;

            /* multiply by 0.1 as described in spec A.4.2.4 */
            gain_scaling_factor = mult16_16_p15(gain_scaling_factor, 3277_i16) as i16;
            /* in Q12, 3277 = 0.1 in Q15*/
        }
        /* Compute the signal according to eq89 (spec 4.2.4 and section A4.2.4) */

        current_adaptative_gain = *previous_adaptative_gain;
        for i in 0..L_SUBFRAME {
            current_adaptative_gain = add16(
                gain_scaling_factor,
                mult16_16_p15(current_adaptative_gain, 29491_i16) as i16,
            ); /* 29492 = 0.9 in Q15, result in Q12 */
            post_filtered_signal[i] = mult16_16_q12(
                current_adaptative_gain,
                short_term_filtered_residual_signal_buffer[NB_LSP_COEFF + i],
            ) as i16;
        }
        *previous_adaptative_gain = current_adaptative_gain;
    }

    if subframe_index > 0 {
        for i in 0..MAXIMUM_INT_PITCH_DELAY {
            residual_signal_buffer[i] = residual_signal_buffer[L_FRAME + i];
            scaled_residual_signal_buffer[i] = scaled_residual_signal_buffer[L_FRAME + i];
        }
    }
}
