use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;
use crate::g729::lp_synthesis_filter::*;

/*****************************************************************************/
/* computeWeightedSpeech: compute wieghted speech according to spec A3.3.3   */

/*      -(i) qLPCoefficients: 20 coefficients(10 for each subframe) in Q12   */

/*           in Q12                                                          */

/*      -(o) LPResidualSignal: 80 values of residual signal in Q0            */

/*****************************************************************************/
#[cfg_attr(target_arch = "xtensa", inline(never))]
pub fn compute_weighted_speech(
    input_signal: &[i16],
    q_lp_coefficients: &[i16],
    weighted_q_lp_coefficients: &[i16],
    weighted_input_signal: &mut [i16],
    lp_residual_signal: &mut [i16],
) {
    /* algo as specified in A3.3.3: */

    let mut weighted_q_lp_low_pass_coefficients = [0; NB_LSP_COEFF]; /* in Q12 */

    /*** compute LPResisualSignal (spec A3.3.3 eqA.3) in Q0 ***/

    for i in 0..L_SUBFRAME {
        let mut acc = sshl(input_signal[NB_LSP_COEFF + i] as i32, 12); /* inputSignal in Q0 is shifted to set acc in Q12 */
        for j in 0..NB_LSP_COEFF {
            acc = mac16_16(
                acc,
                q_lp_coefficients[j],
                input_signal[NB_LSP_COEFF + i - j - 1],
            );
        }
        lp_residual_signal[i] = saturate(pshr(acc, 12), MAX_16 as i32) as i16;
        /* shift back acc to Q0 and saturate it to avoid overflow when going back to 16 bits */
    }

    for i in L_SUBFRAME..L_FRAME {
        let mut acc = sshl(input_signal[NB_LSP_COEFF + i] as i32, 12); /* inputSignal in Q0 is shifted to set acc in Q12 */
        for j in 0..NB_LSP_COEFF {
            acc = mac16_16(
                acc,
                q_lp_coefficients[NB_LSP_COEFF + j],
                input_signal[NB_LSP_COEFF + i - j - 1],
            );
        }
        lp_residual_signal[i] = saturate(pshr(acc, 12), MAX_16 as i32) as i16;
        /* shift back acc to Q0 and saturate it to avoid overflow when going back to 16 bits */
    }

    /*** compute weightedqLPLowPassCoefficients and weightedInputSignal for first subframe ***/

    weighted_q_lp_low_pass_coefficients[0] = sub16(weighted_q_lp_coefficients[0], O7_IN_Q12);
    for i in 1..NB_LSP_COEFF {
        weighted_q_lp_low_pass_coefficients[i] = sub16(
            weighted_q_lp_coefficients[i],
            mult16_16_q12(weighted_q_lp_coefficients[i - 1], O7_IN_Q12) as i16,
        );
    }

    lp_synthesis_filter(
        &lp_residual_signal[0..L_SUBFRAME],
        &weighted_q_lp_low_pass_coefficients,
        &mut weighted_input_signal[0..NB_LSP_COEFF + L_SUBFRAME],
    );

    /*** compute weightedqLPLowPassCoefficients and weightedInputSignal for second subframe ***/

    weighted_q_lp_low_pass_coefficients[0] =
        sub16(weighted_q_lp_coefficients[NB_LSP_COEFF], O7_IN_Q12);
    for i in 1..NB_LSP_COEFF {
        weighted_q_lp_low_pass_coefficients[i] = sub16(
            weighted_q_lp_coefficients[NB_LSP_COEFF + i],
            mult16_16_q12(weighted_q_lp_coefficients[NB_LSP_COEFF + i - 1], O7_IN_Q12) as i16,
        );
    }

    lp_synthesis_filter(
        &lp_residual_signal[L_SUBFRAME..L_FRAME],
        &weighted_q_lp_low_pass_coefficients,
        &mut weighted_input_signal[L_SUBFRAME..L_SUBFRAME + NB_LSP_COEFF + L_SUBFRAME],
    );
}
