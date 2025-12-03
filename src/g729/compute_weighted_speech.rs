use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;
use crate::g729::lp_synthesis_filter::*;

/*****************************************************************************/
/* computeWeightedSpeech: compute wieghted speech according to spec A3.3.3   */
/*    parameters:                                                            */
/*      -(i) inputSignal : 90 values buffer accessed in range [-10, 79] in Q0*/
/*      -(i) qLPCoefficients: 20 coefficients(10 for each subframe) in Q12   */
/*      -(i) weightedqLPCoefficients: 20 coefficients(10 for each subframe)  */
/*           in Q12                                                          */
/*      -(i/o) weightedInputSignal: 90 values in Q0: [-10, -1] as input      */
/*             [0,79] as output in Q0                                        */
/*      -(o) LPResidualSignal: 80 values of residual signal in Q0            */
/*                                                                           */
/*****************************************************************************/
pub fn compute_weighted_speech(
    input_signal: &[Word16],
    q_lp_coefficients: &[Word16],
    weighted_q_lp_coefficients: &[Word16],
    weighted_input_signal: &mut [Word16],
    lp_residual_signal: &mut [Word16],
) {
    /* algo as specified in A3.3.3: */
    /* first compute LPResidualSignal[n] = inputSignal[n] + ∑(i=1..10)qLP[i]*inputSignal[n-i] specA3.3.3 eqA.3 */
    /* then compute qLP' coefficients as qLP'[i] = weightedqLP[i] - 0.7*weightedqLP[i-1] */

    /* finally get the weightedInputSignal[n] = LPResidualSignal[n] - ∑(i=1..10)qLP'[i]*weightedInputSignal[n-i] spec A3.3.3 eqA.2 */

    let mut weighted_q_lp_low_pass_coefficients = [0; NB_LSP_COEFF]; /* in Q12 */

    /*** compute LPResisualSignal (spec A3.3.3 eqA.3) in Q0 ***/
    /* compute residual signal for the first subframe: use the first 10 qLPCoefficients */
    for i in 0..L_SUBFRAME {
        let mut acc = sshl(input_signal[NB_LSP_COEFF + i] as Word32, 12); /* inputSignal in Q0 is shifted to set acc in Q12 */
        for j in 0..NB_LSP_COEFF {
            acc = mac16_16(
                acc,
                q_lp_coefficients[j],
                input_signal[NB_LSP_COEFF + i - j - 1],
            ); /* qLPCoefficients in Q12, inputSignal in Q0 -> acc in Q12 */
        }
        lp_residual_signal[i] = saturate(pshr(acc, 12), MAX_16 as Word32) as Word16;
        /* shift back acc to Q0 and saturate it to avoid overflow when going back to 16 bits */
    }
    /* compute residual signal for the second subframe: use the second part of qLPCoefficients */
    for i in L_SUBFRAME..L_FRAME {
        let mut acc = sshl(input_signal[NB_LSP_COEFF + i] as Word32, 12); /* inputSignal in Q0 is shifted to set acc in Q12 */
        for j in 0..NB_LSP_COEFF {
            acc = mac16_16(
                acc,
                q_lp_coefficients[NB_LSP_COEFF + j],
                input_signal[NB_LSP_COEFF + i - j - 1],
            ); /* qLPCoefficients in Q12, inputSignal in Q0 -> acc in Q12 */
        }
        lp_residual_signal[i] = saturate(pshr(acc, 12), MAX_16 as Word32) as Word16;
        /* shift back acc to Q0 and saturate it to avoid overflow when going back to 16 bits */
    }

    /*** compute weightedqLPLowPassCoefficients and weightedInputSignal for first subframe ***/
    /* spec A3.3.3 a' = weightedqLPLowPassCoefficients[i] =  weightedqLP[i] - 0.7*weightedqLP[i-1] */
    weighted_q_lp_low_pass_coefficients[0] = sub16(weighted_q_lp_coefficients[0], O7_IN_Q12); /* weightedqLP[-1] = 1 -> weightedqLPLowPassCoefficients[0] =  weightedqLPCoefficients[0] - 0.7 */
    for i in 1..NB_LSP_COEFF {
        weighted_q_lp_low_pass_coefficients[i] = sub16(
            weighted_q_lp_coefficients[i],
            mult16_16_q12(weighted_q_lp_coefficients[i - 1], O7_IN_Q12) as Word16,
        );
    }

    /* weightedInputSignal for the first subframe: synthesis filter  1/[A'(z)] */
    // We pass the slice of weighted_input_signal that corresponds to the first subframe + history.
    // The history is at indices [0..NB_LSP_COEFF].
    // The output is at indices [NB_LSP_COEFF..NB_LSP_COEFF+L_SUBFRAME].
    // lp_synthesis_filter expects a slice where index NB_LSP_COEFF is the start of output.
    // So we pass &mut weighted_input_signal[0..NB_LSP_COEFF+L_SUBFRAME].
    // And the excitation is lp_residual_signal[0..L_SUBFRAME].
    lp_synthesis_filter(
        &lp_residual_signal[0..L_SUBFRAME],
        &weighted_q_lp_low_pass_coefficients,
        &mut weighted_input_signal[0..NB_LSP_COEFF + L_SUBFRAME],
    );

    /*** compute weightedqLPLowPassCoefficients and weightedInputSignal for second subframe ***/
    /* spec A3.3.3 a' = weightedqLPLowPassCoefficients[i] =  weightedqLP[i] - 0.7*weightedqLP[i-1] */
    weighted_q_lp_low_pass_coefficients[0] =
        sub16(weighted_q_lp_coefficients[NB_LSP_COEFF], O7_IN_Q12); /* weightedqLP[-1] = 1 -> weightedqLPLowPassCoefficients[0] =  weightedqLPCoefficients[0] - 0.7 */
    for i in 1..NB_LSP_COEFF {
        weighted_q_lp_low_pass_coefficients[i] = sub16(
            weighted_q_lp_coefficients[NB_LSP_COEFF + i],
            mult16_16_q12(weighted_q_lp_coefficients[NB_LSP_COEFF + i - 1], O7_IN_Q12) as Word16,
        );
    }

    /* weightedInputSignal for the second subframe: synthesis filter  1/[A'(z)] */
    // For the second subframe, the history is the previous subframe's output.
    // The previous subframe output ended at index NB_LSP_COEFF + L_SUBFRAME.
    // The history for the second subframe is the last NB_LSP_COEFF samples of the previous output.
    // Wait, `lp_synthesis_filter` expects the slice to start with history.
    // The history for 2nd subframe is at indices [L_SUBFRAME..L_SUBFRAME+NB_LSP_COEFF].
    // No.
    // `weighted_input_signal` has 90 elements.
    // 0..10: History for 1st subframe.
    // 10..50: Output of 1st subframe.
    // 50..90: Output of 2nd subframe.
    // History for 2nd subframe is 40..50.
    // Wait, NB_LSP_COEFF is 10.
    // So history for 2nd subframe is indices [40..50].
    // And output is [50..90].
    // So we pass slice starting at 40.
    // &mut weighted_input_signal[L_SUBFRAME..L_SUBFRAME + NB_LSP_COEFF + L_SUBFRAME]
    // i.e. [40..90].
    lp_synthesis_filter(
        &lp_residual_signal[L_SUBFRAME..L_FRAME],
        &weighted_q_lp_low_pass_coefficients,
        &mut weighted_input_signal[L_SUBFRAME..L_SUBFRAME + NB_LSP_COEFF + L_SUBFRAME],
    );
}
