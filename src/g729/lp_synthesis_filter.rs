use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;

/*****************************************************************************/
/* LPSynthesisFilter : as decribed in spec 4.1.6 eq77                        */
/*    parameters:                                                            */
/*      -(i) excitationVector: u(n), the excitation, 40 values in Q0         */
/*      -(i) LPCoefficients: 10 LP coefficients in Q12                       */
/*      -(i/o) recontructedSpeech: 50 values in Q0                           */
/*             [-NB_LSP_COEFF, -1] of previous values as input               */
/*             [0, L_SUBFRAME[ as output                                     */
/*                                                                           */
/*****************************************************************************/
pub fn lp_synthesis_filter(
    excitation_vector: &[Word16],
    lp_coefficients: &[Word16],
    reconstructed_speech: &mut [Word16],
) {
    /* compute excitationVector[i] - Sum0-9(LPCoefficients[j]*reconstructedSpeech[i-j]) */
    for i in 0..L_SUBFRAME {
        let base = NB_LSP_COEFF + i;
        let mut acc = sshl(excitation_vector[i] as Word32, 12); /* acc get the first term of the sum, in Q12 (excitationVector is in Q0)*/
        // Unrolled 10-tap MAC (same order as j=0..9: highest history index first).
        acc = msu16_16(acc, lp_coefficients[0], reconstructed_speech[base - 1]);
        acc = msu16_16(acc, lp_coefficients[1], reconstructed_speech[base - 2]);
        acc = msu16_16(acc, lp_coefficients[2], reconstructed_speech[base - 3]);
        acc = msu16_16(acc, lp_coefficients[3], reconstructed_speech[base - 4]);
        acc = msu16_16(acc, lp_coefficients[4], reconstructed_speech[base - 5]);
        acc = msu16_16(acc, lp_coefficients[5], reconstructed_speech[base - 6]);
        acc = msu16_16(acc, lp_coefficients[6], reconstructed_speech[base - 7]);
        acc = msu16_16(acc, lp_coefficients[7], reconstructed_speech[base - 8]);
        acc = msu16_16(acc, lp_coefficients[8], reconstructed_speech[base - 9]);
        acc = msu16_16(acc, lp_coefficients[9], reconstructed_speech[base - 10]);
        reconstructed_speech[base] = saturate(pshr(acc, 12), MAX_16 as Word32) as Word16;
        /* shift right acc to get it back in Q0 and check overflow on 16 bits */
    }
}
