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
        let mut acc = sshl(excitation_vector[i] as Word32, 12); /* acc get the first term of the sum, in Q12 (excitationVector is in Q0)*/
        for j in 0..NB_LSP_COEFF {
            // reconstructedSpeech is accessed at i-j-1.
            // Since reconstructedSpeech starts at index 0 in Rust slice, but logically includes history.
            // The C code assumes reconstructedSpeech points to the start of the current subframe,
            // and negative indices access history.
            // In Rust, we pass the slice starting at the current subframe?
            // No, usually we pass the whole buffer and an offset, or we pass a slice that includes history.
            // The C signature is `word16_t *reconstructedSpeech`.
            // And it accesses `reconstructedSpeech[i-j-1]`.
            // If i=0, j=0 -> index -1.
            // So the slice passed to this function MUST include at least NB_LSP_COEFF elements of history BEFORE the current subframe.
            // So `reconstructed_speech` should be a slice where index `NB_LSP_COEFF` corresponds to the start of the subframe?
            // Let's check how it is called in C.
            // `synthesisFilter(LPResidualSignal, weightedqLPLowPassCoefficients, weightedInputSignal);`
            // `weightedInputSignal` is `90 values in Q0: [-10, -1] as input [0,79] as output`.
            // So the pointer passed is likely pointing to index 10.
            // In Rust, I should pass the slice starting at index 0 (which is history start).
            // And the loop should write to `NB_LSP_COEFF + i`.
            // And read from `NB_LSP_COEFF + i - j - 1`.

            // However, the C code says: `reconstructedSpeech[i] = ...`.
            // This implies `reconstructedSpeech` points to the start of the OUTPUT buffer.
            // And `reconstructedSpeech[-1]` is valid.

            // So in Rust, I will assume the slice passed starts at the history.
            // And I will write to `NB_LSP_COEFF + i`.
            // Wait, if I change the signature to take the whole buffer, I need to know where to start.

            // Let's stick to the C logic but adjust indices.
            // I will assume `reconstructed_speech` contains `NB_LSP_COEFF` history + `L_SUBFRAME` output.
            // Total length `NB_LSP_COEFF + L_SUBFRAME`.
            // The loop `i` goes from 0 to `L_SUBFRAME`.
            // The write index is `NB_LSP_COEFF + i`.
            // The read index is `NB_LSP_COEFF + i - j - 1`.

            acc = msu16_16(
                acc,
                lp_coefficients[j],
                reconstructed_speech[NB_LSP_COEFF + i - j - 1],
            );
        }
        reconstructed_speech[NB_LSP_COEFF + i] =
            saturate(pshr(acc, 12), MAX_16 as Word32) as Word16; /* shift right acc to get it back in Q0 and check overflow on 16 bits */
    }
}
