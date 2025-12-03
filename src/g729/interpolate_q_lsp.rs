use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;

/*****************************************************************************/
/* interpolateqLSP : interpolate previous and current qLSP according to      */
/*      spec. 3.2.5 : interpolated = (current+previous)/2                    */
/*    parameters:                                                            */
/*      -(i) previousqLSP : 10 values in Q0.15:  the qLSP of previous frame  */
/*      -(i) currentqLSP : 10 values in Q0.15: the qLSP of current frame     */
/*      -(o) interpolatedqLSP : 10 values in Q0.15 : the interpolated qLSP   */
/*                                                                           */
/*****************************************************************************/
pub fn interpolate_q_lsp(
    previous_q_lsp: &[Word16],
    current_q_lsp: &[Word16],
    interpolated_q_lsp: &mut [Word16],
) {
    /* interpolate previous and current qLSP according to spec. 3.2.5 : interpolated = (current+previous)/2 */
    for i in 0..NB_LSP_COEFF {
        interpolated_q_lsp[i] = pshr(
            add32(previous_q_lsp[i] as Word32, current_q_lsp[i] as Word32),
            1,
        ) as Word16;
    }
}
