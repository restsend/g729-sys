use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;

/*****************************************************************************/
/*                                                                           */
/* Define filter coefficients                                                */
/* Coefficient are given by the filter transfert function :                  */
/*                                                                           */
/*          0.46363718 - 0.92724705z(-1) + 0.46363718z(-2)                   */
/* H(z) = −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−                  */
/*           1 - 1.9059465z(-1) + 0.9114024z(-2)                             */
/*                                                                           */
/* giving:                                                                   */
/*    y[i] = B0*x[i] + B1*x[i-1] + B2*x[i-2]                                 */
/*                   + A1*y[i-1] + A2*y[i-2]                                 */
/*                                                                           */
/*****************************************************************************/

/* coefficients are stored in Q1.13 */
const A1: Word16 = 15836;
const A2: Word16 = -7667;
const B0: Word16 = 7699;
const B1: Word16 = -15398;
const B2: Word16 = 7699;

/* Initialization of context values */
pub fn init_post_processing(
    output_y2: &mut Word32,
    output_y1: &mut Word32,
    input_x0: &mut Word16,
    input_x1: &mut Word16,
) {
    *output_y2 = 0;
    *output_y1 = 0;
    *input_x0 = 0;
    *input_x1 = 0;
}

/*****************************************************************************/
/* postProcessing : high pass filtering and upscaling Spec 4.2.5             */
/*      Algorithm:                                                           */
/*      y[i] = BO*x[i] + B1*x[i-1] + B2*x[i-2] + A1*y[i-1] + A2*y[i-2]       */
/*    parameters:                                                            */
/*      -(i/o) decoderChannelContext : the channel context data              */
/*      -(i/o) signal : 40 values in Q0, reconstructed speech, output        */
/*             replaces the input in buffer                                  */
/*                                                                           */
/*****************************************************************************/
pub fn post_processing(
    output_y2: &mut Word32,
    output_y1: &mut Word32,
    input_x0: &mut Word16,
    input_x1: &mut Word16,
    signal: &mut [Word16],
) {
    let mut input_x2: Word16;
    let mut acc: Word32; /* in Q13 */

    for i in 0..L_SUBFRAME {
        input_x2 = *input_x1;
        *input_x1 = *input_x0;
        *input_x0 = signal[i];

        /* compute with acc and coefficients in Q13 */
        acc = mult16_32_q13(A1, *output_y1); /* Y1 in Q14.13 * A1 in Q1.13 -> acc in Q17.13*/
        acc = mac16_32_q13(acc, A2, *output_y2); /* Y2 in Q14.13 * A2 in Q0.13 -> Q15.13 + acc in Q17.13 -> acc in Q18.12 */
        /* 3*(Xi in Q15.0 * Bi in Q0.13)->Q17.13 + acc in Q18.13 -> acc in 19.13(Overflow??) */
        acc = mac16_16(acc, *input_x0, B0);
        acc = mac16_16(acc, *input_x1, B1);
        acc = saturate(mac16_16(acc, input_x2, B2), MAX_INT29); /* saturate the acc to keep in Q15.13 */

        signal[i] = saturate(pshr(acc, 12), MAX_INT16 as Word32) as Word16; /* acc in Q13 -> *2 and scale back to Q0 */
        *output_y2 = *output_y1;
        *output_y1 = acc;
    }
}
