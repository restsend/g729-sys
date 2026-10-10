use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;

/*****************************************************************************/

/*****************************************************************************/

/* coefficients are stored in Q1.13 */
const A1: i16 = 15836;
const A2: i16 = -7667;
const B0: i16 = 7699;
const B1: i16 = -15398;
const B2: i16 = 7699;

pub fn init_post_processing() -> (i32, i32, i16, i16) {
    (0, 0, 0, 0)
}

/*****************************************************************************/
/* postProcessing : high pass filtering and upscaling Spec 4.2.5             */

/*      -(i/o) signal : 40 values in Q0, reconstructed speech, output        */

/*****************************************************************************/
pub fn post_processing(
    output_y2: &mut i32,
    output_y1: &mut i32,
    input_x0: &mut i16,
    input_x1: &mut i16,
    signal: &mut [i16],
) {
    let mut input_x2: i16;
    let mut acc: i32; /* in Q13 */

    for sample in signal.iter_mut().take(L_SUBFRAME) {
        input_x2 = *input_x1;
        *input_x1 = *input_x0;
        *input_x0 = *sample;

        /* compute with acc and coefficients in Q13 */
        acc = mult16_32_q13(A1, *output_y1);
        acc = mac16_32_q13(acc, A2, *output_y2);

        acc = mac16_16(acc, *input_x0, B0);
        acc = mac16_16(acc, *input_x1, B1);
        acc = saturate(mac16_16(acc, input_x2, B2), MAX_INT29); /* saturate the acc to keep in Q15.13 */

        *sample = saturate(pshr(acc, 12), MAX_INT16 as i32) as i16;
        *output_y2 = *output_y1;
        *output_y1 = acc;
    }
}
