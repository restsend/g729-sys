use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;

pub fn compute_adaptative_codebook_gain(
    target_signal: &[Word16],
    filtered_adaptative_codebook_vector: &[Word16],
    gain_quantization_xy: &mut Word64,
    gain_quantization_yy: &mut Word64,
) -> Word16 {
    *gain_quantization_xy = 0; /* contains the scalar product targetSignal, filteredAdaptativeCodebookVector : numerator */
    *gain_quantization_yy = 0; /* contains the scalar product filteredAdaptativeCodebookVector^2 : denominator */

    for i in 0..L_SUBFRAME {
        *gain_quantization_xy = mac64(
            *gain_quantization_xy,
            target_signal[i] as Word32,
            filtered_adaptative_codebook_vector[i] as Word32,
        );
        *gain_quantization_yy = mac64(
            *gain_quantization_yy,
            filtered_adaptative_codebook_vector[i] as Word32,
            filtered_adaptative_codebook_vector[i] as Word32,
        );
    }

    /* check on values of xx and xy */
    if *gain_quantization_xy <= 0 {
        /* gain would be negative -> return 0 */
        /* this test covers the case of yy(denominator)==0 because if yy==0 then all y==0 and thus xy==0 */
        return 0;
    }

    /* output shall be in Q14 */
    let mut gain = div64(shl64(*gain_quantization_xy, 14), *gain_quantization_yy); /* gain in Q14 */

    /* check if it is not above 1.2 */
    if gain > ONE_POINT_2_IN_Q14 as Word64 {
        gain = ONE_POINT_2_IN_Q14 as Word64;
    }

    gain as Word16
}
