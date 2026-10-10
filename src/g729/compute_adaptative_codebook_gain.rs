use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;

#[cfg_attr(target_arch = "xtensa", inline(never))]
pub fn compute_adaptative_codebook_gain(
    target_signal: &[i16],
    filtered_adaptative_codebook_vector: &[i16],
    gain_quantization_xy: &mut i64,
    gain_quantization_yy: &mut i64,
) -> i16 {
    *gain_quantization_xy = 0;
    *gain_quantization_yy = 0;

    for i in 0..L_SUBFRAME {
        *gain_quantization_xy = mac64(
            *gain_quantization_xy,
            target_signal[i] as i32,
            filtered_adaptative_codebook_vector[i] as i32,
        );
        *gain_quantization_yy = mac64(
            *gain_quantization_yy,
            filtered_adaptative_codebook_vector[i] as i32,
            filtered_adaptative_codebook_vector[i] as i32,
        );
    }

    if *gain_quantization_xy <= 0 {
        return 0;
    }

    /* output shall be in Q14 */
    let mut gain = div64(shl64(*gain_quantization_xy, 14), *gain_quantization_yy); /* gain in Q14 */

    if gain > ONE_POINT_2_IN_Q14 as i64 {
        gain = ONE_POINT_2_IN_Q14 as i64;
    }

    gain as i16
}
