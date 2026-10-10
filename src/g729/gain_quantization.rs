use crate::g729::basic_operations::*;
use crate::g729::codebooks::*;
use crate::g729::fixed_point_math::*;
use crate::g729::ld8k::*;

fn count_leading_zeros(x: i32) -> u16 {
    if x == 0 {
        return 31;
    }
    (x.leading_zeros() - 1) as u16
}

/// Predict the fixed codebook gain using MA prediction.
///
/// # Arguments
///
/// * `previous_gain_prediction_error` - (i16) Previous gain prediction error in Q10.
/// * `fixed_codebook_vector` - (i16) Fixed codebook vector in Q1.13.
///
/// # Returns
///
/// * `predicted_gain` - (i32) Predicted fixed codebook gain in Q16.
#[cfg_attr(target_arch = "xtensa", inline(never))]
pub fn ma_code_gain_prediction(
    previous_gain_prediction_error: &[i16; 4],
    fixed_codebook_vector: &[i16],
) -> i32 {
    let mut acc: i32;
    let mut fixed_codebook_vector_squares_sum: i32 = 0;

    /* compute the sum of squares of fixedCodebookVector in Q26 */
    for &coefficient in fixed_codebook_vector.iter().take(L_SUBFRAME) {
        if coefficient != 0 {
            fixed_codebook_vector_squares_sum =
                mac16_16(fixed_codebook_vector_squares_sum, coefficient, coefficient);
        }
    }

    /* compute E| - E as in eq71, result in Q16 */

    acc = mac16_32_q13(
        8145364,
        -24660,
        g729_log2_q0q16(fixed_codebook_vector_squares_sum),
    ); /* acc in Q16 */

    /* accumulate the MA prediction described in eq69 to the previous Sum, result will be in E~(m) + E| -E as used in eq71 */

    acc = shl(acc, 8); /* acc in Q24 to match the fixed point of next accumulations */
    for i in 0..4 {
        acc = mac16_16(
            acc,
            previous_gain_prediction_error[i],
            MA_PREDICTION_COEFFICIENTS[i],
        );
    }

    /* compute eq71, we already have the exposant in acc so */

    acc = shr(acc, 2);
    acc = mult16_32_q15(5442, acc);
    acc = pshr(acc, 11); /* get acc in Q4.11 */

    g729_exp2_q11q16(acc as i16)
}

/// Update the gain prediction error.
///
/// # Arguments
///
/// * `fixed_codebook_gain_correction_factor` - (i16) Gamma in eq72 in Q3.12.
/// * `previous_gain_prediction_error` - (i16) Previous gain prediction error in Q10.
#[cfg_attr(target_arch = "xtensa", inline(never))]
pub fn compute_gain_prediction_error(
    fixed_codebook_gain_correction_factor: i16,
    previous_gain_prediction_error: &mut [i16; 4],
) {
    /* need to compute eq72: 20log10(fixedCodebookGainCorrectionFactor) */

    let mut current_gain_prediction_error: i32 = sub32(
        g729_log2_q0q16(fixed_codebook_gain_correction_factor as i32),
        786432,
    );
    current_gain_prediction_error = pshr(mult16_32_q12(24660, current_gain_prediction_error), 6);

    previous_gain_prediction_error[3] = previous_gain_prediction_error[2];
    previous_gain_prediction_error[2] = previous_gain_prediction_error[1];
    previous_gain_prediction_error[1] = previous_gain_prediction_error[0];
    previous_gain_prediction_error[0] = current_gain_prediction_error as i16;
}

/// Quantize the adaptive and fixed codebook gains.
///
/// # Arguments
///
/// * `target_signal` - (i16) Target signal in Q0.
/// * `filtered_adaptative_codebook_vector` - (i16) Filtered adaptive codebook vector in Q0.
/// * `convolved_fixed_codebook_vector` - (i16) Convolved fixed codebook vector in Q12.
/// * `fixed_codebook_vector` - (i16) Fixed codebook vector in Q13.
/// * `xy64` - (i64) xy term of eq63 computed previously in Q0.
/// * `yy64` - (i64) yy term of eq63 computed previously in Q0.
/// * `previous_gain_prediction_error` - (i16) Previous gain prediction error in Q10.
/// * `quantized_adaptative_codebook_gain` - (i16) Quantized adaptive codebook gain in Q14.
/// * `quantized_fixed_codebook_gain` - (i16) Quantized fixed codebook gain in Q1.
/// * `gain_codebook_stage1` - (u16) GA parameter value (3 bits).
/// * `gain_codebook_stage2` - (u16) GB parameter value (4 bits).
#[allow(clippy::too_many_arguments)]
#[cfg_attr(target_arch = "xtensa", inline(never))]
pub fn gain_quantization(
    target_signal: &[i16],
    filtered_adaptative_codebook_vector: &[i16],
    convolved_fixed_codebook_vector: &[i16],
    fixed_codebook_vector: &[i16],
    xy64: i64,
    yy64: i64,
    previous_gain_prediction_error: &mut [i16; 4],
    quantized_adaptative_codebook_gain: &mut i16,
    quantized_fixed_codebook_gain: &mut i16,
    gain_codebook_stage1: &mut u16,
    gain_codebook_stage2: &mut u16,
) {
    let mut xz64: i64 = 0;
    let mut yz64: i64 = 0;
    let mut zz64: i64 = 0;
    let mut xy: i32;
    let mut yy: i32;
    let mut xz: i32;
    let mut yz: i32;
    let mut zz: i32;
    let mut min_normalization: u16 = 31;
    let mut current_normalization: u16;
    let best_adaptative_codebook_gain: i32;
    let best_fixed_codebook_gain: i32;

    let mut index_base_ga: usize = 0;
    let mut index_base_gb: usize = 0;
    let mut index_ga: usize = 0;
    let mut index_gb: usize = 0;
    let mut distance_min: i64 = i64::MAX;

    /*** compute spec 3.9 eq63 terms first on 64 bits and then scale them if needed to fit on 32 ***/

    for i in 0..L_SUBFRAME {
        xz64 = mac64(
            xz64,
            target_signal[i] as i32,
            convolved_fixed_codebook_vector[i] as i32,
        ); /* in Q12 */
        yz64 = mac64(
            yz64,
            filtered_adaptative_codebook_vector[i] as i32,
            convolved_fixed_codebook_vector[i] as i32,
        ); /* in Q12 */
        zz64 = mac64(
            zz64,
            convolved_fixed_codebook_vector[i] as i32,
            convolved_fixed_codebook_vector[i] as i32,
        ); /* in Q24 */
    }

    /* now scale this terms to have them fit on 32 bits - terms Xy, Xz and Yz shall fit on 31 bits because used in eq63 with a factor 2 */
    xy = shr64(if xy64 < 0 { -xy64 } else { xy64 }, 30) as i32;
    yy = shr64(yy64, 31) as i32;
    xz = shr64(if xz64 < 0 { -xz64 } else { xz64 }, 30) as i32;
    yz = shr64(if yz64 < 0 { -yz64 } else { yz64 }, 30) as i32;
    zz = shr64(zz64, 31) as i32;

    current_normalization = count_leading_zeros(xy);
    if current_normalization < min_normalization {
        min_normalization = current_normalization;
    }
    current_normalization = count_leading_zeros(xz);
    if current_normalization < min_normalization {
        min_normalization = current_normalization;
    }
    current_normalization = count_leading_zeros(yz);
    if current_normalization < min_normalization {
        min_normalization = current_normalization;
    }
    current_normalization = count_leading_zeros(yy);
    if current_normalization < min_normalization {
        min_normalization = current_normalization;
    }
    current_normalization = count_leading_zeros(zz);
    if current_normalization < min_normalization {
        min_normalization = current_normalization;
    }

    if min_normalization < 31 {
        min_normalization = 31 - min_normalization;
        xy = shr64(xy64, min_normalization as u32) as i32;
        yy = shr64(yy64, min_normalization as u32) as i32;
        xz = shr64(xz64, min_normalization as u32) as i32;
        yz = shr64(yz64, min_normalization as u32) as i32;
        zz = shr64(zz64, min_normalization as u32) as i32;
    } else {
        xy = xy64 as i32; /* in Q0 */
        yy = yy64 as i32; /* in Q0 */
        xz = xz64 as i32; /* in Q12 */
        yz = yz64 as i32; /* in Q12 */
        zz = zz64 as i32; /* in Q24 */
    }

    /*** compute the best gains minimizinq eq63 ***/
    /* Note this bestgain computation is not at all described in the spec, got it from ITU code */

    /* best gain are computed in Q9 and Q2 and fits on 16 bits */
    let denominator: i64 = mac64(mult32_32(yy, zz), -yz, yz); /* (yy*zz) - yz^2) in Q24 (always >= 0)*/

    if denominator == 0 {
        best_adaptative_codebook_gain = shr64(mac64(mult32_32(zz, xy), -xz, yz), 15) as i32;
        best_fixed_codebook_gain = shr64(mac64(mult32_32(yy, xz), -xy, yz), 10) as i32;
    } else {
        /* bestAdaptativeCodebookGain in Q9 */
        let mut numerator_norm: u16;
        let mut numerator: i64 = mac64(mult32_32(zz, xy), -xz, yz); /* in Q24 */
        /* check if we can shift it by 9 without overflow as the bestAdaptativeCodebookGain in computed in Q9 */
        let mut numerator_h: i32 = shr64(numerator, 32) as i32;
        numerator_h = if numerator_h > 0 {
            numerator_h
        } else {
            -numerator_h
        };
        numerator_norm = count_leading_zeros(numerator_h);
        if numerator_norm >= 9 {
            best_adaptative_codebook_gain = div64(sshl64(numerator, 9), denominator) as i32;
        /* bestAdaptativeCodebookGain in Q9 */
        } else {
            let shifted_denominator: i64 = shr64(denominator, (9 - numerator_norm) as u32);
            if shifted_denominator > 0 {
                best_adaptative_codebook_gain =
                    div64(shl64(numerator, numerator_norm as u32), shifted_denominator) as i32;
            /* bestAdaptativeCodebookGain in Q9 */
            } else {
                best_adaptative_codebook_gain = shl32(
                    div64(shl64(numerator, numerator_norm as u32), denominator) as i32,
                    (9 - numerator_norm) as u32,
                ); /* shift left the division result to reach Q9 */
            }
        }

        numerator = mac64(mult32_32(yy, xz), -xy, yz); /* in Q12 */
        /* check if we can shift it by 14(it's in Q12 and denominator in Q24) without overflow as the bestFixedCodebookGain in computed in Q2 */
        numerator_h = shr64(numerator, 32) as i32;
        numerator_h = if numerator_h > 0 {
            numerator_h
        } else {
            -numerator_h
        };
        numerator_norm = count_leading_zeros(numerator_h);

        if numerator_norm >= 14 {
            best_fixed_codebook_gain = div64(sshl64(numerator, 14), denominator) as i32;
        } else {
            let shifted_denominator: i64 = shr64(denominator, (14 - numerator_norm) as u32); /* bestFixedCodebookGain in Q14 */
            if shifted_denominator > 0 {
                best_fixed_codebook_gain =
                    div64(shl64(numerator, numerator_norm as u32), shifted_denominator) as i32;
            /* bestFixedCodebookGain in Q14 */
            } else {
                best_fixed_codebook_gain = shl32(
                    div64(shl64(numerator, numerator_norm as u32), denominator) as i32,
                    (14 - numerator_norm) as u32,
                ); /* shift left the division result to reach Q14 */
            }
        }
    }

    /*** Compute the predicted gain as in spec 3.9.1 eq71 in Q6 ***/
    let predicted_fixed_codebook_gain: i16 = shr32(
        ma_code_gain_prediction(previous_gain_prediction_error, fixed_codebook_vector),
        12,
    ) as i16;

    /***  preselection spec 3.9.2 ***/
    /* Note: spec just says to select the best 50% of each vector, ITU code go through magical constant computation to select the begining of a continuous range */

    while index_base_ga < 6
        && best_fixed_codebook_gain
            > mult16_16_q14(GA_CODEBOOK[index_base_ga][1], predicted_fixed_codebook_gain)
    {
        index_base_ga += 1;
    }
    index_base_ga = index_base_ga.saturating_sub(2);
    while index_base_gb < 12
        && best_adaptative_codebook_gain > shr(GB_CODEBOOK[index_base_gb][0] as i32, 5)
    {
        index_base_gb += 1;
    }
    index_base_gb = index_base_gb.saturating_sub(4);

    /*** test all possibilities of Ga and Gb indexes and select the best one ***/
    xy = -sshl(xy, 1);
    xz = -sshl(xz, 1);
    yz = sshl(yz, 1);

    for i in 0..4 {
        for j in 0..8 {
            let gp: i16 = add16(
                GA_CODEBOOK[i + index_base_ga][0],
                GB_CODEBOOK[j + index_base_gb][0],
            ); /* result in Q14 */
            let gamma: i16 = add16(
                GA_CODEBOOK[i + index_base_ga][1],
                GB_CODEBOOK[j + index_base_gb][1],
            );
            let gc: i32 = mult16_16_q14(gamma, predicted_fixed_codebook_gain);

            /* compute E as in eq63 (first term excluded) */
            let mut acc: i64 = mult32_32(mult16_16(gp, gp), yy);
            acc = mac64(acc, mult16_16(gc as i16, gc as i16), zz);
            acc = mac64(acc, shl32(gp as i32, 14), xy);
            acc = mac64(acc, shl32(gc, 14), xz);
            acc = mac64(acc, mult16_16(gp, gc as i16), yz);

            if acc < distance_min {
                distance_min = acc;
                index_ga = i + index_base_ga;
                index_gb = j + index_base_gb;
                *quantized_adaptative_codebook_gain = gp;
                *quantized_fixed_codebook_gain = shr(gc, 1) as i16;
            }
        }
    }

    compute_gain_prediction_error(
        add16(GA_CODEBOOK[index_ga][1], GB_CODEBOOK[index_gb][1]),
        previous_gain_prediction_error,
    );

    *gain_codebook_stage1 = INDEX_MAPPING_GA[index_ga];
    *gain_codebook_stage2 = INDEX_MAPPING_GB[index_gb];
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gain_quantization() {
        let target_signal: [i16; 40] = [0; 40];
        let filtered_adaptative_codebook_vector: [i16; 40] = [0; 40];
        let convolved_fixed_codebook_vector: [i16; 40] = [0; 40];
        let fixed_codebook_vector: [i16; 40] = [0; 40];
        let xy64: i64 = 0;
        let yy64: i64 = 0;
        let mut previous_gain_prediction_error: [i16; 4] = [-14336, -14336, -14336, -14336];
        let mut quantized_adaptative_codebook_gain: i16 = 0;
        let mut quantized_fixed_codebook_gain: i16 = 0;
        let mut gain_codebook_stage1: u16 = 0;
        let mut gain_codebook_stage2: u16 = 0;

        gain_quantization(
            &target_signal,
            &filtered_adaptative_codebook_vector,
            &convolved_fixed_codebook_vector,
            &fixed_codebook_vector,
            xy64,
            yy64,
            &mut previous_gain_prediction_error,
            &mut quantized_adaptative_codebook_gain,
            &mut quantized_fixed_codebook_gain,
            &mut gain_codebook_stage1,
            &mut gain_codebook_stage2,
        );
    }
}
