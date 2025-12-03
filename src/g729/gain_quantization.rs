use crate::g729::basic_operations::*;
use crate::g729::codebooks::*;
use crate::g729::fixed_point_math::*;
use crate::g729::ld8k::*;

fn count_leading_zeros(x: Word32) -> u16 {
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
pub fn ma_code_gain_prediction(
    previous_gain_prediction_error: &[Word16; 4],
    fixed_codebook_vector: &[Word16],
) -> Word32 {
    let mut acc: Word32;
    let mut fixed_codebook_vector_squares_sum: Word32 = 0;

    /* compute the sum of squares of fixedCodebookVector in Q26 */
    for i in 0..L_SUBFRAME {
        if fixed_codebook_vector[i] != 0 {
            fixed_codebook_vector_squares_sum = mac16_16(
                fixed_codebook_vector_squares_sum,
                fixed_codebook_vector[i],
                fixed_codebook_vector[i],
            );
        }
    }

    /* compute E| - E as in eq71, result in Q16 */
    /* acc = 124.2884 - 3.0103*log2(Sum) */
    /* acc = 8145364[32 bits Q16] - 24660[16 bits Q2.13]*log2(Sum)[32 bits Q16] */
    acc = mac16_32_q13(
        8145364,
        -24660,
        g729_log2_q0q16(fixed_codebook_vector_squares_sum),
    ); /* acc in Q16 */

    /* accumulate the MA prediction described in eq69 to the previous Sum, result will be in E~(m) + E| -E as used in eq71 */
    /* acc in Q16->Q24, previousGainPredictionError in Q10 and MAPredictionCoefficients in Q0.14*/
    acc = shl(acc, 8); /* acc in Q24 to match the fixed point of next accumulations */
    for i in 0..4 {
        acc = mac16_16(
            acc,
            previous_gain_prediction_error[i],
            MA_PREDICTION_COEFFICIENTS[i],
        );
    }
    /* acc range [10, 65] -> Q6.24 */

    /* compute eq71, we already have the exposant in acc so */
    /* g'c = 10^(acc/20)                                    */
    /*     = 2^((acc*ln(10))/(20*ln(2)))                    */
    /*     = 2^(0,1661*acc)                                 */
    acc = shr(acc, 2); /* Q24->Q22 */
    acc = mult16_32_q15(5442, acc); /* 5442 is 0.1661 in Q15 -> acc now in Q4.22 range [1.6, 10.8] */
    acc = pshr(acc, 11); /* get acc in Q4.11 */

    return g729_exp2_q11q16(acc as Word16); /* acc fits on 16 bits, cast it to word16_t to send it to the exp2 function, output in Q16*/
}

/// Update the gain prediction error.
///
/// # Arguments
///
/// * `fixed_codebook_gain_correction_factor` - (i16) Gamma in eq72 in Q3.12.
/// * `previous_gain_prediction_error` - (i16) Previous gain prediction error in Q10.
pub fn compute_gain_prediction_error(
    fixed_codebook_gain_correction_factor: Word16,
    previous_gain_prediction_error: &mut [Word16; 4],
) {
    /* need to compute eq72: 20log10(fixedCodebookGainCorrectionFactor) */
    /*  = (20/log2(10))*log2(fixedCodebookGainCorrectionFactor) */
    /*  = 6.0206*log2(fixedCodebookGainCorrectionFactor) */
    /* log2 input in Q0, output in Q16,fixedCodebookGainCorrectionFactor being in Q12, we shall substract 12 in Q16(786432) to the result of log2 function -> final result in Q2.16 */
    let mut current_gain_prediction_error: Word32 = sub32(
        g729_log2_q0q16(fixed_codebook_gain_correction_factor as Word32),
        786432,
    );
    current_gain_prediction_error = pshr(mult16_32_q12(24660, current_gain_prediction_error), 6); /* 24660 = 6.0206 in Q3.12 -> mult result in Q16, precise shift right to get it in Q4.10 */

    /* shift the array and insert the current Prediction Error */
    previous_gain_prediction_error[3] = previous_gain_prediction_error[2];
    previous_gain_prediction_error[2] = previous_gain_prediction_error[1];
    previous_gain_prediction_error[1] = previous_gain_prediction_error[0];
    previous_gain_prediction_error[0] = current_gain_prediction_error as Word16;
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
pub fn gain_quantization(
    target_signal: &[Word16],
    filtered_adaptative_codebook_vector: &[Word16],
    convolved_fixed_codebook_vector: &[Word16],
    fixed_codebook_vector: &[Word16],
    xy64: Word64,
    yy64: Word64,
    previous_gain_prediction_error: &mut [Word16; 4],
    quantized_adaptative_codebook_gain: &mut Word16,
    quantized_fixed_codebook_gain: &mut Word16,
    gain_codebook_stage1: &mut u16,
    gain_codebook_stage2: &mut u16,
) {
    let mut xz64: Word64 = 0;
    let mut yz64: Word64 = 0;
    let mut zz64: Word64 = 0;
    let mut xy: Word32;
    let mut yy: Word32;
    let mut xz: Word32;
    let mut yz: Word32;
    let mut zz: Word32;
    let mut min_normalization: u16 = 31;
    let mut current_normalization: u16;
    let best_adaptative_codebook_gain: Word32;
    let best_fixed_codebook_gain: Word32;
    let denominator: Word64;
    let predicted_fixed_codebook_gain: Word16;
    let mut index_base_ga: usize = 0;
    let mut index_base_gb: usize = 0;
    let mut index_ga: usize = 0;
    let mut index_gb: usize = 0;
    let mut distance_min: Word64 = i64::MAX;

    /*** compute spec 3.9 eq63 terms first on 64 bits and then scale them if needed to fit on 32 ***/
    /* Xy64 and Yy64 already computed during adaptativeCodebookGain computation */
    for i in 0..L_SUBFRAME {
        xz64 = mac64(
            xz64,
            target_signal[i] as Word32,
            convolved_fixed_codebook_vector[i] as Word32,
        ); /* in Q12 */
        yz64 = mac64(
            yz64,
            filtered_adaptative_codebook_vector[i] as Word32,
            convolved_fixed_codebook_vector[i] as Word32,
        ); /* in Q12 */
        zz64 = mac64(
            zz64,
            convolved_fixed_codebook_vector[i] as Word32,
            convolved_fixed_codebook_vector[i] as Word32,
        ); /* in Q24 */
    }

    /* now scale this terms to have them fit on 32 bits - terms Xy, Xz and Yz shall fit on 31 bits because used in eq63 with a factor 2 */
    xy = shr64(if xy64 < 0 { -xy64 } else { xy64 }, 30) as Word32;
    yy = shr64(yy64, 31) as Word32;
    xz = shr64(if xz64 < 0 { -xz64 } else { xz64 }, 30) as Word32;
    yz = shr64(if yz64 < 0 { -yz64 } else { yz64 }, 30) as Word32;
    zz = shr64(zz64, 31) as Word32;

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
        /* we shall normalise, values are over 32 bits */
        min_normalization = 31 - min_normalization;
        xy = shr64(xy64, min_normalization as u32) as Word32;
        yy = shr64(yy64, min_normalization as u32) as Word32;
        xz = shr64(xz64, min_normalization as u32) as Word32;
        yz = shr64(yz64, min_normalization as u32) as Word32;
        zz = shr64(zz64, min_normalization as u32) as Word32;
    } else {
        /* no need to normalise, values already fit on 32 bits, just cast them */
        xy = xy64 as Word32; /* in Q0 */
        yy = yy64 as Word32; /* in Q0 */
        xz = xz64 as Word32; /* in Q12 */
        yz = yz64 as Word32; /* in Q12 */
        zz = zz64 as Word32; /* in Q24 */
    }

    /*** compute the best gains minimizinq eq63 ***/
    /* Note this bestgain computation is not at all described in the spec, got it from ITU code */
    /* bestAdaptativeCodebookGain = (zz.xy - xz.yz) / (yy*zz) - yz^2) */
    /* bestfixedCodebookGain = (yy*xz - xy*yz) / (yy*zz) - yz^2) */
    /* best gain are computed in Q9 and Q2 and fits on 16 bits */
    denominator = mac64(mult32_32(yy, zz), -yz, yz); /* (yy*zz) - yz^2) in Q24 (always >= 0)*/
    /* avoid division by zero */
    if denominator == 0 {
        /* consider it to be one */
        best_adaptative_codebook_gain = shr64(mac64(mult32_32(zz, xy), -xz, yz), 15) as Word32; /* MAC in Q24 -> Q9 */
        best_fixed_codebook_gain = shr64(mac64(mult32_32(yy, xz), -xy, yz), 10) as Word32;
    /* MAC in Q12 -> Q2 */
    } else {
        /* bestAdaptativeCodebookGain in Q9 */
        let mut numerator_norm: u16;
        let mut numerator: Word64 = mac64(mult32_32(zz, xy), -xz, yz); /* in Q24 */
        /* check if we can shift it by 9 without overflow as the bestAdaptativeCodebookGain in computed in Q9 */
        let mut numerator_h: Word32 = shr64(numerator, 32) as Word32;
        numerator_h = if numerator_h > 0 {
            numerator_h
        } else {
            -numerator_h
        };
        numerator_norm = count_leading_zeros(numerator_h);
        if numerator_norm >= 9 {
            best_adaptative_codebook_gain = div64(sshl64(numerator, 9), denominator) as Word32;
        /* bestAdaptativeCodebookGain in Q9 */
        } else {
            let shifted_denominator: Word64 = shr64(denominator, (9 - numerator_norm) as u32);
            if shifted_denominator > 0 {
                /* can't shift left by 9 the numerator, can we shift right by 9-numeratorNorm the denominator without hiting 0 */
                best_adaptative_codebook_gain =
                    div64(shl64(numerator, numerator_norm as u32), shifted_denominator) as Word32;
            /* bestAdaptativeCodebookGain in Q9 */
            } else {
                best_adaptative_codebook_gain = shl32(
                    div64(shl64(numerator, numerator_norm as u32), denominator) as Word32,
                    (9 - numerator_norm) as u32,
                ); /* shift left the division result to reach Q9 */
            }
        }

        numerator = mac64(mult32_32(yy, xz), -xy, yz); /* in Q12 */
        /* check if we can shift it by 14(it's in Q12 and denominator in Q24) without overflow as the bestFixedCodebookGain in computed in Q2 */
        numerator_h = shr64(numerator, 32) as Word32;
        numerator_h = if numerator_h > 0 {
            numerator_h
        } else {
            -numerator_h
        };
        numerator_norm = count_leading_zeros(numerator_h);

        if numerator_norm >= 14 {
            best_fixed_codebook_gain = div64(sshl64(numerator, 14), denominator) as Word32;
        } else {
            let shifted_denominator: Word64 = shr64(denominator, (14 - numerator_norm) as u32); /* bestFixedCodebookGain in Q14 */
            if shifted_denominator > 0 {
                /* can't shift left by 9 the numerator, can we shift right by 9-numeratorNorm the denominator without hiting 0 */
                best_fixed_codebook_gain =
                    div64(shl64(numerator, numerator_norm as u32), shifted_denominator) as Word32;
            /* bestFixedCodebookGain in Q14 */
            } else {
                best_fixed_codebook_gain = shl32(
                    div64(shl64(numerator, numerator_norm as u32), denominator) as Word32,
                    (14 - numerator_norm) as u32,
                ); /* shift left the division result to reach Q14 */
            }
        }
    }

    /*** Compute the predicted gain as in spec 3.9.1 eq71 in Q6 ***/
    predicted_fixed_codebook_gain = shr32(
        ma_code_gain_prediction(previous_gain_prediction_error, fixed_codebook_vector),
        12,
    ) as Word16; /* in Q16 -> Q4 range [3,1830] */

    /***  preselection spec 3.9.2 ***/
    /* Note: spec just says to select the best 50% of each vector, ITU code go through magical constant computation to select the begining of a continuous range */
    /* much more simple here : vector are ordened in growing order so just select 2 (4 for Gb) indexes before the first value to be superior to the best gain previously computed */
    while index_base_ga < 6
        && best_fixed_codebook_gain
            > mult16_16_q14(GA_CODEBOOK[index_base_ga][1], predicted_fixed_codebook_gain) as Word32
    {
        /* bestFixedCodebookGain> in Q2, GACodebook in Q12 *predictedFixedCodebookGain in Q4 -> Q16-14 */
        index_base_ga += 1;
    }
    if index_base_ga > 0 {
        index_base_ga -= 1;
    }
    if index_base_ga > 0 {
        index_base_ga -= 1;
    }
    while index_base_gb < 12
        && best_adaptative_codebook_gain > shr(GB_CODEBOOK[index_base_gb][0] as Word32, 5)
    {
        index_base_gb += 1;
    }
    if index_base_gb > 0 {
        index_base_gb -= 1;
    }
    if index_base_gb > 0 {
        index_base_gb -= 1;
    }
    if index_base_gb > 0 {
        index_base_gb -= 1;
    }
    if index_base_gb > 0 {
        index_base_gb -= 1;
    }

    /*** test all possibilities of Ga and Gb indexes and select the best one ***/
    xy = -sshl(xy, 1); /* xy term is always used with a -2 factor */
    xz = -sshl(xz, 1); /* xz term is always used with a -2 factor */
    yz = sshl(yz, 1); /* yz term is always used with a 2 factor */

    for i in 0..4 {
        for j in 0..8 {
            /* compute gamma->gc and gp */
            let gp: Word16 = add16(
                GA_CODEBOOK[i + index_base_ga][0],
                GB_CODEBOOK[j + index_base_gb][0],
            ); /* result in Q14 */
            let gamma: Word16 = add16(
                GA_CODEBOOK[i + index_base_ga][1],
                GB_CODEBOOK[j + index_base_gb][1],
            ); /* result in Q3.12 (range [0.185, 5.05])*/
            let gc: Word32 = mult16_16_q14(gamma, predicted_fixed_codebook_gain); /* gamma in Q12, predictedFixedCodebookGain in Q4 -> Q16 -14 -> Q2 */

            /* compute E as in eq63 (first term excluded) */
            let mut acc: Word64 = mult32_32(mult16_16(gp, gp), yy); /* acc = gp^2*yy  gp in Q14, yy in Q0 -> acc in Q28 */
            acc = mac64(acc, mult16_16(gc as Word16, gc as Word16), zz); /* gc in Q2, zz in Q24 -> acc in Q28, note gc is on 32 bits but in a range making gc^2 fitting on 32 bits */
            acc = mac64(acc, shl32(gp as Word32, 14), xy); /* gp in Q14 shifted to Q28, xy in Q0 -> acc in Q28 */
            acc = mac64(acc, shl32(gc, 14), xz); /* gc in Q2 shifted to Q16, xz in Q12 -> acc in Q28 */
            acc = mac64(acc, mult16_16(gp, gc as Word16), yz); /* gp in Q14, gc in Q2 yz in Q12 -> acc in Q28 */

            if acc < distance_min {
                distance_min = acc;
                index_ga = i + index_base_ga;
                index_gb = j + index_base_gb;
                *quantized_adaptative_codebook_gain = gp;
                *quantized_fixed_codebook_gain = shr(gc, 1) as Word16;
            }
        }
    }

    /* update the previous gain prediction error */
    compute_gain_prediction_error(
        add16(GA_CODEBOOK[index_ga][1], GB_CODEBOOK[index_gb][1]),
        previous_gain_prediction_error,
    );

    /* mapping of indexes */
    *gain_codebook_stage1 = INDEX_MAPPING_GA[index_ga];
    *gain_codebook_stage2 = INDEX_MAPPING_GB[index_gb];
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gain_quantization() {
        let target_signal: [Word16; 40] = [0; 40];
        let filtered_adaptative_codebook_vector: [Word16; 40] = [0; 40];
        let convolved_fixed_codebook_vector: [Word16; 40] = [0; 40];
        let fixed_codebook_vector: [Word16; 40] = [0; 40];
        let xy64: Word64 = 0;
        let yy64: Word64 = 0;
        let mut previous_gain_prediction_error: [Word16; 4] = [-14336, -14336, -14336, -14336];
        let mut quantized_adaptative_codebook_gain: Word16 = 0;
        let mut quantized_fixed_codebook_gain: Word16 = 0;
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
