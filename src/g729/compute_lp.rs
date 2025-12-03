use crate::g729::basic_operations::*;
use crate::g729::codebooks::{WLAG, WLP};
use crate::g729::ld8k::*;
use crate::g729::utils::count_leading_zeros;

pub fn auto_correlation_2_lp(
    auto_correlation_coefficients: &[Word32],
    lp_coefficients_q12: &mut [Word16],
    reflection_coefficients: &mut [Word32],
    residual_energy: &mut Word32,
) {
    let mut previous_iteration_lp_coefficients = [0 as Word32; NB_LSP_COEFF + 1];
    let mut lp_coefficients = [0 as Word32; NB_LSP_COEFF + 1];
    let mut e: Word32;
    let mut sum: Word32;

    /* init */
    lp_coefficients[0] = ONE_IN_Q27;
    lp_coefficients[1] = -div32_32_q27(
        auto_correlation_coefficients[1],
        auto_correlation_coefficients[0],
    ) as Word32;
    reflection_coefficients[0] = sshl(lp_coefficients[1], 4); /* k[0] is -r1/r0 in Q31 */

    /* E = r0(1 - a[1]^2) in Q31 */
    e = mult32_32_q31(
        auto_correlation_coefficients[0],
        sub32(
            ONE_IN_Q31,
            mult32_32_q23(lp_coefficients[1], lp_coefficients[1]),
        ),
    );

    for i in 2..=NB_LSP_COEFF {
        /* update the previousIterationLPCoefficients needed for this one */
        for j in 1..i {
            previous_iteration_lp_coefficients[j] = lp_coefficients[j];
        }

        /* sum = r[i] + ∑ a[j]*r[i-j] with j = 1..i-1 (a[0] is always 1) */
        sum = 0;
        for j in 1..i {
            sum = mac32_32_q31(
                sum,
                lp_coefficients[j],
                auto_correlation_coefficients[i - j],
            );
        }
        sum = add32(sshl(sum, 4), auto_correlation_coefficients[i]); /* set sum in Q31 and add r[0] */

        /* a[i] = -sum/E */
        lp_coefficients[i] = -div32_32_q31(sum, e) as Word32;
        reflection_coefficients[i - 1] = lp_coefficients[i];

        /* iterations j = 1..i-1 */
        /* a[j] += a[i]*a[i-j] */
        for j in 1..i {
            lp_coefficients[j] = mac32_32_q31(
                lp_coefficients[j],
                lp_coefficients[i],
                previous_iteration_lp_coefficients[i - j],
            );
        }

        /* E *=(1-a[i]^2) */
        e = mult32_32_q31(
            e,
            sub32(
                ONE_IN_Q31,
                mult32_32_q31(lp_coefficients[i], lp_coefficients[i]),
            ),
        );

        /* set LPCoefficients[i] from Q31 to Q27 */
        lp_coefficients[i] = shr(lp_coefficients[i], 4);
    }
    *residual_energy = e;

    /* convert with rounding the LP Coefficients form Q27 to Q12, ignore first coefficient which is always 1 */
    for i in 0..NB_LSP_COEFF {
        lp_coefficients_q12[i] =
            saturate(pshr(lp_coefficients[i + 1], 15), MAXINT16 as Word32) as Word16;
    }
}

pub fn compute_lp(
    signal: &[Word16],
    lp_coefficients_q12: &mut [Word16],
    reflection_coefficients: &mut [Word32],
    auto_correlation_coefficients: &mut [Word32],
    no_lag_auto_correlation_coefficients: &mut [Word32],
    auto_correlation_coefficients_scale: &mut i8,
    mut auto_correlation_coefficients_number: usize,
) {
    let mut windowed_signal = [0 as Word16; L_LP_ANALYSIS_WINDOW];
    let mut acc64: Word64 = 0;
    let mut right_shift_to_normalise = 0;
    let mut residual_energy: Word32 = 0;

    /* Compute the windowed signal */
    for i in 0..L_LP_ANALYSIS_WINDOW {
        windowed_signal[i] = mult16_16_p15(signal[i], WLP[i]) as Word16;
    }

    /* Compute the autoCorrelation coefficients r[0..10] */
    for i in 0..L_LP_ANALYSIS_WINDOW {
        acc64 = mac64(
            acc64,
            windowed_signal[i] as Word32,
            windowed_signal[i] as Word32,
        );
    }
    if acc64 == 0 {
        acc64 = 1;
    }

    /* normalise the acc64 on 32 bits */
    if acc64 > MAXINT32 as Word64 {
        while acc64 > MAXINT32 as Word64 {
            acc64 = shr64(acc64, 1);
            right_shift_to_normalise += 1;
        }
        auto_correlation_coefficients[0] = acc64 as Word32;
    } else {
        right_shift_to_normalise = -(count_leading_zeros(acc64 as Word32) as i32);
        auto_correlation_coefficients[0] =
            sshl(acc64 as Word32, (-right_shift_to_normalise) as u32);
    }

    *auto_correlation_coefficients_scale = -right_shift_to_normalise as i8;

    if right_shift_to_normalise > 0 {
        for i in 1..auto_correlation_coefficients_number {
            acc64 = 0;
            for j in i..L_LP_ANALYSIS_WINDOW {
                acc64 = add64_32(acc64, mult16_16(windowed_signal[j], windowed_signal[j - i]));
            }
            auto_correlation_coefficients[i] =
                shr64(acc64, right_shift_to_normalise as u32) as Word32;
        }
    } else {
        for i in 1..auto_correlation_coefficients_number {
            let mut acc32: Word32 = 0;
            for j in i..L_LP_ANALYSIS_WINDOW {
                acc32 = mac16_16(acc32, windowed_signal[j], windowed_signal[j - i]);
            }
            auto_correlation_coefficients[i] = sshl(acc32, (-right_shift_to_normalise) as u32);
        }
    }

    /* save autocorrelation before applying lag window */
    for i in 0..auto_correlation_coefficients_number {
        no_lag_auto_correlation_coefficients[i] = auto_correlation_coefficients[i];
    }

    if auto_correlation_coefficients_number > NB_LSP_COEFF + 3 {
        auto_correlation_coefficients_number = NB_LSP_COEFF + 3;
    }

    for i in 1..auto_correlation_coefficients_number {
        auto_correlation_coefficients[i] = mult16_32_p15(WLAG[i], auto_correlation_coefficients[i]);
    }

    auto_correlation_2_lp(
        auto_correlation_coefficients,
        lp_coefficients_q12,
        reflection_coefficients,
        &mut residual_energy,
    );
}
