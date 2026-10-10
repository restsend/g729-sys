use crate::g729::basic_operations::*;
use crate::g729::fixed_point_math::*;
use crate::g729::ld8k::*;
use crate::g729::utils::count_leading_zeros;

/*****************************************************************************/

/*      -(i) inputSignal: 223 values in Q0, buffer accessed in range         */

/*      -the correlation in Q0 on 32 bits                                    */

/*****************************************************************************/
#[cfg_attr(target_arch = "xtensa", inline(never))]
fn get_correlation(buffer: &[i16], current_frame_offset: usize, index: usize) -> i32 {
    let mut correlation: i32 = 0;

    for i in (0..L_FRAME).step_by(2) {
        correlation = mac16_16(
            correlation,
            buffer[current_frame_offset + i],
            buffer[current_frame_offset + i - index],
        );
    }
    correlation
}

/*****************************************************************************/
/* getCorrelation : compute eqA.4 from spec A3.4 on the given range and      */

/*      -(i) inputSignal: signal used to compute the correlation, in Q0      */

/*      - the correlation maximum found on the given range in Q0 on 32 bits  */

/*****************************************************************************/
#[cfg_attr(target_arch = "xtensa", inline(never))]
fn get_correlation_max(
    index: &mut usize,
    buffer: &[i16],
    current_frame_offset: usize,
    range_open: usize,
    range_close: usize,
    step: usize,
) -> i32 {
    let mut correlation_max = MIN_32;

    for i in (range_open..=range_close).step_by(step) {
        let correlation = get_correlation(buffer, current_frame_offset, i);
        if correlation > correlation_max {
            *index = i;
            correlation_max = correlation;
        }
    }

    correlation_max
}

/*****************************************************************************/

/*      -(i) weightedInputSignal: 223 values in Q0, buffer                   */

/*****************************************************************************/
#[cfg_attr(target_arch = "xtensa", inline(never))]
pub fn find_open_loop_pitch_delay(weighted_input_signal: &[i16]) -> u16 {
    let mut scaled_weighted_input_signal_buffer = [0; MAXIMUM_INT_PITCH_DELAY + L_FRAME];
    let mut autocorrelation: i64 = 0;
    let mut index_range1 = 0;
    let mut index_range2 = 0;
    let mut index_range3_even = 0;
    let mut index_range3;

    let mut correlation_max_range3;
    let mut correlation_max_range3_odd;
    let mut auto_correlation_range1;
    let mut auto_correlation_range2;
    let mut auto_correlation_range3;
    let mut normalised_correlation_max_range1;
    let mut normalised_correlation_max_range2;

    let mut index_multiple;

    let current_frame_offset = MAXIMUM_INT_PITCH_DELAY;

    #[allow(clippy::needless_range_loop)] // hot DSP loop
    for i in 0..MAXIMUM_INT_PITCH_DELAY + L_FRAME {
        autocorrelation = mac64(
            autocorrelation,
            weighted_input_signal[i] as i32,
            weighted_input_signal[i] as i32,
        );
    }

    let use_scaled_buffer = if autocorrelation > MAX_32 as i64 {
        let overflow_scale = pshr(
            sub32(
                31,
                count_leading_zeros(shr64(autocorrelation, 31) as i32) as i32,
            ),
            1,
        );
        for i in 0..MAXIMUM_INT_PITCH_DELAY + L_FRAME {
            scaled_weighted_input_signal_buffer[i] =
                shr16(weighted_input_signal[i], overflow_scale as u32);
        }
        true
    } else {
        false
    };

    let buffer = if use_scaled_buffer {
        &scaled_weighted_input_signal_buffer
    } else {
        weighted_input_signal
    };

    /*** compute the correlationMax in the different ranges ***/
    let correlation_max_range1 =
        get_correlation_max(&mut index_range1, buffer, current_frame_offset, 20, 39, 1);
    let correlation_max_range2 =
        get_correlation_max(&mut index_range2, buffer, current_frame_offset, 40, 79, 1);
    correlation_max_range3 = get_correlation_max(
        &mut index_range3_even,
        buffer,
        current_frame_offset,
        80,
        143,
        2,
    );
    index_range3 = index_range3_even;

    /* for the third range, correlationMax shall be computed at +1 and -1 around the maximum found as described in spec A3.4 */
    if index_range3 > 80 {
        correlation_max_range3_odd =
            get_correlation(buffer, current_frame_offset, index_range3 - 1);
        if correlation_max_range3_odd > correlation_max_range3 {
            correlation_max_range3 = correlation_max_range3_odd;
            index_range3 = index_range3_even - 1;
        }
    }
    correlation_max_range3_odd = get_correlation(buffer, current_frame_offset, index_range3 + 1);
    if correlation_max_range3_odd > correlation_max_range3 {
        correlation_max_range3 = correlation_max_range3_odd;
        index_range3 = index_range3_even + 1;
    }

    /*** normalise the correlations ***/

    auto_correlation_range1 = get_correlation(buffer, current_frame_offset - index_range1, 0);
    auto_correlation_range2 = get_correlation(buffer, current_frame_offset - index_range2, 0);
    auto_correlation_range3 = get_correlation(buffer, current_frame_offset - index_range3, 0);

    if auto_correlation_range1 == 0 {
        auto_correlation_range1 = 1;
    }
    if auto_correlation_range2 == 0 {
        auto_correlation_range2 = 1;
    }
    if auto_correlation_range3 == 0 {
        auto_correlation_range3 = 1;
    }

    /* according to ITU code comments, the normalisedCorrelationMax values fit on 16 bits when in Q0, so keep them in Q8 on 32 bits shall not give any overflow */
    normalised_correlation_max_range1 = mult32_32_q23(
        correlation_max_range1,
        g729_inv_sqrt_q0q31(auto_correlation_range1),
    );
    normalised_correlation_max_range2 = mult32_32_q23(
        correlation_max_range2,
        g729_inv_sqrt_q0q31(auto_correlation_range2),
    );
    let normalised_correlation_max_range3 = mult32_32_q23(
        correlation_max_range3,
        g729_inv_sqrt_q0q31(auto_correlation_range3),
    );

    /*** Favouring the delays with the values in the lower range ***/
    /* not clearly documented in spec A3.4, algo from the ITU code */
    index_multiple = shl(index_range2 as i32, 1) as usize;
    if abs(index_multiple as i32 - index_range3 as i32) < 5 {
        normalised_correlation_max_range2 = add32(
            normalised_correlation_max_range2,
            shr(normalised_correlation_max_range3, 2),
        );
    }

    if abs(index_multiple as i32 + index_range2 as i32 - index_range3 as i32) < 7 {
        normalised_correlation_max_range2 = add32(
            normalised_correlation_max_range2,
            shr(normalised_correlation_max_range3, 2),
        );
    }

    index_multiple = shl(index_range1 as i32, 1) as usize;
    if abs(index_multiple as i32 - index_range2 as i32) < 5 {
        normalised_correlation_max_range1 = mac16_32_p15(
            normalised_correlation_max_range1,
            O2_IN_Q15,
            normalised_correlation_max_range2,
        );
    }

    if abs(index_multiple as i32 + index_range1 as i32 - index_range2 as i32) < 7 {
        normalised_correlation_max_range1 = mac16_32_p15(
            normalised_correlation_max_range1,
            O2_IN_Q15,
            normalised_correlation_max_range2,
        );
    }

    /*** return the index corresponding to the greatest normalised Correlation */
    if normalised_correlation_max_range1 < normalised_correlation_max_range2 {
        normalised_correlation_max_range1 = normalised_correlation_max_range2;
        index_range1 = index_range2;
    }
    if normalised_correlation_max_range1 < normalised_correlation_max_range3 {
        index_range1 = index_range3;
    }
    index_range1 as u16
}
