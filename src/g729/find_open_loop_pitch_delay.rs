use crate::g729::basic_operations::*;
use crate::g729::fixed_point_math::*;
use crate::g729::ld8k::*;
use crate::g729::utils::count_leading_zeros;

/*****************************************************************************/
/* getCorrelation : as specified in specA3.4 eqA.4                           */
/*      correlation = ∑(i=0..39)inputSignal[2*i]*inputSignal[2*i-index]      */
/*    paremeters:                                                            */
/*      -(i) inputSignal: 223 values in Q0, buffer accessed in range         */
/*           [-index, L_FRAME[                                               */
/*      -(i) index: integer value in range [20,143]                          */
/*    return value:                                                          */
/*      -the correlation in Q0 on 32 bits                                    */
/*                                                                           */
/*****************************************************************************/
fn get_correlation(buffer: &[Word16], current_frame_offset: usize, index: usize) -> Word32 {
    let mut correlation: Word32 = 0;
    // i goes from 0 to L_FRAME with step 2.
    // j goes from -index to L_FRAME-index with step 2.
    // effectively j = i - index.
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
/*      step(1 compute all the correlation in range, 2 only the even ones)   */
/*      then return the maximum as specified in specA3.4 eqA.4               */
/*    paremeters:                                                            */
/*      -(o) index : the index giving the maximum of correlation on the      */
/*           considered range                                                */
/*      -(i) inputSignal: signal used to compute the correlation, in Q0      */
/*           accessed in range [-rangeClose, L_FRAME[                        */
/*      -(i) rangeOpen and rangeClose : the index range in which looking for */
/*           the correlation max                                             */
/*      -(i) step : incrementing step for the index                          */
/*    return value :                                                         */
/*      - the correlation maximum found on the given range in Q0 on 32 bits  */
/*                                                                           */
/*****************************************************************************/
fn get_correlation_max(
    index: &mut usize,
    buffer: &[Word16],
    current_frame_offset: usize,
    range_open: usize,
    range_close: usize,
    step: usize,
) -> Word32 {
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
/* findOpenLoopPitchDelay : as specified in specA3.4                         */
/*    paremeters:                                                            */
/*      -(i) weightedInputSignal: 223 values in Q0, buffer                   */
/*           accessed in range [-MAXIMUM_INT_PITCH_DELAY(143), L_FRAME(80)[  */
/*    return value:                                                          */
/*      - the openLoopIntegerPitchDelay in Q0 range [20, 143]                */
/*                                                                           */
/*****************************************************************************/
pub fn find_open_loop_pitch_delay(weighted_input_signal: &[Word16]) -> u16 {
    let mut scaled_weighted_input_signal_buffer = [0; MAXIMUM_INT_PITCH_DELAY + L_FRAME];
    let mut autocorrelation: Word64 = 0;
    let mut index_range1 = 0;
    let mut index_range2 = 0;
    let mut index_range3_even = 0;
    let mut index_range3;
    let correlation_max_range1;
    let correlation_max_range2;
    let mut correlation_max_range3;
    let mut correlation_max_range3_odd;
    let mut auto_correlation_range1;
    let mut auto_correlation_range2;
    let mut auto_correlation_range3;
    let mut normalised_correlation_max_range1;
    let mut normalised_correlation_max_range2;
    let normalised_correlation_max_range3;
    let mut index_multiple;

    // weighted_input_signal is expected to be of size MAXIMUM_INT_PITCH_DELAY + L_FRAME (223).
    // The current frame starts at MAXIMUM_INT_PITCH_DELAY.
    let current_frame_offset = MAXIMUM_INT_PITCH_DELAY;

    /* compute on 64 bits the autocorrelation on the input signal and if needed scale to have it on 32 bits */
    for i in 0..MAXIMUM_INT_PITCH_DELAY + L_FRAME {
        // The loop in C is: for (i=-MAXIMUM_INT_PITCH_DELAY; i<L_FRAME; i++)
        // accessing weightedInputSignal[i].
        // In Rust, this corresponds to iterating over the whole buffer.
        autocorrelation = mac64(
            autocorrelation,
            weighted_input_signal[i] as Word32,
            weighted_input_signal[i] as Word32,
        );
    }

    let use_scaled_buffer;
    if autocorrelation > MAX_32 as Word64 {
        let overflow_scale = pshr(
            sub32(
                31,
                count_leading_zeros(shr64(autocorrelation, 31) as Word32) as Word32,
            ),
            1,
        ); /* count number of bits needed over the 31 bits allowed and divide by 2 to get the right scaling for the signal */
        for i in 0..MAXIMUM_INT_PITCH_DELAY + L_FRAME {
            scaled_weighted_input_signal_buffer[i] =
                shr16(weighted_input_signal[i], overflow_scale as u32);
        }
        use_scaled_buffer = true;
    } else {
        use_scaled_buffer = false;
    }

    let buffer = if use_scaled_buffer {
        &scaled_weighted_input_signal_buffer
    } else {
        weighted_input_signal
    };

    /*** compute the correlationMax in the different ranges ***/
    correlation_max_range1 =
        get_correlation_max(&mut index_range1, buffer, current_frame_offset, 20, 39, 1);
    correlation_max_range2 =
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
        /* don't test value out of range [80, 143] */
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
    // getCorrelation(&(scaledWeightedInputSignal[-indexRange1]), 0);
    // In Rust: get_correlation(buffer, current_frame_offset - index_range1, 0)
    auto_correlation_range1 = get_correlation(buffer, current_frame_offset - index_range1, 0);
    auto_correlation_range2 = get_correlation(buffer, current_frame_offset - index_range2, 0);
    auto_correlation_range3 = get_correlation(buffer, current_frame_offset - index_range3, 0);

    if auto_correlation_range1 == 0 {
        auto_correlation_range1 = 1; /* avoid division by 0 */
    }
    if auto_correlation_range2 == 0 {
        auto_correlation_range2 = 1; /* avoid division by 0 */
    }
    if auto_correlation_range3 == 0 {
        auto_correlation_range3 = 1; /* avoid division by 0 */
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
    normalised_correlation_max_range3 = mult32_32_q23(
        correlation_max_range3,
        g729_inv_sqrt_q0q31(auto_correlation_range3),
    );

    /*** Favouring the delays with the values in the lower range ***/
    /* not clearly documented in spec A3.4, algo from the ITU code */
    index_multiple = shl(index_range2 as Word32, 1) as usize; /* indexMultiple = 2*indexRange2 */
    if abs(index_multiple as Word32 - index_range3 as Word32) < 5 {
        /* 2*indexRange2 - indexRange3 < 5 */
        normalised_correlation_max_range2 = add32(
            normalised_correlation_max_range2,
            shr(normalised_correlation_max_range3, 2),
        ); /* Max2 += Max3*0.25 */
    }

    if abs(index_multiple as Word32 + index_range2 as Word32 - index_range3 as Word32) < 7 {
        /* 3*indexRange2 - indexRange3 < 5 */
        normalised_correlation_max_range2 = add32(
            normalised_correlation_max_range2,
            shr(normalised_correlation_max_range3, 2),
        ); /* Max2 += Max3*0.25 */
    }

    index_multiple = shl(index_range1 as Word32, 1) as usize; /* indexMultiple = 2*indexRange1 */
    if abs(index_multiple as Word32 - index_range2 as Word32) < 5 {
        /* 2*indexRange1 - indexRange2 < 5 */
        normalised_correlation_max_range1 = mac16_32_p15(
            normalised_correlation_max_range1,
            O2_IN_Q15,
            normalised_correlation_max_range2,
        ); /* Max1 += Max2*0.2 */
    }

    if abs(index_multiple as Word32 + index_range1 as Word32 - index_range2 as Word32) < 7 {
        /* 3*indexRange1 - indexRange2 < 7 */
        normalised_correlation_max_range1 = mac16_32_p15(
            normalised_correlation_max_range1,
            O2_IN_Q15,
            normalised_correlation_max_range2,
        ); /* Max1 += Max2*0.2 */
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
