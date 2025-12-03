use crate::g729::basic_operations::*;
use crate::g729::codebooks::B30;
use crate::g729::ld8k::*;
use crate::g729::utils::{correlate_vectors, dot_product_16_32_q12};

/// Generates the adaptative codebook vector by interpolation of past excitation
///
/// # Arguments
///
/// * `excitation_vector` - The excitation buffer. The current subframe starts at `current_idx`.
/// * `current_idx` - The index in `excitation_vector` where the current subframe starts.
/// * `int_pitch_delay` - The integer pitch delay.
/// * `frac_pitch_delay` - The fractional pitch delay (-1, 0, or 1).
pub fn generate_adaptative_codebook_vector(
    excitation_vector: &mut [Word16],
    current_idx: usize,
    mut int_pitch_delay: i16,
    mut frac_pitch_delay: i16,
) {
    // fracPitchDelay is in range [-1, 1], convert it to [0,2] needed by eqA.8
    frac_pitch_delay = -frac_pitch_delay;
    if frac_pitch_delay < 0 {
        // if fracPitchDelay is 1 -> pitchDelay of int+(1/3) -> int+1-(2/3)
        int_pitch_delay += 1;
        frac_pitch_delay = 2;
    }

    let delayed_idx = (current_idx as isize - int_pitch_delay as isize) as usize;
    let b30_increased_idx = frac_pitch_delay as usize;
    let b30_decreased_idx = (3 - frac_pitch_delay) as usize;

    for n in 0..L_SUBFRAME {
        let mut acc: Word32 = 0; // acc in Q15
        let mut j = 0;
        for i in 0..10 {
            // j is used as a 3*i index
            // WARNING: spec 3.7.1 and A.8 give an equation leading to delayedExcitationVector[n+i]
            // but ITU code uses delayedExcitationVector[n-i], implemented as code
            // Note: In Rust we need to be careful with usize subtraction.
            // delayedExcitationVector points to excitationVector[-intPitchDelay]
            // So delayedExcitationVector[n-i] is excitationVector[current_idx - intPitchDelay + n - i]

            let idx1 = (delayed_idx as isize + n as isize - i as isize) as usize;
            let idx2 = (delayed_idx as isize + n as isize + 1 + i as isize) as usize;

            acc = mac16_16(acc, excitation_vector[idx1], B30[b30_increased_idx + j]);
            acc = mac16_16(acc, excitation_vector[idx2], B30[b30_decreased_idx + j]);

            j += 3;
        }
        // acc in Q15, shift/round to unscaled value and check overflow on 16 bits
        excitation_vector[current_idx + n] =
            saturate(pshr(acc, 15) as Word32, MAX_16 as Word32) as Word16;
    }
}

/// Compute parameter P1 and P2 as in spec A.3.7
/// Compute also adaptative codebook vector as in spec 3.7.1
///
/// # Arguments
///
/// * `excitation_vector` - The excitation buffer. The current subframe starts at `current_idx`.
/// * `current_idx` - The index in `excitation_vector` where the current subframe starts.
/// * `int_pitch_delay_min` - Low boundary for pitch delay search.
/// * `int_pitch_delay_max` - High boundary for pitch delay search.
/// * `impulse_response` - 40 values as in spec A.3.5 in Q12.
/// * `target_signal` - 40 values as in spec A.3.6 in Q0.
/// * `int_pitch_delay` - Output integer pitch delay.
/// * `frac_pitch_delay` - Output fractional part of pitch delay.
/// * `pitch_delay_codeword` - Output P1 or P2 codeword as in spec 3.7.2.
/// * `sub_frame_index` - 0 for the first subframe, 40 for the second.
pub fn adaptative_codebook_search(
    excitation_vector: &mut [Word16],
    current_idx: usize,
    int_pitch_delay_min: &mut i16,
    int_pitch_delay_max: &mut i16,
    impulse_response: &[Word16],
    target_signal: &[Word16],
    int_pitch_delay: &mut i16,
    frac_pitch_delay: &mut i16,
    pitch_delay_codeword: &mut u16,
    sub_frame_index: u16,
) {
    let mut backward_filtered_target_signal = [0 as Word32; L_SUBFRAME];
    let mut correlation_max: Word32 = Word32::MIN;

    // compute the backward Filtered Target Signal as specified in A.3.7: correlation of target signal and impulse response
    // targetSignal in Q0, impulseResponse in Q12 -> backwardFilteredTargetSignal in Q12
    correlate_vectors(
        target_signal,
        impulse_response,
        &mut backward_filtered_target_signal,
    );

    // maximise the sum as in spec A.3.7, eq A.7
    for i in *int_pitch_delay_min..=*int_pitch_delay_max {
        let idx = (current_idx as isize - i as isize) as usize;
        let correlation = dot_product_16_32_q12(
            &excitation_vector[idx..idx + L_SUBFRAME],
            &backward_filtered_target_signal,
        );

        if correlation > correlation_max {
            correlation_max = correlation;
            *int_pitch_delay = i;
        }
    }

    // compute the adaptativeCodebookVector (with fracPitchDelay at 0)
    // output is in excitationVector[0,L_SUBRAME[ -> excitation_vector[current_idx..current_idx+L_SUBFRAME]
    generate_adaptative_codebook_vector(excitation_vector, current_idx, *int_pitch_delay, 0);

    // if we are at first subframe and intPitchDelay >= 85 -> do not compute fracPitchDelay, set it to 0
    *frac_pitch_delay = 0;
    if !(sub_frame_index == 0 && *int_pitch_delay >= 85) {
        // compute the fracPitchDelay
        let mut adaptative_codebook_vector_backup = [0 as Word16; L_SUBFRAME];

        // search the fractionnal part to get the best correlation
        // we already have in excitationVector for fracPitchDelay = 0 the adaptativeCodebookVector (see specA.3.7)
        correlation_max = dot_product_16_32_q12(
            &excitation_vector[current_idx..current_idx + L_SUBFRAME],
            &backward_filtered_target_signal,
        );
        // backup the adaptativeCodebookVector
        adaptative_codebook_vector_backup
            .copy_from_slice(&excitation_vector[current_idx..current_idx + L_SUBFRAME]);

        // Fractionnal part = -1
        generate_adaptative_codebook_vector(excitation_vector, current_idx, *int_pitch_delay, -1);
        let mut correlation = dot_product_16_32_q12(
            &excitation_vector[current_idx..current_idx + L_SUBFRAME],
            &backward_filtered_target_signal,
        );
        if correlation > correlation_max {
            // fractional part at -1 gives higher correlation
            *frac_pitch_delay = -1;
            correlation_max = correlation;
            // backup the adaptativeCodebookVector
            adaptative_codebook_vector_backup
                .copy_from_slice(&excitation_vector[current_idx..current_idx + L_SUBFRAME]);
        }

        // Fractionnal part = 1
        generate_adaptative_codebook_vector(excitation_vector, current_idx, *int_pitch_delay, 1);
        correlation = dot_product_16_32_q12(
            &excitation_vector[current_idx..current_idx + L_SUBFRAME],
            &backward_filtered_target_signal,
        );
        if correlation > correlation_max {
            // fractional part at 1 gives higher correlation
            *frac_pitch_delay = 1;
        } else {
            // previously computed fractional part gives better result
            // restore the adaptativeCodebookVector
            excitation_vector[current_idx..current_idx + L_SUBFRAME]
                .copy_from_slice(&adaptative_codebook_vector_backup);
        }
    }

    // compute the codeword and intPitchDelayMin/intPitchDelayMax if needed (first subframe only)
    if sub_frame_index == 0 {
        // first subframe
        // compute intPitchDelayMin/intPitchDelayMax as in spec A.3.7
        *int_pitch_delay_min = *int_pitch_delay - 5;
        if *int_pitch_delay_min < 20 {
            *int_pitch_delay_min = 20;
        }
        *int_pitch_delay_max = *int_pitch_delay_min + 9;
        if *int_pitch_delay_max > MAXIMUM_INT_PITCH_DELAY as i16 {
            *int_pitch_delay_max = MAXIMUM_INT_PITCH_DELAY as i16;
            *int_pitch_delay_min = MAXIMUM_INT_PITCH_DELAY as i16 - 9;
        }

        // compute the codeword as in spec 3.7.2
        if *int_pitch_delay <= 85 {
            *pitch_delay_codeword = (3 * (*int_pitch_delay) - 58 + *frac_pitch_delay) as u16;
        } else {
            *pitch_delay_codeword = (*int_pitch_delay + 112) as u16;
        }
    } else {
        // second subframe
        // compute the codeword as in spec 3.7.2
        *pitch_delay_codeword =
            (3 * (*int_pitch_delay - *int_pitch_delay_min) + *frac_pitch_delay + 2) as u16;
    }
}
