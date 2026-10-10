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
    excitation_vector: &mut [i16],
    current_idx: usize,
    mut int_pitch_delay: i16,
    mut frac_pitch_delay: i16,
) {
    frac_pitch_delay = -frac_pitch_delay;
    if frac_pitch_delay < 0 {
        int_pitch_delay += 1;
        frac_pitch_delay = 2;
    }

    let delayed_idx = (current_idx as isize - int_pitch_delay as isize) as usize;
    let b30_increased_idx = frac_pitch_delay as usize;
    let b30_decreased_idx = (3 - frac_pitch_delay) as usize;

    let b30_inc = &B30[b30_increased_idx..b30_increased_idx + 28];
    let b30_dec = &B30[b30_decreased_idx..b30_decreased_idx + 28];

    // 10 taps, stride 3: materialise them as fixed-size arrays so the inner loop
    // is a fully unrollable, bounds-check-free 10-iteration loop.
    let b30_inc_taps: [i16; 10] = core::array::from_fn(|k| b30_inc[3 * k]);
    let b30_dec_taps: [i16; 10] = core::array::from_fn(|k| b30_dec[3 * k]);

    for n in 0..L_SUBFRAME {
        let mut acc: i32 = 0; // acc in Q15

        for (k, (&inc, &dec)) in b30_inc_taps.iter().zip(b30_dec_taps.iter()).enumerate() {
            acc = mac16_16(acc, excitation_vector[delayed_idx + n - k], inc);
            acc = mac16_16(acc, excitation_vector[delayed_idx + n + 1 + k], dec);
        }
        // acc in Q15, shift/round to unscaled value and check overflow on 16 bits
        excitation_vector[current_idx + n] = saturate(pshr(acc, 15), MAX_16 as i32) as i16;
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
#[allow(clippy::too_many_arguments)]
#[cfg_attr(target_arch = "xtensa", inline(never))]
pub fn adaptative_codebook_search(
    excitation_vector: &mut [i16],
    current_idx: usize,
    int_pitch_delay_min: &mut i16,
    int_pitch_delay_max: &mut i16,
    impulse_response: &[i16],
    target_signal: &[i16],
    int_pitch_delay: &mut i16,
    frac_pitch_delay: &mut i16,
    pitch_delay_codeword: &mut u16,
    sub_frame_index: u16,
) {
    let mut backward_filtered_target_signal = [0_i32; L_SUBFRAME];
    let mut correlation_max: i32 = i32::MIN;

    // compute the backward Filtered Target Signal as specified in A.3.7: correlation of target signal and impulse response

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

    generate_adaptative_codebook_vector(excitation_vector, current_idx, *int_pitch_delay, 0);

    *frac_pitch_delay = 0;
    if !(sub_frame_index == 0 && *int_pitch_delay >= 85) {
        let mut adaptative_codebook_vector_backup = [0_i16; L_SUBFRAME];

        correlation_max = dot_product_16_32_q12(
            &excitation_vector[current_idx..current_idx + L_SUBFRAME],
            &backward_filtered_target_signal,
        );

        adaptative_codebook_vector_backup
            .copy_from_slice(&excitation_vector[current_idx..current_idx + L_SUBFRAME]);

        generate_adaptative_codebook_vector(excitation_vector, current_idx, *int_pitch_delay, -1);
        let mut correlation = dot_product_16_32_q12(
            &excitation_vector[current_idx..current_idx + L_SUBFRAME],
            &backward_filtered_target_signal,
        );
        if correlation > correlation_max {
            *frac_pitch_delay = -1;
            correlation_max = correlation;

            adaptative_codebook_vector_backup
                .copy_from_slice(&excitation_vector[current_idx..current_idx + L_SUBFRAME]);
        }

        generate_adaptative_codebook_vector(excitation_vector, current_idx, *int_pitch_delay, 1);
        correlation = dot_product_16_32_q12(
            &excitation_vector[current_idx..current_idx + L_SUBFRAME],
            &backward_filtered_target_signal,
        );
        if correlation > correlation_max {
            *frac_pitch_delay = 1;
        } else {
            excitation_vector[current_idx..current_idx + L_SUBFRAME]
                .copy_from_slice(&adaptative_codebook_vector_backup);
        }
    }

    if sub_frame_index == 0 {
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
        // compute the codeword as in spec 3.7.2
        *pitch_delay_codeword =
            (3 * (*int_pitch_delay - *int_pitch_delay_min) + *frac_pitch_delay + 2) as u16;
    }
}
