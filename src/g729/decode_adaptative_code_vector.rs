use crate::g729::basic_operations::*;
use crate::g729::codebooks::B30;
use crate::g729::ld8k::*;

pub fn init_decode_adaptative_code_vector() -> i16 {
    60
}

/*****************************************************************************/
/* computeAdaptativeCodeVector : as in spec 4.1.3                            */

/*      -(i/o) excitationVector : in Q0 excitation accessed from             */

/*****************************************************************************/
pub fn compute_adaptative_codebook_vector(
    excitation_vector: &mut [i16],
    mut frac_pitch_delay: i16,
    int_pitch_delay: i16,
    current_subframe_offset: usize,
) {
    let excitation_vector_minus_k_idx: usize;

    /* from spec 4.1.3 and 3.7.1 */

    if frac_pitch_delay == 1 {
        excitation_vector_minus_k_idx = current_subframe_offset - (int_pitch_delay + 1) as usize;
        frac_pitch_delay = 2;
    } else {
        frac_pitch_delay = -frac_pitch_delay;

        excitation_vector_minus_k_idx = current_subframe_offset - int_pitch_delay as usize;
    }

    for n in 0..L_SUBFRAME {
        let excitation_vector_n_minus_k_idx = excitation_vector_minus_k_idx + n;
        let excitation_vector_n_minus_k_plus_one_idx = excitation_vector_minus_k_idx + n + 1;

        let b301_idx = frac_pitch_delay as usize;
        let b302_idx = (3 - frac_pitch_delay) as usize;

        let mut acc: i32 = 0; /* in Q15 */
        let mut j = 0;
        for i in 0..10 {
            acc = mac16_16(
                acc,
                excitation_vector[excitation_vector_n_minus_k_idx - i],
                B30[b301_idx + j],
            );

            acc = mac16_16(
                acc,
                excitation_vector[excitation_vector_n_minus_k_plus_one_idx + i],
                B30[b302_idx + j],
            );
            j += 3;
        }

        excitation_vector[current_subframe_offset + n] =
            saturate(pshr(acc, 15), MAX_INT16 as i32) as i16;
    }
}

/*****************************************************************************/
/* decodeAdaptativeCodeVector : as in spec 4.1.3                             */

/*      -(i) adaptativeCodebookIndex : parameter P1 or P2                    */
/*      -(i) parityFlag : based on P1 parity flag : set if parity error      */

/*             P1 on subframe 1. On Subframe 2, contains the intPitchDelay   */

/*      -(i/o) excitationVector : in Q0 excitation accessed from             */

/*****************************************************************************/
pub fn decode_adaptative_code_vector(
    previous_int_pitch_delay: &mut i16,
    sub_frame_index: usize,
    adaptative_codebook_index: u16,
    parity_flag: u8,
    frame_erasure_flag: u8,
    int_pitch_delay: &mut i16,
    excitation_vector: &mut [i16],
) {
    let frac_pitch_delay: i16;

    /*** Compute the Pitch Delay from the Codebook index ***/

    if sub_frame_index == 0 {
        if (parity_flag | frame_erasure_flag) != 0 {
            *int_pitch_delay = *previous_int_pitch_delay; /* set the integer part of Pitch Delay to the last second subframe Pitch Delay computed spec: 4.1.2 */
            /* Note: unable to find anything regarding this part in the spec, just copied it from the ITU source code */
            frac_pitch_delay = 0;
            *previous_int_pitch_delay += 1;
            if *previous_int_pitch_delay > MAXIMUM_INT_PITCH_DELAY as i16 {
                *previous_int_pitch_delay = MAXIMUM_INT_PITCH_DELAY as i16;
            }
        } else {
            /* parity and frameErasure flags are off, do the normal computation (doc 4.1.3) */
            if adaptative_codebook_index < 197 {
                /* *intPitchDelay = (P1 + 2 )/ 3 + 19 */
                *int_pitch_delay = add16(
                    mult16_16_q15(add16(adaptative_codebook_index as i16, 2), 10923_i16) as i16,
                    19,
                ); /* MULT in Q15: 1/3 in Q15: 10923 */
                /* fracPitchDelay = P1 − 3*intPitchDelay  + 58 : fracPitchDelay in -1, 0, 1 */
                frac_pitch_delay = add16(
                    sub16(
                        adaptative_codebook_index as i16,
                        mult16_16(*int_pitch_delay, 3) as i16,
                    ),
                    58,
                );
            } else {
                *int_pitch_delay = sub16(adaptative_codebook_index as i16, 112);
                frac_pitch_delay = 0;
            }

            *previous_int_pitch_delay = *int_pitch_delay;
        }
    } else {
        if frame_erasure_flag != 0 {
            /* unable to find anything regarding this part in the spec, just copied it from the ITU source code */
            *int_pitch_delay = *previous_int_pitch_delay;
            frac_pitch_delay = 0;
            *previous_int_pitch_delay += 1;
            if *previous_int_pitch_delay > MAXIMUM_INT_PITCH_DELAY as i16 {
                *previous_int_pitch_delay = MAXIMUM_INT_PITCH_DELAY as i16;
            }
        } else {
            /* frameErasure flags are off, do the normal computation (doc 4.1.3) */
            let t_min = sub16(*int_pitch_delay, 5).clamp(20, 134);
            /* intPitchDelay = (P2 + 2 )/ 3 − 1 */
            *int_pitch_delay = sub16(
                mult16_16_q15(add16(adaptative_codebook_index as i16, 2), 10923_i16) as i16,
                1,
            );
            /* fracPitchDelay = P2 − 2 − 3((P 2 + 2 )/ 3 − 1) */
            frac_pitch_delay = sub16(
                sub16(
                    adaptative_codebook_index as i16,
                    mult16_16(*int_pitch_delay, 3) as i16,
                ),
                2,
            );
            /* *intPitchDelay = (P2 + 2 )/ 3 − 1 + tMin */
            *int_pitch_delay = add16(*int_pitch_delay, t_min);

            *previous_int_pitch_delay = *int_pitch_delay;
        }
    }

    compute_adaptative_codebook_vector(
        excitation_vector,
        frac_pitch_delay,
        *int_pitch_delay,
        L_PAST_EXCITATION + sub_frame_index,
    );
}
