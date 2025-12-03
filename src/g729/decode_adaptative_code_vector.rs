use crate::g729::basic_operations::*;
use crate::g729::codebooks::B30;
use crate::g729::ld8k::*;

/* init function */
pub fn init_decode_adaptative_code_vector(previous_int_pitch_delay: &mut i16) {
    *previous_int_pitch_delay = 60;
}

/*****************************************************************************/
/* computeAdaptativeCodeVector : as in spec 4.1.3                            */
/*    parameters:                                                            */
/*      -(i/o) excitationVector : in Q0 excitation accessed from             */
/*             [-MAXIMUM_INT_PITCH_DELAY(143), -1] as input                  */
/*             and [0, L_SUBFRAME[ as output to store the adaptative         */
/*             codebook vector                                               */
/*      -(i/o) fracPitchDelay : the fractionnal part of Pitch Delay.         */
/*      -(i/o) intPitchDelay : the integer part of Pitch Delay.              */
/*                                                                           */
/*****************************************************************************/
fn compute_adaptative_codebook_vector(
    excitation_vector: &mut [Word16],
    mut frac_pitch_delay: i16,
    int_pitch_delay: i16,
    current_subframe_offset: usize,
) {
    let excitation_vector_minus_k_idx: usize;

    /* compute the adaptative codebook vector using the pitch delay and the past excitation vector */
    /* from spec 4.1.3 and 3.7.1 */
    /* shall compute v(n ) = ∑ u (n - k + i )b30 (t + 3i ) + ∑ u (n - k + 1 + i )b30 (3 - t + 3i ) for i=0,...,9 and n = 0,...,39 (t in 0, 1, 2) */
    /* with k = intPitchDelay and t = fracPitchDelay wich must be converted from range -1,0,1 to 0,1,2 */
    /* u the past excitation vector */
    /* v the adaptative codebook vector */
    /* b30 an interpolation filter */

    /* scale fracPichDelay from -1,0.1 to 0,1,2 */
    if frac_pitch_delay == 1 {
        // excitationVectorMinusK = &(excitationVector[-(intPitchDelay+1)]); /* fracPitchDelay being positive -> increase by one the integer part and set to 2 the fractional part : -(k+1/3) -> -(k+1)+2/3 */
        excitation_vector_minus_k_idx = current_subframe_offset - (int_pitch_delay + 1) as usize;
        frac_pitch_delay = 2;
    } else {
        frac_pitch_delay = -frac_pitch_delay; /* 0 unchanged, -1 -> +1 */
        // excitationVectorMinusK = &(excitationVector[-intPitchDelay]); /* -(k-1/3) -> -k+1/3  or -(k) -> -k*/
        excitation_vector_minus_k_idx = current_subframe_offset - int_pitch_delay as usize;
    }

    for n in 0..L_SUBFRAME {
        /* loop over the whole subframe */
        // word16_t *excitationVectorNMinusK = &(excitationVectorMinusK[n]); /* point to u(n-k), unscaled value, full range */
        // word16_t *excitationVectorNMinusKPlusOne = &(excitationVectorMinusK[n+1]); /* point to u(n-k+1), unscaled value, full range */
        let excitation_vector_n_minus_k_idx = excitation_vector_minus_k_idx + n;
        let excitation_vector_n_minus_k_plus_one_idx = excitation_vector_minus_k_idx + n + 1;

        // word16_t *b301 = &(b30[fracPitchDelay]); /* point to b30(t) in Q0.15 : sums of all b30 coeffs is < 2, no overflow possible on 32 bits */
        // word16_t *b302 = &(b30[3-fracPitchDelay]); /* point to b30(3-t) in Q0.15*/
        let b301_idx = frac_pitch_delay as usize;
        let b302_idx = (3 - frac_pitch_delay) as usize;

        let mut acc: Word32 = 0; /* in Q15 */
        let mut j = 0;
        for i in 0..10 {
            // acc = MAC16_16(acc, excitationVectorNMinusK[-i], b301[j]); /*  Note : the spec says: u(n−k+i)b30(t+3i) but the ITU code do (and here too) u(n-k-i )b30(t+3i) */
            acc = mac16_16(
                acc,
                excitation_vector[excitation_vector_n_minus_k_idx - i],
                B30[b301_idx + j],
            );
            // acc = MAC16_16(acc, excitationVectorNMinusKPlusOne[i], b302[j]); /* u(n-k+1+i)b30(3-t+3i) */
            acc = mac16_16(
                acc,
                excitation_vector[excitation_vector_n_minus_k_plus_one_idx + i],
                B30[b302_idx + j],
            );
            j += 3;
        }
        // excitationVector[n] = SATURATE(PSHR(acc, 15), MAXINT16); /* acc in Q15, shift/round to unscaled value and check overflow on 16 bits */
        excitation_vector[current_subframe_offset + n] =
            saturate(pshr(acc, 15), MAX_INT16 as Word32) as Word16;
    }
}

/*****************************************************************************/
/* decodeAdaptativeCodeVector : as in spec 4.1.3                             */
/*    parameters:                                                            */
/*      -(i/o) decoderChannelContext : the channel context data              */
/*      -(i) subFrameIndex : 0 or 40 for subframe 1 or subframe 2            */
/*      -(i) adaptativeCodebookIndex : parameter P1 or P2                    */
/*      -(i) parityFlag : based on P1 parity flag : set if parity error      */
/*      -(i) frameErasureFlag : set in case of frame erasure                 */
/*      -(i/o) intPitchDelay : the integer part of Pitch Delay. Computed from*/
/*             P1 on subframe 1. On Subframe 2, contains the intPitchDelay   */
/*             computed on Subframe 1.                                       */
/*      -(i/o) excitationVector : in Q0 excitation accessed from             */
/*             [-MAXIMUM_INT_PITCH_DELAY(143), -1] as input                  */
/*             and [0, L_SUBFRAME[ as output to store the adaptative         */
/*             codebook vector                                               */
/*                                                                           */
/*****************************************************************************/
pub fn decode_adaptative_code_vector(
    previous_int_pitch_delay: &mut i16,
    sub_frame_index: usize,
    adaptative_codebook_index: u16,
    parity_flag: u8,
    frame_erasure_flag: u8,
    int_pitch_delay: &mut i16,
    excitation_vector: &mut [Word16],
) {
    let frac_pitch_delay: i16;

    /*** Compute the Pitch Delay from the Codebook index ***/
    /* fracPitchDelay is computed in the range -1,0,1 */
    if sub_frame_index == 0 {
        /* first subframe */
        if (parity_flag | frame_erasure_flag) != 0 {
            /* there is an error (either parity or frame erased) */
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
                    mult16_16_q15(add16(adaptative_codebook_index as i16, 2), 10923 as Word16)
                        as Word16,
                    19,
                ); /* MULT in Q15: 1/3 in Q15: 10923 */
                /* fracPitchDelay = P1 − 3*intPitchDelay  + 58 : fracPitchDelay in -1, 0, 1 */
                frac_pitch_delay = add16(
                    sub16(
                        adaptative_codebook_index as i16,
                        mult16_16(*int_pitch_delay, 3) as Word16,
                    ),
                    58,
                );
            } else {
                /* adaptativeCodebookIndex>= 197 */
                *int_pitch_delay = sub16(adaptative_codebook_index as i16, 112);
                frac_pitch_delay = 0;
            }

            /* backup the intPitchDelay */
            *previous_int_pitch_delay = *int_pitch_delay;
        }
    } else {
        /* second subframe */
        if frame_erasure_flag != 0 {
            /* there is an error : frame erased, in case of parity error, it has been taken in account at first subframe */
            /* unable to find anything regarding this part in the spec, just copied it from the ITU source code */
            *int_pitch_delay = *previous_int_pitch_delay;
            frac_pitch_delay = 0;
            *previous_int_pitch_delay += 1;
            if *previous_int_pitch_delay > MAXIMUM_INT_PITCH_DELAY as i16 {
                *previous_int_pitch_delay = MAXIMUM_INT_PITCH_DELAY as i16;
            }
        } else {
            /* frameErasure flags are off, do the normal computation (doc 4.1.3) */
            let mut t_min = sub16(*int_pitch_delay, 5); /* intPitchDelay contains the intPitch computed for subframe one */
            if t_min < 20 {
                t_min = 20;
            }
            if t_min > 134 {
                t_min = 134;
            }
            /* intPitchDelay = (P2 + 2 )/ 3 − 1 */
            *int_pitch_delay = sub16(
                mult16_16_q15(add16(adaptative_codebook_index as i16, 2), 10923 as Word16)
                    as Word16,
                1,
            );
            /* fracPitchDelay = P2 − 2 − 3((P 2 + 2 )/ 3 − 1) */
            frac_pitch_delay = sub16(
                sub16(
                    adaptative_codebook_index as i16,
                    mult16_16(*int_pitch_delay, 3) as Word16,
                ),
                2,
            );
            /* *intPitchDelay = (P2 + 2 )/ 3 − 1 + tMin */
            *int_pitch_delay = add16(*int_pitch_delay, t_min);

            /* backup the intPitchDelay */
            *previous_int_pitch_delay = *int_pitch_delay;
        }
    }

    /* compute the adaptative codebook vector using the pitch delay and the past excitation vector */
    compute_adaptative_codebook_vector(
        excitation_vector,
        frac_pitch_delay,
        *int_pitch_delay,
        L_PAST_EXCITATION + sub_frame_index,
    );
}
