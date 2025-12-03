use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;

/*****************************************************************************/
/* decodeFixedCodeVector : compute the fixed codebook vector as in spec 4.1.4*/
/*    parameters:                                                            */
/*      -(i) signs: parameter S(4 signs bit) eq61                            */
/*      -(i) positions: parameter C(4 3bits position and jx bit) eq62        */
/*      -(i) intPitchDelay: integer part of pitch Delay (T)                  */
/*      -(i) boundedPitchGain: Beta in eq47 and eq48, in Q14                 */
/*      -(o) fixedCodebookVector: 40 values in Q13, range:  [-1.8,+1.8]      */
/*                                                                           */
/*****************************************************************************/
pub fn decode_fixed_code_vector(
    mut signs: u16,
    mut positions: u16,
    int_pitch_delay: i16,
    bounded_pitch_gain: Word16,
    fixed_codebook_vector: &mut [Word16],
) {
    let mut positions_array = [0u16; 4];
    let jx: u16;

    /* get the positions into an array: mapping according to eq62 and table7 in spec 3.8 */
    positions_array[0] = (positions & 7) * 5; /* m0 = 5*C, do not use macro here as whatever fixed or floating point computation we use, these are integers */
    positions = shr16(positions as i16, 3) as u16; /* shift right by 3 to get the 3 next bits(m1/5) as LSB */
    positions_array[1] = ((positions & 7) * 5) + 1; /* m1 = 5*C + 1, do not use macro here as whatever fixed or floating point computation we use, these are integers */
    positions = shr16(positions as i16, 3) as u16; /* shift right by 3 to get the 3 next bits(m2/5) as LSB */
    positions_array[2] = ((positions & 7) * 5) + 2; /* m2 = 5*C + 2, do not use macro here as whatever fixed or floating point computation we use, these are integers */
    positions = shr16(positions as i16, 3) as u16; /* shift right by 3 to get the last 4 bits(m3/5 and jx) as LSB */
    jx = positions & 1; /* jx from eq62 is the last bit */
    positions = shr16(positions as i16, 1) as u16; /* shift right by 1 to get the last 3 bits as LSB */
    positions_array[3] = ((positions & 7) * 5) + 3 + jx; /* m3 = 5*C + 3 + jx, do not use macro here as whatever fixed or floating point computation we use, these are integers */

    /* initialise the output Vector */
    for i in 0..L_SUBFRAME {
        fixed_codebook_vector[i] = 0;
    }

    /* get the signs and compute the fixedCodebookVector */
    for i in 0..4 {
        if (signs & 1) != 0 {
            /* sign bit is 1 */
            fixed_codebook_vector[positions_array[i] as usize] = 8192; /* +1 in Q13 */
        } else {
            fixed_codebook_vector[positions_array[i] as usize] = -8192; /* -1 in Q13 */
        }
        signs = shr16(signs as i16, 1) as u16; /* shift signs right by one to get the next sign bit at next loop */
    }

    /* if intPitchDelay is smaller than subframe length, give some correction using boundedPitchGain according to eq48 */
    for i in int_pitch_delay as usize..L_SUBFRAME {
        /* fixedCodebookVector[i] = fixedCodebookVector[i] + fixedCodebookVector[i-intPitchDelay]*boundedPitchGain */
        fixed_codebook_vector[i] = add16(
            fixed_codebook_vector[i],
            mult16_16_p14(
                fixed_codebook_vector[i - int_pitch_delay as usize],
                bounded_pitch_gain,
            ) as Word16,
        );
    }
}
