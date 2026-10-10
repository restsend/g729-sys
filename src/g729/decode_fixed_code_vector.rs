use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;

/*****************************************************************************/
/* decodeFixedCodeVector : compute the fixed codebook vector as in spec 4.1.4*/

/*      -(i) signs: parameter S(4 signs bit) eq61                            */
/*      -(i) positions: parameter C(4 3bits position and jx bit) eq62        */

/*      -(i) boundedPitchGain: Beta in eq47 and eq48, in Q14                 */

/*****************************************************************************/
pub fn decode_fixed_code_vector(
    mut signs: u16,
    mut positions: u16,
    int_pitch_delay: i16,
    bounded_pitch_gain: i16,
    fixed_codebook_vector: &mut [i16],
) {
    let mut positions_array = [0u16; 4];

    /* get the positions into an array: mapping according to eq62 and table7 in spec 3.8 */
    positions_array[0] = (positions & 7) * 5;
    positions = shr16(positions as i16, 3) as u16;
    positions_array[1] = ((positions & 7) * 5) + 1;
    positions = shr16(positions as i16, 3) as u16;
    positions_array[2] = ((positions & 7) * 5) + 2;
    positions = shr16(positions as i16, 3) as u16;
    let jx: u16 = positions & 1; /* jx from eq62 is the last bit */
    positions = shr16(positions as i16, 1) as u16;
    positions_array[3] = ((positions & 7) * 5) + 3 + jx;

    fixed_codebook_vector[..L_SUBFRAME].fill(0);

    for i in 0..4 {
        if (signs & 1) != 0 {
            fixed_codebook_vector[positions_array[i] as usize] = 8192; /* +1 in Q13 */
        } else {
            fixed_codebook_vector[positions_array[i] as usize] = -8192; /* -1 in Q13 */
        }
        signs = shr16(signs as i16, 1) as u16;
    }

    /* if intPitchDelay is smaller than subframe length, give some correction using boundedPitchGain according to eq48 */
    for i in int_pitch_delay as usize..L_SUBFRAME {
        fixed_codebook_vector[i] = add16(
            fixed_codebook_vector[i],
            mult16_16_p14(
                fixed_codebook_vector[i - int_pitch_delay as usize],
                bounded_pitch_gain,
            ) as i16,
        );
    }
}
