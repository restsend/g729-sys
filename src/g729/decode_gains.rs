use crate::g729::basic_operations::*;
use crate::g729::codebooks::*;
use crate::g729::gain_quantization::{compute_gain_prediction_error, ma_code_gain_prediction};

pub fn init_decode_gains() -> [i16; 4] {
    /*init previousGainPredictionError to -14 in Q10 */
    [-14336; 4]
}

/*****************************************************************************/
/* decodeGains : decode adaptive and fixed codebooks gains as in spec 4.1.5  */

/*           in Q1.13.                                                       */

/*             subframe Pitch Gain in Q14                                    */

/*             Codebook Gain in Q1                                           */

/*****************************************************************************/
pub fn decode_gains(
    previous_gain_prediction_error: &mut [i16; 4],
    mut ga: u16,
    mut gb: u16,
    fixed_codebook_vector: &[i16],
    frame_erasure_flag: u8,
    adaptative_codebook_gain: &mut i16,
    fixed_codebook_gain: &mut i16,
) {
    if frame_erasure_flag != 0 {
        /* we have a frame erasure, proceed as described in spec 4.4.2 */
        let mut current_gain_prediction_error: i32 = 0;

        /*  adaptativeCodebookGain as in eq94 */
        if *adaptative_codebook_gain < 16384 {
            /* last subframe gain < 1 in Q14 */
            *adaptative_codebook_gain = mult16_16_q15(*adaptative_codebook_gain, 29491_i16) as i16;
        /* *0.9 in Q15 */
        } else {
            /* bound current subframe gain to 0.9 (14746 in Q14) */
            *adaptative_codebook_gain = 14746;
        }
        /* fixedCodebookGain as in eq93 */
        *fixed_codebook_gain = mult16_16_q15(*fixed_codebook_gain, 32113_i16) as i16; /* *0.98 in Q15 */

        /* And update the previousGainPredictionError according to spec 4.4.3 */
        for &error in previous_gain_prediction_error.iter().take(4) {
            current_gain_prediction_error = add32(current_gain_prediction_error, error as i32);
        }
        current_gain_prediction_error = pshr(current_gain_prediction_error, 2);

        if current_gain_prediction_error < -10240 {
            /* final result is low bounded by -14, so check before doing -4 if it's over -10(-10240 in Q10) or not */
            current_gain_prediction_error = -14336; /* set to -14 in Q10 */
        } else {
            current_gain_prediction_error = sub32(current_gain_prediction_error, 4096);
            /* in Q10 */
        }

        previous_gain_prediction_error[3] = previous_gain_prediction_error[2];
        previous_gain_prediction_error[2] = previous_gain_prediction_error[1];
        previous_gain_prediction_error[1] = previous_gain_prediction_error[0];
        previous_gain_prediction_error[0] = current_gain_prediction_error as i16;

        return;
    }

    /* First recover the GA and GB real index from their mapping tables(spec 3.9.3) */
    ga = REVERSE_INDEX_MAPPING_GA[ga as usize];
    gb = REVERSE_INDEX_MAPPING_GB[gb as usize];

    /* Compute the adaptativeCodebookGain from the tables according to eq73 in spec3.9.2 */

    *adaptative_codebook_gain = add16(GA_CODEBOOK[ga as usize][0], GB_CODEBOOK[gb as usize][0]); /* result in Q1.14 */

    let predicted_fixed_codebook_gain: i32 =
        ma_code_gain_prediction(previous_gain_prediction_error, fixed_codebook_vector); /* predictedFixedCodebookGain on 32 bits in Q11.16 */

    /* get fixed codebook gain correction factor(gama) from the codebooks GA and GB according to eq74 */
    let fixed_codebook_gain_correction_factor: i16 =
        add16(GA_CODEBOOK[ga as usize][1], GB_CODEBOOK[gb as usize][1]);

    /* compute fixedCodebookGain according to eq74 */
    *fixed_codebook_gain = pshr(
        mult16_32_q12(
            fixed_codebook_gain_correction_factor,
            predicted_fixed_codebook_gain,
        ),
        15,
    ) as i16;

    /* use eq72 to compute current prediction error in order to update the previousGainPredictionError array */
    compute_gain_prediction_error(
        fixed_codebook_gain_correction_factor,
        previous_gain_prediction_error,
    );
}
