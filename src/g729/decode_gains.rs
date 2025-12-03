use crate::g729::basic_operations::*;
use crate::g729::codebooks::*;
use crate::g729::gain_quantization::{compute_gain_prediction_error, ma_code_gain_prediction};

/* init function */
pub fn init_decode_gains(previous_gain_prediction_error: &mut [Word16; 4]) {
    /*init previousGainPredictionError to -14 in Q10 */
    previous_gain_prediction_error[0] = -14336;
    previous_gain_prediction_error[1] = -14336;
    previous_gain_prediction_error[2] = -14336;
    previous_gain_prediction_error[3] = -14336;
}

/*****************************************************************************/
/* decodeGains : decode adaptive and fixed codebooks gains as in spec 4.1.5  */
/*    parameters:                                                            */
/*      -(i/o) decoderChannelContext : the channel context data              */
/*      -(i) GA : parameter GA: Gain Codebook Stage 1 (3 bits)               */
/*      -(i) GB : paremeter GB: Gain Codebook Stage 2 (4 bits)               */
/*      -(i) fixedCodebookVector: 40 values current fixed Codebook vector    */
/*           in Q1.13.                                                       */
/*      -(i) frameErasureFlag : set in case of frame erasure                 */
/*      -(i/o) adaptativeCodebookGain : input previous/output current        */
/*             subframe Pitch Gain in Q14                                    */
/*      -(i/o) fixedCodebookGain : input previous/output current Fixed       */
/*             Codebook Gain in Q1                                           */
/*                                                                           */
/*****************************************************************************/
pub fn decode_gains(
    previous_gain_prediction_error: &mut [Word16; 4],
    mut ga: u16,
    mut gb: u16,
    fixed_codebook_vector: &[Word16],
    frame_erasure_flag: u8,
    adaptative_codebook_gain: &mut Word16,
    fixed_codebook_gain: &mut Word16,
) {
    let predicted_fixed_codebook_gain: Word32;
    let fixed_codebook_gain_correction_factor: Word16;

    /* check the erasure flag */
    if frame_erasure_flag != 0 {
        /* we have a frame erasure, proceed as described in spec 4.4.2 */
        let mut current_gain_prediction_error: Word32 = 0;

        /*  adaptativeCodebookGain as in eq94 */
        if *adaptative_codebook_gain < 16384 {
            /* last subframe gain < 1 in Q14 */
            *adaptative_codebook_gain =
                mult16_16_q15(*adaptative_codebook_gain, 29491 as Word16) as Word16;
        /* *0.9 in Q15 */
        } else {
            /* bound current subframe gain to 0.9 (14746 in Q14) */
            *adaptative_codebook_gain = 14746;
        }
        /* fixedCodebookGain as in eq93 */
        *fixed_codebook_gain = mult16_16_q15(*fixed_codebook_gain, 32113 as Word16) as Word16; /* *0.98 in Q15 */

        /* And update the previousGainPredictionError according to spec 4.4.3 */
        for i in 0..4 {
            current_gain_prediction_error = add32(
                current_gain_prediction_error,
                previous_gain_prediction_error[i] as Word32,
            ); /* previousGainPredictionError in Q3.10-> Sum in Q5.10 (on 32 bits) */
        }
        current_gain_prediction_error = pshr(current_gain_prediction_error, 2); /* /4 -> Q3.10 */

        if current_gain_prediction_error < -10240 {
            /* final result is low bounded by -14, so check before doing -4 if it's over -10(-10240 in Q10) or not */
            current_gain_prediction_error = -14336; /* set to -14 in Q10 */
        } else {
            /* substract 4 */
            current_gain_prediction_error = sub32(current_gain_prediction_error, 4096);
            /* in Q10 */
        }

        /* shift the array and insert the current Prediction Error */
        previous_gain_prediction_error[3] = previous_gain_prediction_error[2];
        previous_gain_prediction_error[2] = previous_gain_prediction_error[1];
        previous_gain_prediction_error[1] = previous_gain_prediction_error[0];
        previous_gain_prediction_error[0] = current_gain_prediction_error as Word16;

        return;
    }

    /* First recover the GA and GB real index from their mapping tables(spec 3.9.3) */
    ga = REVERSE_INDEX_MAPPING_GA[ga as usize];
    gb = REVERSE_INDEX_MAPPING_GB[gb as usize];

    /* Compute the adaptativeCodebookGain from the tables according to eq73 in spec3.9.2 */
    /* adaptativeCodebookGain = GACodebook[GA][0] + GBCodebook[GB][0] */
    *adaptative_codebook_gain = add16(GA_CODEBOOK[ga as usize][0], GB_CODEBOOK[gb as usize][0]); /* result in Q1.14 */

    /* Fixed Codebook: MA code-gain prediction */
    predicted_fixed_codebook_gain =
        ma_code_gain_prediction(previous_gain_prediction_error, fixed_codebook_vector); /* predictedFixedCodebookGain on 32 bits in Q11.16 */

    /* get fixed codebook gain correction factor(gama) from the codebooks GA and GB according to eq74 */
    fixed_codebook_gain_correction_factor =
        add16(GA_CODEBOOK[ga as usize][1], GB_CODEBOOK[gb as usize][1]); /* result in Q3.12 (range [0.185, 5.05])*/

    /* compute fixedCodebookGain according to eq74 */
    *fixed_codebook_gain = pshr(
        mult16_32_q12(
            fixed_codebook_gain_correction_factor,
            predicted_fixed_codebook_gain,
        ),
        15,
    ) as Word16; /* Q11.16*Q3.12 -> Q14.16, shift by 15 to get a Q14.1 which fits on 16 bits */

    /* use eq72 to compute current prediction error in order to update the previousGainPredictionError array */
    compute_gain_prediction_error(
        fixed_codebook_gain_correction_factor,
        previous_gain_prediction_error,
    );
}
