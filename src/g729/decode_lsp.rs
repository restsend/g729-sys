use crate::g729::basic_operations::*;
use crate::g729::codebooks::*;
use crate::g729::fixed_point_math::g729_cos_q13q15;
use crate::g729::ld8k::*;
use crate::g729::utils::{insertion_sort, rearrange_coefficients};

/* previous L Code Word initial values (Pi/11 steps) in Q2.13 */
static PREVIOUS_L_CODE_WORD_INIT: [Word16; NB_LSP_COEFF] = [
    2339, 4679, 7018, 9358, 11698, 14037, 16377, 18717, 21056, 23396,
];

pub fn init_decode_lsp(
    previous_l_code_word: &mut [[Word16; NB_LSP_COEFF]; MA_MAX_K],
    last_valid_l0: &mut u16,
    last_q_lsf: &mut [Word16; NB_LSP_COEFF],
) {
    /* init the previousLCodeWord buffer according to doc 3.2.4 -> pi/11 steps */
    for i in 0..MA_MAX_K {
        for j in 0..NB_LSP_COEFF {
            previous_l_code_word[i][j] = PREVIOUS_L_CODE_WORD_INIT[j];
        }
    }

    /* init the last valid values just to avoid problem in case the first frame is a lost one */
    *last_valid_l0 = 0;
    for j in 0..NB_LSP_COEFF {
        last_q_lsf[j] = PREVIOUS_L_CODE_WORD_INIT[j];
    }
}

/*****************************************************************************/
/* computeqLSF : get qLSF extracted from codebooks and process them          */
/*         according to spec 3.2.4                                           */
/*    parameters:                                                            */
/*      -(i/o) codebookqLSF : 10 values i Q2.13 to be updated                */
/*      -(i/o) previousCodeWord : codewords for the last 4 subframes in Q2.13*/
/*                                is updated by this function                */
/*      -(i) L0: the Switched MA predictor retrieved from bitstream          */
/*                                                                           */
/*****************************************************************************/
fn compute_q_lsf(
    codebook_q_lsf: &mut [Word16; NB_LSP_COEFF],
    previous_l_code_word: &mut [[Word16; NB_LSP_COEFF]; MA_MAX_K],
    l0: usize,
    current_ma_predictor: &[[[Word16; NB_LSP_COEFF]; MA_MAX_K]; L0_RANGE],
    current_ma_predictor_sum: &[[Word16; NB_LSP_COEFF]; L0_RANGE],
) {
    let mut acc: Word32; /* Accumulator in Q2.28 */

    /*** rearrange in order to have a minimum distance between two consecutives coefficients ***/
    rearrange_coefficients(codebook_q_lsf, GAP1);
    rearrange_coefficients(codebook_q_lsf, GAP2); /* codebookqLSF still in Q2.13 */

    /*** doc 3.2.4 eq(20) ***/
    /* compute the qLSF as a weighted sum(weighted by MA coefficient selected according to L0 value) of previous and current frame L codewords coefficients */

    /* L0 is the Switched MA predictor of LSP quantizer(1 bit) */
    /* codebookqLSF and previousLCodeWord in Q2.13 */
    /* MAPredictor and MAPredictorSum in Q0.15 with MAPredictorSum[MA switch][i]+Sum[j=0-3](MAPredictor[MA switch][j][i])=1 -> acc will end up being in Q2.28*/
    /* Note : previousLCodeWord array containing the last 4 code words is updated during this phase */

    for i in 0..NB_LSP_COEFF {
        acc = mult16_16(current_ma_predictor_sum[l0][i], codebook_q_lsf[i]);
        for j in (0..MA_MAX_K).rev() {
            acc = mac16_16(
                acc,
                current_ma_predictor[l0][j][i],
                previous_l_code_word[j][i],
            );
            if j > 0 {
                previous_l_code_word[j][i] = previous_l_code_word[j - 1][i];
            } else {
                previous_l_code_word[j][i] = codebook_q_lsf[i];
            }
        }
        /* acc in Q2.28, shift back the acc to a Q2.13 with rounding */
        codebook_q_lsf[i] = pshr(acc, 15) as Word16; /* codebookqLSF in Q2.13 */
    }
    /* Note : codebookqLSF buffer now contains qLSF */

    /*** doc 3.2.4 qLSF stability ***/
    /* qLSF in Q2.13 as are qLSF_MIN and qLSF_MAX and MIN_qLSF_DISTANCE */

    /* sort the codebookqLSF array */
    insertion_sort(codebook_q_lsf, NB_LSP_COEFF);

    /* check for low limit on qLSF[0] */
    if codebook_q_lsf[0] < QLSF_MIN {
        codebook_q_lsf[0] = QLSF_MIN;
    }

    /* check and rectify minimum distance between two consecutive qLSF */
    for i in 0..NB_LSP_COEFF - 1 {
        if sub16(codebook_q_lsf[i + 1], codebook_q_lsf[i]) < MIN_QLSF_DISTANCE {
            codebook_q_lsf[i + 1] = add16(codebook_q_lsf[i], MIN_QLSF_DISTANCE);
        }
    }

    /* check for upper limit on qLSF[NB_LSP_COEFF-1] */
    if codebook_q_lsf[NB_LSP_COEFF - 1] > QLSF_MAX {
        codebook_q_lsf[NB_LSP_COEFF - 1] = QLSF_MAX;
    }
}

/*****************************************************************************/
/* decodeLSP : decode LSP coefficients as in spec 4.1.1/3.2.4                */
/*    parameters:                                                            */
/*      -(i/o) decoderChannelContext : the channel context data              */
/*      -(i) L: 4 elements array containing L[0-3] the first and             */
/*                     second stage vector of LSP quantizer                  */
/*      -(i) frameErased : a boolean, when true, frame has been erased       */
/*      -(o) qLSP: 10 quantized LSP coefficients in Q15 in range [-1,+1[     */
/*                                                                           */
/*****************************************************************************/
pub fn decode_lsp(
    previous_l_code_word: &mut [[Word16; NB_LSP_COEFF]; MA_MAX_K],
    last_valid_l0: &mut u16,
    last_q_lsf: &mut [Word16; NB_LSP_COEFF],
    l: &[u16],
    q_lsp: &mut [Word16; NB_LSP_COEFF],
    frame_erased: u8,
) {
    let mut current_q_lsf = [0 as Word16; NB_LSP_COEFF]; /* buffer to the current qLSF in Q2.13 */

    if frame_erased == 0 {
        /* frame is ok, proceed according to 3.2.4 section of the doc */

        /*** doc 3.2.4 eq(19) ***/
        /* get the L codewords from the codebooks L1, L2 and L3 */
        /* for easier implementation, L2 and L3 5 dimensional codebooks have been stored in one 10 dimensional L2L3 codebook */
        /* get the 5 first coefficient from the L1 and L2 codebooks */
        /* Note : currentqLSF buffer contains L codewords and not qLSF */
        for i in 0..NB_LSP_COEFF / 2 {
            current_q_lsf[i] = add16(L1[l[1] as usize][i], L2L3[l[2] as usize][i]);
            /* codebooks are in Q2.13 for L1 and Q0.13 for L2L3, due to actual values stored in the codebooks, result in Q2.13 */
        }

        /* get the 5 last coefficient from L1 and L3 codebooks */
        for i in NB_LSP_COEFF / 2..NB_LSP_COEFF {
            current_q_lsf[i] = add16(L1[l[1] as usize][i], L2L3[l[3] as usize][i]);
            /* same as previous, output in Q2.13 */
        }

        compute_q_lsf(
            &mut current_q_lsf,
            previous_l_code_word,
            l[0] as usize,
            &MA_PREDICTOR,
            &MA_PREDICTOR_SUM,
        ); /* use regular MAPredictor as this function is not called on SID frame decoding */

        /* backup the qLSF and L0 to restore them in case of frame erased */
        for i in 0..NB_LSP_COEFF {
            last_q_lsf[i] = current_q_lsf[i];
        }
        *last_valid_l0 = l[0];
    } else {
        /* frame erased indicator is set, proceed according to section 4.4 of the specs */
        let mut acc: Word32; /* acc in Q2.28 */

        /* restore the qLSF of last valid frame */
        for i in 0..NB_LSP_COEFF {
            current_q_lsf[i] = last_q_lsf[i];
        }

        /* compute back the codewords from the qLSF and store them in the previousLCodeWord buffer */
        for i in 0..NB_LSP_COEFF {
            /* currentqLSF and previousLCodeWord in Q2.13, MAPredictor in Q0.15 and invMAPredictorSum in Q3.12 */
            acc = shl(last_q_lsf[i] as Word32, 15); /* Q2.13 -> Q2.28 */
            for j in 0..MA_MAX_K {
                acc = msu16_16(
                    acc,
                    MA_PREDICTOR[*last_valid_l0 as usize][j][i],
                    previous_l_code_word[j][i],
                ); /* acc in Q2.28 - MAPredictor in Q0.15 * previousLCodeWord in Q2.13 -> acc in Q2.28 (because 1-Sum(MAPred) < 0.6) */
            }
            acc = mult16_32_q12(INV_MA_PREDICTOR_SUM[*last_valid_l0 as usize][i], acc); /* Q3.12*Q2.28 >>12 -> Q2.28 because invMAPredictor is 1/(1 - Sum(MAPred))*/

            /* update the array of previoux Code Words */
            for j in (0..MA_MAX_K).rev() {
                /* acc in Q2.28, shift back the acc to a Q2.13 with rounding */
                if j > 0 {
                    previous_l_code_word[j][i] = previous_l_code_word[j - 1][i];
                } else {
                    previous_l_code_word[j][i] = pshr(acc, 15) as Word16;
                }
            }
        }
    }

    /* convert qLSF to qLSP: qLSP = cos(qLSF) */
    for i in 0..NB_LSP_COEFF {
        q_lsp[i] = g729_cos_q13q15(current_q_lsf[i]); /* ouput in Q0.15 */
    }
}
