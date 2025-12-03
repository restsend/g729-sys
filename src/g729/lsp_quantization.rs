use crate::g729::basic_operations::*;
use crate::g729::codebooks::*;
use crate::g729::fixed_point_math::*;
use crate::g729::ld8k::*;
use crate::g729::utils::*;

/* static buffers */
const PREVIOUS_QLSF_INIT: [Word16; NB_LSP_COEFF] = [
    2339, 4679, 7018, 9358, 11698, 14037, 16377, 18717, 21056, 23396,
]; /* PI*(float)(j+1)/(float)(M+1) */

/* initialise the stactic buffers */
pub fn init_lsp_quantization(previous_q_lsf: &mut [[Word16; NB_LSP_COEFF]; MA_MAX_K]) {
    for i in 0..MA_MAX_K {
        previous_q_lsf[i].copy_from_slice(&PREVIOUS_QLSF_INIT);
    }
}

/**********************************************************************************/
/* noiseLSPQuantization : Convert LSP to LSF, Quantize LSF and find L parameters, */
/*      qLSF->qLSP as described in spec A3.2.4                                    */
/*    parameters:                                                                 */
/*      -(i/o) previousqLSF : 4 previousqLSF, is updated by this function         */
/*      -(i) LSPCoefficients : 10 LSP coefficients in Q15                         */
/*      -(o) qLSPCoefficients : 10 qLSP coefficients in Q15                       */
/*      -(o) parameters : 3 parameters L0, L1, L2                                 */
/*                                                                                */
/**********************************************************************************/
pub fn noise_lsp_quantization(
    previous_q_lsf: &mut [[Word16; NB_LSP_COEFF]; MA_MAX_K],
    lsp_coefficients: &[Word16],
    q_lsp_coefficients: &mut [Word16],
    parameters: &mut [u8], // Using u8 as in C uint8_t
) {
    let mut lsf = [0; NB_LSP_COEFF]; /* LSF coefficients in Q2.13 range [0, Pi[ */
    let mut weights = [0 as UWord16; NB_LSP_COEFF]; /* weights in Q11 */
    let mut weights_threshold = [0; NB_LSP_COEFF]; /* store in Q13 the threshold used to compute the weights */
    let mut l1_index = [0; L0_RANGE];
    let mut l2_index = [0; L0_RANGE];
    let mut weighted_mean_square_error = [0 as UWord32; L0_RANGE];
    let mut quantizer_output = [0; NB_LSP_COEFF];
    let mut q_lsf = [0; NB_LSP_COEFF];

    /*** compute LSF in Q2.13 : lsf = arcos(lsp) range [0, Pi[ spec 3.2.4 eq18 ***/
    for i in 0..NB_LSP_COEFF {
        lsf[i] = g729_acos_q15q13(lsp_coefficients[i]);
    }

    /*** compute the weights vector as in spec 3.2.4 eq22 ***/
    weights_threshold[0] = sub16(lsf[1], OO4PIPLUS1_IN_Q13);
    for i in 1..NB_LSP_COEFF - 1 {
        weights_threshold[i] = sub16(sub16(lsf[i + 1], lsf[i - 1]), ONE_IN_Q13 as Word16);
    }
    weights_threshold[NB_LSP_COEFF - 1] = sub16(O92PIMINUS1_IN_Q13, lsf[NB_LSP_COEFF - 2]);

    for i in 0..NB_LSP_COEFF {
        if weights_threshold[i] > 0 {
            weights[i] = ONE_IN_Q11 as UWord16;
        } else {
            weights[i] = usaturate(
                add32(
                    pshr(
                        mult16_16(
                            mult16_16_q13(weights_threshold[i], weights_threshold[i]) as Word16,
                            10,
                        ),
                        2,
                    ),
                    ONE_IN_Q11 as Word32,
                ),
                MAX_16 as Word32,
            ) as UWord16;
        }
    }
    weights[4] = mult16_16_q14(weights[4] as Word16, ONE_POINT_2_IN_Q14) as UWord16;
    weights[5] = mult16_16_q14(weights[5] as Word16, ONE_POINT_2_IN_Q14) as UWord16;

    /*** compute the coefficients for the two noise noise MA Predictors ***/
    for l0 in 0..L0_RANGE {
        /* compute the target Vector (l) to be quantized as in spec 3.2.4 eq23 */
        let mut target_vector = [0; NB_LSP_COEFF]; /* vector to be quantized in Q13 */
        let mut mean_square_diff: Word32 = MAX_32;
        let mut quantized_vector = [0; NB_LSP_COEFF]; /* in Q13, the current state of quantized vector */

        for i in 0..NB_LSP_COEFF {
            let mut acc = shl(lsf[i] as Word32, 15); /* acc in Q2.28 */
            for j in 0..MA_MAX_K {
                acc = msu16_16(acc, previous_q_lsf[j][i], NOISE_MA_PREDICTOR[l0][j][i]);
                /* previousqLSF in Q2.13 and MAPredictor in Q0.15-> acc in Q2.28 */
            }
            target_vector[i] =
                mult16_16_q12(pshr(acc, 15) as Word16, INV_NOISE_MA_PREDICTOR_SUM[l0][i]);
            /* acc->Q13 and invMAPredictorSum in Q12 -> targetVector in Q13 */
        }

        /* find closest match for predictionError (minimize mean square diff) in L1 subset codebook: 32 entries from L1 codebook */
        for i in 0..NOISE_L1_RANGE {
            let mut acc: Word32 = 0;
            for j in 0..NB_LSP_COEFF {
                let diff_target_vector_l1 = saturate(
                    sub32(
                        target_vector[j] as Word32,
                        L1[L1_SUBSET_INDEX[i]][j] as Word32,
                    ),
                    MAX_16 as Word32,
                ) as Word16;
                acc = mac16_16(acc, diff_target_vector_l1, diff_target_vector_l1);
            }

            if acc < mean_square_diff {
                mean_square_diff = acc;
                l1_index[l0] = i;
            }
        }

        /* find the closest match in L2 subset wich will minimise the weighted sum of (targetVector - L1 result - L2)^2 */
        /* using eq20, eq21 and eq23 in spec 3.2.4 -> l[i] - l^[i] = (wi - w^[i])/(1-SumMAPred[i]) but ITU code ignores this denominator */
        /* works on the first five coefficients only */
        mean_square_diff = MAX_32;
        for i in 0..NOISE_L2_RANGE {
            let mut acc: Word32 = 0;
            for j in 0..NB_LSP_COEFF / 2 {
                /* commented code : compute in the same way of the ITU code: ignore the denonimator and minimize (wi - w^[i])/(1-SumMAPred[i]) instead of (wi - w^[i]) square sum */
                let diff_target_vector_l1_l2 = saturate(
                    mult16_16_q15(
                        sub32(
                            sub32(
                                target_vector[j] as Word32,
                                L1[L1_SUBSET_INDEX[l1_index[l0]]][j] as Word32,
                            ),
                            L2L3[L2_SUBSET_INDEX[i]][j] as Word32,
                        ) as Word16,
                        NOISE_MA_PREDICTOR_SUM[l0][j],
                    ),
                    MAX_16 as Word32,
                ) as Word16; /* targetVector, L1 and L2L3 in Q13 -> result in Q13 */
                acc = mac16_16(
                    acc,
                    diff_target_vector_l1_l2,
                    mult16_16_q11(diff_target_vector_l1_l2, weights[j] as Word16) as Word16,
                ); /* weights in Q11, diff in Q13 */
            }

            for j in NB_LSP_COEFF / 2..NB_LSP_COEFF {
                let diff_target_vector_l1_l3 = saturate(
                    mult16_16_q15(
                        sub32(
                            sub32(
                                target_vector[j] as Word32,
                                L1[L1_SUBSET_INDEX[l1_index[l0]]][j] as Word32,
                            ),
                            L2L3[L3_SUBSET_INDEX[i]][j] as Word32,
                        ) as Word16,
                        NOISE_MA_PREDICTOR_SUM[l0][j],
                    ),
                    MAX_16 as Word32,
                ) as Word16; /* targetVector, L1 and L2L3 in Q13 -> result in Q13 */
                acc = mac16_16(
                    acc,
                    diff_target_vector_l1_l3,
                    mult16_16_q11(diff_target_vector_l1_l3, weights[j] as Word16) as Word16,
                ); /* weights in Q11, diff in Q13 */
            }

            if acc < mean_square_diff {
                mean_square_diff = acc;
                l2_index[l0] = i;
            }
        }

        /* compute the quantized vector L1+L2/L3 and rearrange it as specified in spec 3.2.4 */
        /* Note: according to the spec, the rearrangement shall be done on each candidate while looking for best match, but the ITU code does it after picking the best match and so we do */
        for i in 0..NB_LSP_COEFF / 2 {
            quantized_vector[i] = add16(
                L1[L1_SUBSET_INDEX[l1_index[l0]]][i],
                L2L3[L2_SUBSET_INDEX[l2_index[l0]]][i],
            );
        }
        for i in NB_LSP_COEFF / 2..NB_LSP_COEFF {
            quantized_vector[i] = add16(
                L1[L1_SUBSET_INDEX[l1_index[l0]]][i],
                L2L3[L3_SUBSET_INDEX[l2_index[l0]]][i],
            );
        }

        /* rearrange with a minimum distance of 0.0012 */
        for i in 1..NB_LSP_COEFF / 2 {
            if quantized_vector[i - 1] > sub16(quantized_vector[i], GAP1) {
                quantized_vector[i - 1] = pshr(
                    sub16(add16(quantized_vector[i], quantized_vector[i - 1]), GAP1) as Word32,
                    1,
                ) as Word16;
                quantized_vector[i] = pshr(
                    add16(add16(quantized_vector[i], quantized_vector[i - 1]), GAP1) as Word32,
                    1,
                ) as Word16;
            }
        }
        for i in NB_LSP_COEFF / 2 + 1..NB_LSP_COEFF {
            if quantized_vector[i - 1] > sub16(quantized_vector[i], GAP1) {
                quantized_vector[i - 1] = pshr(
                    sub16(add16(quantized_vector[i], quantized_vector[i - 1]), GAP1) as Word32,
                    1,
                ) as Word16;
                quantized_vector[i] = pshr(
                    add16(add16(quantized_vector[i], quantized_vector[i - 1]), GAP1) as Word32,
                    1,
                ) as Word16;
            }
        }

        /* rearrange the whole quantizedVector with a distance of 0.0006 */
        for i in 1..NB_LSP_COEFF {
            if quantized_vector[i - 1] > sub16(quantized_vector[i], GAP2) {
                quantized_vector[i - 1] = pshr(
                    sub16(add16(quantized_vector[i], quantized_vector[i - 1]), GAP2) as Word32,
                    1,
                ) as Word16;
                quantized_vector[i] = pshr(
                    add16(add16(quantized_vector[i], quantized_vector[i - 1]), GAP2) as Word32,
                    1,
                ) as Word16;
            }
        }

        /* compute the weighted mean square distance using the final quantized vector according to eq21 */
        weighted_mean_square_error[l0] = 0;
        for i in 0..NB_LSP_COEFF {
            let diff_target_vector_quantized_vector = usaturate(
                abs(mult16_32_q15(
                    NOISE_MA_PREDICTOR_SUM[l0][i],
                    sub32(target_vector[i] as Word32, quantized_vector[i] as Word32),
                )) as Word32,
                MAX_U16 as Word32,
            ) as UWord16; /* targetVector and quantizedVector in Q13 -> result in Q13 */
            weighted_mean_square_error[l0] = umac16_16(
                weighted_mean_square_error[l0],
                diff_target_vector_quantized_vector,
                mult16_16_q11(
                    diff_target_vector_quantized_vector as Word16,
                    weights[i] as Word16,
                ) as UWord16,
            ); /* weights in Q11, diff in Q13 */
        }
    }

    /* now select L0 and copy the selected coefficients to the output buffer */
    if weighted_mean_square_error[0] < weighted_mean_square_error[1] {
        parameters[0] = 0;
        parameters[1] = l1_index[0] as u8;
        parameters[2] = l2_index[0] as u8;
    } else {
        parameters[0] = 1;
        parameters[1] = l1_index[1] as u8;
        parameters[2] = l2_index[1] as u8;
    }

    /*** Compute the quantized LSF from the L coefficients ***/
    /* reconstruct vector from the codebooks using the selected parameters spec 3.2.4 eq19 */
    for i in 0..NB_LSP_COEFF / 2 {
        quantizer_output[i] = add16(
            L1[L1_SUBSET_INDEX[parameters[1] as usize]][i],
            L2L3[L2_SUBSET_INDEX[parameters[2] as usize]][i],
        ); /* codebooks are in Q2.13 for L1 and Q0.13 for L2L3, due to actual values stored in the codebooks, result in Q2.13 */
    }
    for i in NB_LSP_COEFF / 2..NB_LSP_COEFF {
        quantizer_output[i] = add16(
            L1[L1_SUBSET_INDEX[parameters[1] as usize]][i],
            L2L3[L3_SUBSET_INDEX[parameters[2] as usize]][i],
        ); /* same as previous, output in Q2.13 */
    }
    /* rearrange in order to have a minimum distance between two consecutives coefficients spec 3.2.4 */
    rearrange_coefficients(&mut quantizer_output, GAP1);
    rearrange_coefficients(&mut quantizer_output, GAP2); /* currentqLSF still in Q2.13 */

    /* compute qLSF spec 3.2.4 eq20 */
    for i in 0..NB_LSP_COEFF {
        let mut acc = mult16_16(
            NOISE_MA_PREDICTOR_SUM[parameters[0] as usize][i],
            quantizer_output[i],
        ); /* (1 - ∑Pi,k)*lˆi(m) Q15 * Q13 -> Q28 */
        for j in 0..MA_MAX_K {
            acc = mac16_16(
                acc,
                NOISE_MA_PREDICTOR[parameters[0] as usize][j][i],
                previous_q_lsf[j][i],
            );
        }
        /* acc in Q2.28, shift back the acc to a Q2.13 with rounding */
        q_lsf[i] = pshr(acc, 15) as Word16; /* qLSF in Q2.13 */
    }

    /* update the previousqLSF buffer with current quantizer output */
    for i in (1..MA_MAX_K).rev() {
        previous_q_lsf[i] = previous_q_lsf[i - 1];
    }
    previous_q_lsf[0].copy_from_slice(&quantizer_output);

    /*** qLSF stability check ***/
    insertion_sort(&mut q_lsf, NB_LSP_COEFF);

    /* check for low limit on qLSF[0] */
    if q_lsf[1] < QLSF_MIN {
        q_lsf[1] = QLSF_MIN;
    }

    /* check and rectify minimum distance between two consecutive qLSF */
    for i in 0..NB_LSP_COEFF - 1 {
        if sub16(q_lsf[i + 1], q_lsf[i]) < MIN_QLSF_DISTANCE {
            q_lsf[i + 1] = q_lsf[i] + MIN_QLSF_DISTANCE;
        }
    }

    /* check for upper limit on qLSF[NB_LSP_COEFF-1] */
    if q_lsf[NB_LSP_COEFF - 1] > QLSF_MAX {
        q_lsf[NB_LSP_COEFF - 1] = QLSF_MAX;
    }

    /* convert qLSF to qLSP: qLSP = cos(qLSF) */
    for i in 0..NB_LSP_COEFF {
        q_lsp_coefficients[i] = g729_cos_q13q15(q_lsf[i]); /* ouput in Q0.15 */
    }
}

/*****************************************************************************/
/* LSPQuantization : Convert LSP to LSF, Quantize LSF and find L parameters, */
/*      qLSF->qLSP as described in spec A3.2.4                               */
/*    parameters:                                                            */
/*      -(i/o) encoderChannelContext : the channel context data              */
/*      -(i) LSPCoefficients : 10 LSP coefficients in Q15                    */
/*      -(o) qLSPCoefficients : 10 qLSP coefficients in Q15                  */
/*      -(o) parameters : 4 parameters L0, L1, L2, L3                        */
/*                                                                           */
/*****************************************************************************/
pub fn lsp_quantization(
    previous_q_lsf: &mut [[Word16; NB_LSP_COEFF]; MA_MAX_K],
    lsp_coefficients: &[Word16],
    q_lsp_coefficients: &mut [Word16],
    parameters: &mut [u16], // Using u16 as in C uint16_t
) {
    let mut lsf = [0; NB_LSP_COEFF]; /* LSF coefficients in Q2.13 range [0, Pi[ */
    let mut weights = [0 as UWord16; NB_LSP_COEFF]; /* weights in Q11 */
    let mut weights_threshold = [0; NB_LSP_COEFF]; /* store in Q13 the threshold used to compute the weights */
    let mut l1_index = [0; L0_RANGE];
    let mut l2_index = [0; L0_RANGE];
    let mut l3_index = [0; L0_RANGE];
    let mut weighted_mean_square_error = [0 as UWord32; L0_RANGE];
    let mut quantizer_output = [0; NB_LSP_COEFF];
    let mut q_lsf = [0; NB_LSP_COEFF];

    /*** compute LSF in Q2.13 : lsf = arcos(lsp) range [0, Pi[ spec 3.2.4 eq18 ***/
    for i in 0..NB_LSP_COEFF {
        lsf[i] = g729_acos_q15q13(lsp_coefficients[i]);
    }

    /*** compute the weights vector as in spec 3.2.4 eq22 ***/
    weights_threshold[0] = sub16(lsf[1], OO4PIPLUS1_IN_Q13);
    for i in 1..NB_LSP_COEFF - 1 {
        weights_threshold[i] = sub16(sub16(lsf[i + 1], lsf[i - 1]), ONE_IN_Q13 as Word16);
    }
    weights_threshold[NB_LSP_COEFF - 1] = sub16(O92PIMINUS1_IN_Q13, lsf[NB_LSP_COEFF - 2]);

    for i in 0..NB_LSP_COEFF {
        if weights_threshold[i] > 0 {
            weights[i] = ONE_IN_Q11 as UWord16;
        } else {
            weights[i] = usaturate(
                add32(
                    pshr(
                        mult16_16(
                            mult16_16_q13(weights_threshold[i], weights_threshold[i]) as Word16,
                            10,
                        ),
                        2,
                    ),
                    ONE_IN_Q11 as Word32,
                ),
                MAX_16 as Word32,
            ) as UWord16;
        }
    }
    weights[4] = mult16_16_q14(weights[4] as Word16, ONE_POINT_2_IN_Q14) as UWord16;
    weights[5] = mult16_16_q14(weights[5] as Word16, ONE_POINT_2_IN_Q14) as UWord16;

    /*** compute the coefficients for the two MA Predictors ***/
    for l0 in 0..L0_RANGE {
        /* compute the target Vector (l) to be quantized as in spec 3.2.4 eq23 */
        let mut target_vector = [0; NB_LSP_COEFF]; /* vector to be quantized in Q13 */
        let mut mean_square_diff: Word32 = MAX_32;
        let mut quantized_vector = [0; NB_LSP_COEFF]; /* in Q13, the current state of quantized vector */

        for i in 0..NB_LSP_COEFF {
            let mut acc = shl(lsf[i] as Word32, 15); /* acc in Q2.28 */
            for j in 0..MA_MAX_K {
                acc = msu16_16(acc, previous_q_lsf[j][i], MA_PREDICTOR[l0][j][i]);
                /* previousqLSF in Q2.13 and MAPredictor in Q0.15-> acc in Q2.28 */
            }
            target_vector[i] = mult16_16_q12(pshr(acc, 15) as Word16, INV_MA_PREDICTOR_SUM[l0][i]);
            /* acc->Q13 and invMAPredictorSum in Q12 -> targetVector in Q13 */
        }

        /* find closest match for predictionError (minimize mean square diff) in L1 subset codebook: 32 entries from L1 codebook */
        for i in 0..L1_RANGE {
            let mut acc: Word32 = 0;
            for j in 0..NB_LSP_COEFF {
                let diff_target_vector_l1 = saturate(
                    sub32(target_vector[j] as Word32, L1[i][j] as Word32),
                    MAX_16 as Word32,
                ) as Word16;
                acc = mac16_16(acc, diff_target_vector_l1, diff_target_vector_l1);
            }

            if acc < mean_square_diff {
                mean_square_diff = acc;
                l1_index[l0] = i;
            }
        }

        /* find the closest match in L2 wich will minimise the weighted sum of (targetVector - L1 result - L2)^2 */
        /* using eq20, eq21 and eq23 in spec 3.2.4 -> l[i] - l^[i] = (wi - w^[i])/(1-SumMAPred[i]) but ITU code ignores this denominator */
        /* works on the first five coefficients only */
        mean_square_diff = MAX_32;
        for i in 0..L2_RANGE {
            let mut acc: Word32 = 0;
            for j in 0..NB_LSP_COEFF / 2 {
                /* commented code : compute in the same way of the ITU code: ignore the denonimator and minimize (wi - w^[i])/(1-SumMAPred[i]) instead of (wi - w^[i]) square sum */
                let diff_target_vector_l1_l2 = saturate(
                    mult16_16_q15(
                        sub32(
                            sub32(target_vector[j] as Word32, L1[l1_index[l0]][j] as Word32),
                            L2L3[i][j] as Word32,
                        ) as Word16,
                        MA_PREDICTOR_SUM[l0][j],
                    ),
                    MAX_16 as Word32,
                ) as Word16; /* targetVector, L1 and L2L3 in Q13 -> result in Q13 */
                acc = mac16_16(
                    acc,
                    diff_target_vector_l1_l2,
                    mult16_16_q11(diff_target_vector_l1_l2, weights[j] as Word16) as Word16,
                ); /* weights in Q11, diff in Q13 */
            }

            if acc < mean_square_diff {
                mean_square_diff = acc;
                l2_index[l0] = i;
            }
        }

        /* find the closest match in L3 wich will minimise the weighted sum of (targetVector - L1 result - L3)^2 */
        /* using eq20, eq21 and eq23 in spec 3.2.4 -> l[i] - l^[i] = (wi - w^[i])/(1-SumMAPred[i]) but ITU code ignores this denominator */
        /* works on the first five coefficients only */
        mean_square_diff = MAX_32;
        for i in 0..L2_RANGE {
            let mut acc: Word32 = 0;
            for j in NB_LSP_COEFF / 2..NB_LSP_COEFF {
                /* commented code : compute in the same way of the ITU code: ignore the denonimator and minimize (wi - w^[i])/(1-SumMAPred[i]) instead of (wi - w^[i]) square sum */
                let diff_target_vector_l1_l3 = saturate(
                    mult16_16_q15(
                        sub32(
                            sub32(target_vector[j] as Word32, L1[l1_index[l0]][j] as Word32),
                            L2L3[i][j] as Word32,
                        ) as Word16,
                        MA_PREDICTOR_SUM[l0][j],
                    ),
                    MAX_16 as Word32,
                ) as Word16; /* targetVector, L1 and L2L3 in Q13 -> result in Q13 */
                acc = mac16_16(
                    acc,
                    diff_target_vector_l1_l3,
                    mult16_16_q11(diff_target_vector_l1_l3, weights[j] as Word16) as Word16,
                ); /* weights in Q11, diff in Q13 */
            }

            if acc < mean_square_diff {
                mean_square_diff = acc;
                l3_index[l0] = i;
            }
        }

        /* compute the quantized vector L1+L2/L3 and rearrange it as specified in spec 3.2.4(first the higher part (L2) and then the lower part (L3)) */
        /* Note: according to the spec, the rearrangement shall be done on each candidate while looking for best match, but the ITU code does it after picking the best match and so we do */
        for i in 0..NB_LSP_COEFF / 2 {
            quantized_vector[i] = add16(L1[l1_index[l0]][i], L2L3[l2_index[l0]][i]);
        }
        for i in NB_LSP_COEFF / 2..NB_LSP_COEFF {
            quantized_vector[i] = add16(L1[l1_index[l0]][i], L2L3[l3_index[l0]][i]);
        }

        /* rearrange with a minimum distance of 0.0012 */
        for i in 1..NB_LSP_COEFF / 2 {
            if quantized_vector[i - 1] > sub16(quantized_vector[i], GAP1) {
                quantized_vector[i - 1] = pshr(
                    sub16(add16(quantized_vector[i], quantized_vector[i - 1]), GAP1) as Word32,
                    1,
                ) as Word16;
                quantized_vector[i] = pshr(
                    add16(add16(quantized_vector[i], quantized_vector[i - 1]), GAP1) as Word32,
                    1,
                ) as Word16;
            }
        }
        for i in NB_LSP_COEFF / 2 + 1..NB_LSP_COEFF {
            if quantized_vector[i - 1] > sub16(quantized_vector[i], GAP1) {
                quantized_vector[i - 1] = pshr(
                    sub16(add16(quantized_vector[i], quantized_vector[i - 1]), GAP1) as Word32,
                    1,
                ) as Word16;
                quantized_vector[i] = pshr(
                    add16(add16(quantized_vector[i], quantized_vector[i - 1]), GAP1) as Word32,
                    1,
                ) as Word16;
            }
        }

        /* rearrange the whole quantizedVector with a distance of 0.0006 */
        for i in 1..NB_LSP_COEFF {
            if quantized_vector[i - 1] > sub16(quantized_vector[i], GAP2) {
                quantized_vector[i - 1] = pshr(
                    sub16(add16(quantized_vector[i], quantized_vector[i - 1]), GAP2) as Word32,
                    1,
                ) as Word16;
                quantized_vector[i] = pshr(
                    add16(add16(quantized_vector[i], quantized_vector[i - 1]), GAP2) as Word32,
                    1,
                ) as Word16;
            }
        }

        /* compute the weighted mean square distance using the final quantized vector according to eq21 */
        weighted_mean_square_error[l0] = 0;
        for i in 0..NB_LSP_COEFF {
            let diff_target_vector_quantized_vector = usaturate(
                abs(mult16_32_q15(
                    MA_PREDICTOR_SUM[l0][i],
                    sub32(target_vector[i] as Word32, quantized_vector[i] as Word32),
                )) as Word32,
                MAX_U16 as Word32,
            ) as UWord16; /* targetVector and quantizedVector in Q13 -> result in Q13 */
            weighted_mean_square_error[l0] = umac16_16(
                weighted_mean_square_error[l0],
                diff_target_vector_quantized_vector,
                mult16_16_q11(
                    diff_target_vector_quantized_vector as Word16,
                    weights[i] as Word16,
                ) as UWord16,
            ); /* weights in Q11, diff in Q13 */
        }
    }

    /* now select L0 and copy the selected coefficients to the output buffer */
    if weighted_mean_square_error[0] < weighted_mean_square_error[1] {
        parameters[0] = 0;
        parameters[1] = l1_index[0] as u16;
        parameters[2] = l2_index[0] as u16;
        parameters[3] = l3_index[0] as u16;
    } else {
        parameters[0] = 1;
        parameters[1] = l1_index[1] as u16;
        parameters[2] = l2_index[1] as u16;
        parameters[3] = l3_index[1] as u16;
    }

    /*** Compute the quantized LSF from the L coefficients ***/
    /* reconstruct vector from the codebooks using the selected parameters spec 3.2.4 eq19 */
    for i in 0..NB_LSP_COEFF / 2 {
        quantizer_output[i] = add16(
            L1[parameters[1] as usize][i],
            L2L3[parameters[2] as usize][i],
        ); /* codebooks are in Q2.13 for L1 and Q0.13 for L2L3, due to actual values stored in the codebooks, result in Q2.13 */
    }
    for i in NB_LSP_COEFF / 2..NB_LSP_COEFF {
        quantizer_output[i] = add16(
            L1[parameters[1] as usize][i],
            L2L3[parameters[3] as usize][i],
        ); /* same as previous, output in Q2.13 */
    }
    /* rearrange in order to have a minimum distance between two consecutives coefficients spec 3.2.4 */
    rearrange_coefficients(&mut quantizer_output, GAP1);
    rearrange_coefficients(&mut quantizer_output, GAP2); /* currentqLSF still in Q2.13 */

    /* compute qLSF spec 3.2.4 eq20 */
    for i in 0..NB_LSP_COEFF {
        let mut acc = mult16_16(
            MA_PREDICTOR_SUM[parameters[0] as usize][i],
            quantizer_output[i],
        ); /* (1 - ∑Pi,k)*lˆi(m) Q15 * Q13 -> Q28 */
        for j in 0..MA_MAX_K {
            acc = mac16_16(
                acc,
                MA_PREDICTOR[parameters[0] as usize][j][i],
                previous_q_lsf[j][i],
            );
        }
        /* acc in Q2.28, shift back the acc to a Q2.13 with rounding */
        q_lsf[i] = pshr(acc, 15) as Word16; /* qLSF in Q2.13 */
    }

    /* update the previousqLSF buffer with current quantizer output */
    for i in (1..MA_MAX_K).rev() {
        previous_q_lsf[i] = previous_q_lsf[i - 1];
    }
    previous_q_lsf[0].copy_from_slice(&quantizer_output);

    /*** qLSF stability check ***/
    insertion_sort(&mut q_lsf, NB_LSP_COEFF);

    /* check for low limit on qLSF[0] */
    if q_lsf[1] < QLSF_MIN {
        q_lsf[1] = QLSF_MIN;
    }

    /* check and rectify minimum distance between two consecutive qLSF */
    for i in 0..NB_LSP_COEFF - 1 {
        if sub16(q_lsf[i + 1], q_lsf[i]) < MIN_QLSF_DISTANCE {
            q_lsf[i + 1] = q_lsf[i] + MIN_QLSF_DISTANCE;
        }
    }

    /* check for upper limit on qLSF[NB_LSP_COEFF-1] */
    if q_lsf[NB_LSP_COEFF - 1] > QLSF_MAX {
        q_lsf[NB_LSP_COEFF - 1] = QLSF_MAX;
    }

    /* convert qLSF to qLSP: qLSP = cos(qLSF) */
    for i in 0..NB_LSP_COEFF {
        q_lsp_coefficients[i] = g729_cos_q13q15(q_lsf[i]); /* ouput in Q0.15 */
    }
}
