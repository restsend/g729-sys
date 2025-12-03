use crate::g729::basic_operations::*;
use crate::g729::fixed_point_math::g729_sqrt_q0q7;
use crate::g729::ld8k::*;
use crate::g729::lp_synthesis_filter::lp_synthesis_filter;
use crate::g729::utils::count_leading_zeros;

/* init function */
pub fn init_post_filter(
    residual_signal_buffer: &mut [Word16],
    scaled_residual_signal_buffer: &mut [Word16],
    long_term_filtered_residual_signal_buffer: &mut [Word16],
    short_term_filtered_residual_signal_buffer: &mut [Word16],
    previous_adaptative_gain: &mut Word16,
) {
    /* set to zero the residual signal memory */
    for x in residual_signal_buffer.iter_mut() {
        *x = 0;
    }
    for x in scaled_residual_signal_buffer.iter_mut() {
        *x = 0;
    }

    /* set to zero the one word of longTermFilteredResidualSignal needed as memory for tilt compensation filter */
    long_term_filtered_residual_signal_buffer[0] = 0;

    /* intialise the shortTermFilteredResidualSignal filter memory and pointer*/
    for x in short_term_filtered_residual_signal_buffer.iter_mut() {
        *x = 0;
    }

    /* initialise the previous Gain for adaptative gain control */
    *previous_adaptative_gain = 4096; /* 1 in Q12 */
}

/*****************************************************************************/
/* postFilter: filter the reconstructed speech according to spec A.4.2       */
/*    parameters:                                                            */
/*      -(i/o) decoderChannelContext : the channel context data              */
/*      -(i) LPCoefficients: 10 LP coeff for current subframe in Q12         */
/*      -(i) reconstructedSpeech: output of LP Synthesis, 50 values in Q0    */
/*           10 values of previous subframe, accessed in range [-10, 39]     */
/*      -(i) intPitchDelay: the integer part of Pitch Delay in Q0            */
/*      -(i) subframeIndex: 0 or L_SUBFRAME for subframe 0 or 1              */
/*      -(o) postFilteredSignal: 40 values in Q0                             */
/*                                                                           */
/*****************************************************************************/
pub fn post_filter(
    residual_signal_buffer: &mut [Word16],
    scaled_residual_signal_buffer: &mut [Word16],
    long_term_filtered_residual_signal_buffer: &mut [Word16],
    short_term_filtered_residual_signal_buffer: &mut [Word16],
    previous_adaptative_gain: &mut Word16,
    lp_coefficients: &[Word16],
    reconstructed_speech: &[Word16], /* Expecting slice starting at history, length NB_LSP_COEFF + L_SUBFRAME */
    mut int_pitch_delay: i16,
    subframe_index: usize,
    post_filtered_signal: &mut [Word16],
) {
    let mut lp_gamma_n_coefficients = [0 as Word16; NB_LSP_COEFF]; /* in Q12 */
    let mut correlation_max: Word32 = MININT32;
    let mut best_int_pitch_delay: i16 = 0;
    let mut residual_signal_energy: Word32 = 0; /* in Q-4 */
    let mut delayed_residual_signal_energy: Word32 = 0; /* in Q-4 */
    let mut maximum_three: Word32;
    let mut correlation_max_word16: Word16 = 0;
    let mut residual_signal_energy_word16: Word16 = 0;
    let mut delayed_residual_signal_energy_word16: Word16 = 0;
    let mut lp_gamma_d_coefficients = [0 as Word16; NB_LSP_COEFF]; /* in Q12 */
    let mut hf = [0 as Word16; 22]; /* the truncated impulse response to short term filter Hf in Q12 */
    let mut rh1: Word32;
    let mut tilt_compensated_signal = [0 as Word16; L_SUBFRAME]; /* in Q0 */
    let mut gain_scaling_factor: Word16; /* in Q12 */
    let mut short_term_filtered_residual_signal_square_sum: u32 = 0;

    /********************************************************************/
    /* Long Term Post Filter                                            */
    /********************************************************************/
    /*** Compute LPGammaN and LPGammaD coefficients : LPGamma[0] = LP[0]*Gamma^(i+1) (i=0..9) ***/
    /* GAMMA_XX constants are in Q15 */
    lp_gamma_n_coefficients[0] = mult16_16_p15(lp_coefficients[0], GAMMA_N1) as Word16;
    lp_gamma_n_coefficients[1] = mult16_16_p15(lp_coefficients[1], GAMMA_N2) as Word16;
    lp_gamma_n_coefficients[2] = mult16_16_p15(lp_coefficients[2], GAMMA_N3) as Word16;
    lp_gamma_n_coefficients[3] = mult16_16_p15(lp_coefficients[3], GAMMA_N4) as Word16;
    lp_gamma_n_coefficients[4] = mult16_16_p15(lp_coefficients[4], GAMMA_N5) as Word16;
    lp_gamma_n_coefficients[5] = mult16_16_p15(lp_coefficients[5], GAMMA_N6) as Word16;
    lp_gamma_n_coefficients[6] = mult16_16_p15(lp_coefficients[6], GAMMA_N7) as Word16;
    lp_gamma_n_coefficients[7] = mult16_16_p15(lp_coefficients[7], GAMMA_N8) as Word16;
    lp_gamma_n_coefficients[8] = mult16_16_p15(lp_coefficients[8], GAMMA_N9) as Word16;
    lp_gamma_n_coefficients[9] = mult16_16_p15(lp_coefficients[9], GAMMA_N10) as Word16;

    /*** Compute the residual signal as described in spec 4.2.1 eq79 ***/
    /* Compute also a scaled residual signal: shift right by 2 to avoid overflows on 32 bits when computing correlation and energy */

    /* pointers to current subframe beginning */
    // residualSignal = &(decoderChannelContext->residualSignalBuffer[MAXIMUM_INT_PITCH_DELAY+subframeIndex]);
    // scaledResidualSignal = &(decoderChannelContext->scaledResidualSignalBuffer[MAXIMUM_INT_PITCH_DELAY+subframeIndex]);
    let residual_signal_offset = MAXIMUM_INT_PITCH_DELAY + subframe_index;

    for i in 0..L_SUBFRAME {
        let mut acc = sshl(reconstructed_speech[NB_LSP_COEFF + i] as Word32, 12); /* reconstructedSpeech in Q0 shifted to set acc in Q12 */
        for j in 0..NB_LSP_COEFF {
            acc = mac16_16(
                acc,
                lp_gamma_n_coefficients[j],
                reconstructed_speech[NB_LSP_COEFF + i - j - 1],
            ); /* LPGammaNCoefficients in Q12, reconstructedSpeech in Q0 -> acc in Q12 */
        }
        residual_signal_buffer[residual_signal_offset + i] =
            saturate(pshr(acc, 12), MAX_INT16 as Word32) as Word16; /* shift back acc to Q0 and saturate it to avoid overflow when going back to 16 bits */
        scaled_residual_signal_buffer[residual_signal_offset + i] = pshr(
            residual_signal_buffer[residual_signal_offset + i] as Word32,
            2,
        ) as Word16; /* shift acc to Q-2 and saturate it to get the scaled version of the signal */
    }

    /*** Compute the maximum correlation on scaledResidualSignal delayed by intPitchDelay +/- 3 to get the best delay. Spec 4.2.1 eq80 ***/
    /* using a scaled(Q-2) signals gives correlation in Q-4. */
    if int_pitch_delay > (MAXIMUM_INT_PITCH_DELAY - 3) as i16 {
        /* intPitchDelay shall be < MAXIMUM_INT_PITCH_DELAY-3 (140) */
        int_pitch_delay = (MAXIMUM_INT_PITCH_DELAY - 3) as i16;
    }

    for i in (int_pitch_delay - 3)..=(int_pitch_delay + 3) {
        let mut correlation: Word32 = 0;
        // delayedResidualSignal = &(scaledResidualSignal[-i]); /* delayedResidualSignal points to scaledResidualSignal[-i] */
        let delayed_residual_signal_offset = residual_signal_offset as isize - i as isize;

        /* compute correlation: ∑r(n)*rk(n) */
        for j in 0..L_SUBFRAME {
            correlation = mac16_16(
                correlation,
                scaled_residual_signal_buffer
                    [(delayed_residual_signal_offset + j as isize) as usize],
                scaled_residual_signal_buffer[residual_signal_offset + j],
            );
        }
        /* if we have a maximum correlation */
        if correlation > correlation_max {
            correlation_max = correlation;
            best_int_pitch_delay = i; /* get the intPitchDelay */
        }
    }

    /* saturate correlation to a positive integer */
    if correlation_max < 0 {
        correlation_max = 0;
    }

    /*** Compute the signal energy ∑r(n)*r(n) and delayed signal energy ∑rk(n)*rk(n) which shall be used to compute gl spec 4.2.1 eq81, eq 82 and eq83 ***/
    // delayedResidualSignal = &(scaledResidualSignal[-bestIntPitchDelay]); /* in Q-2, points to the residual signal delayed to give the higher correlation: rk(n) */
    let delayed_residual_signal_offset =
        residual_signal_offset as isize - best_int_pitch_delay as isize;

    for i in 0..L_SUBFRAME {
        residual_signal_energy = mac16_16(
            residual_signal_energy,
            scaled_residual_signal_buffer[residual_signal_offset + i],
            scaled_residual_signal_buffer[residual_signal_offset + i],
        );
        delayed_residual_signal_energy = mac16_16(
            delayed_residual_signal_energy,
            scaled_residual_signal_buffer[(delayed_residual_signal_offset + i as isize) as usize],
            scaled_residual_signal_buffer[(delayed_residual_signal_offset + i as isize) as usize],
        );
    }

    /*** Scale correlationMax, residualSignalEnergy and delayedResidualSignalEnergy to the best fit on 16 bits ***/
    /* these variables must fit on 16bits for the following computation, to avoid loosing information, scale them */
    /* at best fit: scale the higher of three to get the value over 2^14 and shift the other two from the same amount */
    /* Note: all three value are >= 0 */
    maximum_three = correlation_max;
    if maximum_three < residual_signal_energy {
        maximum_three = residual_signal_energy;
    }
    if maximum_three < delayed_residual_signal_energy {
        maximum_three = delayed_residual_signal_energy;
    }

    if maximum_three > 0 {
        /* if all of them a null, just do nothing otherwise shift right to get the max number in range [0x4000,0x8000[ */
        let leading_zeros = count_leading_zeros(maximum_three) as i16;
        if leading_zeros < 16 {
            correlation_max_word16 = shr32(correlation_max, (16 - leading_zeros) as u32) as Word16;
            residual_signal_energy_word16 =
                shr32(residual_signal_energy, (16 - leading_zeros) as u32) as Word16;
            delayed_residual_signal_energy_word16 =
                shr32(delayed_residual_signal_energy, (16 - leading_zeros) as u32) as Word16;
        } else {
            /* if the values already fit on 16 bits, no need to shift */
            correlation_max_word16 = correlation_max as Word16;
            residual_signal_energy_word16 = residual_signal_energy as Word16;
            delayed_residual_signal_energy_word16 = delayed_residual_signal_energy as Word16;
        }
    }

    /* eq78: Hp(z)=(1 + γp*gl*z(−T))/(1 + γp*gl) -> (with g=γp*gl) Hp(z)=1/(1+g) + (g/(1+g))*z(-T) = g0 + g1*z(-T) */
    /* g = gl/2 (as γp=0.5)= (eq83) correlationMax/(2*delayedResidualSignalEnergy) */
    /* compute g0 = 1/(1+g)=  delayedResidualSignalEnergy/(delayedResidualSignalEnergy+correlationMax/2) = 1-g1*/
    /* compute g1 = g/(1+g) = correlationMax/(2*delayedResidualSignalEnergy+correlationMax) = 1-g0 */

    /*** eq82 -> (correlationMax^2)/(residualSignalEnergy*delayedResidualSignalEnergy)<0.5 ***/
    /* (correlationMax^2) < (residualSignalEnergy*delayedResidualSignalEnergy)*0.5 */
    if (mult16_16(correlation_max_word16, correlation_max_word16) < shr(mult16_16(residual_signal_energy_word16, delayed_residual_signal_energy_word16), 1)) /* eq82 */
        || ((correlation_max_word16 == 0) && (delayed_residual_signal_energy_word16 == 0))
    {
        /* correlationMax and delayedResidualSignalEnergy values are 0 -> unable to compute g0 and g1 -> disable filter */
        /* long term post filter disabled */
        for i in 0..L_SUBFRAME {
            // decoderChannelContext->longTermFilteredResidualSignal[i] = residualSignal[i];
            long_term_filtered_residual_signal_buffer[1 + i] =
                residual_signal_buffer[residual_signal_offset + i];
        }
    } else {
        /* eq82 gives long term filter enabled, */
        let g0: Word16;
        let g1: Word16;
        /* eq83: gl = correlationMax/delayedResidualSignalEnergy bounded in ]0,1] */
        /* check if gl > 1 -> gl=1 -> g=1/2 -> g0=2/3 and g1=1/3 */
        if correlation_max > delayed_residual_signal_energy {
            g0 = 21845; /* 2/3 in Q15 */
            g1 = 10923; /* 1/3 in Q15 */
        } else {
            /* g1 = correlationMax/(2*delayedResidualSignalEnergy+correlationMax) */
            g1 = div32(
                shl32(correlation_max_word16 as Word32, 15),
                add32(
                    shl32(delayed_residual_signal_energy_word16 as Word32, 1),
                    correlation_max_word16 as Word32,
                ),
            ) as Word16; /* g1 in Q15 */
            g0 = sub16(32767, g1); /* g0 = 1 - g1 in Q15 */
        }

        /* longTermFilteredResidualSignal[i] = g0*residualSignal[i] + g1*delayedResidualSignal[i]*/
        // delayedResidualSignal = &(residualSignal[-bestIntPitchDelay]);
        let delayed_residual_signal_offset =
            residual_signal_offset as isize - best_int_pitch_delay as isize;
        for i in 0..L_SUBFRAME {
            long_term_filtered_residual_signal_buffer[1 + i] = saturate(
                pshr(
                    add32(
                        mult16_16(g0, residual_signal_buffer[residual_signal_offset + i]),
                        mult16_16(
                            g1,
                            residual_signal_buffer
                                [(delayed_residual_signal_offset + i as isize) as usize],
                        ),
                    ),
                    15,
                ),
                MAX_INT16 as Word32,
            ) as Word16;
        }
    }

    /********************************************************************/
    /* Tilt Compensation Filter                                         */
    /********************************************************************/

    /* compute hf the truncated (to 22 coefficients) impulse response of the filter A(z/γn)/A(z/γd) described in spec 4.2.2 eq84 */
    /* hf(i) = LPGammaNCoeff[i] - ∑[j:0..9]LPGammaDCoeff[j]*hf[i-j-1]) */
    /* GAMMA_XX constants are in Q15 */
    lp_gamma_d_coefficients[0] = mult16_16_p15(lp_coefficients[0], GAMMA_D1) as Word16;
    lp_gamma_d_coefficients[1] = mult16_16_p15(lp_coefficients[1], GAMMA_D2) as Word16;
    lp_gamma_d_coefficients[2] = mult16_16_p15(lp_coefficients[2], GAMMA_D3) as Word16;
    lp_gamma_d_coefficients[3] = mult16_16_p15(lp_coefficients[3], GAMMA_D4) as Word16;
    lp_gamma_d_coefficients[4] = mult16_16_p15(lp_coefficients[4], GAMMA_D5) as Word16;
    lp_gamma_d_coefficients[5] = mult16_16_p15(lp_coefficients[5], GAMMA_D6) as Word16;
    lp_gamma_d_coefficients[6] = mult16_16_p15(lp_coefficients[6], GAMMA_D7) as Word16;
    lp_gamma_d_coefficients[7] = mult16_16_p15(lp_coefficients[7], GAMMA_D8) as Word16;
    lp_gamma_d_coefficients[8] = mult16_16_p15(lp_coefficients[8], GAMMA_D9) as Word16;
    lp_gamma_d_coefficients[9] = mult16_16_p15(lp_coefficients[9], GAMMA_D10) as Word16;

    hf[0] = 4096; /* 1 in Q12 as LPGammaNCoefficients and LPGammaDCoefficient doesn't contain the first element which is 1 and past values of hf are 0 */
    for i in 1..11 {
        let mut acc = sshl(lp_gamma_n_coefficients[i - 1] as Word32, 12); /* LPGammaNCoefficients in Q12 -> acc in Q24 */
        for j in 0..NB_LSP_COEFF {
            if j < i {
                /* j<i to avoid access to negative index of hf(past values are 0 anyway) */
                acc = msu16_16(acc, lp_gamma_d_coefficients[j], hf[i - j - 1]); /* LPGammaDCoefficient in Q12, hf in Q12 -> Q24 TODO: Possible overflow?? */
            }
        }
        hf[i] = saturate(pshr(acc, 12), MAX_INT16 as Word32) as Word16; /* get result back in Q12 and saturate on 16 bits */
    }
    for i in 11..22 {
        let mut acc: Word32 = 0;
        for j in 0..NB_LSP_COEFF {
            /* j<i to avoid access to negative index of hf(past values are 0 anyway) */
            acc = msu16_16(acc, lp_gamma_d_coefficients[j], hf[i - j - 1]); /* LPGammaDCoefficient in Q12, hf in Q12 -> Q24 TODO: Possible overflow?? */
        }
        hf[i] = saturate(pshr(acc, 12), MAX_INT16 as Word32) as Word16; /* get result back in Q12 and saturate on 16 bits */
    }

    /* hf is then used to compute k'1 spec 4.2.3 eq87: k'1 = -rh1/rh0 */
    /* rh0 = ∑[i:0..21]hf[i]*hf[i] */
    /* rh1 = ∑[i:0..20]hf[i]*hf[i+1] */
    rh1 = mult16_16(hf[0], hf[1]);
    for i in 1..21 {
        rh1 = mac16_16(rh1, hf[i], hf[i + 1]); /* rh1 in Q24 */
    }

    /* tiltCompensationGain is set to 0 if k'1>0 -> rh1<0 (as rh0 is always>0) */
    if rh1 < 0 {
        /* tiltCompensationGain = 0 -> no gain filter is off, just copy the input */
        // memcpy(tiltCompensatedSignal, decoderChannelContext->longTermFilteredResidualSignal, L_SUBFRAME*sizeof(word16_t));
        for i in 0..L_SUBFRAME {
            tilt_compensated_signal[i] = long_term_filtered_residual_signal_buffer[1 + i];
        }
    } else {
        /*compute tiltCompensationGain = k'1*γt */
        let tilt_compensation_gain: Word16;
        let mut rh0 = mult16_16(hf[0], hf[0]);
        for i in 1..22 {
            rh0 = mac16_16(rh0, hf[i], hf[i]); /* rh0 in Q24 */
        }
        rh1 = mult16_32_q15(GAMMA_T, rh1); /* GAMMA_T in Q15, rh1 in Q24*/
        tilt_compensation_gain = saturate(div32(rh1, pshr(rh0, 12)), MAX_INT16 as Word32) as Word16; /* rh1 in Q24, PSHR(rh0,12) in Q12 -> tiltCompensationGain in Q12 */

        /* compute filter Ht (spec A.4.2.3 eqA14) = 1 + gain*z(-1) */
        for i in 0..L_SUBFRAME {
            tilt_compensated_signal[i] = msu16_16_q12(
                long_term_filtered_residual_signal_buffer[1 + i] as Word32,
                tilt_compensation_gain,
                long_term_filtered_residual_signal_buffer[1 + i - 1],
            ) as Word16;
        }
    }
    /* update memory word of longTermFilteredResidualSignal for next subframe */
    // decoderChannelContext->longTermFilteredResidualSignal[-1] = decoderChannelContext->longTermFilteredResidualSignal[L_SUBFRAME-1];
    long_term_filtered_residual_signal_buffer[0] =
        long_term_filtered_residual_signal_buffer[1 + L_SUBFRAME - 1];

    /********************************************************************/
    /* synthesis filter 1/[Â(z /γd)] spec A.4.2.2                       */
    /*                                                                  */
    /*   Note: Â(z/γn) was done before when computing residual signal   */
    /********************************************************************/
    /* shortTermFilteredResidualSignal is accessed in range [-NB_LSP_COEFF,L_SUBFRAME[ */
    // synthesisFilter(tiltCompensatedSignal, LPGammaDCoefficients, decoderChannelContext->shortTermFilteredResidualSignal);
    lp_synthesis_filter(
        &tilt_compensated_signal,
        &lp_gamma_d_coefficients,
        short_term_filtered_residual_signal_buffer,
    );

    /* get the last NB_LSP_COEFF of shortTermFilteredResidualSignal and set them as memory for next subframe(they do not overlap so use memcpy) */
    // memcpy(decoderChannelContext->shortTermFilteredResidualSignalBuffer, &(decoderChannelContext->shortTermFilteredResidualSignalBuffer[L_SUBFRAME]), NB_LSP_COEFF*sizeof(word16_t));
    for i in 0..NB_LSP_COEFF {
        short_term_filtered_residual_signal_buffer[i] =
            short_term_filtered_residual_signal_buffer[L_SUBFRAME + i];
    }

    /********************************************************************/
    /* Adaptive Gain Control spec A.4.2.4                               */
    /*                                                                  */
    /********************************************************************/

    /*** compute G(gain scaling factor) according to eqA15 : G = Sqrt((∑s(n)^2)/∑sf(n)^2 ) ***/
    /* compute ∑sf(n)^2, scale the signal shifting right by 4 to avoid possible overflow on 32 bits sum */
    for i in 0..L_SUBFRAME {
        // short_term_filtered_residual_signal_square_sum = umac16_16_q4(short_term_filtered_residual_signal_square_sum, short_term_filtered_residual_signal_buffer[NB_LSP_COEFF + i] as u16, short_term_filtered_residual_signal_buffer[NB_LSP_COEFF + i] as u16); /* inputs are both in Q0, output is in Q-4 */
        let val = short_term_filtered_residual_signal_buffer[NB_LSP_COEFF + i];
        let prod = mult16_16(val, val);
        let prod_q4 = shr(prod, 4);
        short_term_filtered_residual_signal_square_sum = uadd32(
            short_term_filtered_residual_signal_square_sum,
            prod_q4 as u32,
        );
    }

    /* if the sum is null we can't compute gain -> output of postfiltering is the output of shortTermFilter and previousAdaptativeGain is set to 0 */
    /* the reset of previousAdaptativeGain is not mentionned in the spec but in ITU code only */
    if short_term_filtered_residual_signal_square_sum == 0 {
        *previous_adaptative_gain = 0;
        for i in 0..L_SUBFRAME {
            post_filtered_signal[i] = short_term_filtered_residual_signal_buffer[NB_LSP_COEFF + i];
        }
    } else {
        /* we can compute adaptativeGain and output signal */
        let mut current_adaptative_gain: Word16;
        /* compute ∑s(n)^2 scale the signal shifting right by 4 to avoid possible overflow on 32 bits sum, same shift was applied at denominator */
        let mut reconstructed_speech_square_sum: u32 = 0;
        for i in 0..L_SUBFRAME {
            // reconstructed_speech_square_sum = umac16_16_q4(reconstructed_speech_square_sum, reconstructed_speech[NB_LSP_COEFF + i] as u16, reconstructed_speech[NB_LSP_COEFF + i] as u16); /* inputs are both in Q0, output is in Q-4 */
            let val = reconstructed_speech[NB_LSP_COEFF + i];
            let prod = mult16_16(val, val);
            let prod_q4 = shr(prod, 4);
            reconstructed_speech_square_sum =
                uadd32(reconstructed_speech_square_sum, prod_q4 as u32);
        }

        if reconstructed_speech_square_sum == 0 {
            /* numerator is null -> current gain is null */
            gain_scaling_factor = 0;
        } else {
            let mut fraction_result: u32; /* stores  ∑s(n)^2)/∑sf(n)^2 in Q10 on a 32 bit unsigned */
            let scaled_short_term_filtered_residual_signal_square_sum: u32;
            /* Compute ∑s(n)^2)/∑sf(n)^2  result shall be in Q10 */
            /* normalise the numerator on 32 bits */
            let numerator_shift = reconstructed_speech_square_sum.leading_zeros() as i16;
            reconstructed_speech_square_sum =
                ushl(reconstructed_speech_square_sum, numerator_shift as u32); /* reconstructedSpeechSquareSum*2^numeratorShift */

            /* normalise denominator to get the result directly in Q10 if possible */
            // scaledShortTermFilteredResidualSignalSquareSum = VSHR32(shortTermFilteredResidualSignalSquareSum, 10-numeratorShift); /* shortTermFilteredResidualSignalSquareSum*2^(numeratorShift-10)*/
            if 10 - numerator_shift >= 0 {
                scaled_short_term_filtered_residual_signal_square_sum =
                    short_term_filtered_residual_signal_square_sum >> (10 - numerator_shift);
            } else {
                scaled_short_term_filtered_residual_signal_square_sum =
                    short_term_filtered_residual_signal_square_sum << (numerator_shift - 10);
            }

            if scaled_short_term_filtered_residual_signal_square_sum == 0 {
                /* shift might have sent to zero the denominator */
                fraction_result = udiv32(
                    reconstructed_speech_square_sum,
                    short_term_filtered_residual_signal_square_sum,
                ); /* result in QnumeratorShift */
                // fractionResult = VSHR32(fractionResult, numeratorShift-10); /* result in Q10 */
                if numerator_shift - 10 >= 0 {
                    fraction_result = fraction_result >> (numerator_shift - 10);
                } else {
                    fraction_result = fraction_result << (10 - numerator_shift);
                }
            } else {
                /* ok denominator is still > 0 */
                fraction_result = udiv32(
                    reconstructed_speech_square_sum,
                    scaled_short_term_filtered_residual_signal_square_sum,
                ); /* result in Q10 */
            }
            /* now compute current Gain =  Sqrt((∑s(n)^2)/∑sf(n)^2 ) */
            /* g729Sqrt_Q0Q7(Q0)->Q7, by giving a Q10 as input, output is in Q12 */
            gain_scaling_factor =
                saturate(g729_sqrt_q0q7(fraction_result), MAX_INT16 as Word32) as Word16;

            /* multiply by 0.1 as described in spec A.4.2.4 */
            gain_scaling_factor = mult16_16_p15(gain_scaling_factor, 3277 as Word16) as Word16;
            /* in Q12, 3277 = 0.1 in Q15*/
        }
        /* Compute the signal according to eq89 (spec 4.2.4 and section A4.2.4) */
        /* currentGain = 0.9*previousGain + 0.1*gainScalingFactor the 0.1 factor has already been integrated in the variable gainScalingFactor */
        /* outputsignal = currentGain*shortTermFilteredResidualSignal */
        current_adaptative_gain = *previous_adaptative_gain;
        for i in 0..L_SUBFRAME {
            current_adaptative_gain = add16(
                gain_scaling_factor,
                mult16_16_p15(current_adaptative_gain, 29491 as Word16) as Word16,
            ); /* 29492 = 0.9 in Q15, result in Q12 */
            post_filtered_signal[i] = mult16_16_q12(
                current_adaptative_gain,
                short_term_filtered_residual_signal_buffer[NB_LSP_COEFF + i],
            ) as Word16;
        }
        *previous_adaptative_gain = current_adaptative_gain;
    }

    /* shift buffers if needed */
    if subframe_index > 0 {
        /* only after 2nd subframe treatment */
        /* shift left by L_FRAME the residualSignal and scaledResidualSignal buffers */
        // memmove(decoderChannelContext->residualSignalBuffer, &(decoderChannelContext->residualSignalBuffer[L_FRAME]), MAXIMUM_INT_PITCH_DELAY*sizeof(word16_t));
        // memmove(decoderChannelContext->scaledResidualSignalBuffer, &(decoderChannelContext->scaledResidualSignalBuffer[L_FRAME]), MAXIMUM_INT_PITCH_DELAY*sizeof(word16_t));
        for i in 0..MAXIMUM_INT_PITCH_DELAY {
            residual_signal_buffer[i] = residual_signal_buffer[L_FRAME + i];
            scaled_residual_signal_buffer[i] = scaled_residual_signal_buffer[L_FRAME + i];
        }
    }
}
