use crate::g729::basic_operations::*;
use crate::g729::cng::*;
use crate::g729::decode_adaptative_code_vector::*;
use crate::g729::decode_fixed_code_vector::*;
use crate::g729::decode_gains::*;
use crate::g729::decode_lsp::*;
use crate::g729::interpolate_q_lsp::interpolate_q_lsp;
use crate::g729::ld8k::*;
use crate::g729::lp_synthesis_filter::lp_synthesis_filter;
use crate::g729::post_filter::*;
use crate::g729::post_processing::*;
use crate::g729::q_lsp_2_lp::q_lsp_2_lp;
use crate::g729::utils::*;

/* buffers allocation */
static PREVIOUS_Q_LSP_INITIAL_VALUES: [Word16; NB_LSP_COEFF] = [
    30000, 26000, 21000, 15000, 8000, 0, -8000, -15000, -21000, -26000,
]; /* in Q0.15 the initials values for the previous qLSP buffer */

pub struct DecoderChannelContext {
    /*** buffers used in decoder bloc ***/
    pub previous_q_lsp: [Word16; NB_LSP_COEFF], /* previous quantised LSP in Q0.15 */
    pub excitation_vector: [Word16; L_PAST_EXCITATION + L_FRAME], /* in Q0 this vector contains:
                                                0->153 : the past excitation vector.(length is Max Pitch Delay: 144 + interpolation window size : 10)
                                                154-> 154+L_FRAME-1 : the current frame adaptative Code Vector first used to compute then the excitation vector */
    pub bounded_adaptative_codebook_gain: Word16, /* the pitch gain from last subframe bounded in range [0.2,0.8] in Q0.14 */
    pub adaptative_codebook_gain: Word16, /* the gains needs to be stored in case of frame erasure in Q14 */
    pub fixed_codebook_gain: Word16,      /* in Q14.1 */
    pub reconstructed_speech: [Word16; NB_LSP_COEFF + L_FRAME], /* in Q0, output of the LP synthesis filter, the first 10 words store the previous frame output */
    pub pseudo_random_seed: u16, /* seed used in the pseudo random number generator */
    pub cng_pseudo_random_seed: u16, /* seed used in the pseudo random number generator for CNG */

    /*** buffers used in decodeLSP bloc ***/
    pub last_q_lsf: [Word16; NB_LSP_COEFF], /* this buffer stores the last qLSF to be used in case of frame lost in Q2.13 */
    /* buffer to store the last 4 frames codewords, used to compute the current qLSF */
    pub previous_l_code_word: [[Word16; NB_LSP_COEFF]; MA_MAX_K], /* in Q2.13, buffer to store the last 4 frames codewords, used to compute the current qLSF */
    /* the values stored are the codewords computed from the codebooks and rearranged */
    pub last_valid_l0: u16, /* this one store the L0 of last valid frame to be used in case of frame erased */

    /*** buffer used in decodeAdaptativeCodeVector bloc ***/
    pub previous_int_pitch_delay: i16, /* store the last valid Integer Pitch Delay computed, used in case of parity error or frame erased */

    /*** buffer used in decodeGains bloc ***/
    pub previous_gain_prediction_error: [Word16; 4], /* the last four gain prediction error U(m) eq69 and eq72, spec3.9.1 in Q10*/

    /*** buffers used in postFilter bloc ***/
    pub residual_signal_buffer: [Word16; MAXIMUM_INT_PITCH_DELAY + L_FRAME], /* store the residual signal (current subframe and MAXIMUM_INT_PITCH_DELAY of previous values) in Q0 */
    pub scaled_residual_signal_buffer: [Word16; MAXIMUM_INT_PITCH_DELAY + L_FRAME], /* same as previous but in Q-2 */
    pub long_term_filtered_residual_signal_buffer: [Word16; 1 + L_SUBFRAME], /* the output of long term filter in Q0, need 1 word from previous subframe for tilt compensation filter */
    pub short_term_filtered_residual_signal_buffer: [Word16; NB_LSP_COEFF + L_SUBFRAME], /* the output of short term filter(synthesis filter) in Q0, need NB_LSP_COEFF word from previous subframe as filter memory */
    pub previous_adaptative_gain: Word16, /* previous gain for adaptative gain control */

    /*** buffers used in postProcessing bloc ***/
    pub input_x0: Word16,
    pub input_x1: Word16,
    pub output_y2: Word32,
    pub output_y1: Word32,

    /* SID frame management */
    pub cng_channel_context: CngChannelContext, /* store informations specific to CNG */
    pub previous_frame_is_active_flag: u8,      /* store last processed frame type */
}

pub fn init_bcg729_decoder_channel() -> DecoderChannelContext {
    let mut ctx = DecoderChannelContext {
        previous_q_lsp: PREVIOUS_Q_LSP_INITIAL_VALUES,
        excitation_vector: [0; L_PAST_EXCITATION + L_FRAME],
        bounded_adaptative_codebook_gain: BOUNDED_PITCH_GAIN_MIN,
        adaptative_codebook_gain: 0,
        fixed_codebook_gain: 0,
        reconstructed_speech: [0; NB_LSP_COEFF + L_FRAME],
        pseudo_random_seed: 21845,
        cng_pseudo_random_seed: CNG_DTX_RANDOM_SEED_INIT,
        last_q_lsf: [0; NB_LSP_COEFF],
        previous_l_code_word: [[0; NB_LSP_COEFF]; MA_MAX_K],
        last_valid_l0: 0,
        previous_int_pitch_delay: 0,
        previous_gain_prediction_error: [0; 4],
        residual_signal_buffer: [0; MAXIMUM_INT_PITCH_DELAY + L_FRAME],
        scaled_residual_signal_buffer: [0; MAXIMUM_INT_PITCH_DELAY + L_FRAME],
        long_term_filtered_residual_signal_buffer: [0; 1 + L_SUBFRAME],
        short_term_filtered_residual_signal_buffer: [0; NB_LSP_COEFF + L_SUBFRAME],
        previous_adaptative_gain: 0,
        input_x0: 0,
        input_x1: 0,
        output_y2: 0,
        output_y1: 0,
        cng_channel_context: init_bcg729_cng_channel(),
        previous_frame_is_active_flag: 1,
    };

    init_decode_lsp(
        &mut ctx.previous_l_code_word,
        &mut ctx.last_valid_l0,
        &mut ctx.last_q_lsf,
    );
    init_decode_adaptative_code_vector(&mut ctx.previous_int_pitch_delay);
    init_decode_gains(&mut ctx.previous_gain_prediction_error);
    init_post_filter(
        &mut ctx.residual_signal_buffer,
        &mut ctx.scaled_residual_signal_buffer,
        &mut ctx.long_term_filtered_residual_signal_buffer,
        &mut ctx.short_term_filtered_residual_signal_buffer,
        &mut ctx.previous_adaptative_gain,
    );
    init_post_processing(
        &mut ctx.output_y2,
        &mut ctx.output_y1,
        &mut ctx.input_x0,
        &mut ctx.input_x1,
    );

    ctx
}

pub fn bcg729_decoder(
    decoder_channel_context: &mut DecoderChannelContext,
    bit_stream: Option<&[u8]>,
    bit_stream_length: u8,
    frame_erasure_flag: u8,
    mut sid_frame_flag: u8,
    rfc3389_payload_flag: u8,
    signal: &mut [Word16],
) {
    let mut parameters = [0 as u16; NB_PARAMETERS];
    /* internal buffers which we do not need to keep between calls */
    let mut q_lsp = [0 as Word16; NB_LSP_COEFF]; /* store the qLSP coefficients in Q0.15 */
    let mut interpolated_q_lsp = [0 as Word16; NB_LSP_COEFF]; /* store the interpolated qLSP coefficient in Q0.15 */
    let mut lp = [0 as Word16; 2 * NB_LSP_COEFF]; /* store the 2 sets of LP coefficients in Q12 */
    let mut int_pitch_delay: i16 = 0; /* store the Pitch Delay in and out of decodeAdaptativeCodeVector, in for decodeFixedCodeVector */
    let mut fixed_codebook_vector = [0 as Word16; L_SUBFRAME]; /* the fixed Codebook Vector in Q1.13*/
    let mut post_filtered_signal = [0 as Word16; L_SUBFRAME]; /* store the postfiltered signal in Q0 */

    let parity_error_flag: u8;
    let mut parameters_index = 4; /* this is used to select the right parameter according to the subframe currently computed, start pointing to P1 */
    let mut lp_coefficients_index = 0; /* this is used to select the right LP Coefficients according to the subframe currently computed */

    /*** parse the bitstream and get all parameter into an array as in spec 4 - Table 8 ***/
    if let Some(bs) = bit_stream {
        if sid_frame_flag == 0 {
            parameters_bit_stream_2_array(bs, &mut parameters);
        }
    } else {
        for i in 0..NB_PARAMETERS {
            parameters[i] = 0;
        }
    }

    /* manage frameErasure and CNG as specified in B.27 */
    if frame_erasure_flag != 0 {
        if decoder_channel_context.previous_frame_is_active_flag != 0 {
            sid_frame_flag = 0;
        } else {
            sid_frame_flag = 1;
        }
    }

    /* this is a SID frame, process it using the dedicated function */
    if sid_frame_flag == 1 {
        decode_sid_frame(
            &mut decoder_channel_context.cng_channel_context,
            decoder_channel_context.previous_frame_is_active_flag,
            bit_stream,
            bit_stream_length,
            &mut decoder_channel_context.excitation_vector[L_PAST_EXCITATION..],
            &mut decoder_channel_context.previous_q_lsp,
            &mut lp,
            &mut decoder_channel_context.cng_pseudo_random_seed,
            &mut decoder_channel_context.previous_l_code_word,
            rfc3389_payload_flag,
        );
        decoder_channel_context.previous_frame_is_active_flag = 0;

        /* loop over the two subframes */
        for subframe_index in (0..L_FRAME).step_by(L_SUBFRAME) {
            /* reconstruct speech using LP synthesis filter spec 4.1.6 eq77 */
            /* excitationVector in Q0, LP in Q12, recontructedSpeech in Q0 -> +NB_LSP_COEFF on the index of this one because the first NB_LSP_COEFF elements store the previous frame filter output */
            // LPSynthesisFilter(&(decoderChannelContext->excitationVector[L_PAST_EXCITATION + subframeIndex]), &(LP[LPCoefficientsIndex]), &(decoderChannelContext->reconstructedSpeech[NB_LSP_COEFF+subframeIndex]) );
            lp_synthesis_filter(
                &decoder_channel_context.excitation_vector[L_PAST_EXCITATION + subframe_index..],
                &lp[lp_coefficients_index..],
                &mut decoder_channel_context.reconstructed_speech[subframe_index..],
            );

            /* NOTE: ITU code check for overflow after LP Synthesis Filter computation and if it happened, divide excitation buffer by 2 and recompute the LP Synthesis Filter */
            /*	here, possible overflows are managed directly inside the Filter by saturation at MAXINT16 on each result */

            /* postFilter */
            post_filter(
                &mut decoder_channel_context.residual_signal_buffer,
                &mut decoder_channel_context.scaled_residual_signal_buffer,
                &mut decoder_channel_context.long_term_filtered_residual_signal_buffer,
                &mut decoder_channel_context.short_term_filtered_residual_signal_buffer,
                &mut decoder_channel_context.previous_adaptative_gain,
                &lp[lp_coefficients_index..],
                &decoder_channel_context.reconstructed_speech[subframe_index..],
                decoder_channel_context.previous_int_pitch_delay,
                subframe_index,
                &mut post_filtered_signal,
            );

            /* postProcessing */
            post_processing(
                &mut decoder_channel_context.output_y2,
                &mut decoder_channel_context.output_y1,
                &mut decoder_channel_context.input_x0,
                &mut decoder_channel_context.input_x1,
                &mut post_filtered_signal,
            );

            /* copy postProcessing Output to the signal output buffer */
            for i in 0..L_SUBFRAME {
                signal[subframe_index + i] = post_filtered_signal[i];
            }

            /* increase LPCoefficient Indexes */
            lp_coefficients_index += NB_LSP_COEFF;
        }

        decoder_channel_context.bounded_adaptative_codebook_gain = BOUNDED_PITCH_GAIN_MIN;

        /* Shift Excitation Vector by L_FRAME left */
        // memmove(decoderChannelContext->excitationVector, &(decoderChannelContext->excitationVector[L_FRAME]), L_PAST_EXCITATION*sizeof(word16_t));
        for i in 0..L_PAST_EXCITATION {
            decoder_channel_context.excitation_vector[i] =
                decoder_channel_context.excitation_vector[L_FRAME + i];
        }
        /* Copy the last 10 words of reconstructed Speech to the begining of the array for next frame computation */
        // memcpy(decoderChannelContext->reconstructedSpeech, &(decoderChannelContext->reconstructedSpeech[L_FRAME]), NB_LSP_COEFF*sizeof(word16_t));
        for i in 0..NB_LSP_COEFF {
            decoder_channel_context.reconstructed_speech[i] =
                decoder_channel_context.reconstructed_speech[L_FRAME + i];
        }

        return;
    }

    decoder_channel_context.previous_frame_is_active_flag = 1;
    /* re-init the CNG pseudo random seed at each active frame spec B.4 */
    decoder_channel_context.cng_pseudo_random_seed = CNG_DTX_RANDOM_SEED_INIT; /* re-initialise CNG pseudo Random seed to 11111 according to ITU code */

    /*****************************************************************************************/
    /*** on frame basis : decodeLSP, interpolate them with previous ones and convert to LP ***/
    decode_lsp(
        &mut decoder_channel_context.previous_l_code_word,
        &mut decoder_channel_context.last_valid_l0,
        &mut decoder_channel_context.last_q_lsf,
        &parameters[0..4],
        &mut q_lsp,
        frame_erasure_flag,
    ); /* decodeLSP need the first 4 parameters: L0-L3 */

    interpolate_q_lsp(
        &decoder_channel_context.previous_q_lsp,
        &q_lsp,
        &mut interpolated_q_lsp,
    );
    /* copy the currentqLSP to previousqLSP buffer */
    for i in 0..NB_LSP_COEFF {
        decoder_channel_context.previous_q_lsp[i] = q_lsp[i];
    }

    /* call the qLSP2LP function for first subframe */
    q_lsp_2_lp(&interpolated_q_lsp, &mut lp[0..NB_LSP_COEFF]);
    /* call the qLSP2LP function for second subframe */
    q_lsp_2_lp(&q_lsp, &mut lp[NB_LSP_COEFF..]);

    /* check the parity on the adaptativeCodebookIndexSubframe1(P1) with the received one (P0)*/
    parity_error_flag = (compute_parity(parameters[4]) ^ parameters[5]) as u8;

    /* loop over the two subframes */
    for subframe_index in (0..L_FRAME).step_by(L_SUBFRAME) {
        /* decode the adaptative Code Vector */
        decode_adaptative_code_vector(
            &mut decoder_channel_context.previous_int_pitch_delay,
            subframe_index,
            parameters[parameters_index],
            parity_error_flag,
            frame_erasure_flag,
            &mut int_pitch_delay,
            &mut decoder_channel_context.excitation_vector,
        );
        if subframe_index == 0 {
            /* at first subframe we have P0 between P1 and C1 */
            parameters_index += 2;
        } else {
            parameters_index += 1;
        }

        /* in case of frame erasure we shall generate pseudoRandom signs and index for fixed code vector decoding according to spec 4.4.4 */
        if frame_erasure_flag != 0 {
            parameters[parameters_index] =
                pseudo_random(&mut decoder_channel_context.pseudo_random_seed) & 0x1fff; /* signs are set to the 13 LSB of the first pseudoRandom number */
            parameters[parameters_index + 1] =
                pseudo_random(&mut decoder_channel_context.pseudo_random_seed) & 0x000f;
            /* signs are set to the 4 LSB of the second pseudoRandom number */
        }

        /* decode the fixed Code Vector */
        decode_fixed_code_vector(
            parameters[parameters_index + 1],
            parameters[parameters_index],
            int_pitch_delay,
            decoder_channel_context.bounded_adaptative_codebook_gain,
            &mut fixed_codebook_vector,
        );
        parameters_index += 2;

        /* decode gains */
        decode_gains(
            &mut decoder_channel_context.previous_gain_prediction_error,
            parameters[parameters_index],
            parameters[parameters_index + 1],
            &fixed_codebook_vector,
            frame_erasure_flag,
            &mut decoder_channel_context.adaptative_codebook_gain,
            &mut decoder_channel_context.fixed_codebook_gain,
        );

        parameters_index += 2;

        /* update bounded Adaptative Codebook Gain (in Q14) according to eq47 */
        decoder_channel_context.bounded_adaptative_codebook_gain =
            decoder_channel_context.adaptative_codebook_gain;
        if decoder_channel_context.bounded_adaptative_codebook_gain > BOUNDED_PITCH_GAIN_MAX {
            decoder_channel_context.bounded_adaptative_codebook_gain = BOUNDED_PITCH_GAIN_MAX;
        }
        if decoder_channel_context.bounded_adaptative_codebook_gain < BOUNDED_PITCH_GAIN_MIN {
            decoder_channel_context.bounded_adaptative_codebook_gain = BOUNDED_PITCH_GAIN_MIN;
        }

        /* compute excitation vector according to eq75 */
        /* excitationVector = adaptative Codebook Vector * adaptativeCodebookGain + fixed Codebook Vector * fixedCodebookGain */
        /* the adaptative Codebook Vector is in the excitationVector buffer [L_PAST_EXCITATION + subframeIndex] */
        /* with adaptative Codebook Vector in Q0, adaptativeCodebookGain in Q14, fixed Codebook Vector in Q1.13 and fixedCodebookGain in Q14.1 -> result in Q14 on 32 bits */
        /* -> shift right 14 bits and store the value in Q0 in a 16 bits type */
        for i in 0..L_SUBFRAME {
            decoder_channel_context.excitation_vector[L_PAST_EXCITATION + subframe_index + i] =
                saturate(
                    pshr(
                        add32(
                            mult16_16(
                                decoder_channel_context.excitation_vector
                                    [L_PAST_EXCITATION + subframe_index + i],
                                decoder_channel_context.adaptative_codebook_gain,
                            ),
                            mult16_16(
                                fixed_codebook_vector[i],
                                decoder_channel_context.fixed_codebook_gain,
                            ),
                        ),
                        14,
                    ),
                    MAX_INT16 as Word32,
                ) as Word16;
        }

        /* reconstruct speech using LP synthesis filter spec 4.1.6 eq77 */
        /* excitationVector in Q0, LP in Q12, recontructedSpeech in Q0 -> +NB_LSP_COEFF on the index of this one because the first NB_LSP_COEFF elements store the previous frame filter output */
        // LPSynthesisFilter(&(decoderChannelContext->excitationVector[L_PAST_EXCITATION + subframeIndex]), &(LP[LPCoefficientsIndex]), &(decoderChannelContext->reconstructedSpeech[NB_LSP_COEFF+subframeIndex]) );
        lp_synthesis_filter(
            &decoder_channel_context.excitation_vector[L_PAST_EXCITATION + subframe_index..],
            &lp[lp_coefficients_index..],
            &mut decoder_channel_context.reconstructed_speech[subframe_index..],
        );

        /* NOTE: ITU code check for overflow after LP Synthesis Filter computation and if it happened, divide excitation buffer by 2 and recompute the LP Synthesis Filter */
        /*	here, possible overflows are managed directly inside the Filter by saturation at MAXINT16 on each result */

        /* postFilter */
        post_filter(
            &mut decoder_channel_context.residual_signal_buffer,
            &mut decoder_channel_context.scaled_residual_signal_buffer,
            &mut decoder_channel_context.long_term_filtered_residual_signal_buffer,
            &mut decoder_channel_context.short_term_filtered_residual_signal_buffer,
            &mut decoder_channel_context.previous_adaptative_gain,
            &lp[lp_coefficients_index..],
            &decoder_channel_context.reconstructed_speech[subframe_index..],
            int_pitch_delay,
            subframe_index,
            &mut post_filtered_signal,
        ); /* postProcessing */

        post_processing(
            &mut decoder_channel_context.output_y2,
            &mut decoder_channel_context.output_y1,
            &mut decoder_channel_context.input_x0,
            &mut decoder_channel_context.input_x1,
            &mut post_filtered_signal,
        );

        /* copy postProcessing Output to the signal output buffer */
        for i in 0..L_SUBFRAME {
            signal[subframe_index + i] = post_filtered_signal[i];
        }

        /* increase LPCoefficient Indexes */
        lp_coefficients_index += NB_LSP_COEFF;
    }

    /* Shift Excitation Vector by L_FRAME left */
    // memmove(decoderChannelContext->excitationVector, &(decoderChannelContext->excitationVector[L_FRAME]), L_PAST_EXCITATION*sizeof(word16_t));
    for i in 0..L_PAST_EXCITATION {
        decoder_channel_context.excitation_vector[i] =
            decoder_channel_context.excitation_vector[L_FRAME + i];
    }
    /* Copy the last 10 words of reconstructed Speech to the begining of the array for next frame computation */
    // memcpy(decoderChannelContext->reconstructedSpeech, &(decoderChannelContext->reconstructedSpeech[L_FRAME]), NB_LSP_COEFF*sizeof(word16_t));
    for i in 0..NB_LSP_COEFF {
        decoder_channel_context.reconstructed_speech[i] =
            decoder_channel_context.reconstructed_speech[L_FRAME + i];
    }
}

pub struct Decoder {
    context: DecoderChannelContext,
}

impl Decoder {
    pub fn new() -> Self {
        Decoder {
            context: init_bcg729_decoder_channel(),
        }
    }

    pub fn decode(
        &mut self,
        bit_stream: Option<&[u8]>,
        bit_stream_length: u8,
        frame_erasure_flag: u8,
        sid_frame_flag: u8,
        rfc3389_payload_flag: u8,
        signal: &mut [i16],
    ) {
        bcg729_decoder(
            &mut self.context,
            bit_stream,
            bit_stream_length,
            frame_erasure_flag,
            sid_frame_flag,
            rfc3389_payload_flag,
            signal,
        );
    }
}
