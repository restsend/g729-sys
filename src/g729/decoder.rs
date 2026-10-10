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

static PREVIOUS_Q_LSP_INITIAL_VALUES: [i16; NB_LSP_COEFF] = [
    30000, 26000, 21000, 15000, 8000, 0, -8000, -15000, -21000, -26000,
]; /* in Q0.15 the initials values for the previous qLSP buffer */

pub struct DecoderChannelContext {
    /*** buffers used in decoder bloc ***/
    pub previous_q_lsp: [i16; NB_LSP_COEFF], /* previous quantised LSP in Q0.15 */
    pub excitation_vector: [i16; L_PAST_EXCITATION + L_FRAME],
    pub bounded_adaptative_codebook_gain: i16,
    pub adaptative_codebook_gain: i16, /* the gains needs to be stored in case of frame erasure in Q14 */
    pub fixed_codebook_gain: i16,      /* in Q14.1 */
    pub reconstructed_speech: [i16; NB_LSP_COEFF + L_FRAME], /* in Q0, output of the LP synthesis filter, the first 10 words store the previous frame output */
    pub pseudo_random_seed: u16,
    pub cng_pseudo_random_seed: u16,

    /*** buffers used in decodeLSP bloc ***/
    pub last_q_lsf: [i16; NB_LSP_COEFF], /* this buffer stores the last qLSF to be used in case of frame lost in Q2.13 */

    pub previous_l_code_word: [[i16; NB_LSP_COEFF]; MA_MAX_K], /* in Q2.13, buffer to store the last 4 frames codewords, used to compute the current qLSF */

    pub last_valid_l0: u16, /* this one store the L0 of last valid frame to be used in case of frame erased */

    /*** buffer used in decodeAdaptativeCodeVector bloc ***/
    pub previous_int_pitch_delay: i16,

    /*** buffer used in decodeGains bloc ***/
    pub previous_gain_prediction_error: [i16; 4], /* the last four gain prediction error U(m) eq69 and eq72, spec3.9.1 in Q10*/

    /*** buffers used in postFilter bloc ***/
    pub residual_signal_buffer: [i16; MAXIMUM_INT_PITCH_DELAY + L_FRAME], /* store the residual signal (current subframe and MAXIMUM_INT_PITCH_DELAY of previous values) in Q0 */
    pub scaled_residual_signal_buffer: [i16; MAXIMUM_INT_PITCH_DELAY + L_FRAME], /* same as previous but in Q-2 */
    pub long_term_filtered_residual_signal_buffer: [i16; 1 + L_SUBFRAME], /* the output of long term filter in Q0, need 1 word from previous subframe for tilt compensation filter */
    pub short_term_filtered_residual_signal_buffer: [i16; NB_LSP_COEFF + L_SUBFRAME], /* the output of short term filter(synthesis filter) in Q0, need NB_LSP_COEFF word from previous subframe as filter memory */
    pub previous_adaptative_gain: i16,

    /*** buffers used in postProcessing bloc ***/
    pub input_x0: i16,
    pub input_x1: i16,
    pub output_y2: i32,
    pub output_y1: i32,

    pub cng_channel_context: CngChannelContext,
    pub previous_frame_is_active_flag: u8,
}

pub fn init_bcg729_decoder_channel() -> DecoderChannelContext {
    let (previous_l_code_word, last_valid_l0, last_q_lsf) = init_decode_lsp();
    let (
        residual_signal_buffer,
        scaled_residual_signal_buffer,
        long_term_filtered_residual_signal_buffer,
        short_term_filtered_residual_signal_buffer,
        previous_adaptative_gain,
    ) = init_post_filter();
    let (output_y2, output_y1, input_x0, input_x1) = init_post_processing();

    DecoderChannelContext {
        previous_q_lsp: PREVIOUS_Q_LSP_INITIAL_VALUES,
        excitation_vector: [0; L_PAST_EXCITATION + L_FRAME],
        bounded_adaptative_codebook_gain: BOUNDED_PITCH_GAIN_MIN,
        adaptative_codebook_gain: 0,
        fixed_codebook_gain: 0,
        reconstructed_speech: [0; NB_LSP_COEFF + L_FRAME],
        pseudo_random_seed: 21845,
        cng_pseudo_random_seed: CNG_DTX_RANDOM_SEED_INIT,
        last_q_lsf,
        previous_l_code_word,
        last_valid_l0,
        previous_int_pitch_delay: init_decode_adaptative_code_vector(),
        previous_gain_prediction_error: init_decode_gains(),
        residual_signal_buffer,
        scaled_residual_signal_buffer,
        long_term_filtered_residual_signal_buffer,
        short_term_filtered_residual_signal_buffer,
        previous_adaptative_gain,
        input_x0,
        input_x1,
        output_y2,
        output_y1,
        cng_channel_context: init_bcg729_cng_channel(),
        previous_frame_is_active_flag: 1,
    }
}

pub fn bcg729_decoder(
    decoder_channel_context: &mut DecoderChannelContext,
    bit_stream: Option<&[u8]>,
    bit_stream_length: u8,
    frame_erasure_flag: u8,
    mut sid_frame_flag: u8,
    rfc3389_payload_flag: u8,
    signal: &mut [i16],
) {
    let mut parameters = [0_u16; NB_PARAMETERS];

    let mut q_lsp = [0_i16; NB_LSP_COEFF]; /* store the qLSP coefficients in Q0.15 */
    let mut interpolated_q_lsp = [0_i16; NB_LSP_COEFF]; /* store the interpolated qLSP coefficient in Q0.15 */
    let mut lp = [0_i16; 2 * NB_LSP_COEFF]; /* store the 2 sets of LP coefficients in Q12 */
    let mut int_pitch_delay: i16 = 0;
    let mut fixed_codebook_vector = [0_i16; L_SUBFRAME]; /* the fixed Codebook Vector in Q1.13*/
    let mut post_filtered_signal = [0_i16; L_SUBFRAME]; /* store the postfiltered signal in Q0 */

    let mut parameters_index = 4; /* this is used to select the right parameter according to the subframe currently computed, start pointing to P1 */
    let mut lp_coefficients_index = 0;

    /*** parse the bitstream and get all parameter into an array as in spec 4 - Table 8 ***/
    if let Some(bs) = bit_stream {
        if sid_frame_flag == 0 {
            parameters_bit_stream_2_array(bs, &mut parameters);
        }
    } else {
        parameters[..NB_PARAMETERS].fill(0);
    }

    /* manage frameErasure and CNG as specified in B.27 */
    if frame_erasure_flag != 0 {
        if decoder_channel_context.previous_frame_is_active_flag != 0 {
            sid_frame_flag = 0;
        } else {
            sid_frame_flag = 1;
        }
    }

    if sid_frame_flag == 1 {
        decode_sid_frame(
            &mut decoder_channel_context.cng_channel_context,
            decoder_channel_context.previous_frame_is_active_flag,
            bit_stream,
            bit_stream_length,
            &mut decoder_channel_context.excitation_vector,
            &mut decoder_channel_context.previous_q_lsp,
            &mut lp,
            &mut decoder_channel_context.cng_pseudo_random_seed,
            &mut decoder_channel_context.previous_l_code_word,
            rfc3389_payload_flag,
        );
        decoder_channel_context.previous_frame_is_active_flag = 0;

        for subframe_index in (0..L_FRAME).step_by(L_SUBFRAME) {
            /* reconstruct speech using LP synthesis filter spec 4.1.6 eq77 */

            lp_synthesis_filter(
                &decoder_channel_context.excitation_vector[L_PAST_EXCITATION + subframe_index..],
                &lp[lp_coefficients_index..],
                &mut decoder_channel_context.reconstructed_speech[subframe_index..],
            );

            /* NOTE: ITU code check for overflow after LP Synthesis Filter computation and if it happened, divide excitation buffer by 2 and recompute the LP Synthesis Filter */

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

            post_processing(
                &mut decoder_channel_context.output_y2,
                &mut decoder_channel_context.output_y1,
                &mut decoder_channel_context.input_x0,
                &mut decoder_channel_context.input_x1,
                &mut post_filtered_signal,
            );

            signal[subframe_index..subframe_index + L_SUBFRAME]
                .copy_from_slice(&post_filtered_signal[..L_SUBFRAME]);

            lp_coefficients_index += NB_LSP_COEFF;
        }

        decoder_channel_context.bounded_adaptative_codebook_gain = BOUNDED_PITCH_GAIN_MIN;

        for i in 0..L_PAST_EXCITATION {
            decoder_channel_context.excitation_vector[i] =
                decoder_channel_context.excitation_vector[L_FRAME + i];
        }

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

    decoder_channel_context
        .previous_q_lsp
        .copy_from_slice(&q_lsp);

    /* call the qLSP2LP function for first subframe */
    q_lsp_2_lp(&interpolated_q_lsp, &mut lp[0..NB_LSP_COEFF]);
    /* call the qLSP2LP function for second subframe */
    q_lsp_2_lp(&q_lsp, &mut lp[NB_LSP_COEFF..]);

    /* check the parity on the adaptativeCodebookIndexSubframe1(P1) with the received one (P0)*/
    let parity_error_flag: u8 = (compute_parity(parameters[4]) ^ parameters[5]) as u8;

    for subframe_index in (0..L_FRAME).step_by(L_SUBFRAME) {
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
                pseudo_random(&mut decoder_channel_context.pseudo_random_seed) & 0x1fff;
            parameters[parameters_index + 1] =
                pseudo_random(&mut decoder_channel_context.pseudo_random_seed) & 0x000f;
        }

        decode_fixed_code_vector(
            parameters[parameters_index + 1],
            parameters[parameters_index],
            int_pitch_delay,
            decoder_channel_context.bounded_adaptative_codebook_gain,
            &mut fixed_codebook_vector,
        );
        parameters_index += 2;

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
        decoder_channel_context.bounded_adaptative_codebook_gain = decoder_channel_context
            .adaptative_codebook_gain
            .clamp(BOUNDED_PITCH_GAIN_MIN, BOUNDED_PITCH_GAIN_MAX);

        /* compute excitation vector according to eq75 */

        let adaptative_gain = decoder_channel_context.adaptative_codebook_gain;
        let fixed_gain = decoder_channel_context.fixed_codebook_gain;
        for (i, sample) in decoder_channel_context.excitation_vector
            [L_PAST_EXCITATION + subframe_index..L_PAST_EXCITATION + subframe_index + L_SUBFRAME]
            .iter_mut()
            .enumerate()
        {
            *sample = saturate(
                pshr(
                    add32(
                        mult16_16(*sample, adaptative_gain),
                        mult16_16(fixed_codebook_vector[i], fixed_gain),
                    ),
                    14,
                ),
                MAX_INT16 as i32,
            ) as i16;
        }

        /* reconstruct speech using LP synthesis filter spec 4.1.6 eq77 */

        lp_synthesis_filter(
            &decoder_channel_context.excitation_vector[L_PAST_EXCITATION + subframe_index..],
            &lp[lp_coefficients_index..],
            &mut decoder_channel_context.reconstructed_speech[subframe_index..],
        );

        /* NOTE: ITU code check for overflow after LP Synthesis Filter computation and if it happened, divide excitation buffer by 2 and recompute the LP Synthesis Filter */

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
        );

        post_processing(
            &mut decoder_channel_context.output_y2,
            &mut decoder_channel_context.output_y1,
            &mut decoder_channel_context.input_x0,
            &mut decoder_channel_context.input_x1,
            &mut post_filtered_signal,
        );

        signal[subframe_index..subframe_index + L_SUBFRAME]
            .copy_from_slice(&post_filtered_signal[..L_SUBFRAME]);

        lp_coefficients_index += NB_LSP_COEFF;
    }

    for i in 0..L_PAST_EXCITATION {
        decoder_channel_context.excitation_vector[i] =
            decoder_channel_context.excitation_vector[L_FRAME + i];
    }

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
