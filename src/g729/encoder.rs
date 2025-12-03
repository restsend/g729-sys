use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;
// use crate::g729::fixed_point_math::*;
use crate::g729::adaptative_codebook_search::*;
use crate::g729::compute_adaptative_codebook_gain::*;
use crate::g729::compute_lp::*;
use crate::g729::compute_weighted_speech::*;
use crate::g729::find_open_loop_pitch_delay::*;
use crate::g729::fixed_codebook_search::*;
use crate::g729::gain_quantization::*;
use crate::g729::interpolate_q_lsp::*;
use crate::g729::lp2lsp_conversion::*;
use crate::g729::lp_synthesis_filter::lp_synthesis_filter;
use crate::g729::lsp_quantization::*;
use crate::g729::pre_processing::*;
use crate::g729::q_lsp_2_lp::*;
use crate::g729::utils::*;

pub struct EncoderChannelContext {
    /* buffers used in decoder bloc */
    /* Signal buffer mapping : 240 word16_t length */
    pub signal_buffer: [Word16; L_LP_ANALYSIS_WINDOW],

    // Indices for signal buffer
    // signalLastInputFrame index = L_LP_ANALYSIS_WINDOW - L_FRAME
    // signalCurrentFrame index = L_LP_ANALYSIS_WINDOW - L_SUBFRAME - L_FRAME
    pub previous_lsp_coefficients: [Word16; NB_LSP_COEFF],
    pub previous_q_lsp_coefficients: [Word16; NB_LSP_COEFF],

    pub weighted_input_signal: [Word16; MAXIMUM_INT_PITCH_DELAY + L_FRAME],
    pub excitation_vector: [Word16; L_PAST_EXCITATION + L_FRAME],

    pub target_signal: [Word16; NB_LSP_COEFF + L_SUBFRAME],

    pub last_quantized_adaptative_codebook_gain: Word16,

    /*** buffer used in preProcessing ***/
    pub pre_processing_state: PreProcessingState,

    /*** buffer used in LSPQuantization ***/
    pub previous_q_lsf: [[Word16; NB_LSP_COEFF]; MA_MAX_K],

    /*** buffer used in gainQuantization ***/
    pub previous_gain_prediction_error: [Word16; 4],
    // VAD/DTX context placeholders
    // pub vad_channel_context: Option<VADChannelContext>,
    // pub dtx_channel_context: Option<DTXChannelContext>,
}

const PREVIOUS_LSP_INITIAL_VALUES: [Word16; NB_LSP_COEFF] = [
    30000, 26000, 21000, 15000, 8000, 0, -8000, -15000, -21000, -26000,
];

// Indices
const SIGNAL_LAST_INPUT_FRAME_IDX: usize = L_LP_ANALYSIS_WINDOW - L_FRAME;
const SIGNAL_CURRENT_FRAME_IDX: usize = L_LP_ANALYSIS_WINDOW - L_SUBFRAME - L_FRAME;

impl EncoderChannelContext {
    pub fn new(_enable_vad: bool) -> Self {
        let mut ctx = EncoderChannelContext {
            signal_buffer: [0; L_LP_ANALYSIS_WINDOW],
            previous_lsp_coefficients: PREVIOUS_LSP_INITIAL_VALUES,
            previous_q_lsp_coefficients: PREVIOUS_LSP_INITIAL_VALUES,
            weighted_input_signal: [0; MAXIMUM_INT_PITCH_DELAY + L_FRAME],
            excitation_vector: [0; L_PAST_EXCITATION + L_FRAME],
            target_signal: [0; NB_LSP_COEFF + L_SUBFRAME],
            last_quantized_adaptative_codebook_gain: O2_IN_Q14,
            pre_processing_state: PreProcessingState::new(),
            previous_q_lsf: [[0; NB_LSP_COEFF]; MA_MAX_K],
            previous_gain_prediction_error: [0; 4],
        };

        // initPreProcessing
        // Already done by PreProcessingState::new()

        // initLSPQuantization
        init_lsp_quantization(&mut ctx.previous_q_lsf);

        // initGainQuantization
        // In C: initGainQuantization sets previousGainPredictionError to -14dB in Q10.
        // -14dB in Q10 is -14336? No.
        // Let's check initGainQuantization in gainQuantization.c
        // It sets it to -14336.
        for i in 0..4 {
            ctx.previous_gain_prediction_error[i] = -14336;
        }

        ctx
    }

    pub fn encode(
        &mut self,
        input_frame: &[Word16],
        bit_stream: &mut [u8],
        bit_stream_length: &mut u8,
    ) {
        let mut parameters = [0u16; NB_PARAMETERS];

        // internal buffers
        let mut lp_coefficients = [0i16; NB_LSP_COEFF];
        let _lsf_coefficients = [0i16; NB_LSP_COEFF];
        let mut q_lp_coefficients = [0i16; 2 * NB_LSP_COEFF];
        let mut weighted_q_lp_coefficients = [0i16; 2 * NB_LSP_COEFF];
        let mut lsp_coefficients = [0i16; NB_LSP_COEFF];
        let mut q_lsp_coefficients = [0i16; NB_LSP_COEFF];
        let mut interpolated_q_lsp = [0i16; NB_LSP_COEFF];

        let mut parameters_index = 4;
        let mut impulse_response_input = [0i16; L_SUBFRAME];

        // VAD placeholders
        let mut reflection_coefficients = [0i32; NB_LSP_COEFF];
        let mut auto_correlation_coefficients = [0i32; NB_LSP_COEFF + 3];
        let mut no_lag_auto_correlation_coefficients = [0i32; NB_LSP_COEFF + 3];
        let mut auto_correlation_coefficients_scale = 0i8;

        // Pre-processing
        // preProcessing(encoderChannelContext, inputFrame, encoderChannelContext->signalLastInputFrame);

        self.pre_processing_state.pre_processing(
            input_frame,
            &mut self.signal_buffer[SIGNAL_LAST_INPUT_FRAME_IDX..],
        );

        // Compute LP
        // computeLP(encoderChannelContext->signalBuffer, LPCoefficients, reflectionCoefficients, ...);
        // signalBuffer is used as input.
        compute_lp(
            &self.signal_buffer,
            &mut lp_coefficients,
            &mut reflection_coefficients,
            &mut auto_correlation_coefficients,
            &mut no_lag_auto_correlation_coefficients,
            &mut auto_correlation_coefficients_scale,
            NB_LSP_COEFF + 1, // VAD disabled for now
        );

        // LP to LSP
        if !lp2lsp_conversion(&lp_coefficients, &mut lsp_coefficients) {
            lsp_coefficients.copy_from_slice(&self.previous_lsp_coefficients);
        }

        // VAD check would go here.
        *bit_stream_length = 10;

        // LSP Quantization
        lsp_quantization(
            &mut self.previous_q_lsf,
            &mut lsp_coefficients,
            &mut q_lsp_coefficients,
            &mut parameters,
        );
        // Interpolate qLSP
        interpolate_q_lsp(
            &self.previous_q_lsp_coefficients,
            &q_lsp_coefficients,
            &mut interpolated_q_lsp,
        );

        // Update previous qLSP
        self.previous_q_lsp_coefficients
            .copy_from_slice(&q_lsp_coefficients);

        // qLSP to LP
        // first subframe
        q_lsp_2_lp(&interpolated_q_lsp, &mut q_lp_coefficients[0..NB_LSP_COEFF]);
        // second subframe
        q_lsp_2_lp(&q_lsp_coefficients, &mut q_lp_coefficients[NB_LSP_COEFF..]);

        // Compute weighted qLP
        for i in 0..2 * NB_LSP_COEFF {
            // weightedqLPCoefficients[i] = qLPCoefficients[i]*Gamma^(i%10+1)
            // We can use a helper or just loop.
            // The C code unrolls it.
            let gamma_idx = i % NB_LSP_COEFF;
            let gamma = match gamma_idx {
                0 => GAMMA_E1,
                1 => GAMMA_E2,
                2 => GAMMA_E3,
                3 => GAMMA_E4,
                4 => GAMMA_E5,
                5 => GAMMA_E6,
                6 => GAMMA_E7,
                7 => GAMMA_E8,
                8 => GAMMA_E9,
                9 => GAMMA_E10,
                _ => 0,
            };
            weighted_q_lp_coefficients[i] = mult16_16_p15(q_lp_coefficients[i], gamma) as Word16;
        }

        // Compute weighted speech
        // computeWeightedSpeech(encoderChannelContext->signalCurrentFrame, qLPCoefficients, weightedqLPCoefficients, &(encoderChannelContext->weightedInputSignal[MAXIMUM_INT_PITCH_DELAY]), &(encoderChannelContext->excitationVector[L_PAST_EXCITATION]));
        // signalCurrentFrame is at SIGNAL_CURRENT_FRAME_IDX.
        // weightedInputSignal output starts at MAXIMUM_INT_PITCH_DELAY.
        // excitationVector output starts at L_PAST_EXCITATION.

        compute_weighted_speech(
            &self.signal_buffer[SIGNAL_CURRENT_FRAME_IDX - NB_LSP_COEFF..],
            &q_lp_coefficients,
            &weighted_q_lp_coefficients,
            &mut self.weighted_input_signal[MAXIMUM_INT_PITCH_DELAY - NB_LSP_COEFF..],
            &mut self.excitation_vector[L_PAST_EXCITATION..],
        );

        // Open loop pitch search
        let open_loop_pitch_delay = find_open_loop_pitch_delay(&self.weighted_input_signal);
        let mut int_pitch_delay_min = open_loop_pitch_delay as i16 - 3;
        if int_pitch_delay_min < 20 {
            int_pitch_delay_min = 20;
        }
        let mut int_pitch_delay_max = int_pitch_delay_min + 6;
        if int_pitch_delay_max > MAXIMUM_INT_PITCH_DELAY as i16 {
            int_pitch_delay_max = MAXIMUM_INT_PITCH_DELAY as i16;
            int_pitch_delay_min = MAXIMUM_INT_PITCH_DELAY as i16 - 6;
        }

        // Subframe loop
        impulse_response_input[0] = ONE_IN_Q12 as Word16;
        for i in 1..L_SUBFRAME {
            impulse_response_input[i] = 0;
        }

        let mut lp_coefficients_index = 0;

        for subframe_index in (0..L_FRAME).step_by(L_SUBFRAME) {
            let mut int_pitch_delay: i16 = 0;
            let mut frac_pitch_delay: i16 = 0;
            let adaptative_codebook_gain: Word16;

            let mut impulse_response_buffer = [0i16; NB_LSP_COEFF + L_SUBFRAME];
            let mut filtered_adaptative_codebook_vector = [0i16; NB_LSP_COEFF + L_SUBFRAME];
            let mut gain_quantization_xy: Word64 = 0;
            let mut gain_quantization_yy: Word64 = 0;
            let mut fixed_codebook_vector = [0i16; L_SUBFRAME];
            let mut convolved_fixed_codebook_vector = [0i16; L_SUBFRAME];
            let mut quantized_adaptative_codebook_gain: Word16 = 0;
            let mut quantized_fixed_codebook_gain: Word16 = 0;

            // Compute impulse response
            // synthesisFilter(impulseResponseInput, &(weightedqLPCoefficients[LPCoefficientsIndex]), &(impulseResponseBuffer[NB_LSP_COEFF]));
            // In Rust, synthesis_filter takes the full buffer and writes to the end.
            lp_synthesis_filter(
                &impulse_response_input,
                &weighted_q_lp_coefficients[lp_coefficients_index..],
                &mut impulse_response_buffer,
            );

            // Compute target signal
            // synthesisFilter( &(encoderChannelContext->excitationVector[L_PAST_EXCITATION+subframeIndex]), &(weightedqLPCoefficients[LPCoefficientsIndex]), &(encoderChannelContext->targetSignal[NB_LSP_COEFF]));
            lp_synthesis_filter(
                &self.excitation_vector[L_PAST_EXCITATION + subframe_index
                    ..L_PAST_EXCITATION + subframe_index + L_SUBFRAME],
                &weighted_q_lp_coefficients[lp_coefficients_index..],
                &mut self.target_signal,
            );

            // Adaptative Codebook Search
            let mut param_pitch_delay: u16 = 0;
            adaptative_codebook_search(
                &mut self.excitation_vector,
                L_PAST_EXCITATION + subframe_index,
                &mut int_pitch_delay_min,
                &mut int_pitch_delay_max,
                &impulse_response_buffer[NB_LSP_COEFF..],
                &self.target_signal[NB_LSP_COEFF..],
                &mut int_pitch_delay,
                &mut frac_pitch_delay,
                &mut param_pitch_delay,
                subframe_index as u16,
            );
            parameters[parameters_index] = param_pitch_delay;

            // Compute adaptative codebook gain
            // synthesisFilter(&(encoderChannelContext->excitationVector[L_PAST_EXCITATION + subframeIndex]), &(weightedqLPCoefficients[LPCoefficientsIndex]), &(filteredAdaptativeCodebookVector[NB_LSP_COEFF]));
            lp_synthesis_filter(
                &self.excitation_vector[L_PAST_EXCITATION + subframe_index
                    ..L_PAST_EXCITATION + subframe_index + L_SUBFRAME],
                &weighted_q_lp_coefficients[lp_coefficients_index..],
                &mut filtered_adaptative_codebook_vector,
            );

            adaptative_codebook_gain = compute_adaptative_codebook_gain(
                &self.target_signal[NB_LSP_COEFF..],
                &filtered_adaptative_codebook_vector[NB_LSP_COEFF..],
                &mut gain_quantization_xy,
                &mut gain_quantization_yy,
            );

            parameters_index += 1;
            if subframe_index == 0 {
                parameters[parameters_index] = compute_parity(parameters[parameters_index - 1]);
                parameters_index += 1;
            }

            // Fixed Codebook Search
            let mut param_fixed_codebook_idx: u16 = 0;
            let mut param_fixed_codebook_sign: u16 = 0;

            fixed_codebook_search(
                &self.target_signal[NB_LSP_COEFF..],
                &mut impulse_response_buffer[NB_LSP_COEFF..],
                int_pitch_delay,
                self.last_quantized_adaptative_codebook_gain,
                &filtered_adaptative_codebook_vector[NB_LSP_COEFF..],
                adaptative_codebook_gain,
                &mut param_fixed_codebook_idx,
                &mut param_fixed_codebook_sign,
                &mut fixed_codebook_vector,
                &mut convolved_fixed_codebook_vector,
            );
            parameters[parameters_index] = param_fixed_codebook_idx;
            parameters[parameters_index + 1] = param_fixed_codebook_sign;
            parameters_index += 2;

            // Gain Quantization
            let mut param_gain_stage1: u16 = 0;
            let mut param_gain_stage2: u16 = 0;

            gain_quantization(
                &self.target_signal[NB_LSP_COEFF..],
                &filtered_adaptative_codebook_vector[NB_LSP_COEFF..],
                &convolved_fixed_codebook_vector,
                &fixed_codebook_vector,
                gain_quantization_xy,
                gain_quantization_yy,
                &mut self.previous_gain_prediction_error,
                &mut quantized_adaptative_codebook_gain,
                &mut quantized_fixed_codebook_gain,
                &mut param_gain_stage1,
                &mut param_gain_stage2,
            );
            parameters[parameters_index] = param_gain_stage1;
            parameters[parameters_index + 1] = param_gain_stage2;
            parameters_index += 2;

            // Memory updates
            lp_coefficients_index += NB_LSP_COEFF;
            self.last_quantized_adaptative_codebook_gain = quantized_adaptative_codebook_gain;
            if self.last_quantized_adaptative_codebook_gain > ONE_POINT_2_IN_Q14 {
                self.last_quantized_adaptative_codebook_gain = ONE_POINT_2_IN_Q14;
            }
            if self.last_quantized_adaptative_codebook_gain < O2_IN_Q14 {
                self.last_quantized_adaptative_codebook_gain = O2_IN_Q14;
            }

            // Compute excitation
            for i in 0..L_SUBFRAME {
                self.excitation_vector[L_PAST_EXCITATION + subframe_index + i] = saturate(
                    pshr(
                        add32(
                            mult16_16(
                                self.excitation_vector[L_PAST_EXCITATION + subframe_index + i],
                                quantized_adaptative_codebook_gain,
                            ),
                            mult16_16(fixed_codebook_vector[i], quantized_fixed_codebook_gain),
                        ),
                        14,
                    ),
                    MAX_INT16 as Word32,
                )
                    as Word16;
            }

            // Update targetSignal memory
            let quantized_adaptative_codebook_gain_q13 =
                pshr(quantized_adaptative_codebook_gain as Word32, 1) as Word16;
            for i in 0..NB_LSP_COEFF {
                let acc = mac16_16(
                    mult16_16(
                        quantized_adaptative_codebook_gain_q13,
                        filtered_adaptative_codebook_vector[L_SUBFRAME + i],
                    ),
                    quantized_fixed_codebook_gain,
                    convolved_fixed_codebook_vector[L_SUBFRAME - NB_LSP_COEFF + i],
                );
                self.target_signal[i] = saturate(
                    sub32(self.target_signal[L_SUBFRAME + i] as Word32, pshr(acc, 13)),
                    MAX_INT16 as Word32,
                ) as Word16;
            }
        }

        // Frame basis memory updates
        // shift left by L_FRAME the signal buffer
        self.signal_buffer.copy_within(L_FRAME.., 0);

        // update previousLSP coefficient buffer
        self.previous_lsp_coefficients
            .copy_from_slice(&lsp_coefficients);
        self.previous_q_lsp_coefficients
            .copy_from_slice(&q_lsp_coefficients);

        // shift left by L_FRAME the weightedInputSignal buffer
        self.weighted_input_signal.copy_within(L_FRAME.., 0);

        // shift left by L_FRAME the excitationVector
        self.excitation_vector.copy_within(L_FRAME.., 0);

        // Convert array of parameters into bitStream
        parameters_array_2_bit_stream(&parameters, bit_stream);
    }
}

pub struct Encoder {
    context: EncoderChannelContext,
}

impl Encoder {
    pub fn new(enable_vad: bool) -> Self {
        Encoder {
            context: EncoderChannelContext::new(enable_vad),
        }
    }

    pub fn encode(
        &mut self,
        input_frame: &[i16],
        bit_stream: &mut [u8],
        bit_stream_length: &mut u8,
    ) {
        self.context
            .encode(input_frame, bit_stream, bit_stream_length);
    }
}
