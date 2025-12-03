use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;

/* coefficients a stored in Q1.12 for A1 (the only one having a float value > 1) and Q0.12 for all the others */
const A1: Word16 = 7807;
const A2: Word16 = -3733;
const B0: Word16 = 1899;
const B1: Word16 = -3798;
const B2: Word16 = 1899;

pub struct PreProcessingState {
    output_y2: Word32,
    output_y1: Word32,
    input_x0: Word16,
    input_x1: Word16,
}

impl PreProcessingState {
    pub fn new() -> Self {
        Self {
            output_y2: 0,
            output_y1: 0,
            input_x0: 0,
            input_x1: 0,
        }
    }

    pub fn pre_processing(&mut self, signal: &[Word16], pre_processed_signal: &mut [Word16]) {
        let mut input_x2: Word16;
        let mut acc: Word32;

        for i in 0..L_FRAME {
            input_x2 = self.input_x1;
            self.input_x1 = self.input_x0;
            self.input_x0 = signal[i];

            /* compute with acc and coefficients in Q12 */
            acc = mult16_32_q12(A1, self.output_y1);
            acc = mac16_32_q12(acc, A2, self.output_y2);
            acc = mac16_16(acc, self.input_x0, B0);
            acc = mac16_16(acc, self.input_x1, B1);
            acc = mac16_16(acc, input_x2, B2);

            acc = saturate(acc, MAXINT28);

            pre_processed_signal[i] = pshr(acc, 12) as Word16;
            self.output_y2 = self.output_y1;
            self.output_y1 = acc;
        }
    }
}
