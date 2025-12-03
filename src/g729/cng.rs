use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;

pub struct CngChannelContext {
    pub received_sid_gain: Word16,     /* gain in Q3 */
    pub smoothed_sid_gain: Word16,     /* gain in Q3 */
    pub q_lsp: [Word16; NB_LSP_COEFF], /* qLSP in Q0.15 */
    pub last_frame_energy: Word64,     /* in Q0 */
}

static SID_Q_LSP_INITIAL_VALUES: [Word16; NB_LSP_COEFF] = [
    31441, 27566, 21458, 13612, 4663, -4663, -13612, -21458, -27566, -31441,
];

pub fn init_bcg729_cng_channel() -> CngChannelContext {
    CngChannelContext {
        received_sid_gain: 0,
        smoothed_sid_gain: 0,
        q_lsp: SID_Q_LSP_INITIAL_VALUES,
        last_frame_energy: 0,
    }
}
pub fn decode_sid_frame(
    _cng_channel_context: &mut CngChannelContext,
    _previous_frame_is_active_flag: u8,
    _bit_stream: Option<&[u8]>,
    _bit_stream_length: u8,
    excitation_vector: &mut [Word16],
    _previous_q_lsp: &mut [Word16; NB_LSP_COEFF],
    _lp: &mut [Word16],
    _pseudo_random_seed: &mut u16,
    _previous_l_code_word: &mut [[Word16; NB_LSP_COEFF]; MA_MAX_K],
    _rfc3389_payload_flag: u8,
) {
    // TODO: Implement CNG/SID frame decoding
    // For now, just fill excitation with zeros to avoid noise if called
    for x in excitation_vector.iter_mut() {
        *x = 0;
    }
}
