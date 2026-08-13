#![cfg_attr(not(feature = "std"), no_std)]

pub mod g729;

/// One G.729 frame contains 80 16-bit PCM samples at 8 kHz.
pub const FRAME_SAMPLES: usize = 80;
/// Voice frame payload length in bytes (10 bytes = 80 bits).
pub const VOICE_FRAME_BYTES: usize = 10;

/// Encoder wrapper.
pub struct Encoder {
    inner: g729::encoder::Encoder,
}

impl Encoder {
    /// Create a new G.729 encoder.
    ///
    /// `enable_vad` toggles Annex B (VAD/DTX).
    pub fn new(enable_vad: bool) -> Self {
        Self {
            inner: g729::encoder::Encoder::new(enable_vad),
        }
    }

    /// Encode one 80-sample frame into the caller-provided buffer.
    ///
    /// Returns the number of bytes written to `out` (always `<= VOICE_FRAME_BYTES`).
    /// This method is always available, including in `no_std`.
    pub fn encode_into(
        &mut self,
        input_80_samples: &[i16; FRAME_SAMPLES],
        out: &mut [u8; VOICE_FRAME_BYTES],
    ) -> u8 {
        let mut len: u8 = 0;
        self.inner.encode(input_80_samples, out, &mut len);
        len
    }

    /// Encode one 80-sample frame into a fresh `Vec<u8>`.
    ///
    /// Only available with the `std` feature (enabled by default).
    #[cfg(feature = "std")]
    pub fn encode(&mut self, input_80_samples: &[i16; FRAME_SAMPLES]) -> Vec<u8> {
        let mut out = [0u8; VOICE_FRAME_BYTES];
        let len = self.encode_into(input_80_samples, &mut out);
        out[..len as usize].to_vec()
    }

    pub fn rfc3389_payload(&mut self) -> [u8; 11] {
        // Not implemented in Rust backend yet, return zeros or implement if needed
        [0u8; 11]
    }
}

/// Decoder wrapper.
pub struct Decoder {
    inner: g729::decoder::Decoder,
}

impl Default for Decoder {
    fn default() -> Self {
        Self::new()
    }
}

impl Decoder {
    /// Create a new G.729 decoder.
    pub fn new() -> Self {
        Self {
            inner: g729::decoder::Decoder::new(),
        }
    }

    pub fn decode(
        &mut self,
        payload: &[u8],
        frame_erased: bool,
        is_sid: bool,
        rfc3389: bool,
    ) -> [i16; FRAME_SAMPLES] {
        let mut out = [0i16; FRAME_SAMPLES];
        let len = payload.len() as u8;
        let payload_opt = if len > 0 { Some(payload) } else { None };

        self.inner.decode(
            payload_opt,
            len,
            frame_erased as u8,
            is_sid as u8,
            rfc3389 as u8,
            &mut out,
        );
        out
    }
}
