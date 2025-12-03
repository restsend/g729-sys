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
    pub fn new(enable_vad: bool) -> anyhow::Result<Self> {
        Ok(Self {
            inner: g729::encoder::Encoder::new(enable_vad),
        })
    }

    pub fn encode(&mut self, input_80_samples: &[i16; FRAME_SAMPLES]) -> Vec<u8> {
        let mut out = [0u8; VOICE_FRAME_BYTES];
        let mut len: u8 = 0;
        self.inner.encode(input_80_samples, &mut out, &mut len);
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

impl Decoder {
    pub fn new() -> anyhow::Result<Self> {
        Ok(Self {
            inner: g729::decoder::Decoder::new(),
        })
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
