#![deny(missing_docs)]
//! Safe wrappers around the bcg729 C library bundled in this crate.

use std::ffi::c_void;

#[allow(non_camel_case_types)]
type int16_t = i16;
#[allow(non_camel_case_types)]
type uint8_t = u8;

#[allow(non_camel_case_types)]
#[repr(C)]
struct bcg729EncoderChannelContextStruct(c_void);
#[allow(non_camel_case_types)]
#[repr(C)]
struct bcg729DecoderChannelContextStruct(c_void);

#[link(name = "bcg729", kind = "static")]
unsafe extern "C" {
    fn initBcg729EncoderChannel(enable_vad: uint8_t) -> *mut bcg729EncoderChannelContextStruct;
    fn closeBcg729EncoderChannel(ctx: *mut bcg729EncoderChannelContextStruct);
    fn bcg729Encoder(
        ctx: *mut bcg729EncoderChannelContextStruct,
        input_frame: *const int16_t,
        bit_stream: *mut uint8_t,
        bit_stream_length: *mut uint8_t,
    );
    fn bcg729GetRFC3389Payload(ctx: *mut bcg729EncoderChannelContextStruct, payload: *mut uint8_t);

    fn initBcg729DecoderChannel() -> *mut bcg729DecoderChannelContextStruct;
    fn closeBcg729DecoderChannel(ctx: *mut bcg729DecoderChannelContextStruct);
    fn bcg729Decoder(
        ctx: *mut bcg729DecoderChannelContextStruct,
        bit_stream: *const uint8_t,
        bit_stream_length: uint8_t,
        frame_erasure_flag: uint8_t,
        sid_frame_flag: uint8_t,
        rfc3389_payload_flag: uint8_t,
        signal: *mut int16_t,
    );
}

/// One G.729 frame contains 80 16-bit PCM samples at 8 kHz.
pub const FRAME_SAMPLES: usize = 80;
/// Voice frame payload length in bytes (10 bytes = 80 bits).
pub const VOICE_FRAME_BYTES: usize = 10;

/// Encoder wrapper.
pub struct Encoder {
    ctx: *mut bcg729EncoderChannelContextStruct,
}

impl Encoder {
    /// Create an encoder. Set enable_vad to true to enable VAD/DTX.
    pub fn new(enable_vad: bool) -> anyhow::Result<Self> {
        let ctx = unsafe { initBcg729EncoderChannel(if enable_vad { 1 } else { 0 }) };
        if ctx.is_null() {
            anyhow::bail!("initBcg729EncoderChannel returned null");
        }
        Ok(Self { ctx })
    }

    /// Encode one 10ms frame of 80 PCM samples. Returns the produced payload bytes.
    /// With VAD, this may return 0 (no frame), 2 (SID), or 10 (voice) bytes.
    pub fn encode(&mut self, input_80_samples: &[i16; FRAME_SAMPLES]) -> Vec<u8> {
        let mut out = [0u8; VOICE_FRAME_BYTES];
        let mut len: u8 = 0;
        unsafe {
            bcg729Encoder(
                self.ctx,
                input_80_samples.as_ptr(),
                out.as_mut_ptr(),
                &mut len,
            );
        }
        out[..len as usize].to_vec()
    }

    /// If last frame was SID, get RFC3389 CN payload (11 bytes) for comfort noise.
    pub fn rfc3389_payload(&mut self) -> [u8; 11] {
        let mut p = [0u8; 11];
        unsafe { bcg729GetRFC3389Payload(self.ctx, p.as_mut_ptr()) };
        p
    }
}

impl Drop for Encoder {
    fn drop(&mut self) {
        if !self.ctx.is_null() {
            unsafe { closeBcg729EncoderChannel(self.ctx) };
            self.ctx = std::ptr::null_mut();
        }
    }
}

/// Decoder wrapper.
pub struct Decoder {
    ctx: *mut bcg729DecoderChannelContextStruct,
}

impl Decoder {
    /// Create a decoder.
    pub fn new() -> anyhow::Result<Self> {
        let ctx = unsafe { initBcg729DecoderChannel() };
        if ctx.is_null() {
            anyhow::bail!("initBcg729DecoderChannel returned null");
        }
        Ok(Self { ctx })
    }

    /// Decode a payload into one 80-sample frame.
    /// Provide flags for erasure/SID/RFC3389 as appropriate.
    pub fn decode(
        &mut self,
        payload: &[u8],
        frame_erased: bool,
        is_sid: bool,
        rfc3389: bool,
    ) -> [i16; FRAME_SAMPLES] {
        let mut out = [0i16; FRAME_SAMPLES];
        let len = payload.len() as u8;
        unsafe {
            bcg729Decoder(
                self.ctx,
                payload.as_ptr(),
                len,
                frame_erased as u8,
                is_sid as u8,
                rfc3389 as u8,
                out.as_mut_ptr(),
            );
        }
        out
    }
}

impl Drop for Decoder {
    fn drop(&mut self) {
        if !self.ctx.is_null() {
            unsafe { closeBcg729DecoderChannel(self.ctx) };
            self.ctx = std::ptr::null_mut();
        }
    }
}
