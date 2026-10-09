#![no_main]

//! Fuzz the G.729 encoder with arbitrary 16-bit PCM.

use g729_sys::{Encoder, FRAME_SAMPLES, VOICE_FRAME_BYTES};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut enc_vad = Encoder::new(true);
    let mut enc = Encoder::new(false);

    // Build an 80-sample frame from the input (little-endian), padding with 0.
    let mut pcm = [0i16; FRAME_SAMPLES];
    for (i, s) in pcm.iter_mut().enumerate() {
        let lo = data.get(i * 2).copied().unwrap_or(0) as u16;
        let hi = data.get(i * 2 + 1).copied().unwrap_or(0) as u16;
        *s = (lo | (hi << 8)) as i16;
    }

    let mut out = [0u8; VOICE_FRAME_BYTES];
    let _ = enc.encode_into(&pcm, &mut out);
    let _ = enc_vad.encode_into(&pcm, &mut out);

    // Feed several frames to exercise inter-frame state.
    for f in 0..3 {
        for s in pcm.iter_mut() {
            *s = s.wrapping_add(f * 997);
        }
        let _ = enc.encode_into(&pcm, &mut out);
    }
});
