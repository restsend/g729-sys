#![no_main]

//! Fuzz encode → decode round-trips, and decode of mutated packets.

use g729_sys::{Decoder, Encoder, FRAME_SAMPLES, VOICE_FRAME_BYTES};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut enc = Encoder::new(false);
    let mut dec = Decoder::new();

    let mut pcm = [0i16; FRAME_SAMPLES];
    for (i, s) in pcm.iter_mut().enumerate() {
        let lo = data.get(i * 2).copied().unwrap_or(0) as u16;
        let hi = data.get(i * 2 + 1).copied().unwrap_or(0) as u16;
        *s = (lo | (hi << 8)) as i16;
    }

    let mut pkt = [0u8; VOICE_FRAME_BYTES];
    let len = enc.encode_into(&pcm, &mut pkt) as usize;
    let _ = dec.decode(&pkt[..len], false, false, false);

    // Decode the (possibly mutated) tail bytes as an arbitrary packet too.
    let _ = dec.decode(data, false, false, false);
    let _ = dec.decode(data, true, false, false);
});
