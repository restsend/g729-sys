#![no_main]

//! Fuzz the G.729 decoder with arbitrary payloads.
//!
//! Layout: `data[0]` packs the three flags (bit0=erased, bit1=sid, bit2=rfc3389),
//! `data[1..]` is the RTP payload handed to the decoder. Payloads of arbitrary
//! length (0 bytes, 2-byte SID, 10-byte voice, truncated, oversized, ...) must
//! never panic.

use g729_sys::Decoder;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if data.is_empty() {
        return;
    }
    let flags = data[0];
    let payload = &data[1..];

    let mut dec = Decoder::new();
    let _ = dec.decode(
        payload,
        flags & 0b001 != 0,
        flags & 0b010 != 0,
        flags & 0b100 != 0,
    );

    // Also drive a few more packets on the same decoder to exercise stateful
    // paths (concealment / CNG) after a random first packet.
    for chunk in payload.chunks(10).take(4) {
        let _ = dec.decode(chunk, false, false, false);
        let _ = dec.decode(chunk, true, false, false);
    }
});
