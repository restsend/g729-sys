//! Regression tests for fuzzing-found panics.
//!
//! The decoder must not panic on malformed/truncated RTP payloads (it used to
//! index past the end of the bit-stream buffer in
//! `parameters_bit_stream_2_array`).

use g729_sys::Decoder;

#[test]
fn decode_short_payloads_do_not_panic() {
    let mut dec = Decoder::new();
    let flags = [
        (false, false, false),
        (true, false, false),
        (false, true, false),
        (false, false, true),
        (true, true, true),
    ];
    for len in 0..=16usize {
        let payload = vec![0xAAu8; len];
        for (erased, sid, rfc) in flags {
            // Must not panic for any length (0-byte, truncated voice frame, ...).
            let _ = dec.decode(&payload, erased, sid, rfc);
        }
    }
}

#[test]
fn decode_arbitrary_payloads_do_not_panic() {
    // A spread of lengths around the interesting boundaries: 0, 1, 2 (SID),
    // 10 (voice), and just below/above.
    let mut dec = Decoder::new();
    let mut seed: u32 = 0x1234_5678;
    for len in [0usize, 1, 2, 3, 9, 10, 11, 20, 40, 100, 255, 256, 300] {
        for _ in 0..64 {
            let payload: Vec<u8> = (0..len)
                .map(|_| {
                    seed = seed.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    (seed >> 24) as u8
                })
                .collect();
            let _ = dec.decode(&payload, false, false, false);
            let _ = dec.decode(&payload, true, false, false);
        }
    }
}
