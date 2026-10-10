#![no_main]

//! Fuzz the G.729 encoder with arbitrary 16-bit PCM.
//!
//! Runs a longer session (several dozen frames) on a single stateful encoder so
//! that Annex B VAD initialisation/smoothing and DTX (SID/untransmitted
//! scheduling) are actually exercised.

use g729_sys::{Encoder, FRAME_SAMPLES, VOICE_FRAME_BYTES};
use libfuzzer_sys::fuzz_target;

struct Lcg(u64);
impl Lcg {
    #[inline]
    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as u32
    }
}

fuzz_target!(|data: &[u8]| {
    // Seed from the input so different inputs explore different PCM / sessions.
    let seed = {
        let mut s: u64 = 0x9e3779b97f4a7c15;
        for &b in data {
            s = s.wrapping_mul(1099511628211).wrapping_add(b as u64 + 1);
        }
        s
    };
    let mut lcg = Lcg(seed);

    let mut out = [0u8; VOICE_FRAME_BYTES];

    for vad in [false, true] {
        let mut enc = Encoder::new(vad);
        for frame in 0..64u32 {
            let mut pcm = [0i16; FRAME_SAMPLES];
            // occasional low-level (silence-like) frames to drive VAD.
            let silent = (frame % 5) == 4;
            for s in pcm.iter_mut() {
                if silent {
                    *s = (lcg.next_u32() % 401) as i16 - 200;
                } else {
                    *s = (lcg.next_u32() % 65536) as u16 as i16;
                }
            }
            let len = enc.encode_into(&pcm, &mut out) as usize;
            assert!(matches!(len, 0 | 2 | VOICE_FRAME_BYTES));
            if vad {
                let _ = enc.rfc3389_payload();
            }
        }
    }
});
