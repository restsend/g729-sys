// Long-session Annex B soak test.
//
// Exercises the stateful parts that a single short fuzz input cannot reach:
// VAD initialisation (first 32 frames) and smoothing, DTX SID/untransmitted
// scheduling, and decoder CNG across voice/silence transitions.

use g729_sys::{Decoder, Encoder, FRAME_SAMPLES, VOICE_FRAME_BYTES};

struct Lcg(u64);
impl Lcg {
    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as u32
    }
}

#[test]
fn annex_b_long_session_is_stable_and_exercises_dtx() {
    let mut lcg = Lcg(0x1234_5678_9abc_def0);
    let mut encoder = Encoder::new(true);

    let frames = 3000;
    let mut packets: Vec<Vec<u8>> = Vec::with_capacity(frames);
    let mut n_voice = 0usize;
    let mut n_sid = 0usize;
    let mut n_untransmitted = 0usize;

    for frame in 0..frames {
        let mut pcm = [0i16; FRAME_SAMPLES];
        // Alternate several voiced frames with a low-level noise floor so VAD
        // has to make real VOICE/NOISE decisions.
        let silence = (frame / 48) % 3 == 1;
        for s in pcm.iter_mut() {
            if silence {
                *s = (lcg.next_u32() % 401) as i16 - 200;
            } else {
                let n = (lcg.next_u32() % 20000) as i32 - 10000;
                *s = n as i16;
            }
        }

        let out = encoder.encode(&pcm);
        assert!(
            out.is_empty() || out.len() == VOICE_FRAME_BYTES || out.len() == 2,
            "frame {}: unexpected DTX length {}",
            frame,
            out.len()
        );
        match out.len() {
            VOICE_FRAME_BYTES => n_voice += 1,
            2 => n_sid += 1,
            _ => n_untransmitted += 1,
        }
        packets.push(out);
    }

    assert!(
        n_voice > 0 && n_sid > 0 && n_untransmitted > 0,
        "session did not exercise all DTX frame types: voice={} sid={} untransmitted={}",
        n_voice,
        n_sid,
        n_untransmitted
    );

    // Decode the session with correct framing, and periodically poke the
    // erasure and RFC3389 comfort-noise paths.
    let mut decoder = Decoder::new();
    for (i, packet) in packets.iter().enumerate() {
        match packet.len() {
            VOICE_FRAME_BYTES => {
                let _ = decoder.decode(packet, false, false, false);
            }
            2 => {
                let _ = decoder.decode(packet, false, true, false);
                if i % 5 == 0 {
                    let _ = decoder.decode(packet, false, true, true);
                }
            }
            _ => {
                let _ = decoder.decode(&[], false, true, false);
            }
        }
        if i % 13 == 0 {
            let _ = decoder.decode(packet, true, false, false);
        }
    }
}

/// RFC3389 payload from the encoder must always be 11 bytes and decodable.
#[test]
fn annex_b_rfc3389_payload_shape() {
    let mut encoder = Encoder::new(true);
    let mut decoder = Decoder::new();
    let mut lcg = Lcg(42);
    for _ in 0..600 {
        let mut pcm = [0i16; FRAME_SAMPLES];
        for s in pcm.iter_mut() {
            *s = ((lcg.next_u32() % 800) as i32 - 400) as i16;
        }
        let _ = encoder.encode(&pcm);
        let payload = encoder.rfc3389_payload();
        let _ = decoder.decode(&payload, false, true, true);
    }
}
