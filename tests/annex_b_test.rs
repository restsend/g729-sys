// Annex B (VAD/DTX/CNG) bit-exact regression test.
//
// The reference dumps under fixtures/ were produced by bcg729 (commit
// fb195811fe630415c9bf433c0c0bebde14314ff2) using its own encoder/decoder with
// VAD enabled on fixtures/en_8k_16bit.pcm. This test asserts that the Rust
// implementation produces byte-identical bitstreams (including SID/DTX frame
// lengths and the RFC3389 payload) and sample-identical decoded output.

use g729_sys::{Decoder, Encoder, FRAME_SAMPLES};

fn read_input_pcm() -> Vec<i16> {
    let bytes = std::fs::read("fixtures/en_8k_16bit.pcm").expect("read input pcm");
    bytes
        .as_chunks::<2>()
        .0
        .iter()
        .map(|c| i16::from_le_bytes(*c))
        .collect()
}

struct RefFrame {
    len: usize,
    bytes: Vec<u8>,
    rfc3389: [u8; 11],
}

fn read_ref_enc() -> Vec<RefFrame> {
    let text = std::fs::read_to_string("fixtures/annex_b_ref_enc.txt").expect("read ref enc");
    let mut frames = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let nums: Vec<i64> = line
            .split_whitespace()
            .map(|t| t.parse().expect("parse int"))
            .collect();
        let len = nums[1] as usize;
        assert!(
            len == 0 || len == 2 || len == 10,
            "unexpected ref len {}",
            len
        );
        let bytes = nums[2..2 + len].iter().map(|&v| v as u8).collect();
        let mut rfc3389 = [0u8; 11];
        for i in 0..11 {
            rfc3389[i] = nums[2 + len + i] as u8;
        }
        frames.push(RefFrame {
            len,
            bytes,
            rfc3389,
        });
    }
    frames
}

#[test]
fn annex_b_encoder_bit_exact() {
    let pcm = read_input_pcm();
    let refs = read_ref_enc();
    let mut encoder = Encoder::new(true);

    let mut n_voice = 0;
    let mut n_sid = 0;
    let mut n_untransmitted = 0;

    for (idx, reference) in refs.iter().enumerate() {
        let start = idx * FRAME_SAMPLES;
        if start + FRAME_SAMPLES > pcm.len() {
            break;
        }
        let mut frame = [0i16; FRAME_SAMPLES];
        frame.copy_from_slice(&pcm[start..start + FRAME_SAMPLES]);

        let out = encoder.encode(&frame);
        assert_eq!(
            out.len(),
            reference.len,
            "frame {}: DTX length mismatch (got {}, want {})",
            idx,
            out.len(),
            reference.len
        );
        assert_eq!(
            out, reference.bytes,
            "frame {}: bitstream mismatch (got {:?}, want {:?})",
            idx, out, reference.bytes
        );

        let rfc = encoder.rfc3389_payload();
        assert_eq!(
            rfc, reference.rfc3389,
            "frame {}: RFC3389 payload mismatch (got {:?}, want {:?})",
            idx, rfc, reference.rfc3389
        );

        match reference.len {
            10 => n_voice += 1,
            2 => n_sid += 1,
            _ => n_untransmitted += 1,
        }
    }

    // Sanity: the fixture actually exercises all three frame types.
    assert!(
        n_voice > 0 && n_sid > 0 && n_untransmitted > 0,
        "fixture did not exercise all DTX frame types: voice={} sid={} untransmitted={}",
        n_voice,
        n_sid,
        n_untransmitted
    );
}

#[test]
fn annex_b_decoder_bit_exact() {
    let refs = read_ref_enc();
    let ref_bytes = std::fs::read("fixtures/annex_b_ref_dec.pcm").expect("read ref dec");
    let ref_pcm: Vec<i16> = ref_bytes
        .as_chunks::<2>()
        .0
        .iter()
        .map(|c| i16::from_le_bytes(*c))
        .collect();

    let mut decoder = Decoder::new();

    for (idx, reference) in refs.iter().enumerate() {
        let (payload, is_sid): (&[u8], bool) = match reference.len {
            10 => (&reference.bytes, false),
            2 => (&reference.bytes, true),
            _ => (&[], true),
        };
        let out = decoder.decode(payload, false, is_sid, false);
        let want = &ref_pcm[idx * FRAME_SAMPLES..(idx + 1) * FRAME_SAMPLES];
        assert_eq!(
            out, want,
            "frame {}: decoded PCM mismatch (len {})",
            idx, reference.len
        );
    }
}
