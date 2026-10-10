use criterion::{black_box, criterion_group, criterion_main, Criterion};
use g729_sys::{Decoder, Encoder, FRAME_SAMPLES, VOICE_FRAME_BYTES};

fn benchmark_encoder(c: &mut Criterion) {
    let mut group = c.benchmark_group("Encoder");
    let input = [0i16; FRAME_SAMPLES];

    group.bench_function("Rust Encoder", |b| {
        let mut encoder = Encoder::new(false);
        b.iter(|| {
            let mut out = [0u8; VOICE_FRAME_BYTES];
            encoder.encode_into(black_box(&input), &mut out);
        })
    });
    group.finish();
}

fn benchmark_decoder(c: &mut Criterion) {
    let mut group = c.benchmark_group("Decoder");
    // 10 bytes of silence payload (approximate)
    let payload = [0u8; 10];

    group.bench_function("Rust Decoder", |b| {
        let mut decoder = Decoder::new();
        b.iter(|| {
            black_box(decoder.decode(black_box(&payload), false, false, false));
        })
    });
    group.finish();
}

criterion_group!(benches, benchmark_encoder, benchmark_decoder);
criterion_main!(benches);
