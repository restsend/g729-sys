use criterion::{black_box, criterion_group, criterion_main, Criterion};
use g729_sys::{g729, FRAME_SAMPLES};

fn benchmark_encoder(c: &mut Criterion) {
    let mut group = c.benchmark_group("Encoder");
    let input = [0i16; FRAME_SAMPLES];

    group.bench_function("Rust Encoder", |b| {
        let mut encoder = g729::encoder::Encoder::new(false);
        b.iter(|| {
            let mut out = [0u8; 10];
            let mut len = 0;
            encoder.encode(black_box(&input), &mut out, &mut len);
        })
    });
    group.finish();
}

fn benchmark_decoder(c: &mut Criterion) {
    let mut group = c.benchmark_group("Decoder");
    // 10 bytes of silence payload (approximate)
    let payload = [0u8; 10];

    group.bench_function("Rust Decoder", |b| {
        let mut decoder = g729::decoder::Decoder::new();
        b.iter(|| {
            let mut out = [0i16; FRAME_SAMPLES];
            decoder.decode(Some(black_box(&payload)), 10, 0, 0, 0, &mut out);
        })
    });
    group.finish();
}

criterion_group!(benches, benchmark_encoder, benchmark_decoder);
criterion_main!(benches);
