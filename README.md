## g729-sys

Rust implementation of the G.729 codec.

### Features
- Pure Rust implementation
- `no_std` compatible (disable the default `std` feature; the codec and Annex B use only `core`)
- SIMD optimized (`aarch64` / `x86_64` for the hot correlation kernels)
- G.729 Annex B: VAD (voice activity detection), DTX (SID / untransmitted frames) and CNG (comfort noise generation), including RFC3389 comfort noise payloads

### `no_std`

```toml
[dependencies]
g729-sys = { version = "0.2", default-features = false }
```

The `std` feature is enabled by default and only gates the convenience
`Encoder::encode` / `Vec`-returning API. The encoder, decoder and all Annex B
state (VAD/DTX/CNG) are heap-free and work without `std`.

### Performance

On-device numbers from the ESP32-S3 (QFN56, dual-core 240 MHz, **no FPU**)
`wifi_ua` firmware baseline in [`../rtcembed`](../rtcembed)
(`scripts/baseline/BASELINE.md`), `opt-level = 3`, per 20 ms frame:

| path | µs / 20 ms frame | % of one 240 MHz core |
|---|---|---|
| G.729 encode (DC / noise-like signal) | 4425 / 4476 | ≈ 22 % |
| G.729 decode | 848 / 868 | ≈ 3.6 % |
| full-duplex call (encode + decode) | — | ≈ 26.4 % |

Decode is ~5× cheaper than encode; a full-duplex G.729 call is a **real-time
factor of ≈ 0.34**. The board has no FPU, so the integer CELP/QMF math
dominates.

Bit-exact encode micro-optimisations (20 000 frames, byte-identical bitstream),
on-device:

| build | DC | noise-like signal |
|---|---|---|
| original | 5867–6047 µs | 5748–5948 µs |
| + Φ-fold | 5169–5184 µs | 5486–5505 µs |
| + synthesis unroll / index hoist | 4995–5036 µs | 5185–5265 µs |
| + L1 clamp/zip | **4924–4935 µs** | **5089–5106 µs** |

≈ −17 % (DC) / −13 % (signal); decode unchanged. A per-package
`opt-level = 3` override (vs the firmware's global `"z"`) additionally took
encode 6753 → 5986 µs (−11 %) and decode 1249 → 867 µs (−31 %).

### Testing

The encoder bitstream (including SID/DTX frame lengths and the RFC3389 payload)
and the decoder output are verified **bit-exact against the bcg729 reference**
on the regression corpus (`tests/annex_b_test.rs`); `tests/annex_b_soak_test.rs`
additionally exercises long VAD/DTX/CNG sessions.

### License
MIT
