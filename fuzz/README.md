# g729-sys fuzz targets

`cargo-fuzz` targets for the G.729 encoder/decoder. They must never
panic/abort on arbitrary input — malformed or truncated RTP payloads are
expected in the field.

## Targets

| target | what it feeds |
|---|---|
| `fuzz_decoder` | arbitrary bytes as the RTP payload (`data[0]` = erased/sid/rfc3389 flags), incl. 0-byte, 2-byte SID, 10-byte voice, truncated and oversized packets |
| `fuzz_encoder` | arbitrary 16-bit PCM as one frame, then a few more frames |
| `fuzz_roundtrip` | encode a frame, decode it, then decode mutated bytes |

## Run

```bash
# one-off (installs cargo-fuzz if missing)
./fuzz.sh fuzz_decoder 60

# or directly
cargo +nightly fuzz run fuzz_decoder -- -max_total_time=60
cargo +nightly fuzz run fuzz_encoder  -- -max_total_time=60
cargo +nightly fuzz run fuzz_roundtrip -- -max_total_time=60
```

Requires a nightly toolchain (`cargo-fuzz` builds with the sanitizer runtime).

## Deterministic smoke harness

`fuzz_decoder` exercises the same surface as the `decode_short_payloads_do_not_panic`
regression test in `tests/`, which runs under a plain `cargo test`.
