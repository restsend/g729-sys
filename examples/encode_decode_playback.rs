use g729_sys::{Decoder, Encoder, FRAME_SAMPLES, VOICE_FRAME_BYTES};
use std::fs::File;
use std::io::{Read, Write};

fn main() -> std::io::Result<()> {
    let input_path = "fixtures/en_8k_16bit.pcm";
    let output_path = "fixtures/output.pcm";

    println!("Opening input file: {}", input_path);
    let mut input_file = File::open(input_path)?;
    let mut output_file = File::create(output_path)?;

    let mut encoder = Encoder::new(false); // VAD disabled
    let mut decoder = Decoder::new();

    let mut input_buffer = [0u8; FRAME_SAMPLES * 2]; // 80 samples * 2 bytes
    let mut bit_stream = [0u8; VOICE_FRAME_BYTES]; // 10 bytes for a G.729 frame
    let mut output_buffer = [0u8; FRAME_SAMPLES * 2];

    let mut frame_count = 0;

    loop {
        let bytes_read = input_file.read(&mut input_buffer)?;
        if bytes_read < input_buffer.len() {
            break;
        }

        // Convert bytes to i16 (Little Endian).
        let mut pcm_buffer = [0i16; FRAME_SAMPLES];
        for (sample, chunk) in pcm_buffer.iter_mut().zip(input_buffer.as_chunks::<2>().0) {
            *sample = i16::from_le_bytes(*chunk);
        }

        // Encode, then decode the produced bitstream.
        let len = encoder.encode_into(&pcm_buffer, &mut bit_stream) as usize;
        let decoded_pcm = decoder.decode(&bit_stream[..len], false, false, false);

        // Convert i16 to bytes (Little Endian).
        for (bytes, sample) in output_buffer
            .as_chunks_mut::<2>()
            .0
            .iter_mut()
            .zip(decoded_pcm.iter())
        {
            bytes.copy_from_slice(&sample.to_le_bytes());
        }

        output_file.write_all(&output_buffer)?;
        frame_count += 1;
    }

    println!("Processed {} frames.", frame_count);
    println!("Output saved to: {}", output_path);
    println!("To play the output, run:");
    println!("ffplay -f s16le -ar 8000 -ac 1 {}", output_path);

    Ok(())
}
