use g729_sys::g729::decoder::Decoder;
use g729_sys::g729::encoder::Encoder;
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

    let mut input_buffer = [0u8; 160]; // 80 samples * 2 bytes
    let mut pcm_buffer = [0i16; 80];
    let mut bit_stream = [0u8; 10]; // 10 bytes for G.729 frame
    let mut bit_stream_length = 0u8;
    let mut decoded_pcm = [0i16; 80];
    let mut output_buffer = [0u8; 160];

    let mut frame_count = 0;

    loop {
        let bytes_read = input_file.read(&mut input_buffer)?;
        if bytes_read < 160 {
            break;
        }

        // Convert bytes to i16 (Little Endian)
        for i in 0..80 {
            pcm_buffer[i] = i16::from_le_bytes([input_buffer[2 * i], input_buffer[2 * i + 1]]);
        }

        // Encode
        encoder.encode(&pcm_buffer, &mut bit_stream, &mut bit_stream_length);

        // Decode
        // bit_stream is Option<&[u8]>
        // frame_erasure_flag = 0
        // sid_frame_flag = 0 (since VAD is disabled)
        // rfc3389_payload_flag = 0
        decoder.decode(
            Some(&bit_stream),
            bit_stream_length,
            0,
            0,
            0,
            &mut decoded_pcm,
        );

        // Convert i16 to bytes (Little Endian)
        for i in 0..80 {
            let bytes = decoded_pcm[i].to_le_bytes();
            output_buffer[2 * i] = bytes[0];
            output_buffer[2 * i + 1] = bytes[1];
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
