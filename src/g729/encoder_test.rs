#[cfg(test)]
mod tests {
    use crate::g729::encoder::EncoderChannelContext;
    use crate::g729::ld8k::L_FRAME;

    #[test]
    fn test_encoder_init() {
        let mut encoder = EncoderChannelContext::new(false);
        let input_frame = [0i16; L_FRAME];
        let mut bit_stream = [0u8; 10]; // 80 bits = 10 bytes
        let mut bit_stream_length = 0;

        encoder.encode(&input_frame, &mut bit_stream, &mut bit_stream_length);

        assert_eq!(bit_stream_length, 10);
    }

    #[test]
    fn test_encoder_sine_wave() {
        let mut encoder = EncoderChannelContext::new(false);
        let mut bit_stream = [0u8; 10];
        let mut bit_stream_length = 0;

        // Generate a sine wave
        let mut input_frame = [0i16; L_FRAME];
        for i in 0..L_FRAME {
            let angle = i as f32 * 2.0 * 3.14159 * 440.0 / 8000.0;
            input_frame[i] = (angle.sin() * 10000.0) as i16;
        }

        // Encode multiple frames
        for _ in 0..10 {
            encoder.encode(&input_frame, &mut bit_stream, &mut bit_stream_length);
            assert_eq!(bit_stream_length, 10);
        }

        let sum: u32 = bit_stream.iter().map(|&x| x as u32).sum();
        assert!(
            sum > 0,
            "Bitstream should not be all zeros for sine wave input"
        );
    }
}
