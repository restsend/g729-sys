#[cfg(target_arch = "aarch64")]
use core::arch::aarch64::*;
#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::*;

use crate::g729::basic_operations::*;
use crate::g729::codebooks::MA_PREDICTION_COEFFICIENTS;
use crate::g729::fixed_point_math::*;
use crate::g729::ld8k::*;

pub fn insertion_sort(x: &mut [Word16], length: usize) {
    for i in 1..length {
        let current_value = x[i];
        let mut j = (i as i32) - 1;
        while j >= 0 && x[j as usize] > current_value {
            x[(j + 1) as usize] = x[j as usize];
            j -= 1;
        }
        x[(j + 1) as usize] = current_value;
    }
}

pub fn get_min_in_array(x: &[Word16], length: usize) -> Word16 {
    let mut min = MAX_16;
    for i in 0..length {
        if x[i] < min {
            min = x[i];
        }
    }
    min
}

pub fn compute_parity(mut adaptative_codebook_index: UWord16) -> UWord16 {
    let mut parity = 1;
    adaptative_codebook_index = adaptative_codebook_index >> 2; /* ignore the two LSB */

    for _ in 0..6 {
        parity ^= adaptative_codebook_index & 1; /* XOR with the LSB */
        adaptative_codebook_index = adaptative_codebook_index >> 1;
    }
    parity
}

pub fn rearrange_coefficients(q_lsp: &mut [Word16], j: Word16) {
    /* qLSP in Q2.13 and J in Q0.13(fitting on 4 bits: possible values 10 and 5) */
    for i in 1..NB_LSP_COEFF {
        let delta = (add16(sub16(q_lsp[i - 1], q_lsp[i]), j)) / 2; /* delta = (l[i-1] - l[i] +J)/2 */
        if delta > 0 {
            q_lsp[i - 1] = sub16(q_lsp[i - 1], delta); /* qLSP still in Q2.13 */
            q_lsp[i] = add16(q_lsp[i], delta);
        }
    }
}

pub fn synthesis_filter(
    input_signal: &[Word16],
    filter_coefficients: &[Word16],
    filtered_signal: &mut [Word16],
) {
    // filtered_signal is accessed in ranges [-10, -1] as input and [0, 39] as output.
    // In Rust, we pass the slice starting at index 0 of the output buffer.
    // But we need access to previous 10 elements.
    // So the caller should pass a slice that includes the memory.
    // The C function signature: void synthesisFilter(word16_t inputSignal[], word16_t filterCoefficients[], word16_t filteredSignal[]);
    // In C, filteredSignal points to the start of the output buffer (index 0).
    // But the loop accesses filteredSignal[i-j-1]. When i=0, j=0, it accesses filteredSignal[-1].
    // So in Rust, we must pass the slice such that index 10 corresponds to "0".
    // Or we pass the full buffer and an offset.
    // Let's assume the passed `filtered_signal` slice starts at the memory (index -10 in C terms).
    // So `filtered_signal[10]` corresponds to `filteredSignal[0]` in C.
    // Wait, the C comment says: "filteredSignal: 50 values in Q0 accessed in ranges [-10,-1] as input and [0, 39] as output."
    // So the pointer passed to C function is likely pointing to the 11th element of the allocated array.
    // In Rust, let's change the signature to take the full buffer and the start index?
    // Or just assume the slice passed starts at "0" and we can't access negative indices.
    // If I follow C pointer arithmetic, I should pass the slice starting at "0". But Rust slices don't support negative indexing.
    // So I will change the signature to take the whole buffer and assume the "current" frame starts at index `NB_LSP_COEFF`.
    // But `inputSignal` is just the input.

    // Let's look at how it's called in C.
    // synthesisFilter(&(encoderChannelContext->targetSignal[NB_LSP_COEFF]), ...)
    // So it passes a pointer to the "current" part.

    // In Rust, I'll implement it taking `filtered_signal` as the slice starting at the "memory" (index -NB_LSP_COEFF).
    // So `filtered_signal` has length `NB_LSP_COEFF + L_SUBFRAME`.
    // The output is written to `filtered_signal[NB_LSP_COEFF..]`.

    for i in 0..L_SUBFRAME {
        let mut acc = sshl(input_signal[i] as Word32, 12); /* acc get the first term of the sum, in Q12 (inputSignal is in Q0)*/
        for j in 0..NB_LSP_COEFF {
            // filteredSignal[i-j-1] in C.
            // In Rust slice starting at -10:
            // index 0 corresponds to C index -10.
            // index 10 corresponds to C index 0.
            // C index k maps to Rust index k + 10.
            // C: i-j-1. Rust: i-j-1 + 10 = i - j + 9.
            acc = msu16_16(
                acc,
                filter_coefficients[j],
                filtered_signal[i + NB_LSP_COEFF - j - 1],
            );
        }
        filtered_signal[i + NB_LSP_COEFF] = saturate(pshr(acc, 12), MAX_16 as Word32) as Word16;
    }
}

pub fn dot_product(x: &[Word16], y: &[Word16]) -> Word32 {
    let len = x.len().min(y.len());
    #[allow(unused)]
    let mut sum: Word32 = 0;
    let mut i = 0;

    #[cfg(target_arch = "aarch64")]
    unsafe {
        let mut sum_vec = vdupq_n_s32(0);
        while i + 8 <= len {
            let a = vld1q_s16(x.as_ptr().add(i));
            let b = vld1q_s16(y.as_ptr().add(i));
            let a_low = vget_low_s16(a);
            let b_low = vget_low_s16(b);
            let a_high = vget_high_s16(a);
            let b_high = vget_high_s16(b);

            sum_vec = vmlal_s16(sum_vec, a_low, b_low);
            sum_vec = vmlal_s16(sum_vec, a_high, b_high);

            i += 8;
        }
        sum = vaddvq_s32(sum_vec);
    }

    #[cfg(target_arch = "x86_64")]
    unsafe {
        let mut sum_vec = _mm_setzero_si128();
        while i + 8 <= len {
            let a = _mm_loadu_si128(x.as_ptr().add(i) as *const _);
            let b = _mm_loadu_si128(y.as_ptr().add(i) as *const _);
            let prod = _mm_madd_epi16(a, b);
            sum_vec = _mm_add_epi32(sum_vec, prod);
            i += 8;
        }
        let high = _mm_unpackhi_epi64(sum_vec, sum_vec);
        let sum_vec = _mm_add_epi32(sum_vec, high);
        let high = _mm_shuffle_epi32(sum_vec, 1);
        let sum_vec = _mm_add_epi32(sum_vec, high);
        sum = _mm_cvtsi128_si32(sum_vec);
    }

    while i < len {
        sum = mac16_16(sum, x[i], y[i]);
        i += 1;
    }
    sum
}

pub fn dot_product_16_32_q12(x: &[Word16], y: &[Word32]) -> Word32 {
    let len = x.len().min(y.len());
    #[allow(unused)]
    let mut sum: Word32 = 0;
    let mut i = 0;

    #[cfg(target_arch = "aarch64")]
    unsafe {
        let mut sum_vec = vdupq_n_s32(0);
        while i + 4 <= len {
            let a = vld1_s16(x.as_ptr().add(i)); // Load 4 i16
            let b = vld1q_s32(y.as_ptr().add(i)); // Load 4 i32

            let a_32 = vmovl_s16(a); // Extend to 4 i32

            // Multiply -> 4 i64 (low and high)
            let prod_low = vmull_s32(vget_low_s32(a_32), vget_low_s32(b));
            let prod_high = vmull_high_s32(a_32, b);

            // Shift right by 12
            let prod_low_shifted = vshrq_n_s64(prod_low, 12);
            let prod_high_shifted = vshrq_n_s64(prod_high, 12);

            // Narrow back to i32
            let res_low = vmovn_s64(prod_low_shifted); // 2 i32
            let res_high = vmovn_s64(prod_high_shifted); // 2 i32

            // Combine
            let res = vcombine_s32(res_low, res_high);

            // Accumulate
            sum_vec = vaddq_s32(sum_vec, res);

            i += 4;
        }
        sum = vaddvq_s32(sum_vec);
    }

    #[cfg(target_arch = "x86_64")]
    unsafe {
        let mut sum_vec = _mm_setzero_si128();
        while i + 4 <= len {
            let a_16 = _mm_loadl_epi64(x.as_ptr().add(i) as *const _); // Load 4 i16 (64 bits) into low part of XMM
                                                                       // Sign extend to 32-bit using SSE2 trick
            let a_32 = _mm_srai_epi32(_mm_unpacklo_epi16(a_16, a_16), 16);

            let b = _mm_loadu_si128(y.as_ptr().add(i) as *const _); // Load 4 i32

            // Low 2 elements
            let prod_even = _mm_mul_epi32(a_32, b);

            // High 2 elements (shuffle to move 1,3 to 0,2)
            let a_high = _mm_shuffle_epi32(a_32, 0xF5); // 3,3,1,1 -> 11 11 01 01
            let b_high = _mm_shuffle_epi32(b, 0xF5);
            let prod_odd = _mm_mul_epi32(a_high, b_high);

            // Shift right by 12
            let prod_even_shifted = _mm_srai_epi64(prod_even, 12);
            let prod_odd_shifted = _mm_srai_epi64(prod_odd, 12);

            // Shuffle to get low 32 bits of each 64 bit
            let res_even = _mm_shuffle_epi32(prod_even_shifted, 0x08); // 00 00 10 00 -> 2,0
            let res_odd = _mm_shuffle_epi32(prod_odd_shifted, 0x08);

            // Unpack
            let res = _mm_unpacklo_epi32(res_even, res_odd); // 0e, 0o, 2e, 2o -> 0, 1, 2, 3

            sum_vec = _mm_add_epi32(sum_vec, res);

            i += 4;
        }
        // Horizontal sum
        let high = _mm_unpackhi_epi64(sum_vec, sum_vec);
        let sum_vec = _mm_add_epi32(sum_vec, high);
        let high = _mm_shuffle_epi32(sum_vec, 1);
        let sum_vec = _mm_add_epi32(sum_vec, high);
        sum = _mm_cvtsi128_si32(sum_vec);
    }

    while i < len {
        sum = mac16_32_q12(sum, x[i], y[i]);
        i += 1;
    }
    sum
}

pub fn vec_mult_16_16(x: &[Word16], y: &[Word16], out: &mut [Word32]) {
    let len = x.len().min(y.len()).min(out.len());
    let mut i = 0;

    #[cfg(target_arch = "aarch64")]
    unsafe {
        while i + 8 <= len {
            let a = vld1q_s16(x.as_ptr().add(i));
            let b = vld1q_s16(y.as_ptr().add(i));

            let prod_low = vmull_s16(vget_low_s16(a), vget_low_s16(b));
            let prod_high = vmull_high_s16(a, b);

            vst1q_s32(out.as_mut_ptr().add(i), prod_low);
            vst1q_s32(out.as_mut_ptr().add(i + 4), prod_high);

            i += 8;
        }
    }

    #[cfg(target_arch = "x86_64")]
    unsafe {
        while i + 8 <= len {
            let a = _mm_loadu_si128(x.as_ptr().add(i) as *const _);
            let b = _mm_loadu_si128(y.as_ptr().add(i) as *const _);

            let lo = _mm_mullo_epi16(a, b);
            let hi = _mm_mulhi_epi16(a, b);

            let res_lo = _mm_unpacklo_epi16(lo, hi);
            let res_hi = _mm_unpackhi_epi16(lo, hi);

            _mm_storeu_si128(out.as_mut_ptr().add(i) as *mut _, res_lo);
            _mm_storeu_si128(out.as_mut_ptr().add(i + 4) as *mut _, res_hi);

            i += 8;
        }
    }

    while i < len {
        out[i] = mult16_16(x[i], y[i]);
        i += 1;
    }
}

pub fn correlate_vectors(x: &[Word16], y: &[Word16], c: &mut [Word32]) {
    for i in 0..L_SUBFRAME {
        c[i] = dot_product(&x[i..L_SUBFRAME], &y[0..L_SUBFRAME - i]);
    }
}

#[inline]
pub fn count_leading_zeros(x: Word32) -> UWord16 {
    if x == 0 {
        31
    } else {
        (x as u32).leading_zeros() as UWord16 - 1
    }
}

#[inline]
pub fn unsigned_count_leading_zeros(x: UWord32) -> UWord16 {
    x.leading_zeros() as UWord16
}

pub fn ma_code_gain_prediction(
    previous_gain_prediction_error: &[Word16],
    fixed_codebook_vector: &[Word16],
) -> Word32 {
    /* compute the sum of squares of fixedCodebookVector in Q26 */
    let mut fixed_codebook_vector_squares_sum: Word32 = 0;

    for i in 0..L_SUBFRAME {
        if fixed_codebook_vector[i] != 0 {
            fixed_codebook_vector_squares_sum = mac16_16(
                fixed_codebook_vector_squares_sum,
                fixed_codebook_vector[i],
                fixed_codebook_vector[i],
            );
        }
    }

    /* compute E| - E as in eq71, result in Q16 */
    let mut acc = mac16_32_q13(
        8145364,
        -24660,
        g729_log2_q0q16(fixed_codebook_vector_squares_sum),
    ); /* acc in Q16 */

    /* accumulate the MA prediction described in eq69 to the previous Sum, result will be in E~(m) + E| -E as used in eq71 */
    acc = shl(acc, 8); /* acc in Q24 to match the fixed point of next accumulations */
    for i in 0..4 {
        acc = mac16_16(
            acc,
            previous_gain_prediction_error[i],
            MA_PREDICTION_COEFFICIENTS[i],
        );
    }

    /* compute eq71, we already have the exposant in acc so */
    /* g'c = 10^(acc/20)                                    */
    /*     = 2^((acc*ln(10))/(20*ln(2)))                    */
    /*     = 2^(0,1661*acc)                                 */
    acc = shr(acc, 2); /* Q24->Q22 */
    acc = mult16_32_q15(5442, acc); /* 5442 is 0.1661 in Q15 -> acc now in Q4.22 range [1.6, 10.8] */
    acc = pshr(acc, 11); /* get acc in Q4.11 */

    g729_exp2_q11q16(acc as Word16)
}

pub fn compute_gain_prediction_error(
    fixed_codebook_gain_correction_factor: Word16,
    previous_gain_prediction_error: &mut [Word16],
) {
    /* need to compute eq72: 20log10(fixedCodebookGainCorrectionFactor) */
    /*  = (20/log2(10))*log2(fixedCodebookGainCorrectionFactor) */
    /*  = 6.0206*log2(fixedCodebookGainCorrectionFactor) */
    /* log2 input in Q0, output in Q16,fixedCodebookGainCorrectionFactor being in Q12, we shall substract 12 in Q16(786432) to the result of log2 function -> final result in Q2.16 */
    let mut current_gain_prediction_error = sub32(
        g729_log2_q0q16(fixed_codebook_gain_correction_factor as Word32),
        786432,
    );
    current_gain_prediction_error = pshr(mult16_32_q12(24660, current_gain_prediction_error), 6); /* 24660 = 6.0206 in Q3.12 -> mult result in Q16, precise shift right to get it in Q4.10 */

    /* shift the array and insert the current Prediction Error */
    previous_gain_prediction_error[3] = previous_gain_prediction_error[2];
    previous_gain_prediction_error[2] = previous_gain_prediction_error[1];
    previous_gain_prediction_error[1] = previous_gain_prediction_error[0];
    previous_gain_prediction_error[0] = current_gain_prediction_error as Word16;
}

pub fn parameters_array_2_bit_stream(parameters: &[UWord16], bit_stream: &mut [u8]) {
    bit_stream[0] = (((parameters[0] & 0x1) << 7) | (parameters[1] & 0x7f)) as u8;

    bit_stream[1] = (((parameters[2] & 0x1f) << 3) | ((parameters[3] >> 2) & 0x7)) as u8;

    bit_stream[2] = (((parameters[3] & 0x3) << 6) | ((parameters[4] >> 2) & 0x3f)) as u8;

    bit_stream[3] = (((parameters[4] & 0x3) << 6)
        | ((parameters[5] & 0x1) << 5)
        | ((parameters[6] >> 8) & 0x1f)) as u8;

    bit_stream[4] = (parameters[6] & 0xff) as u8;

    bit_stream[5] = (((parameters[7] & 0xf) << 4)
        | ((parameters[8] & 0x7) << 1)
        | ((parameters[9] >> 3) & 0x1)) as u8;

    bit_stream[6] = (((parameters[9] & 0x7) << 5) | (parameters[10] & 0x1f)) as u8;

    bit_stream[7] = ((parameters[11] >> 5) & 0xff) as u8;

    bit_stream[8] = (((parameters[11] & 0x1f) << 3) | ((parameters[12] >> 1) & 0x7)) as u8;

    bit_stream[9] = (((parameters[12] & 0x1) << 7)
        | ((parameters[13] & 0x7) << 4)
        | (parameters[14] & 0xf)) as u8;
}

pub fn cng_parameters_array_2_bit_stream(parameters: &[UWord16], bit_stream: &mut [u8]) {
    bit_stream[0] = (((parameters[0] & 0x1) << 7)
        | ((parameters[1] & 0x1f) << 2)
        | ((parameters[2] >> 2) & 0x3)) as u8;

    bit_stream[1] = (((parameters[2] & 0x03) << 6) | ((parameters[3] & 0x1f) << 1)) as u8;
}

pub fn parameters_bit_stream_2_array(bit_stream: &[u8], parameters: &mut [UWord16]) {
    parameters[0] = ((bit_stream[0] >> 7) & 0x1) as UWord16;
    parameters[1] = (bit_stream[0] & 0x7f) as UWord16;
    parameters[2] = ((bit_stream[1] >> 3) & 0x1f) as UWord16;
    parameters[3] =
        (((bit_stream[1] & 0x7) as UWord16) << 2) | ((bit_stream[2] >> 6) & 0x3) as UWord16;
    parameters[4] =
        (((bit_stream[2] & 0x3f) as UWord16) << 2) | ((bit_stream[3] >> 6) & 0x3) as UWord16;
    parameters[5] = ((bit_stream[3] >> 5) & 0x1) as UWord16;
    parameters[6] = (((bit_stream[3] & 0x1f) as UWord16) << 8) | bit_stream[4] as UWord16;
    parameters[7] = ((bit_stream[5] >> 4) & 0xf) as UWord16;
    parameters[8] = ((bit_stream[5] >> 1) & 0x7) as UWord16;
    parameters[9] =
        (((bit_stream[5] & 0x1) as UWord16) << 3) | ((bit_stream[6] >> 5) & 0x7) as UWord16;
    parameters[10] = (bit_stream[6] & 0x1f) as UWord16;
    parameters[11] = ((bit_stream[7] as UWord16) << 5) | ((bit_stream[8] >> 3) & 0x1f) as UWord16;
    parameters[12] =
        (((bit_stream[8] & 0x7) as UWord16) << 1) | ((bit_stream[9] >> 7) & 0x1) as UWord16;
    parameters[13] = ((bit_stream[9] >> 4) & 0x7) as UWord16;
    parameters[14] = (bit_stream[9] & 0xf) as UWord16;
}

pub fn pseudo_random(random_generator_seed: &mut UWord16) -> UWord16 {
    /* pseudoRandomSeed is stored in an uint16_t var, we shall not worry about overflow here */
    /* pseudoRandomSeed*31821 + 13849; */
    // MAC16_16(13849, (*randomGeneratorSeed), 31821);
    // MAC16_16 returns Word32.
    // In C: *randomGeneratorSeed = MAC16_16(...)
    // It implicitly casts Word32 to uint16_t (truncates).
    let res = mac16_16(13849, *random_generator_seed as Word16, 31821);
    *random_generator_seed = res as UWord16;
    *random_generator_seed
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dot_product() {
        let x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let y = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        // Sum = 55
        assert_eq!(dot_product(&x, &y), 55);

        let x = [MAX_16, 1];
        let y = [1, 1];
        assert_eq!(dot_product(&x, &y), MAX_16 as Word32 + 1);

        // Test negative
        let x = [-1, -2];
        let y = [1, 1];
        assert_eq!(dot_product(&x, &y), -3);
    }
}
