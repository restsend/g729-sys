#[cfg(target_arch = "aarch64")]
use core::arch::aarch64::*;
#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::*;

use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;

pub fn insertion_sort(x: &mut [i16]) {
    let length = x.len();
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

pub fn get_min_in_array(x: &[i16]) -> i16 {
    x.iter().copied().min().unwrap_or(MAX_16)
}

pub fn compute_parity(mut adaptative_codebook_index: u16) -> u16 {
    let mut parity = 1;
    adaptative_codebook_index >>= 2;

    for _ in 0..6 {
        parity ^= adaptative_codebook_index & 1;
        adaptative_codebook_index >>= 1;
    }
    parity
}

pub fn rearrange_coefficients(q_lsp: &mut [i16], j: i16) {
    /* qLSP in Q2.13 and J in Q0.13(fitting on 4 bits: possible values 10 and 5) */
    for i in 1..NB_LSP_COEFF {
        let delta = (add16(sub16(q_lsp[i - 1], q_lsp[i]), j)) / 2;
        if delta > 0 {
            q_lsp[i - 1] = sub16(q_lsp[i - 1], delta); /* qLSP still in Q2.13 */
            q_lsp[i] = add16(q_lsp[i], delta);
        }
    }
}

/// Vectorized prefix of [`dot_product`]. Returns the number of leading elements
/// handled by SIMD; `*sum` holds their (partial) result.
#[cfg(target_arch = "aarch64")]
fn dot_product_simd(x: &[i16], y: &[i16], len: usize, sum: &mut i32) -> usize {
    let mut i = 0;
    // SAFETY: the loop only runs while `i + 8 <= len` and `len <= x.len()/y.len()`,

    unsafe {
        let mut sum_vec = vdupq_n_s32(0);
        while i + 8 <= len {
            let a = vld1q_s16(x.as_ptr().add(i));
            let b = vld1q_s16(y.as_ptr().add(i));
            sum_vec = vmlal_s16(sum_vec, vget_low_s16(a), vget_low_s16(b));
            sum_vec = vmlal_s16(sum_vec, vget_high_s16(a), vget_high_s16(b));
            i += 8;
        }
        *sum = vaddvq_s32(sum_vec);
    }
    i
}

#[cfg(target_arch = "x86_64")]
fn dot_product_simd(x: &[i16], y: &[i16], len: usize, sum: &mut i32) -> usize {
    let mut i = 0;
    // SAFETY: the loop only runs while `i + 8 <= len` and `len <= x.len()/y.len()`,

    unsafe {
        let mut sum_vec = _mm_setzero_si128();
        while i + 8 <= len {
            let a = _mm_loadu_si128(x.as_ptr().add(i) as *const _);
            let b = _mm_loadu_si128(y.as_ptr().add(i) as *const _);
            sum_vec = _mm_add_epi32(sum_vec, _mm_madd_epi16(a, b));
            i += 8;
        }
        let high = _mm_unpackhi_epi64(sum_vec, sum_vec);
        let sum_vec = _mm_add_epi32(sum_vec, high);
        let high = _mm_shuffle_epi32(sum_vec, 1);
        let sum_vec = _mm_add_epi32(sum_vec, high);
        *sum = _mm_cvtsi128_si32(sum_vec);
    }
    i
}

pub fn dot_product(x: &[i16], y: &[i16]) -> i32 {
    let len = x.len().min(y.len());
    let mut sum: i32 = 0;

    #[cfg(target_arch = "aarch64")]
    let mut i = dot_product_simd(x, y, len, &mut sum);
    #[cfg(target_arch = "x86_64")]
    let mut i = dot_product_simd(x, y, len, &mut sum);
    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    let mut i = 0;

    while i < len {
        sum = mac16_16(sum, x[i], y[i]);
        i += 1;
    }
    sum
}

/// Vectorized prefix of [`dot_product_16_32_q12`].
#[cfg(target_arch = "aarch64")]
fn dot_product_16_32_q12_simd(x: &[i16], y: &[i32], len: usize, sum: &mut i32) -> usize {
    let mut i = 0;
    // SAFETY: the loop only runs while `i + 4 <= len` and `len <= x.len()/y.len()`,

    unsafe {
        let mut sum_vec = vdupq_n_s32(0);
        while i + 4 <= len {
            let a = vmovl_s16(vld1_s16(x.as_ptr().add(i)));
            let b = vld1q_s32(y.as_ptr().add(i));
            let prod_low = vmull_s32(vget_low_s32(a), vget_low_s32(b));
            let prod_high = vmull_high_s32(a, b);
            let res_low = vmovn_s64(vshrq_n_s64(prod_low, 12));
            let res_high = vmovn_s64(vshrq_n_s64(prod_high, 12));
            sum_vec = vaddq_s32(sum_vec, vcombine_s32(res_low, res_high));
            i += 4;
        }
        *sum = vaddvq_s32(sum_vec);
    }
    i
}

pub fn dot_product_16_32_q12(x: &[i16], y: &[i32]) -> i32 {
    let len = x.len().min(y.len());
    let mut sum: i32 = 0;

    #[cfg(target_arch = "aarch64")]
    let mut i = dot_product_16_32_q12_simd(x, y, len, &mut sum);
    #[cfg(not(target_arch = "aarch64"))]
    let mut i = 0;

    while i < len {
        sum = mac16_32_q12(sum, x[i], y[i]);
        i += 1;
    }
    sum
}

/// Vectorized prefix of [`vec_mult_16_16`]; writes the products into `out`.
#[cfg(target_arch = "aarch64")]
fn vec_mult_16_16_simd(x: &[i16], y: &[i16], out: &mut [i32], len: usize) -> usize {
    let mut i = 0;
    // SAFETY: the loop only runs while `i + 8 <= len` and `len` is the minimum of

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
    i
}

#[cfg(target_arch = "x86_64")]
fn vec_mult_16_16_simd(x: &[i16], y: &[i16], out: &mut [i32], len: usize) -> usize {
    let mut i = 0;
    // SAFETY: the loop only runs while `i + 8 <= len` and `len` is the minimum of

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
    i
}

pub fn vec_mult_16_16(x: &[i16], y: &[i16], out: &mut [i32]) {
    let len = x.len().min(y.len()).min(out.len());

    #[cfg(target_arch = "aarch64")]
    let mut i = vec_mult_16_16_simd(x, y, out, len);
    #[cfg(target_arch = "x86_64")]
    let mut i = vec_mult_16_16_simd(x, y, out, len);
    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    let mut i = 0;

    while i < len {
        out[i] = mult16_16(x[i], y[i]);
        i += 1;
    }
}

pub fn correlate_vectors(x: &[i16], y: &[i16], c: &mut [i32]) {
    for i in 0..L_SUBFRAME {
        c[i] = dot_product(&x[i..L_SUBFRAME], &y[0..L_SUBFRAME - i]);
    }
}

#[inline]
pub fn count_leading_zeros(x: i32) -> u16 {
    if x == 0 {
        31
    } else {
        (x as u32).leading_zeros() as u16 - 1
    }
}

#[inline]
pub fn unsigned_count_leading_zeros(x: u32) -> u16 {
    x.leading_zeros() as u16
}

pub fn parameters_array_2_bit_stream(parameters: &[u16], bit_stream: &mut [u8]) {
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

pub fn parameters_bit_stream_2_array(bit_stream: &[u8], parameters: &mut [u16]) {
    // A well-formed G.729 voice frame is 10 bytes, but callers may hand us a

    #[inline]
    fn byte(bit_stream: &[u8], i: usize) -> u8 {
        bit_stream.get(i).copied().unwrap_or(0)
    }

    parameters[0] = ((byte(bit_stream, 0) >> 7) & 0x1) as u16;
    parameters[1] = (byte(bit_stream, 0) & 0x7f) as u16;
    parameters[2] = ((byte(bit_stream, 1) >> 3) & 0x1f) as u16;
    parameters[3] =
        (((byte(bit_stream, 1) & 0x7) as u16) << 2) | ((byte(bit_stream, 2) >> 6) & 0x3) as u16;
    parameters[4] =
        (((byte(bit_stream, 2) & 0x3f) as u16) << 2) | ((byte(bit_stream, 3) >> 6) & 0x3) as u16;
    parameters[5] = ((byte(bit_stream, 3) >> 5) & 0x1) as u16;
    parameters[6] = (((byte(bit_stream, 3) & 0x1f) as u16) << 8) | byte(bit_stream, 4) as u16;
    parameters[7] = ((byte(bit_stream, 5) >> 4) & 0xf) as u16;
    parameters[8] = ((byte(bit_stream, 5) >> 1) & 0x7) as u16;
    parameters[9] =
        (((byte(bit_stream, 5) & 0x1) as u16) << 3) | ((byte(bit_stream, 6) >> 5) & 0x7) as u16;
    parameters[10] = (byte(bit_stream, 6) & 0x1f) as u16;
    parameters[11] =
        ((byte(bit_stream, 7) as u16) << 5) | ((byte(bit_stream, 8) >> 3) & 0x1f) as u16;
    parameters[12] =
        (((byte(bit_stream, 8) & 0x7) as u16) << 1) | ((byte(bit_stream, 9) >> 7) & 0x1) as u16;
    parameters[13] = ((byte(bit_stream, 9) >> 4) & 0x7) as u16;
    parameters[14] = (byte(bit_stream, 9) & 0xf) as u16;
}

pub fn pseudo_random(random_generator_seed: &mut u16) -> u16 {
    let res = mac16_16(13849, *random_generator_seed as i16, 31821);
    *random_generator_seed = res as u16;
    *random_generator_seed
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dot_product() {
        let x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let y = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

        assert_eq!(dot_product(&x, &y), 55);

        let x = [MAX_16, 1];
        let y = [1, 1];
        assert_eq!(dot_product(&x, &y), MAX_16 as i32 + 1);

        let x = [-1, -2];
        let y = [1, 1];
        assert_eq!(dot_product(&x, &y), -3);
    }
}
