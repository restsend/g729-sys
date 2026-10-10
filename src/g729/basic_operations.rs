// Basic operations and macros ported from basicOperationsMacros.h and fixedPointMacros.h

pub const MAX_16: i16 = 0x7fff;
#[allow(dead_code)] // used by the fixed-point unit tests
pub const MIN_16: i16 = -0x8000;
pub const MAX_U16: u16 = 0xffff;
pub const MAX_32: i32 = 0x7fffffff;
pub const MIN_32: i32 = -0x80000000;

#[inline]
pub fn extend32(x: i16) -> i32 {
    x as i32
}

#[inline]
pub fn neg16(x: i16) -> i16 {
    x.wrapping_neg()
}

#[inline]
pub fn neg32(x: i32) -> i32 {
    x.wrapping_neg()
}

/*** shifts ***/
#[inline]
pub fn shr(a: i32, shift: u32) -> i32 {
    a >> shift
}

#[inline]
pub fn shl(a: i32, shift: u32) -> i32 {
    a << shift
}

#[inline]
pub fn sshl(a: i32, shift: u32) -> i32 {
    // Rust's left shift on signed integers is defined as arithmetic shift (preserves sign bit behavior in 2's complement)
    // which matches the intent of the C macro workaround for UB.
    a << shift
}

#[inline]
pub fn ushl(a: u32, shift: u32) -> u32 {
    a << shift
}

/* shift right with rounding: used to extract the integer value of a Qa number */
#[inline]
pub fn pshr(a: i32, shift: u32) -> i32 {
    if shift == 0 {
        return a;
    }
    let half = 1 << (shift - 1);
    (a.wrapping_add(half)) >> shift
}

/* shift right with checking on sign of shift value */
#[inline]
pub fn vshr32(a: i32, shift: i32) -> i32 {
    if shift > 0 {
        shr32(a, shift as u32)
    } else {
        shl32(a, (-shift) as u32)
    }
}

#[inline]
pub fn svshr32(a: i32, shift: i32) -> i32 {
    if shift > 0 {
        shr32(a, shift as u32)
    } else {
        sshl(a, (-shift) as u32)
    }
}

#[inline]
pub fn shr16(a: i16, shift: u32) -> i16 {
    a >> shift
}

#[inline]
pub fn shl16(a: i16, shift: u32) -> i16 {
    a << shift
}

#[inline]
pub fn shr32(a: i32, shift: u32) -> i32 {
    a >> shift
}

#[inline]
pub fn shl32(a: i32, shift: u32) -> i32 {
    a << shift
}

#[inline]
pub fn shr64(a: i64, shift: u32) -> i64 {
    a >> shift
}

#[inline]
pub fn shl64(a: i64, shift: u32) -> i64 {
    a << shift
}

#[inline]
pub fn sshl64(a: i64, shift: u32) -> i64 {
    a << shift
}

/* avoid overflows: a+1 is used to check on negative value because range of a 2n signed bits int is -2pow(n) - 2pow(n)-1 */
/* SATURATE Macro shall be called with MAXINT(nbits). Ex: SATURATE(x,MAXINT16) with MAXINT16  defined to 2pow(16) - 1 */
#[inline]
pub fn saturate(x: i32, a: i32) -> i32 {
    // (((x)>(a) ? (a) : (x)<-(a+1) ? -(a+1) : (x)))
    if x > a {
        a
    } else if x < -(a.wrapping_add(1)) {
        -(a.wrapping_add(1))
    } else {
        x
    }
}

#[inline]
pub fn usaturate(x: i32, a: i32) -> i32 {
    if x > a {
        a
    } else {
        x
    }
}

/* absolute value */
#[inline]
pub fn abs(a: i32) -> i32 {
    a.abs()
}

#[inline]
pub fn abs16(a: i16) -> i16 {
    a.abs()
}

/*** add and sub ***/
#[inline]
pub fn add16(a: i16, b: i16) -> i16 {
    (a as i32 + b as i32) as i16
}

#[inline]
pub fn sub16(a: i16, b: i16) -> i16 {
    (a as i32 - b as i32) as i16
}

#[inline]
pub fn add32(a: i32, b: i32) -> i32 {
    a.wrapping_add(b)
}

#[inline]
pub fn uadd32(a: u32, b: u32) -> u32 {
    a.wrapping_add(b)
}

#[inline]
pub fn sub32(a: i32, b: i32) -> i32 {
    a.wrapping_sub(b)
}

/*** Multiplications/Accumulations ***/
#[inline]
pub fn mult16_16(a: i16, b: i16) -> i32 {
    (a as i32) * (b as i32)
}

#[inline]
pub fn mult16_32(a: i16, b: i32) -> i32 {
    (a as i32).wrapping_mul(b)
}

#[inline]
pub fn umult16_16(a: u16, b: u16) -> u32 {
    (a as u32).wrapping_mul(b as u32)
}

#[inline]
pub fn mac16_16(c: i32, a: i16, b: i16) -> i32 {
    add32(c, mult16_16(a, b))
}

#[inline]
pub fn umac16_16(c: u32, a: u16, b: u16) -> u32 {
    uadd32(c, umult16_16(a, b))
}

#[inline]
pub fn msu16_16(c: i32, a: i16, b: i16) -> i32 {
    sub32(c, mult16_16(a, b))
}

#[inline]
pub fn div32(a: i32, b: i32) -> i32 {
    a / b
}

#[inline]
pub fn udiv32(a: u32, b: u32) -> u32 {
    a / b
}

/* Q3 operations */
#[inline]
pub fn mult16_16_q3(a: i16, b: i16) -> i32 {
    shr(mult16_16(a, b), 3)
}

#[inline]
pub fn mac16_16_q3(c: i32, a: i16, b: i16) -> i32 {
    add32(c, mult16_16_q3(a, b))
}

/* Q4 operations */
#[inline]
pub fn mult16_16_q4(a: i16, b: i16) -> i32 {
    shr(mult16_16(a, b), 4)
}

#[inline]
pub fn mac16_16_q4(c: i32, a: i16, b: i16) -> i32 {
    add32(c, mult16_16_q4(a, b))
}

/* Q11 operations */
#[inline]
pub fn mult16_16_q11(a: i16, b: i16) -> i32 {
    shr(mult16_16(a, b), 11)
}

#[inline]
pub fn mult16_16_p11(a: i16, b: i16) -> i32 {
    shr(add32(1024, mult16_16(a, b)), 11)
}

/* Q12 operations */
#[inline]
pub fn mult16_32_q12(a: i16, b: i32) -> i32 {
    ((a as i64 * b as i64) >> 12) as i32
}

#[inline]
pub fn mac16_32_q12(c: i32, a: i16, b: i32) -> i32 {
    add32(c, mult16_32_q12(a, b))
}

#[inline]
pub fn mult16_16_q12(a: i16, b: i16) -> i32 {
    shr(mult16_16(a, b), 12)
}

#[inline]
pub fn msu16_16_q12(c: i32, a: i16, b: i16) -> i32 {
    sub32(c, mult16_16_q12(a, b))
}

/* Q13 operations */
#[inline]
pub fn mult16_16_q13(a: i16, b: i16) -> i32 {
    shr(mult16_16(a, b), 13)
}

#[inline]
pub fn mult16_32_q13(a: i16, b: i32) -> i32 {
    ((a as i64 * b as i64) >> 13) as i32
}

#[inline]
pub fn mac16_16_q13(c: i32, a: i16, b: i16) -> i32 {
    add32(c, mult16_16_q13(a, b))
}

#[inline]
pub fn mac16_32_q13(c: i32, a: i16, b: i32) -> i32 {
    add32(c, mult16_32_q13(a, b))
}

/* Q14 operations */
#[inline]
pub fn mult16_32_p14(a: i16, b: i32) -> i32 {
    ((a as i64 * b as i64 + 8192) >> 14) as i32
}

#[inline]
pub fn mult16_32_q14(a: i16, b: i32) -> i32 {
    ((a as i64 * b as i64) >> 14) as i32
}

#[inline]
pub fn mult16_16_p14(a: i16, b: i16) -> i32 {
    shr(add32(8192, mult16_16(a, b)), 14)
}

#[inline]
pub fn mult16_16_q14(a: i16, b: i16) -> i32 {
    shr(mult16_16(a, b), 14)
}

#[inline]
pub fn mac16_16_q14(c: i32, a: i16, b: i16) -> i32 {
    add32(c, mult16_16_q14(a, b))
}

#[inline]
pub fn msu16_16_q14(c: i32, a: i16, b: i16) -> i32 {
    sub32(c, mult16_16_q14(a, b))
}

#[inline]
pub fn mac16_32_q14(c: i32, a: i16, b: i32) -> i32 {
    add32(c, mult16_32_q14(a, b))
}

/* Q15 operations */
#[inline]
pub fn mult16_16_q15(a: i16, b: i16) -> i32 {
    shr(mult16_16(a, b), 15)
}

#[inline]
pub fn mult16_16_p15(a: i16, b: i16) -> i32 {
    shr(add32(16384, mult16_16(a, b)), 15)
}

#[inline]
pub fn mult16_32_p15(a: i16, b: i32) -> i32 {
    ((a as i64 * b as i64 + 16384) >> 15) as i32
}

#[inline]
pub fn mult16_32_q15(a: i16, b: i32) -> i32 {
    ((a as i64 * b as i64) >> 15) as i32
}

/// Mirror of bcg729's `MULT16_16_Q15` when the first operand is a full 32-bit
/// value (the macro multiplies 32-bit operands, it does not truncate to 16 bits).
#[inline]
pub fn mult32_16_q15(a: i32, b: i16) -> i32 {
    ((a as i64 * b as i64) >> 15) as i32
}

/// Mirror of bcg729's `MULT16_16` when the first operand is a full 32-bit value.
#[inline]
pub fn mult32_16(a: i32, b: i16) -> i32 {
    (a as i64 * b as i64) as i32
}

/// Mirror of bcg729's `MULT16_32_Q15` with both operands 32-bit.
#[inline]
pub fn mult32_32_q15(a: i32, b: i32) -> i32 {
    ((a as i64 * b as i64) >> 15) as i32
}

#[inline]
pub fn mac16_32_p15(c: i32, a: i16, b: i32) -> i32 {
    add32(c, mult16_32_p15(a, b))
}

/* 64 bits operations */
#[inline]
pub fn add64(a: i64, b: i64) -> i64 {
    a.wrapping_add(b)
}

#[inline]
pub fn add64_32(a: i64, b: i32) -> i64 {
    a.wrapping_add(b as i64)
}

#[inline]
pub fn mult32_32(a: i32, b: i32) -> i64 {
    (a as i64).wrapping_mul(b as i64)
}

#[inline]
pub fn div64(a: i64, b: i64) -> i64 {
    a / b
}

#[inline]
pub fn mac64(c: i64, a: i32, b: i32) -> i64 {
    c.wrapping_add((a as i64).wrapping_mul(b as i64))
}

/* Divisions */
#[inline]
pub fn div32_32_q24(a: i32, b: i32) -> i64 {
    ((a as i64) << 24) / (b as i64)
}

#[inline]
pub fn div32_32_q27(a: i32, b: i32) -> i64 {
    sshl64(a as i64, 27) / (b as i64)
}

#[inline]
pub fn div32_32_q31(a: i32, b: i32) -> i64 {
    sshl64(a as i64, 31) / (b as i64)
}

#[inline]
pub fn mult32_32_q23(a: i32, b: i32) -> i32 {
    shr64((a as i64).wrapping_mul(b as i64), 23) as i32
}

#[inline]
pub fn mult32_32_q31(a: i32, b: i32) -> i32 {
    shr64((a as i64).wrapping_mul(b as i64), 31) as i32
}

#[inline]
pub fn mac32_32_q31(c: i32, a: i32, b: i32) -> i32 {
    add32(c, mult32_32_q31(a, b))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_saturate() {
        assert_eq!(saturate(40000, MAX_16 as i32), MAX_16 as i32);
        assert_eq!(saturate(-40000, MAX_16 as i32), MIN_16 as i32);
        assert_eq!(saturate(100, MAX_16 as i32), 100);
    }

    #[test]
    fn test_add16() {
        assert_eq!(add16(10, 20), 30);
        // Overflow check? The C macro ADD16 is just cast: ((word16_t)((word16_t)(a)+(word16_t)(b)))
        // It does NOT saturate. It wraps.
        assert_eq!(add16(MAX_16, 1), MIN_16);
    }

    #[test]
    fn test_mult16_16() {
        assert_eq!(mult16_16(10, 20), 200);
        assert_eq!(mult16_16(MAX_16, 2), 65534);
    }

    #[test]
    fn test_pshr() {
        // PSHR(a,shift) (SHR((a)+((EXTEND32(1)<<((shift))>>1)),shift))
        // shift=1: a + (1<<1>>1) = a+1 >> 1. (Rounding)
        assert_eq!(pshr(3, 1), 2); // (3+1)>>1 = 2
        assert_eq!(pshr(2, 1), 1); // (2+1)>>1 = 1
        assert_eq!(pshr(5, 2), 1); // (5 + (1<<2>>1))>>2 = (5+2)>>2 = 7>>2 = 1.
                                   // Wait. 1<<2 = 4. 4>>1 = 2. 5+2=7. 7>>2 = 1. Correct.
                                   // 6, 2 -> (6+2)>>2 = 2.
    }
}
