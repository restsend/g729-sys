// Basic operations and macros ported from basicOperationsMacros.h and fixedPointMacros.h

pub type Word16 = i16;
pub type UWord16 = u16;
pub type Word32 = i32;
pub type UWord32 = u32;
pub type Word64 = i64;

pub const MAX_16: Word16 = 0x7fff;
pub const MIN_16: Word16 = -0x8000;
pub const MAX_U16: UWord16 = 0xffff;
pub const MAX_32: Word32 = 0x7fffffff;
pub const MIN_32: Word32 = -0x80000000;

#[inline]
pub fn extend32(x: Word16) -> Word32 {
    x as Word32
}

#[inline]
pub fn neg16(x: Word16) -> Word16 {
    x.wrapping_neg()
}

#[inline]
pub fn neg32(x: Word32) -> Word32 {
    x.wrapping_neg()
}

#[inline]
pub fn neg64(x: Word64) -> Word64 {
    x.wrapping_neg()
}

/*** shifts ***/
#[inline]
pub fn shr(a: Word32, shift: u32) -> Word32 {
    a >> shift
}

#[inline]
pub fn shl(a: Word32, shift: u32) -> Word32 {
    a << shift
}

#[inline]
pub fn sshl(a: Word32, shift: u32) -> Word32 {
    // Rust's left shift on signed integers is defined as arithmetic shift (preserves sign bit behavior in 2's complement)
    // which matches the intent of the C macro workaround for UB.
    a << shift
}

#[inline]
pub fn ushl(a: UWord32, shift: u32) -> UWord32 {
    a << shift
}

/* shift right with rounding: used to extract the integer value of a Qa number */
#[inline]
pub fn pshr(a: Word32, shift: u32) -> Word32 {
    if shift == 0 {
        return a;
    }
    let half = 1 << (shift - 1);
    (a.wrapping_add(half)) >> shift
}

/* shift right with checking on sign of shift value */
#[inline]
pub fn vshr32(a: Word32, shift: i32) -> Word32 {
    if shift > 0 {
        shr32(a, shift as u32)
    } else {
        shl32(a, (-shift) as u32)
    }
}

#[inline]
pub fn svshr32(a: Word32, shift: i32) -> Word32 {
    if shift > 0 {
        shr32(a, shift as u32)
    } else {
        sshl(a, (-shift) as u32)
    }
}

#[inline]
pub fn shr16(a: Word16, shift: u32) -> Word16 {
    a >> shift
}

#[inline]
pub fn shl16(a: Word16, shift: u32) -> Word16 {
    a << shift
}

#[inline]
pub fn shr32(a: Word32, shift: u32) -> Word32 {
    a >> shift
}

#[inline]
pub fn shl32(a: Word32, shift: u32) -> Word32 {
    a << shift
}

#[inline]
pub fn shr64(a: Word64, shift: u32) -> Word64 {
    a >> shift
}

#[inline]
pub fn shl64(a: Word64, shift: u32) -> Word64 {
    a << shift
}

#[inline]
pub fn sshl64(a: Word64, shift: u32) -> Word64 {
    a << shift
}

/* avoid overflows: a+1 is used to check on negative value because range of a 2n signed bits int is -2pow(n) - 2pow(n)-1 */
/* SATURATE Macro shall be called with MAXINT(nbits). Ex: SATURATE(x,MAXINT16) with MAXINT16  defined to 2pow(16) - 1 */
#[inline]
pub fn saturate(x: Word32, a: Word32) -> Word32 {
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
pub fn usaturate(x: Word32, a: Word32) -> Word32 {
    if x > a {
        a
    } else {
        x
    }
}

/* absolute value */
#[inline]
pub fn abs(a: Word32) -> Word32 {
    a.abs()
}

#[inline]
pub fn abs16(a: Word16) -> Word16 {
    a.abs()
}

#[inline]
pub fn abs_s(a: Word16) -> Word16 {
    if a == MIN_16 {
        MAX_16
    } else {
        a.abs()
    }
}

/*** add and sub ***/
#[inline]
pub fn add16(a: Word16, b: Word16) -> Word16 {
    (a as Word32 + b as Word32) as Word16
}

#[inline]
pub fn sub16(a: Word16, b: Word16) -> Word16 {
    (a as Word32 - b as Word32) as Word16
}

#[inline]
pub fn add32(a: Word32, b: Word32) -> Word32 {
    a.wrapping_add(b)
}

#[inline]
pub fn uadd32(a: UWord32, b: UWord32) -> UWord32 {
    a.wrapping_add(b)
}

#[inline]
pub fn sub32(a: Word32, b: Word32) -> Word32 {
    a.wrapping_sub(b)
}

/*** Multiplications/Accumulations ***/
#[inline]
pub fn mult16_16(a: Word16, b: Word16) -> Word32 {
    (a as Word32) * (b as Word32)
}

#[inline]
pub fn mult16_32(a: Word16, b: Word32) -> Word32 {
    (a as Word32) * b
}

#[inline]
pub fn umult16_16(a: UWord16, b: UWord16) -> UWord32 {
    (a as UWord32).wrapping_mul(b as UWord32)
}

#[inline]
pub fn mac16_16(c: Word32, a: Word16, b: Word16) -> Word32 {
    add32(c, mult16_16(a, b))
}

#[inline]
pub fn umac16_16(c: UWord32, a: UWord16, b: UWord16) -> UWord32 {
    uadd32(c, umult16_16(a, b))
}

#[inline]
pub fn msu16_16(c: Word32, a: Word16, b: Word16) -> Word32 {
    sub32(c, mult16_16(a, b))
}

#[inline]
pub fn div32(a: Word32, b: Word32) -> Word32 {
    a / b
}

#[inline]
pub fn udiv32(a: UWord32, b: UWord32) -> UWord32 {
    a / b
}

/* Q3 operations */
#[inline]
pub fn mult16_16_q3(a: Word16, b: Word16) -> Word32 {
    shr(mult16_16(a, b), 3)
}

#[inline]
pub fn mult16_32_q3(a: Word16, b: Word32) -> Word32 {
    add32(
        (a as Word32).wrapping_mul(shr(b, 3)),
        shr(mult16_16(a, (b & 0x00000007) as Word16), 3),
    )
}

#[inline]
pub fn mac16_16_q3(c: Word32, a: Word16, b: Word16) -> Word32 {
    add32(c, mult16_16_q3(a, b))
}

/* Q4 operations */
#[inline]
pub fn mult16_16_q4(a: Word16, b: Word16) -> Word32 {
    shr(mult16_16(a, b), 4)
}

#[inline]
pub fn umult16_16_q4(a: UWord16, b: UWord16) -> UWord32 {
    // C macro says SHR. #define UMULT16_16_Q4(a,b) (SHR(UMULT16_16((a),(b)),4))
    (umult16_16(a, b)) >> 4
}

#[inline]
pub fn umac16_16_q4(c: UWord32, a: UWord16, b: UWord16) -> UWord32 {
    uadd32(c, umult16_16_q4(a, b))
}

#[inline]
pub fn mac16_16_q4(c: Word32, a: Word16, b: Word16) -> Word32 {
    add32(c, mult16_16_q4(a, b))
}

/* Q11 operations */
#[inline]
pub fn mult16_16_q11(a: Word16, b: Word16) -> Word32 {
    shr(mult16_16(a, b), 11)
}

#[inline]
pub fn mult16_16_p11(a: Word16, b: Word16) -> Word32 {
    shr(add32(1024, mult16_16(a, b)), 11)
}

/* Q12 operations */
#[inline]
pub fn mult16_32_q12(a: Word16, b: Word32) -> Word32 {
    ((a as i64 * b as i64) >> 12) as Word32
}

#[inline]
pub fn mac16_32_q12(c: Word32, a: Word16, b: Word32) -> Word32 {
    add32(c, mult16_32_q12(a, b))
}

#[inline]
pub fn mult16_16_q12(a: Word16, b: Word16) -> Word32 {
    shr(mult16_16(a, b), 12)
}

#[inline]
pub fn mac16_16_q12(c: Word32, a: Word16, b: Word16) -> Word32 {
    add32(c, mult16_16_q12(a, b))
}

#[inline]
pub fn msu16_16_q12(c: Word32, a: Word16, b: Word16) -> Word32 {
    sub32(c, mult16_16_q12(a, b))
}

/* Q13 operations */
#[inline]
pub fn mult16_16_q13(a: Word16, b: Word16) -> Word32 {
    shr(mult16_16(a, b), 13)
}

#[inline]
pub fn mult16_16_p13(a: Word16, b: Word16) -> Word32 {
    shr(add32(4096, mult16_16(a, b)), 13)
}

#[inline]
pub fn mult16_32_q13(a: Word16, b: Word32) -> Word32 {
    ((a as i64 * b as i64) >> 13) as Word32
}

#[inline]
pub fn mac16_16_q13(c: Word32, a: Word16, b: Word16) -> Word32 {
    add32(c, mult16_16_q13(a, b))
}

#[inline]
pub fn mac16_32_q13(c: Word32, a: Word16, b: Word32) -> Word32 {
    add32(c, mult16_32_q13(a, b))
}

/* Q14 operations */
#[inline]
pub fn mult16_32_p14(a: Word16, b: Word32) -> Word32 {
    ((a as i64 * b as i64 + 8192) >> 14) as Word32
}

#[inline]
pub fn mult16_32_q14(a: Word16, b: Word32) -> Word32 {
    ((a as i64 * b as i64) >> 14) as Word32
}

#[inline]
pub fn mult16_16_p14(a: Word16, b: Word16) -> Word32 {
    shr(add32(8192, mult16_16(a, b)), 14)
}

#[inline]
pub fn mult16_16_q14(a: Word16, b: Word16) -> Word32 {
    shr(mult16_16(a, b), 14)
}

#[inline]
pub fn mac16_16_q14(c: Word32, a: Word16, b: Word16) -> Word32 {
    add32(c, mult16_16_q14(a, b))
}

#[inline]
pub fn msu16_16_q14(c: Word32, a: Word16, b: Word16) -> Word32 {
    sub32(c, mult16_16_q14(a, b))
}

#[inline]
pub fn mac16_32_q14(c: Word32, a: Word16, b: Word32) -> Word32 {
    add32(c, mult16_32_q14(a, b))
}

/* Q15 operations */
#[inline]
pub fn mult16_16_q15(a: Word16, b: Word16) -> Word32 {
    shr(mult16_16(a, b), 15)
}

#[inline]
pub fn mult16_16_p15(a: Word16, b: Word16) -> Word32 {
    shr(add32(16384, mult16_16(a, b)), 15)
}

#[inline]
pub fn mult16_32_p15(a: Word16, b: Word32) -> Word32 {
    ((a as i64 * b as i64 + 16384) >> 15) as Word32
}

#[inline]
pub fn mult16_32_q15(a: Word16, b: Word32) -> Word32 {
    ((a as i64 * b as i64) >> 15) as Word32
}

#[inline]
pub fn mac16_32_p15(c: Word32, a: Word16, b: Word32) -> Word32 {
    add32(c, mult16_32_p15(a, b))
}

/* 64 bits operations */
#[inline]
pub fn add64(a: Word64, b: Word64) -> Word64 {
    a.wrapping_add(b)
}

#[inline]
pub fn sub64(a: Word64, b: Word32) -> Word64 {
    a.wrapping_sub(b as Word64)
}

#[inline]
pub fn add64_32(a: Word64, b: Word32) -> Word64 {
    a.wrapping_add(b as Word64)
}

#[inline]
pub fn mult32_32(a: Word32, b: Word32) -> Word64 {
    (a as Word64).wrapping_mul(b as Word64)
}

#[inline]
pub fn div64(a: Word64, b: Word64) -> Word64 {
    a / b
}

#[inline]
pub fn mac64(c: Word64, a: Word32, b: Word32) -> Word64 {
    c.wrapping_add((a as Word64).wrapping_mul(b as Word64))
}

/* Divisions */
#[inline]
pub fn div32_32_q24(a: Word32, b: Word32) -> Word64 {
    ((a as Word64) << 24) / (b as Word64)
}

#[inline]
pub fn div32_32_q27(a: Word32, b: Word32) -> Word64 {
    sshl64(a as Word64, 27) / (b as Word64)
}

#[inline]
pub fn div32_32_q31(a: Word32, b: Word32) -> Word64 {
    sshl64(a as Word64, 31) / (b as Word64)
}

#[inline]
pub fn mult32_32_q23(a: Word32, b: Word32) -> Word32 {
    shr64((a as Word64).wrapping_mul(b as Word64), 23) as Word32
}

#[inline]
pub fn mult32_32_q31(a: Word32, b: Word32) -> Word32 {
    shr64((a as Word64).wrapping_mul(b as Word64), 31) as Word32
}

#[inline]
pub fn mac32_32_q31(c: Word32, a: Word32, b: Word32) -> Word32 {
    add32(c, mult32_32_q31(a, b))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_saturate() {
        assert_eq!(saturate(40000, MAX_16 as Word32), MAX_16 as Word32);
        assert_eq!(saturate(-40000, MAX_16 as Word32), MIN_16 as Word32);
        assert_eq!(saturate(100, MAX_16 as Word32), 100);
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
