use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;
use crate::g729::utils::{count_leading_zeros, unsigned_count_leading_zeros};

/* constants defined in Q16: actual values:
 KL0 = -2.059978
 KL1 = 5.770780
 KL2 = -3.847187
 KL3 = 1.139907
*/
const KL0: Word32 = -135003;
const KL1: Word32 = 378194;
const KL2: Word32 = -252129;
const KL3: Word32 = 74705;

pub fn g729_log2_q0q16(x: Word32) -> Word32 {
    /* first get the integer part and put it in the 16 MSB of return value (in Q16) */
    let leading_zeros = count_leading_zeros(x); /* note: MSB is excluded as considered as sign bit */
    let ret_value = shl32(sub16(30, leading_zeros as Word16) as Word32, 16);

    /* now shift the number to have it on this form 01XX XXXX XXXX XXXX, and keep only 16 bits -> consider it as a number in range [0.5, 1[ in Q0.15 */
    let acc = vshr32(x, 16 - leading_zeros as i32) as Word16;

    /* log2(x) ~= -3.059978 + 5.770780*x - 3.847187*x^2 + 1.139907*x^3 (for .5 < x < 1) Taylor Serie log2(x) at x near 0.75 */
    /* log2(acc) +1 = -135003 +acc*(378194 +acc*(-252129 + acc*74705)) acc in Q15 and constants in Q16 -> final result will be log2(x) -int in Q16(on 32 bits) */
    let acc32 = add32(
        KL0,
        mult16_32_q15(
            acc,
            add32(KL1, mult16_32_q15(acc, add32(KL2, mult16_32_q15(acc, KL3)))),
        ),
    );

    add32(ret_value, acc32)
}

/* constants defined in Q15: actual values:
 E0 = 1
 E1 = log(2)
 E2 = 3-4*log(2)
 E3 = 3*log(2) - 2
*/
const E0: Word16 = 16384;
const E1: Word16 = 11356;
const E2: Word16 = 3726;
const E3: Word16 = 1301;

pub fn g729_exp2_q11q16(x: Word16) -> Word32 {
    let integer = shr16(x, 11);
    if integer > 14 {
        return 0x7fffffff;
    } else {
        if integer < -15 {
            return 0;
        }
    }
    let mut frac = shl16(x - shl16(integer, 11), 3);
    frac = add16(
        E0,
        mult16_16_q14(
            frac,
            add16(
                E1,
                mult16_16_q14(frac, add16(E2, mult16_16_q14(E3, frac) as Word16) as Word16)
                    as Word16,
            ) as Word16,
        ) as Word16,
    );
    vshr32(extend32(frac), -(integer as i32) - 2)
}

/* constants in Q14 */
const C0: Word16 = 3634;
const C1: Word16 = 21173;
const C2: Word16 = -12627;
const C3: Word16 = 4204;

pub fn g729_sqrt_q0q7(x: UWord32) -> Word32 {
    if x == 0 {
        return 0;
    }
    /* set x in Q14 in range [0.25,1[ */
    let k = (19 - unsigned_count_leading_zeros(x) as i32) >> 1;
    let x_scaled = vshr32(x as Word32, k * 2); /* x = x.2^-2k */

    /* sqrt(x) ~= 0.22178 + 1.29227*x - 0.77070*x^2 + 0.25659*x^3 (for .25 < x < 1) */
    /* consider x as in Q14: y = x.2^(-2k-14) -> and give sqrt(y).2^14 = sqrt(x).2^(-k-7).2^14 */
    // Note: x_scaled is Word32, but used as Word16 in MULT16_16_Q14.
    // In C: MULT16_16_Q14(x, ...) where x is word32_t. It casts to word16_t.
    let x16 = x_scaled as Word16;
    let mut rt = add16(
        C0,
        mult16_16_q14(
            x16,
            add16(
                C1,
                mult16_16_q14(x16, add16(C2, mult16_16_q14(x16, C3) as Word16) as Word16) as Word16,
            ) as Word16,
        ) as Word16,
    ) as Word32; /* rt = sqrt(x).2^(7-k)*/
    rt = vshr32(rt, -k); /* rt = sqrt(x).2^7 */
    rt
}

pub fn g729_inv_sqrt_q0q31(x: Word32) -> Word32 {
    if x == 1 {
        return MAXINT32;
    }
    div32_32_q24(g729_sqrt_q0q7(x as UWord32), x) as Word32
}

/* constants Q0.15 */
const KCOS1: Word32 = 32768;
const KCOS2: Word32 = -16384;
const KCOS3: Word32 = 1365;
const KCOS4: Word32 = -46;

const KSIN1: Word32 = 32768;
const KSIN2: Word32 = -5461;
const KSIN3: Word32 = 273;
const KSIN4: Word32 = -7;

pub fn g729_cos_q13q15(mut x: Word16) -> Word16 {
    /* input var x in Q2.13 and in ]0, Pi[ */
    let x2: Word16;
    let mut x_scaled: Word16;
    if x < 12868 {
        if x < 6434 {
            /* x in ]0, Pi/4[ */
            x2 = mult16_16_p11(x, x) as Word16; /* in Q0.15 */
            saturate(
                add32(
                    KCOS1,
                    mult16_16_p15(
                        x2,
                        add32(
                            KCOS2,
                            mult16_16_p15(
                                x2,
                                add32(KCOS3, mult16_16_p15(KCOS4 as Word16, x2)) as Word16,
                            ),
                        ) as Word16,
                    ),
                ),
                MAXINT16 as Word32,
            ) as Word16
        } else {
            /* x in [Pi/4, Pi/2[ */
            x = sub16(12868, x); /* x = pi/2 -x, x in [0, Pi/4] in Q0.13 */
            x2 = mult16_16_p11(x, x) as Word16; /* in Q0.15 */
            mult16_16_p13(
                x,
                add32(
                    KSIN1,
                    mult16_16_p15(
                        x2,
                        add32(
                            KSIN2,
                            mult16_16_p15(
                                x2,
                                add32(KSIN3, mult16_16_p15(KSIN4 as Word16, x2)) as Word16,
                            ),
                        ) as Word16,
                    ),
                ) as Word16,
            ) as Word16
        }
    } else {
        /* x in [Pi/2, Pi[ */
        x_scaled = sub16(25736, x); /* xScaled = Pi - x -> in [0,Pi/2] with cos(Pi-x) = -cos(x) and sin(Pi-x) =  */
        if x < 19302 {
            /* x in [Pi/2, 3Pi/4], xScaled in [Pi/4, Pi/2] */
            x_scaled = sub16(12868, x_scaled); /* xScaled = pi/2 - xScaled = x - Pi/2, xScaled in [0, Pi/4] in Q0.13 */
            x2 = mult16_16_p11(x_scaled, x_scaled) as Word16; /* in Q0.15 */
            mult16_16_p13(
                neg16(x_scaled),
                add32(
                    KSIN1,
                    mult16_16_p15(
                        x2,
                        add32(
                            KSIN2,
                            mult16_16_p15(
                                x2,
                                add32(KSIN3, mult16_16_p15(KSIN4 as Word16, x2)) as Word16,
                            ),
                        ) as Word16,
                    ),
                ) as Word16,
            ) as Word16
        } else {
            /* x in [3Pi/4, Pi[ -> xScaled in [0, Pi/4], cos(xScaled) = -cos(x) */
            x2 = mult16_16_p11(x_scaled, x_scaled) as Word16; /* in Q0.15 */
            sub32(
                -KCOS1,
                mult16_16_p15(
                    x2,
                    add32(
                        KCOS2,
                        mult16_16_p15(
                            x2,
                            add32(KCOS3, mult16_16_p15(KCOS4 as Word16, x2)) as Word16,
                        ),
                    ) as Word16,
                ),
            ) as Word16
        }
    }
}

/* KPI6 = pi/6 in Q15 */
const KPI6: Word32 = 17157;
/* KtanPI6 = tan(pi/6) in Q15 */
const KTAN_PI6: Word32 = 18919;
/* KtanPI12 = tan(pi/12) in Q15 */
const KTAN_PI12: Word32 = 8780;

/* B = 0.257977658811405 in Q15 */
const ATAN_B: Word16 = 8453;
/* C = 0.59120450521312 in Q15 */
const ATAN_C: Word16 = 19373;

pub fn g729_atan_q15q13(mut x: Word32) -> Word16 {
    /* constants for rational polynomial */
    let mut angle: Word32;
    let x2: Word16;
    let mut high_segment = false;
    let mut sign = false;
    let mut complement = false;

    /* make argument positive */
    if x < 0 {
        x = neg32(x);
        sign = true;
    }

    /* limit argument to 0..1 */
    if x > ONE_IN_Q15 {
        complement = true;
        x = div32(ONE_IN_Q30, x); /* 1/x in Q15 */
    }

    /* determine segmentation */
    if x > KTAN_PI12 {
        high_segment = true;
        /* x = (x - k)/(1 + k*x); */
        x = div32(
            sshl(sub32(x, KTAN_PI6), 15),
            add32(mult16_32_q15(KTAN_PI6 as Word16, x), ONE_IN_Q15),
        );
    }

    /* argument is now < tan(15 degrees) */
    /* approximate the function */
    x2 = mult16_16_q15(x as Word16, x as Word16) as Word16;
    // angle = div32(mult16_16(x as Word16, add32(ONE_IN_Q15, mult16_16_q15(ATAN_B, x2)) as Word16), add32(ONE_IN_Q15, mult16_16_q15(ATAN_C, x2)));  /* ang = x*(1 + B*x2)/(1 + C*x2) */
    let num = (x as Word32) * add32(ONE_IN_Q15, mult16_16_q15(ATAN_B, x2));
    let den = add32(ONE_IN_Q15, mult16_16_q15(ATAN_C, x2));
    angle = div32(num, den);

    /* now restore offset if needed */
    if high_segment {
        angle += KPI6;
    }

    /* restore complement if needed */
    if complement {
        angle = sub32(HALF_PI_Q15_32, angle);
    }

    /* set result in Q13 */
    angle = pshr(angle, 2);

    /* restore sign if needed */
    if sign {
        neg16(angle as Word16)
    } else {
        angle as Word16
    }
}

pub fn g729_asin_q15q13(x: Word16) -> Word16 {
    let xx = mult16_16(x, x);
    let sub = sub32(ONE_IN_Q30, xx);
    let sqrt = g729_sqrt_q0q7(sub as UWord32);
    let denom = pshr(sqrt, 7);
    let num = sshl(x as Word32, 15);
    let div = div32(num, denom);
    let atan = g729_atan_q15q13(div);
    atan
}

pub fn g729_acos_q15q13(x: Word16) -> Word16 {
    let asin = g729_asin_q15q13(x);
    sub16(HALF_PI_Q13, asin)
}
