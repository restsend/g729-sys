use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;

pub const NB_COMPUTED_VALUES_CHEBYSHEV_POLYNOMIAL: usize = 51;

/* x = cos(w) with w in [0,Pi] in 50 steps */
/* in Q15 */
static COS_W0_PI: [Word16; NB_COMPUTED_VALUES_CHEBYSHEV_POLYNOMIAL] = [
    32760, 32703, 32509, 32187, 31738, 31164, 30466, 29649, 28714, 27666, 26509, 25248, 23886,
    22431, 20887, 19260, 17557, 15786, 13951, 12062, 10125, 8149, 6140, 4106, 2057, 0, -2057,
    -4106, -6140, -8149, -10125, -12062, -13951, -15786, -17557, -19260, -20887, -22431, -23886,
    -25248, -26509, -27666, -28714, -29649, -30466, -31164, -31738, -32187, -32509, -32703, -32760,
];

/*****************************************************************************/
/* ChebyshevPolynomial : Compute the Chebyshev polynomial, spec 3.2.3 eq17   */
/*    parameters:                                                            */
/*      -(i) x : input value of polynomial function in Q15                   */
/*      -(i) f : the polynome coefficients, 6 values in Q15 on 32 bits       */
/*           f[0] is not used                                                */
/*    return value :                                                         */
/*      - result of polynomial function in Q15                               */
/*                                                                           */
/*****************************************************************************/
fn chebyshev_polynomial(x: Word16, f: &[Word32]) -> Word32 {
    /* bk in Q15*/
    let mut bk: Word32;
    let mut bk1 = add32(shl(x as Word32, 1), f[1]); /* init: b4=2x+f1 */
    let mut bk2 = ONE_IN_Q15 as Word32; /* init: b5=1 */

    for k in (1..=3).rev() {
        /* at the end of loop execution we have b1 in bk1 and b2 in bk2 */
        bk = sub32(add32(shl(mult16_32_q15(x, bk1), 1), f[5 - k]), bk2); /* bk = 2*x*bk1 − bk2 + f(5-k) all in Q15*/
        bk2 = bk1;
        bk1 = bk;
    }

    sub32(add32(mult16_32_q15(x, bk1), shr(f[5], 1)), bk2) /* C(x) = x*b1 - b2 + f(5)/2 */
}

/*****************************************************************************/
/* LP2LSPConversion : Compute polynomials, find their roots as in spec A3.2.3*/
/*    parameters:                                                            */
/*      -(i) LPCoefficients[] : 10 coefficients in Q12                       */
/*      -(o) LSPCoefficients[] : 10 coefficients in Q15                      */
/*                                                                           */
/*    return value :                                                         */
/*      - boolean: 1 if all roots found, 0 if unable to compute 10 roots     */
/*                                                                           */
/*****************************************************************************/
pub fn lp2lsp_conversion(lp_coefficients: &[Word16], lsp_coefficients: &mut [Word16]) -> bool {
    let mut f1 = [0 as Word32; 6];
    let mut f2 = [0 as Word32; 6]; /* coefficients for polynomials F1 anf F2 in Q12 for computation, then converted in Q15 for the Chebyshev Polynomial function */
    let mut number_of_root_found = 0; /* used to check the final number of roots found and exit the loop on each polynomial computation when we have 10 roots */
    let mut previous_cx: Word32;
    let mut cx: Word32; /* value of Chebyshev Polynomial at current point in Q15 */

    /*** Compute the polynomials coefficients according to spec 3.2.3 eq15 ***/
    f1[0] = ONE_IN_Q12 as Word32; /* values 0 are not part of the output, they are just used for computation purpose */
    f2[0] = ONE_IN_Q12 as Word32;

    for i in 0..5 {
        f1[i + 1] = add32(
            lp_coefficients[i] as Word32,
            sub32(lp_coefficients[9 - i] as Word32, f1[i]),
        ); /* note: index on LPCoefficients are -1 respect to spec because the unused value 0 is not stored */
        f2[i + 1] = add32(
            f2[i],
            sub32(
                lp_coefficients[i] as Word32,
                lp_coefficients[9 - i] as Word32,
            ),
        ); /* note: index on LPCoefficients are -1 respect to spec because the unused value 0 is not stored */
    }
    /* convert the coefficients from Q12 to Q15 to be used by the Chebyshev Polynomial function (f1/2[0] aren't used so they are not converted) */
    for i in 1..6 {
        f1[i] = shl(f1[i], 3);
        f2[i] = shl(f2[i], 3);
    }

    /*** Compute at each step(50 steps for the AnnexA version) the Chebyshev polynomial to find the 10 roots ***/
    /* start using f1 polynomials coefficients and altern with f2 after founding each root (spec 3.2.3 eq13 and eq14) */
    let mut use_f1 = true; /* start with f1 coefficients */
    previous_cx = chebyshev_polynomial(COS_W0_PI[0], &f1); /* compute the first point and store it as the previous value for polynomial */

    for i in 1..NB_COMPUTED_VALUES_CHEBYSHEV_POLYNOMIAL {
        cx = chebyshev_polynomial(COS_W0_PI[i], if use_f1 { &f1 } else { &f2 });
        if ((previous_cx ^ cx) & 0x10000000) != 0 {
            /* check signe change by XOR on the value of first bit */
            /* divide 2 times the interval to find a more accurate root */
            let mut x_low = COS_W0_PI[i - 1];
            let mut x_high = COS_W0_PI[i];
            let mut x_mean: Word16;

            for _j in 0..2 {
                let middle_cx: Word32;
                x_mean = shr(add32(x_low as Word32, x_high as Word32), 1) as Word16;
                middle_cx = chebyshev_polynomial(x_mean, if use_f1 { &f1 } else { &f2 }); /* compute the polynome for the value in the middle of current interval */

                if ((previous_cx ^ middle_cx) & 0x10000000) != 0 {
                    /* check signe change by XOR on the value of first bit */
                    x_high = x_mean;
                    cx = middle_cx; /* used for linear interpolation on root */
                } else {
                    x_low = x_mean;
                    previous_cx = middle_cx;
                }
            }

            /* toggle the polynomial coefficients in use between f1 and f2 */
            use_f1 = !use_f1;

            /* linear interpolation for better root accuracy */
            /* xMean = xLow - (xHigh-xLow)* previousCx/(Cx-previousCx); */
            x_mean = sub32(
                x_low as Word32,
                mult16_32_q15(
                    sub32(x_high as Word32, x_low as Word32) as Word16,
                    div32(
                        shl(saturate(previous_cx, MAXINT17), 14),
                        shr(sub32(cx, previous_cx), 1),
                    ),
                ),
            ) as Word16; /* Cx are in Q2.15 so we can shift them left 14 bits, the denominator is shifted righ by 1 so the division result is in Q15 */

            /* recompute previousCx with the new coefficients */
            previous_cx = chebyshev_polynomial(x_mean, if use_f1 { &f1 } else { &f2 });

            lsp_coefficients[number_of_root_found] = x_mean;

            number_of_root_found += 1;
            if number_of_root_found == NB_LSP_COEFF {
                break;
            } /* exit the for loop as soon as we habe all the LSP*/
        }
    }
    if number_of_root_found != NB_LSP_COEFF {
        return false;
    } /* we were not able to find the 10 roots */

    true
}
