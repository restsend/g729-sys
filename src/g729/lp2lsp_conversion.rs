use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;

pub const NB_COMPUTED_VALUES_CHEBYSHEV_POLYNOMIAL: usize = 51;

/* in Q15 */
static COS_W0_PI: [i16; NB_COMPUTED_VALUES_CHEBYSHEV_POLYNOMIAL] = [
    32760, 32703, 32509, 32187, 31738, 31164, 30466, 29649, 28714, 27666, 26509, 25248, 23886,
    22431, 20887, 19260, 17557, 15786, 13951, 12062, 10125, 8149, 6140, 4106, 2057, 0, -2057,
    -4106, -6140, -8149, -10125, -12062, -13951, -15786, -17557, -19260, -20887, -22431, -23886,
    -25248, -26509, -27666, -28714, -29649, -30466, -31164, -31738, -32187, -32509, -32703, -32760,
];

/*****************************************************************************/
/* ChebyshevPolynomial : Compute the Chebyshev polynomial, spec 3.2.3 eq17   */

/*      -(i) x : input value of polynomial function in Q15                   */
/*      -(i) f : the polynome coefficients, 6 values in Q15 on 32 bits       */

/*      - result of polynomial function in Q15                               */

/*****************************************************************************/
fn chebyshev_polynomial(x: i16, f: &[i32]) -> i32 {
    /* bk in Q15*/
    let mut bk: i32;
    let mut bk1 = add32(shl(x as i32, 1), f[1]);
    let mut bk2 = ONE_IN_Q15;

    for k in (1..=3).rev() {
        bk = sub32(add32(shl(mult16_32_q15(x, bk1), 1), f[5 - k]), bk2); /* bk = 2*x*bk1 − bk2 + f(5-k) all in Q15*/
        bk2 = bk1;
        bk1 = bk;
    }

    sub32(add32(mult16_32_q15(x, bk1), shr(f[5], 1)), bk2)
}

/*****************************************************************************/
/* LP2LSPConversion : Compute polynomials, find their roots as in spec A3.2.3*/

/*      -(i) LPCoefficients[] : 10 coefficients in Q12                       */
/*      -(o) LSPCoefficients[] : 10 coefficients in Q15                      */

/*****************************************************************************/
pub fn lp2lsp_conversion(lp_coefficients: &[i16], lsp_coefficients: &mut [i16]) -> bool {
    let mut f1 = [0_i32; 6];
    let mut f2 = [0_i32; 6]; /* coefficients for polynomials F1 anf F2 in Q12 for computation, then converted in Q15 for the Chebyshev Polynomial function */
    let mut number_of_root_found = 0;
    let mut previous_cx: i32;
    let mut cx: i32; /* value of Chebyshev Polynomial at current point in Q15 */

    /*** Compute the polynomials coefficients according to spec 3.2.3 eq15 ***/
    f1[0] = ONE_IN_Q12;
    f2[0] = ONE_IN_Q12;

    for i in 0..5 {
        f1[i + 1] = add32(
            lp_coefficients[i] as i32,
            sub32(lp_coefficients[9 - i] as i32, f1[i]),
        ); /* note: index on LPCoefficients are -1 respect to spec because the unused value 0 is not stored */
        f2[i + 1] = add32(
            f2[i],
            sub32(lp_coefficients[i] as i32, lp_coefficients[9 - i] as i32),
        ); /* note: index on LPCoefficients are -1 respect to spec because the unused value 0 is not stored */
    }

    for i in 1..6 {
        f1[i] = shl(f1[i], 3);
        f2[i] = shl(f2[i], 3);
    }

    /*** Compute at each step(50 steps for the AnnexA version) the Chebyshev polynomial to find the 10 roots ***/
    /* start using f1 polynomials coefficients and altern with f2 after founding each root (spec 3.2.3 eq13 and eq14) */
    let mut use_f1 = true;
    previous_cx = chebyshev_polynomial(COS_W0_PI[0], &f1);

    for i in 1..NB_COMPUTED_VALUES_CHEBYSHEV_POLYNOMIAL {
        cx = chebyshev_polynomial(COS_W0_PI[i], if use_f1 { &f1 } else { &f2 });
        if ((previous_cx ^ cx) & 0x10000000) != 0 {
            let mut x_low = COS_W0_PI[i - 1];
            let mut x_high = COS_W0_PI[i];
            let mut x_mean: i16;

            for _j in 0..2 {
                x_mean = shr(add32(x_low as i32, x_high as i32), 1) as i16;
                let middle_cx: i32 = chebyshev_polynomial(x_mean, if use_f1 { &f1 } else { &f2 });

                if ((previous_cx ^ middle_cx) & 0x10000000) != 0 {
                    x_high = x_mean;
                    cx = middle_cx;
                } else {
                    x_low = x_mean;
                    previous_cx = middle_cx;
                }
            }

            use_f1 = !use_f1;

            let delta_x = sub32(x_high as i32, x_low as i32);
            let interp = if previous_cx == cx {
                mult32_32_q15(delta_x, if previous_cx > 0 { MAXINT32 } else { MININT32 })
            } else {
                mult32_32_q15(
                    delta_x,
                    shl(
                        div32(
                            shl(saturate(previous_cx, MAXINT17), 14),
                            sub32(cx, previous_cx),
                        ),
                        1,
                    ),
                )
            };
            x_mean = sub32(x_low as i32, interp) as i16;

            previous_cx = chebyshev_polynomial(x_mean, if use_f1 { &f1 } else { &f2 });

            lsp_coefficients[number_of_root_found] = x_mean;

            number_of_root_found += 1;
            if number_of_root_found == NB_LSP_COEFF {
                break;
            }
        }
    }
    if number_of_root_found != NB_LSP_COEFF {
        return false;
    }

    true
}
