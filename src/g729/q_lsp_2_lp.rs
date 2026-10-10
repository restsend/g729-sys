use crate::g729::basic_operations::*;

/*****************************************************************************/
/* computePolynomialCoefficients : according to spec. 3.2.6                  */

/*      -(o) f : 6 values in Q24 : polynomial coefficients on 32 bits        */

/*****************************************************************************/
#[cfg_attr(target_arch = "xtensa", inline(never))]
fn compute_polynomial_coefficients(q_lsp: &[i16], f: &mut [i32]) {
    f[0] = 16777216;

    /* correspont to i=1 in the algorithm description in spec. 3.2.6 */
    f[1] = mult16_16(q_lsp[0], -1024);

    /* Note : index of qLSP are -1 respect of what is in the spec because qLSP array is indexed from 0-9 and not 1-10 */
    for i in 2..6 {
        f[i] = sshl(
            sub32(f[i - 2], mult16_32_p15(q_lsp[2 * i - 2], f[i - 1])),
            1,
        ); /* with qLSP in Q0.15 and f in Q24 */
        for j in (2..i).rev() {
            f[j] = add32(
                f[j],
                sub32(f[j - 2], mult16_32_p14(q_lsp[2 * i - 2], f[j - 1])),
            ); /* qLPS in Q0.15 and f in Q24, using MULT16_32_P14 instead of P15 does the *2 on qLSP. Result in Q24 */
        }

        f[1] = sub32(f[1], sshl(q_lsp[2 * i - 2] as i32, 10)); /* qLSP in Q0.15, must be shift by 9 to get in Q24 and one more to be *2 */
    }
}

/*****************************************************************************/
/* qLSP2LP : convert qLSP into LP parameters according to spec. 3.2.6        */

/*      -(o) LP : 10 LP coefficients in Q12                                  */

/*****************************************************************************/
#[cfg_attr(target_arch = "xtensa", inline(never))]
pub fn q_lsp_2_lp(q_lsp: &[i16], lp: &mut [i16]) {
    let mut f1 = [0_i32; 6];
    let mut f2 = [0_i32; 6];

    compute_polynomial_coefficients(q_lsp, &mut f1);

    compute_polynomial_coefficients(&q_lsp[1..], &mut f2);

    for i in (1..6).rev() {
        f1[i] = add32(f1[i], f1[i - 1]); /* f1 is still in Q24 */
        f2[i] = sub32(f2[i], f2[i - 1]); /* f1 is still in Q24 */
    }

    for i in 0..5 {
        lp[i] = pshr(add32(f1[i + 1], f2[i + 1]), 13) as i16; /* f1 and f2 in Q24, LP in Q12 */

        lp[9 - i] = pshr(sub32(f1[i + 1], f2[i + 1]), 13) as i16; /* f1 and f2 in Q24, LP in Q12 */
    }
}
