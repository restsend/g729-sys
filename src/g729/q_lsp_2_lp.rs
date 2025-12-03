use crate::g729::basic_operations::*;
// use crate::g729::ld8k::*;

/*****************************************************************************/
/* computePolynomialCoefficients : according to spec. 3.2.6                  */
/*    parameters:                                                            */
/*      -(i) qLSP : 10 LSP to be converted in Q0.15 range [-1, +1[           */
/*      -(o) f : 6 values in Q24 : polynomial coefficients on 32 bits        */
/*                                                                           */
/*****************************************************************************/
fn compute_polynomial_coefficients(q_lsp: &[Word16], f: &mut [Word32]) {
    /* init values */
    /* f[-1] which is not available is to be egal at 0, it is directly removed in the following when used as a factor */
    f[0] = 16777216; /* 1 in Q24, f[0] is used only for the other coefficient computation */

    /* correspont to i=1 in the algorithm description in spec. 3.2.6 */
    f[1] = mult16_16(q_lsp[0], -1024); /* f[1] = -2*qLSP[0], with qLSP in Q0.15 * -2*2^9(-1024) -> f[1] in Q1.24 */

    /* so we start with i=2 */
    /* Note : index of qLSP are -1 respect of what is in the spec because qLSP array is indexed from 0-9 and not 1-10 */
    for i in 2..6 {
        /* spec: f[i] = 2(f[i-2] - qLSP[2i-1]*f[i-1]) */
        f[i] = sshl(
            sub32(f[i - 2], mult16_32_p15(q_lsp[2 * i - 2], f[i - 1])),
            1,
        ); /* with qLSP in Q0.15 and f in Q24 */
        for j in (2..i).rev() {
            /* case of j=1 is just doing f[1] -= 2*qLSP[2i-1], done after the loop in order to avoid reference to f[-1] */
            /* spec: f[j] = f[j] -2*qLSP[2i-i]*f[j-1] + f[j-2] (all right terms refer to the value at the previous iteration on i indexed loop) */
            f[j] = add32(
                f[j],
                sub32(f[j - 2], mult16_32_p14(q_lsp[2 * i - 2], f[j - 1])),
            ); /* qLPS in Q0.15 and f in Q24, using MULT16_32_P14 instead of P15 does the *2 on qLSP. Result in Q24 */
        }
        /* f[1] -=  2*qLSP[2i-1] */
        f[1] = sub32(f[1], sshl(q_lsp[2 * i - 2] as Word32, 10)); /* qLSP in Q0.15, must be shift by 9 to get in Q24 and one more to be *2 */
    }
}

/*****************************************************************************/
/* qLSP2LP : convert qLSP into LP parameters according to spec. 3.2.6        */
/*    parameters:                                                            */
/*      -(i) qLSP : 10 LSP to be converted in Q0.15 range [-1, +1[           */
/*      -(o) LP : 10 LP coefficients in Q12                                  */
/*                                                                           */
/*****************************************************************************/
pub fn q_lsp_2_lp(q_lsp: &[Word16], lp: &mut [Word16]) {
    let mut f1 = [0 as Word32; 6];
    let mut f2 = [0 as Word32; 6]; /* define two buffer to store the polynomials coefficients (size is 6 for 5 coefficient because fx[0] is used during computation as a buffer) */

    compute_polynomial_coefficients(q_lsp, &mut f1);

    compute_polynomial_coefficients(&q_lsp[1..], &mut f2);

    /* convert f coefficient in f' coefficient f'1[i] = f1[i]+f[i-1] and f'2[i] = f2[i] - f2[i-1] */
    for i in (1..6).rev() {
        f1[i] = add32(f1[i], f1[i - 1]); /* f1 is still in Q24 */
        f2[i] = sub32(f2[i], f2[i - 1]); /* f1 is still in Q24 */
    }

    /* compute the LP coefficients */
    /* fx[0] is not needed anymore, acces the f1 and f2 buffers via the fp1 and fp2 pointers */
    /* pointers to access directly the elements 1 of f1 and f2 */
    /* In Rust we just index from 1 */
    for i in 0..5 {
        /* i in [0,5[ : LP[i] = (f1[i] + f2[i])/2 */
        lp[i] = pshr(add32(f1[i + 1], f2[i + 1]), 13) as Word16; /* f1 and f2 in Q24, LP in Q12 */
        /* i in [5,9[ : LP[i] = (f1[9-i] - f2[9-i])/2 */
        lp[9 - i] = pshr(sub32(f1[i + 1], f2[i + 1]), 13) as Word16; /* f1 and f2 in Q24, LP in Q12 */
    }
}
