use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;
use crate::g729::utils::{count_leading_zeros, dot_product, vec_mult_16_16};

/// Compute a diagonal of Phi values: start from Phi(39,j) and step Phi(38, j-1) down to Phi(39-j, 0)
/// Phi(i,j) = Phi(i+1,j+1) + h(39-i)*h(39-j)
#[cfg_attr(target_arch = "xtensa", inline(never))]
fn compute_phi_diagonal(
    j: isize,
    impulse_response: &[i16],
    phi: &mut [[i32; L_SUBFRAME]; L_SUBFRAME],
    phi_scaling: u16,
    sign: &[i16],
) {
    let j_orig = j as usize;
    let len = j_orig + 1;
    let delta = L_SUBFRAME - 1 - j_orig;

    let mut prod = [0i32; L_SUBFRAME];
    vec_mult_16_16(
        &impulse_response[0..len],
        &impulse_response[delta..delta + len],
        &mut prod[0..len],
    );

    let mut acc: i32 = 0;
    for (k, &p) in prod[0..len].iter().enumerate() {
        acc = add32(acc, p);
        let row = L_SUBFRAME - 1 - k;
        let col = j_orig - k;
        let mut s = if phi_scaling == 0 {
            acc
        } else {
            shr(acc, phi_scaling as u32)
        };

        s *= sign[row] as i32;
        s *= sign[col] as i32;
        phi[row][col] = s;
        phi[col][row] = s;
    }
}

/// computeImpulseResponseCorrelationMatrix: as in spec 3.8.1 eq51, eq56, eq57
///
/// # Arguments
///
/// * `impulse_response` - 40 values in Q12
/// * `correlation_signal` - 40 values in Q12 get absolute value of input as output as specified in spec 3.8.1
/// * `correlation_signal_sign` - 40 values of -1 or 1 : the sign of the input correlationSignal elements
/// * `phi` - a triangular matrix composed of Phi(i,j) in Q24
#[cfg_attr(target_arch = "xtensa", inline(never))]
fn compute_impulse_response_correlation_matrix(
    impulse_response: &[i16],
    correlation_signal: &mut [i16],
    correlation_signal_sign: &mut [i16],
    phi: &mut [[i32; L_SUBFRAME]; L_SUBFRAME],
) {
    let mut acc: i32 = 0;
    let mut phi_scaling: u16 = 0;

    // this diagonal must be divided by 2 according to spec 3.8.1 eq57
    let mut i_comp = L_SUBFRAME - 1;
    for &h in impulse_response.iter().take(L_SUBFRAME) {
        acc = mac16_16(acc, h, h);
        phi[i_comp][i_comp] = shr(acc, 1); // divide by 2: eq57
        i_comp = i_comp.saturating_sub(1);
    }

    if phi[0][0] > 0x6666666 {
        // bcg729 countLeadingZeros excludes the sign bit.
        let scaled = (phi[0][0].wrapping_shl(1)).wrapping_add(0x3333333);
        phi_scaling = (3 - count_leading_zeros(scaled) as i32) as u16;
        for (i, row) in phi.iter_mut().enumerate() {
            row[i] = shr(row[i], phi_scaling as u32);
        }
    }

    for i in 0..L_SUBFRAME {
        if correlation_signal[i] >= 0 {
            correlation_signal_sign[i] = 1;
        } else {
            correlation_signal_sign[i] = -1;
            correlation_signal[i] = -correlation_signal[i];
        }
    }

    for i in 0..8 {
        for j in 0..4 {
            compute_phi_diagonal(
                (5 * i + j) as isize,
                impulse_response,
                phi,
                phi_scaling,
                correlation_signal_sign,
            );
        }
    }
}

/// fixedCodebookSearch: compute fixed codebook parameters (codeword and sign)
///      compute also fixed codebook vector as in spec 3.8.1
///
/// # Arguments
///
/// * `target_signal` - 40 values as in spec A.3.6 in Q0
/// * `impulse_response` - 40 values as in spec A.3.5 in Q12
/// * `int_pitch_delay` - current integer pitch delay
/// * `last_quantized_adaptative_codebook_gain` - previous subframe pitch gain quantized in Q14
/// * `filtered_adaptative_codebook_vector` - 40 values in Q0
/// * `adaptative_codebook_gain` - in Q14
/// * `fixed_codebook_parameter` - Output fixed codebook parameter
/// * `fixed_codebook_pulses_signs` - Output fixed codebook pulses signs
/// * `fixed_codebook_vector` - Output 40 values as in spec 3.8, eq45 in Q13
/// * `fixed_codebook_vector_convolved` - Output 40 values as in spec 3.9, eq64 in Q12
#[allow(clippy::too_many_arguments)]
#[cfg_attr(target_arch = "xtensa", inline(never))]
pub fn fixed_codebook_search(
    target_signal: &[i16],
    impulse_response: &mut [i16],
    int_pitch_delay: i16,
    last_quantized_adaptative_codebook_gain: i16,
    filtered_adaptative_codebook_vector: &[i16],
    adaptative_codebook_gain: i16,
    fixed_codebook_parameter: &mut u16,
    fixed_codebook_pulses_signs: &mut u16,
    fixed_codebook_vector: &mut [i16],
    fixed_codebook_vector_convolved: &mut [i16],
) {
    let mut fixed_codebook_target_signal = [0_i16; L_SUBFRAME];
    let mut correlation_signal_32 = [0_i32; L_SUBFRAME]; // on 32 bits in Q12
    let mut correlation_signal = [0_i16; L_SUBFRAME];
    let mut correlation_signal_max: i32 = 0;

    let mut correlation_signal_sign_i16 = [0_i16; L_SUBFRAME];

    let mut phi = [[0_i32; L_SUBFRAME]; L_SUBFRAME];
    let mut i0 = 0;
    let mut i1 = 0;
    let mut i2 = 0;
    let mut i3 = 0;
    let mut correlation_square_max: i32 = -1;
    let mut energy_max: i32 = 1;
    let mut m0 = 0;
    let mut m1 = 0;
    let mut m2 = 0;
    let mut m3 = 0;
    let mut m_switch = [[2, 3, 0, 1], [3, 0, 1, 2]];
    let mut jx = 0;

    for i in 0..L_SUBFRAME {
        fixed_codebook_target_signal[i] = msu16_16_q14(
            target_signal[i] as i32,
            filtered_adaptative_codebook_vector[i],
            adaptative_codebook_gain,
        ) as i16; // adaptativeCodebookGain in Q14, other values in Q0
    }

    // update impulse vector as in spec 3.8 eq49
    for i in int_pitch_delay as usize..L_SUBFRAME {
        impulse_response[i] = mac16_16_q14(
            impulse_response[i] as i32,
            impulse_response[i - int_pitch_delay as usize],
            last_quantized_adaptative_codebook_gain,
        ) as i16;
    }

    // compute the correlation signal as in spec 3.8.1 eq52

    for n in 0..L_SUBFRAME {
        correlation_signal_32[n] = dot_product(
            &fixed_codebook_target_signal[n..L_SUBFRAME],
            &impulse_response[0..L_SUBFRAME - n],
        );
        let absc_correlation_signal_32 = if correlation_signal_32[n] >= 0 {
            correlation_signal_32[n]
        } else {
            -correlation_signal_32[n]
        };
        if absc_correlation_signal_32 > correlation_signal_max {
            correlation_signal_max = absc_correlation_signal_32;
        }
    }

    let correlation_signal_max_norm = count_leading_zeros(correlation_signal_max) as u32;

    if correlation_signal_max_norm < 18 {
        for i in 0..L_SUBFRAME {
            correlation_signal[i] =
                shr(correlation_signal_32[i], 18 - correlation_signal_max_norm) as i16;
        }
    } else {
        for i in 0..L_SUBFRAME {
            correlation_signal[i] = correlation_signal_32[i] as i16;
        }
    }

    compute_impulse_response_correlation_matrix(
        impulse_response,
        &mut correlation_signal,
        &mut correlation_signal_sign_i16,
        &mut phi,
    );

    // search for impulses leading to a max in C^2/E : spec 3.8.1 eq53

    let mut m3_base = 3;
    while m3_base < 5 {
        for (m_index, track) in m_switch.iter().enumerate() {
            // define for this loop on m3 track the Correlation and Energy giving the maximum of eq53
            let mut m3_track_correlation_square: i32 = -1;
            let mut m3_track_energy: i32 = 1;

            let mut first_m2 = 0;
            let mut correlation_m2_m3_max: i16 = 0;

            for _ in 0..2 {
                let mut correlation_m2: i16 = -1;
                let mut current_m2 = 0;

                let mut j = track[0];
                while j < L_SUBFRAME {
                    if correlation_signal[j] > correlation_m2 && j != first_m2 {
                        current_m2 = j;
                        correlation_m2 = correlation_signal[j];
                    }
                    j += 5;
                }
                first_m2 = current_m2;

                let energy_m2: i32 = phi[current_m2][current_m2]; // compute the energy with terms of eq55 using m2 only: Phi'(m2,m2)

                let mut j = track[1];
                while j < L_SUBFRAME {
                    let correlation_m2_m3 = add16(correlation_m2, correlation_signal[j]);
                    let energy_m2_m3 = add32(energy_m2, add32(phi[current_m2][j], phi[j][j])); // compute the energy if eq55 using term including m2 and m3: Phi'(m2,m2) is already in energyM2 + Phi'(m2,m3) + Phi'(m3,m3)
                    let correlation_m2_m3_square = mult16_16(correlation_m2_m3, correlation_m2_m3);

                    if mult32_32(m3_track_energy, correlation_m2_m3_square)
                        > mult32_32(energy_m2_m3, m3_track_correlation_square)
                    {
                        m3_track_correlation_square = correlation_m2_m3_square;
                        m3_track_energy = energy_m2_m3;
                        correlation_m2_m3_max = correlation_m2_m3;
                        m3 = j;
                        m2 = current_m2;
                    }
                    j += 5;
                }
            }
            let energy_m2_m3_max: i32 = m3_track_energy;

            m3_track_correlation_square = -1;
            m3_track_energy = 1;

            let mut i = track[2];
            while i < L_SUBFRAME {
                let correlation_m2_m3_m0 = add16(correlation_m2_m3_max, correlation_signal[i]);
                let energy_m2_m3_m0 = add32(
                    energy_m2_m3_max,
                    add32(phi[i][i], add32(phi[i][m2], phi[i][m3])),
                ); // add to the previously computed energy the terms of eq59 we can compute with the selected m0: Phi'(m0,m0) + Phi'(m0,m2) + Phi'(m0,m3)

                let mut j = track[3];
                while j < L_SUBFRAME {
                    let correlation_m2_m3_m0_m1 =
                        add16(correlation_m2_m3_m0, correlation_signal[j]);
                    let energy_m2_m3_m0_m1 = add32(
                        energy_m2_m3_m0,
                        add32(phi[j][i], add32(phi[j][j], add32(phi[j][m2], phi[j][m3]))),
                    ); // add to the previously computed energy the terms of eq59 we can compute with the selected m1: Phi'(m1,m0) + Phi'(m1,m1) + Phi'(m1,m2) + Phi'(m1,m3)
                    let correlation_m2_m3_m0_m1_square =
                        mult16_16(correlation_m2_m3_m0_m1, correlation_m2_m3_m0_m1);

                    if mult32_32(m3_track_energy, correlation_m2_m3_m0_m1_square)
                        > mult32_32(energy_m2_m3_m0_m1, m3_track_correlation_square)
                    {
                        m3_track_correlation_square = correlation_m2_m3_m0_m1_square;
                        m3_track_energy = energy_m2_m3_m0_m1;
                        m1 = j;
                        m0 = i;
                    }
                    j += 5;
                }
                i += 5;
            }

            if mult32_32(energy_max, m3_track_correlation_square)
                > mult32_32(m3_track_energy, correlation_square_max)
            {
                correlation_square_max = m3_track_correlation_square;
                energy_max = m3_track_energy;
                if m_index == 0 {
                    i0 = m0;
                    i1 = m1;
                    i2 = m2;
                    i3 = m3;
                } else {
                    i0 = m3;
                    i1 = m0;
                    i2 = m1;
                    i3 = m2;
                }
                jx = m3_base - 3; // needed for parameter computation apec 3.8.2 eq62
            }
        }
        m_switch[0][1] += 1;
        m_switch[1][0] += 1;
        m3_base += 1;
    }

    fixed_codebook_vector[..L_SUBFRAME].fill(0);

    // set the four pulses, in Q13
    fixed_codebook_vector[i0] = sshl(correlation_signal_sign_i16[i0] as i32, 13) as i16;
    fixed_codebook_vector[i1] = sshl(correlation_signal_sign_i16[i1] as i32, 13) as i16;
    fixed_codebook_vector[i2] = sshl(correlation_signal_sign_i16[i2] as i32, 13) as i16;
    fixed_codebook_vector[i3] = sshl(correlation_signal_sign_i16[i3] as i32, 13) as i16;

    // adapt it according to eq48
    for i in int_pitch_delay as usize..L_SUBFRAME {
        fixed_codebook_vector[i] = mac16_16_q14(
            fixed_codebook_vector[i] as i32,
            fixed_codebook_vector[i - int_pitch_delay as usize],
            last_quantized_adaptative_codebook_gain,
        ) as i16;
    }

    *fixed_codebook_parameter = (mult16_16_q15(i0 as i16, O2_IN_Q15)
        + ((mult16_16_q15(i1 as i16, O2_IN_Q15)) << 3)
        + ((mult16_16_q15(i2 as i16, O2_IN_Q15)) << 6)
        + ((((mult16_16_q15(i3 as i16, O2_IN_Q15)) << 1) + jx) << 9))
        as u16;

    *fixed_codebook_pulses_signs = (((correlation_signal_sign_i16[i0] + 1) >> 1) as u16)
        | ((((correlation_signal_sign_i16[i1] + 1) >> 1) << 1) as u16)
        | ((((correlation_signal_sign_i16[i2] + 1) >> 1) << 2) as u16)
        | ((((correlation_signal_sign_i16[i3] + 1) >> 1) << 3) as u16);

    // compute the fixedCodebook vector convolved with impulse response spec 3.9 eq64

    // eq64 make use of fixedCodebook vector adapted by eq48, using the impulse position(and thus fixed codebook vector before the adaptation)  but
    // the impulse response adapted as in eq49 gives the same output

    fixed_codebook_vector_convolved[..i0].fill(0);

    if correlation_signal_sign_i16[i0] > 0 {
        for (i, j) in (i0..L_SUBFRAME).zip(0..L_SUBFRAME) {
            fixed_codebook_vector_convolved[i] = impulse_response[j];
        }
    } else {
        for (i, j) in (i0..L_SUBFRAME).zip(0..L_SUBFRAME) {
            fixed_codebook_vector_convolved[i] = neg16(impulse_response[j]);
        }
    }

    if correlation_signal_sign_i16[i1] > 0 {
        for (i, j) in (i1..L_SUBFRAME).zip(0..L_SUBFRAME) {
            fixed_codebook_vector_convolved[i] =
                add16(fixed_codebook_vector_convolved[i], impulse_response[j]);
        }
    } else {
        for (i, j) in (i1..L_SUBFRAME).zip(0..L_SUBFRAME) {
            fixed_codebook_vector_convolved[i] =
                sub16(fixed_codebook_vector_convolved[i], impulse_response[j]);
        }
    }

    if correlation_signal_sign_i16[i2] > 0 {
        for (i, j) in (i2..L_SUBFRAME).zip(0..L_SUBFRAME) {
            fixed_codebook_vector_convolved[i] =
                add16(fixed_codebook_vector_convolved[i], impulse_response[j]);
        }
    } else {
        for (i, j) in (i2..L_SUBFRAME).zip(0..L_SUBFRAME) {
            fixed_codebook_vector_convolved[i] =
                sub16(fixed_codebook_vector_convolved[i], impulse_response[j]);
        }
    }

    if correlation_signal_sign_i16[i3] > 0 {
        for (i, j) in (i3..L_SUBFRAME).zip(0..L_SUBFRAME) {
            fixed_codebook_vector_convolved[i] =
                add16(fixed_codebook_vector_convolved[i], impulse_response[j]);
        }
    } else {
        for (i, j) in (i3..L_SUBFRAME).zip(0..L_SUBFRAME) {
            fixed_codebook_vector_convolved[i] =
                sub16(fixed_codebook_vector_convolved[i], impulse_response[j]);
        }
    }
}
