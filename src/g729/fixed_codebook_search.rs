use crate::g729::basic_operations::*;
use crate::g729::ld8k::*;
use crate::g729::utils::{dot_product, vec_mult_16_16};

/// Compute a diagonal of Phi values: start from Phi(39,j) and step Phi(38, j-1) down to Phi(39-j, 0)
/// Phi(i,j) = Phi(i+1,j+1) + h(39-i)*h(39-j)
fn compute_phi_diagonal(
    j: isize,
    impulse_response: &[Word16],
    phi: &mut [[Word32; L_SUBFRAME]; L_SUBFRAME],
    phi_scaling: u16,
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

    let mut acc: Word32 = 0;
    if phi_scaling == 0 {
        for k in 0..len {
            acc = add32(acc, prod[k]);
            let row = L_SUBFRAME - 1 - k;
            let col = j_orig - k;
            phi[row][col] = acc;
        }
    } else {
        for k in 0..len {
            acc = add32(acc, prod[k]);
            let row = L_SUBFRAME - 1 - k;
            let col = j_orig - k;
            phi[row][col] = shr(acc, phi_scaling as u32);
        }
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
fn compute_impulse_response_correlation_matrix(
    impulse_response: &[Word16],
    correlation_signal: &mut [Word16],
    correlation_signal_sign: &mut [i16],
    phi: &mut [[Word32; L_SUBFRAME]; L_SUBFRAME],
) {
    let mut acc: Word32 = 0;
    let mut phi_scaling: u16 = 0;
    let mut correlation_signal_sign_inv = [0i16; L_SUBFRAME];

    // first compute the diagonal Phi(x,x) : Phi(39,39) = h[0]^2 # Phi(38,38) = Phi(39,39)+h[1]^2
    // this diagonal must be divided by 2 according to spec 3.8.1 eq57
    let mut i_comp = L_SUBFRAME - 1;
    for i in 0..L_SUBFRAME {
        acc = mac16_16(acc, impulse_response[i], impulse_response[i]); // impulseResponse in Q12 -> acc in Q24
        phi[i_comp][i_comp] = shr(acc, 1); // divide by 2: eq57
        if i_comp > 0 {
            i_comp -= 1;
        }
    }

    // check for possible overflow: Phi will be summed 10 times, so max Phi (by construction Phi[0][0]*2 is the max of Phi-> 2*Phi[0][0]*10 must be < 0x7fff ffff -> Phi[0][0]< 0x06666666 - otherwise scale Phi)
    if phi[0][0] > 0x6666666 {
        // complement 0xccccccc adding 0x3333333 to shift by one when max(2*Phi[0][0]) is in 0x0fffffff < max < 0xcccccc
        phi_scaling = (3 - ((phi[0][0] << 1) + 0x3333333).leading_zeros()) as u16;
        for i in 0..L_SUBFRAME {
            phi[i][i] = shr(phi[i][i], phi_scaling as u32);
        }
    }

    // Compute all diagonals but the 34, 29, 24, 19, 14, 9 and 4
    for i in 0..8 {
        for j in 0..4 {
            compute_phi_diagonal((5 * i + j) as isize, impulse_response, phi, phi_scaling);
        }
    }

    // correlationSignal -> absolute value and get sign (and his inverse in an array)
    for i in 0..L_SUBFRAME {
        if correlation_signal[i] >= 0 {
            correlation_signal_sign[i] = 1;
            correlation_signal_sign_inv[i] = -1;
        } else {
            // correlationSignal < 0
            correlation_signal_sign[i] = -1;
            correlation_signal_sign_inv[i] = 1;
            correlation_signal[i] = -correlation_signal[i];
        }
    }

    // modify the signs according to eq56
    for i in 0..L_SUBFRAME {
        let sign_of_correlation_signal_j: &[i16] = if correlation_signal_sign[i] > 0 {
            correlation_signal_sign
        } else {
            &correlation_signal_sign_inv
        };

        for j in 0..=i {
            // multiply by the selected sign the matrix element
            // Note : even the not needed and thus not computed elements are multiplicated... might found other way to do this sign stuff to be more efficient
            phi[i][j] = phi[i][j] * (sign_of_correlation_signal_j[j] as Word32);
        }
    }

    // duplicate the usefull values to their symetric part to get easier acces to the matrix elements
    for i in 0..8 {
        for j in 0..4 {
            let start_index = 5 * i + j;
            for k in 0..=start_index {
                phi[start_index - k][L_SUBFRAME - 1 - k] = phi[L_SUBFRAME - 1 - k][start_index - k];
            }
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
pub fn fixed_codebook_search(
    target_signal: &[Word16],
    impulse_response: &mut [Word16],
    int_pitch_delay: i16,
    last_quantized_adaptative_codebook_gain: Word16,
    filtered_adaptative_codebook_vector: &[Word16],
    adaptative_codebook_gain: Word16,
    fixed_codebook_parameter: &mut u16,
    fixed_codebook_pulses_signs: &mut u16,
    fixed_codebook_vector: &mut [Word16],
    fixed_codebook_vector_convolved: &mut [Word16],
) {
    let mut fixed_codebook_target_signal = [0 as Word16; L_SUBFRAME];
    let mut correlation_signal_32 = [0 as Word32; L_SUBFRAME]; // on 32 bits in Q12
    let mut correlation_signal = [0 as Word16; L_SUBFRAME]; // normalised to fit on 13 bits
    let mut correlation_signal_max: Word32 = 0;
    // Wait, compute_impulse_response_correlation_matrix expects i16 for sign.
    let mut correlation_signal_sign_i16 = [0 as i16; L_SUBFRAME];

    let mut phi = [[0 as Word32; L_SUBFRAME]; L_SUBFRAME];
    let mut i0 = 0;
    let mut i1 = 0;
    let mut i2 = 0;
    let mut i3 = 0;
    let mut correlation_square_max: Word32 = -1;
    let mut energy_max: Word32 = 1;
    let mut m0 = 0;
    let mut m1 = 0;
    let mut m2 = 0;
    let mut m3 = 0;
    let mut m_switch = [[2, 3, 0, 1], [3, 0, 1, 2]];
    let mut jx = 0;

    // compute the target signal for fixed codebook spec 3.8.1 eq50 : fixedCodebookTargetSignal[i] = targetSignal[i] - (adaptativeCodebookGain * filteredAdaptativeCodebookVector[i])
    for i in 0..L_SUBFRAME {
        fixed_codebook_target_signal[i] = msu16_16_q14(
            target_signal[i] as Word32,
            filtered_adaptative_codebook_vector[i],
            adaptative_codebook_gain,
        ) as Word16; // adaptativeCodebookGain in Q14, other values in Q0
    }

    // update impulse vector as in spec 3.8 eq49
    for i in int_pitch_delay as usize..L_SUBFRAME {
        impulse_response[i] = mac16_16_q14(
            impulse_response[i] as Word32,
            impulse_response[i - int_pitch_delay as usize],
            last_quantized_adaptative_codebook_gain,
        ) as Word16; // h[n] = h[n] + β*h[n-T], impulseResponse in Q12, lastQuantizedAdaptativeCodebookGain in Q14
    }

    // compute the correlation signal as in spec 3.8.1 eq52
    // compute on 32 bits and get the maximum
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

    // normalise on 13 bits
    // C uses a custom countLeadingZeros which excludes the sign bit.
    // Rust leading_zeros() includes the sign bit (which is 0 here since correlation_signal_max is abs).
    // So C_norm = Rust_norm - 1.
    let correlation_signal_max_norm = correlation_signal_max.leading_zeros().saturating_sub(1);

    if correlation_signal_max_norm < 18 {
        // if it doesn't already fit on 13 bits
        for i in 0..L_SUBFRAME {
            correlation_signal[i] =
                shr(correlation_signal_32[i], 18 - correlation_signal_max_norm) as Word16;
        }
    } else {
        // it fits on 13 bits, just copy it to the 16 bits buffer
        for i in 0..L_SUBFRAME {
            correlation_signal[i] = correlation_signal_32[i] as Word16;
        }
    }

    compute_impulse_response_correlation_matrix(
        impulse_response,
        &mut correlation_signal,
        &mut correlation_signal_sign_i16,
        &mut phi,
    );

    // search for impulses leading to a max in C^2/E : spec 3.8.1 eq53
    // ... (algorithm description omitted) ...

    let mut m3_base = 3;
    while m3_base < 5 {
        for m_index in 0..2 {
            // define for this loop on m3 track the Correlation and Energy giving the maximum of eq53
            let mut m3_track_correlation_square: Word32 = -1;
            let mut m3_track_energy: Word32 = 1;

            // Loop on the two maxima of correlation in the m2 index
            let mut first_m2 = 0; // save the first maximum index to not select it again
            let mut correlation_m2_m3_max: Word16 = 0; // stores the contribution of m2 and m3 impulses to the correlation for the maximum selected
            let energy_m2_m3_max: Word32; // same thing but for the energy

            for _ in 0..2 {
                let mut correlation_m2: Word16 = -1;
                let mut current_m2 = 0;
                let energy_m2: Word32;

                let mut j = m_switch[m_index][0];
                while j < L_SUBFRAME {
                    // in the m2 range, find the correlation Max -> select m2
                    if correlation_signal[j] > correlation_m2 && j != first_m2 {
                        current_m2 = j;
                        correlation_m2 = correlation_signal[j];
                    }
                    j += 5;
                }
                first_m2 = current_m2; // to avoid selecting the same maximum at next iteration

                energy_m2 = phi[current_m2][current_m2]; // compute the energy with terms of eq55 using m2 only: Phi'(m2,m2)

                // with selected m2, test the 8 m3 possibilities for the current m3 track
                let mut j = m_switch[m_index][1];
                while j < L_SUBFRAME {
                    let correlation_m2_m3 = add16(correlation_m2, correlation_signal[j]); // compute the correlation sum due to m2 and m3 pulses
                    let energy_m2_m3 = add32(energy_m2, add32(phi[current_m2][j], phi[j][j])); // compute the energy if eq55 using term including m2 and m3: Phi'(m2,m2) is already in energyM2 + Phi'(m2,m3) + Phi'(m3,m3)
                    let correlation_m2_m3_square = mult16_16(correlation_m2_m3, correlation_m2_m3);
                    // check if the current correlation/energy couple gives better results than the stored one : maximise C^2/E -> C^2/E > C^2max/Emax => Emax*C^2 > C^2max*E
                    if mult32_32(m3_track_energy, correlation_m2_m3_square) as i64
                        > mult32_32(energy_m2_m3, m3_track_correlation_square) as i64
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
            energy_m2_m3_max = m3_track_energy;

            // reset the current m3 track correlationSquare and energy
            m3_track_correlation_square = -1;
            m3_track_energy = 1;

            let mut i = m_switch[m_index][2];
            while i < L_SUBFRAME {
                // test the 8 possibilities for m0 track
                let correlation_m2_m3_m0 = add16(correlation_m2_m3_max, correlation_signal[i]); // compute correlation with current m0 taking in account the previously selected m2 and m3
                let energy_m2_m3_m0 = add32(
                    energy_m2_m3_max,
                    add32(phi[i][i], add32(phi[i][m2], phi[i][m3])),
                ); // add to the previously computed energy the terms of eq59 we can compute with the selected m0: Phi'(m0,m0) + Phi'(m0,m2) + Phi'(m0,m3)

                let mut j = m_switch[m_index][3];
                while j < L_SUBFRAME {
                    // test the 8 possibilities for m1 track
                    let correlation_m2_m3_m0_m1 =
                        add16(correlation_m2_m3_m0, correlation_signal[j]); // compute correlation with current m1 taking in account the previously selected m2, m3 and m0
                    let energy_m2_m3_m0_m1 = add32(
                        energy_m2_m3_m0,
                        add32(phi[j][i], add32(phi[j][j], add32(phi[j][m2], phi[j][m3]))),
                    ); // add to the previously computed energy the terms of eq59 we can compute with the selected m1: Phi'(m1,m0) + Phi'(m1,m1) + Phi'(m1,m2) + Phi'(m1,m3)
                    let correlation_m2_m3_m0_m1_square =
                        mult16_16(correlation_m2_m3_m0_m1, correlation_m2_m3_m0_m1);
                    // check if the current correlation/energy couple gives better results than the stored one : maximise C^2/E -> C^2/E > C^2max/Emax => Emax*C^2 > C^2max*E
                    if mult32_32(m3_track_energy, correlation_m2_m3_m0_m1_square) as i64
                        > mult32_32(energy_m2_m3_m0_m1, m3_track_correlation_square) as i64
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

            // check with currently selected indexes if this one is better
            if mult32_32(energy_max, m3_track_correlation_square) as i64
                > mult32_32(m3_track_energy, correlation_square_max) as i64
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
        m_switch[1][0] += 1; // increment the m3Base into the mSwitch
        m3_base += 1;
    }

    // compute the fixedCodebookVector
    for i in 0..L_SUBFRAME {
        fixed_codebook_vector[i] = 0; // reset the vector
    }

    // set the four pulses, in Q13
    fixed_codebook_vector[i0] = sshl(correlation_signal_sign_i16[i0] as Word32, 13) as Word16;
    fixed_codebook_vector[i1] = sshl(correlation_signal_sign_i16[i1] as Word32, 13) as Word16;
    fixed_codebook_vector[i2] = sshl(correlation_signal_sign_i16[i2] as Word32, 13) as Word16;
    fixed_codebook_vector[i3] = sshl(correlation_signal_sign_i16[i3] as Word32, 13) as Word16;

    // adapt it according to eq48
    for i in int_pitch_delay as usize..L_SUBFRAME {
        fixed_codebook_vector[i] = mac16_16_q14(
            fixed_codebook_vector[i] as Word32,
            fixed_codebook_vector[i - int_pitch_delay as usize],
            last_quantized_adaptative_codebook_gain,
        ) as Word16; // h[n] = h[n] + β*h[n-T], fixedCodebookVector in Q13, lastQuantizedAdaptativeCodebookGain in Q14
    }

    // compute the parameters
    *fixed_codebook_parameter = (mult16_16_q15(i0 as Word16, O2_IN_Q15)
        + ((mult16_16_q15(i1 as Word16, O2_IN_Q15)) << 3)
        + ((mult16_16_q15(i2 as Word16, O2_IN_Q15)) << 6)
        + ((((mult16_16_q15(i3 as Word16, O2_IN_Q15)) << 1) + jx as Word32) << 9))
        as u16;

    *fixed_codebook_pulses_signs = (((correlation_signal_sign_i16[i0] + 1) >> 1) as u16)
        | ((((correlation_signal_sign_i16[i1] + 1) >> 1) << 1) as u16)
        | ((((correlation_signal_sign_i16[i2] + 1) >> 1) << 2) as u16)
        | ((((correlation_signal_sign_i16[i3] + 1) >> 1) << 3) as u16);

    // compute the fixedCodebook vector convolved with impulse response spec 3.9 eq64
    // this vector is used in gain quantization but computed here because it's faster doing it having directly the impulses positions
    // eq64 make use of fixedCodebook vector adapted by eq48, using the impulse position(and thus fixed codebook vector before the adaptation)  but
    // the impulse response adapted as in eq49 gives the same output
    // reset the vector
    for i in 0..i0 {
        fixed_codebook_vector_convolved[i] = 0;
    }

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
