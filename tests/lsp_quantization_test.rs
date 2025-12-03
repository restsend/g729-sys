#[cfg(test)]
mod tests {
    use g729_sys::g729::ld8k::*;
    use g729_sys::g729::lsp_quantization::*;

    #[test]
    fn test_init_lsp_quantization() {
        let mut previous_q_lsf = [[0; NB_LSP_COEFF]; MA_MAX_K];
        init_lsp_quantization(&mut previous_q_lsf);

        // Check if initialized correctly (checking first few values based on PREVIOUS_QLSF_INIT)
        // 2339, 4679, 7018...
        assert_eq!(previous_q_lsf[0][0], 2339);
        assert_eq!(previous_q_lsf[0][1], 4679);
        assert_eq!(previous_q_lsf[3][9], 23396);
    }

    #[test]
    fn test_lsp_quantization_runs() {
        let mut previous_q_lsf = [[0; NB_LSP_COEFF]; MA_MAX_K];
        init_lsp_quantization(&mut previous_q_lsf);

        let lsp_coefficients = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000];
        let mut q_lsp_coefficients = [0; NB_LSP_COEFF];
        let mut parameters = [0u16; 4];

        lsp_coefficients.iter().for_each(|&x| assert!(x >= 0)); // Basic sanity check

        lsp_quantization(
            &mut previous_q_lsf,
            &lsp_coefficients,
            &mut q_lsp_coefficients,
            &mut parameters,
        );

        // Just check that we got some output and didn't crash
        // The exact values would depend on the complex logic, but we can check ranges if we knew them.
        // For now, just ensuring it runs is a good step.
    }
}
