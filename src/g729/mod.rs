// G.729 Rust implementation

pub mod adaptative_codebook_search;
pub mod basic_operations;
pub mod codebooks;
pub mod compute_adaptative_codebook_gain;
pub mod compute_lp;
pub mod compute_weighted_speech;
pub mod encoder;
pub mod find_open_loop_pitch_delay;
pub mod fixed_codebook_search;
pub mod fixed_point_math;
pub mod interpolate_q_lsp;
pub mod ld8k;
pub mod lp2lsp_conversion;
pub mod lp_synthesis_filter;
pub mod lsp_quantization;
pub mod pre_processing;
pub mod q_lsp_2_lp;
pub mod utils;

#[cfg(test)]
mod encoder_test;

pub mod cng;
pub mod decode_adaptative_code_vector;
pub mod decode_fixed_code_vector;
pub mod decode_gains;
pub mod decode_lsp;
pub mod decoder;
pub mod gain_quantization;
pub mod post_filter;
pub mod post_processing;
