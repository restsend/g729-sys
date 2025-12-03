// Constants and definitions
use crate::g729::basic_operations::{Word16, Word32};

pub const L_FRAME: usize = 80;
pub const L_SUBFRAME: usize = 40;
pub const NB_LSP_COEFF: usize = 10;
pub const L_LP_ANALYSIS_WINDOW: usize = 240;
pub const MAXIMUM_INT_PITCH_DELAY: usize = 143;
pub const L_PAST_EXCITATION: usize = 154;
pub const NB_PARAMETERS: usize = 15;

pub const L0_RANGE: usize = 2;
pub const L1_RANGE: usize = 128;
pub const L2_RANGE: usize = 32;
pub const L3_RANGE: usize = 32;
pub const MA_MAX_K: usize = 4;
pub const NOISE_L1_RANGE: usize = 32;
pub const NOISE_L2_RANGE: usize = 16;

pub const GAP1: Word16 = 10;
pub const GAP2: Word16 = 5;
pub const QLSF_MIN: Word16 = 40;
pub const QLSF_MAX: Word16 = 25681;
pub const MIN_QLSF_DISTANCE: Word16 = 321;

/* MAXINTXX define the maximum signed integer value on XX bits(2^(XX-1) - 1) */
/* used to check on overflows in fixed point mode */
pub const MAXINT16: Word16 = 0x7fff;
pub const MAX_INT16: Word16 = MAXINT16; // Alias
pub const MAXUINT16: u16 = 0xffff;
pub const MAXINT17: i32 = 0xffff;
pub const MAXINT28: i32 = 0x7ffffff;
pub const MAXINT29: i32 = 0xfffffff;
pub const MAX_INT29: i32 = MAXINT29; // Alias
pub const MININT32: Word32 = -0x80000000;
pub const MAXINT32: Word32 = 0x7fffffff;
pub const MAXINT64: i64 = 0x7fffffffffffffff;

/* several values used for inits */
pub const ONE_IN_Q31: Word32 = 0x7FFFFFFF;
pub const ONE_IN_Q30: Word32 = 0x40000000;
pub const ONE_IN_Q27: Word32 = 0x08000000;
pub const ONE_IN_Q15: Word32 = 0x00008000;
pub const ONE_IN_Q13: Word32 = 0x00002000;
pub const ONE_IN_Q12: Word32 = 0x00001000;
pub const ONE_IN_Q11: Word32 = 0x00000800;

pub const HALF_PI_Q13: Word16 = 12868;
pub const HALF_PI_Q15_32: Word32 = 51472;

/* 0.04*Pi + 1 and 0.92*Pi - 1 used by LSPQuantization */
pub const OO4PIPLUS1_IN_Q13: Word16 = 9221;
pub const O92PIMINUS1_IN_Q13: Word16 = 15485;
/* 1.2 in Q14 */
pub const ONE_POINT_2_IN_Q14: Word16 = 19661;
/* 0.7 in Q12 */
pub const O7_IN_Q12: Word16 = 2867;
/* 0.2 in Q14 */
pub const O2_IN_Q14: Word16 = 3277;
/* 0.2 in Q15 */
pub const O2_IN_Q15: Word16 = 6554;

/* weighted speech for open-loop pitch delay (spec A3.3.3) in Q15 0.75^(1..10)*/
pub const GAMMA_E1: Word16 = 24756;
pub const GAMMA_E2: Word16 = 18432;
pub const GAMMA_E3: Word16 = 13824;
pub const GAMMA_E4: Word16 = 10368;
pub const GAMMA_E5: Word16 = 7776;
pub const GAMMA_E6: Word16 = 5832;
pub const GAMMA_E7: Word16 = 4374;
pub const GAMMA_E8: Word16 = 3280;
pub const GAMMA_E9: Word16 = 2460;
pub const GAMMA_E10: Word16 = 1845;

/* pitch gain boundaries in Q14 */
pub const BOUNDED_PITCH_GAIN_MIN: Word16 = 3277;
pub const BOUNDED_PITCH_GAIN_MAX: Word16 = 13107;

/* post filters values defined in 4.2.2 in Q15 pow 1 to 10 */
pub const GAMMA_N1: Word16 = 18022;
pub const GAMMA_N2: Word16 = 9912;
pub const GAMMA_N3: Word16 = 5452;
pub const GAMMA_N4: Word16 = 2998;
pub const GAMMA_N5: Word16 = 1649;
pub const GAMMA_N6: Word16 = 907;
pub const GAMMA_N7: Word16 = 499;
pub const GAMMA_N8: Word16 = 274;
pub const GAMMA_N9: Word16 = 151;
pub const GAMMA_N10: Word16 = 83;

pub const GAMMA_D1: Word16 = 22938;
pub const GAMMA_D2: Word16 = 16056;
pub const GAMMA_D3: Word16 = 11239;
pub const GAMMA_D4: Word16 = 7868;
pub const GAMMA_D5: Word16 = 5507;
pub const GAMMA_D6: Word16 = 3855;
pub const GAMMA_D7: Word16 = 2699;
pub const GAMMA_D8: Word16 = 1889;
pub const GAMMA_D9: Word16 = 1322;
pub const GAMMA_D10: Word16 = 926;

/* post filter value GAMMA_T 0.8 in Q15 (spec A.4.2.3)*/
pub const GAMMA_T: Word16 = 26214;

/*** CNG ***/
/* 1/2sqrt*(40) in Q1.13 */
pub const GAUSSIAN_EXCITATION_COEFF_FACTOR: Word16 = 25905;
/* 0.75 in Q15 */
pub const COEFF_K: Word16 = 24576;

/*** VAD ***/
pub const LOG2_240_Q16: Word32 = 518186;
pub const INV_LOG2_10_Q15: Word16 = 9864;

/* buffer size for history on different values (defined in spec table B.1)*/
pub const NI: usize = 32;
pub const N0: usize = 128;

/*** DTX ***/
/* 1.12202 in Q20 */
pub const THRESHOLD3_IN_Q20: Word32 = 1176553;
/* 1.20226 in Q20 */
pub const THRESHOLD1_IN_Q20: Word32 = 1260661;
pub const CNG_DTX_RANDOM_SEED_INIT: u16 = 11111;
