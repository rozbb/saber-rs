/// The modulus q of our base ring Z/2^13 Z
pub(crate) const MODULUS_Q: u16 = 1 << 13;
/// The bitlength of the modulus q
pub(crate) const MODULUS_Q_BITS: usize = 13;
/// The bitlength of the modulus p
pub(crate) const MODULUS_P_BITS: usize = 10;
/// The degree of the polynomial ring over Z/qZ
pub(crate) const RING_DEG: usize = 256;

// Key for the following constants:
//   * L is the dimension of the vectors/matrices we use.
//     E.g., a secret in Lightsaber is an element of R_q^2
//   * MODULUS_T_BITS is the modulus for some ciphertext values.
//     E.g., the c values in Lightsaber are in R_3
//   * MU is roughly the width of the binomial distribution we sample from.
//     E.g., in Lightsaber, we sample from a binomial distribution with width 10

pub(crate) const LIGHTSABER_L: usize = 2;
pub(crate) const LIGHTSABER_MODULUS_T_BITS: usize = 3;
pub(crate) const LIGHTSABER_MU: usize = 10;

pub(crate) const SABER_L: usize = 3;
pub(crate) const SABER_MODULUS_T_BITS: usize = 4;
pub(crate) const SABER_MU: usize = 8;

pub(crate) const FIRESABER_L: usize = 4;
pub(crate) const FIRESABER_MODULUS_T_BITS: usize = 6;
pub(crate) const FIRESABER_MU: usize = 6;

// We need to store maximum values because we can't do const arithmetic on some buffer sizes

/// The maximum possible μ value is μ=10, set by Lightsaber
pub(crate) const MAX_MU: usize = 10;
/// The maximum possible l value is l=4, set by FireSaber
pub(crate) const MAX_L: usize = 4;
// The maximum possible log(T) value is T=6, set by Firesaber
pub(crate) const MAX_MODULUS_T_BITS: usize = 6;
