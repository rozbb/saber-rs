//! Contains modules for ring and matrix arithmetic

mod matrix_arith;
#[cfg(all(feature = "neon", target_arch = "aarch64"))]
mod neon;
mod ring_arith;

// Export all the underlying types
pub(crate) use matrix_arith::*;
pub(crate) use ring_arith::*;
