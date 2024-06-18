#[cfg(test)]
#[macro_use]
extern crate std;

mod consts;
mod gen;
pub mod kem;
mod matrix_arith;
mod pke;
mod ring_arith;
mod util;

// Known-answer tests
#[cfg(test)]
mod kat;

pub use pke::{ciphertext_len, max_ciphertext_len};
