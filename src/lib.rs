#![no_std]
#![forbid(unsafe_code)]
#![warn(
    clippy::unwrap_used,
    missing_docs,
    rust_2018_idioms,
    unused_lifetimes,
    unused_qualifications
)]
#![doc = include_str!("../README.md")]

mod arithmetic;
mod consts;
mod gen;
mod impls;
mod kem;
mod pke;
mod ser;

pub use impls::*;

use turboshake::{
    digest::{ExtendableOutput, Update, XofReader},
    CTurboShake256,
};

/// Helper function that computes the 32-bytes digest the concatenation of the given
/// inputs using TurboSHAKE256 with the given domain separator `DS`
pub(crate) fn turboshake256_hash<const DS: u8>(input: &[&[u8]]) -> [u8; 32] {
    let mut hasher = CTurboShake256::<DS>::default();
    for chunk in input {
        hasher.update(chunk);
    }

    let mut out = [0u8; 32];
    let mut reader = hasher.finalize_xof();
    reader.read(&mut out);

    out
}
