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

/// Helper function that computes the 32-bytes digest of the concatenation of the given inputs
/// using TurboSHAKE256 with the given domain separator `DS` Pass an empty slice for `input1` to
/// hash a single input.
// Note: we cannot take a `&[&[u8]]` because that's a nested borrow, which aeneas doesn't support
// yet.
pub(crate) fn turboshake256_hash<const DS: u8>(input0: &[u8], input1: &[u8]) -> [u8; 32] {
    let mut hasher = CTurboShake256::<DS>::default();
    hasher.update(input0);
    hasher.update(input1);

    let mut out = [0u8; 32];
    let mut reader = hasher.finalize_xof();
    reader.read(&mut out);

    out
}
