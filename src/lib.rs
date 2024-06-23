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

mod consts;
mod gen;
mod impls;
mod kem;
mod matrix_arith;
mod pke;
mod ring_arith;
mod ser;

pub extern crate kem as kem_traits;
pub use impls::*;
