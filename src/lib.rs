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
