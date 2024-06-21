#[cfg(test)]
#[macro_use]
extern crate std;

mod consts;
mod gen;
mod impls;
pub mod kem;
mod matrix_arith;
mod pke;
mod ring_arith;
mod util;

// Known-answer tests
#[cfg(test)]
mod kat;

pub extern crate kem as kem_traits;
