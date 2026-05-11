//! NEON-accelerated polynomial arithmetic for aarch64.
//!
//! This module provides SIMD-accelerated versions of the core arithmetic operations
//! using ARM NEON intrinsics. All functions operate on `u16` coefficient arrays that
//! represent polynomials in the ring ℤ\[X\]/(X^256 + 1).
//!
//! NEON is part of the base AArch64 ISA, so no runtime feature detection is needed.

#![allow(unsafe_code)]

use core::arch::aarch64::*;

use super::RingElem;
use crate::consts::RING_DEG;

/// Half the ring degree — we split 256-coeff polys into two 128-coeff halves for Karatsuba.
const HALF: usize = RING_DEG / 2;

// ---------------------------------------------------------------------------
// Core polynomial multiplication
// ---------------------------------------------------------------------------

/// NEON-accelerated schoolbook multiplication of two 128-coefficient polynomials.
///
/// Computes `out += a * b` where `a` and `b` are degree-127 polynomials.
/// `out` must be zeroed on entry and have room for 256 coefficients.
///
/// Strategy: broadcast each `a[i]` to all 8 lanes, then sweep across `b` in
/// chunks of 8, using `vmlaq_u16` (fused multiply-accumulate) to process 8
/// products per iteration.
#[inline(never)]
fn schoolbook_128_neon(out: &mut [u16; RING_DEG], a: &[u16; HALF], b: &[u16; HALF]) {
    // Safety: guarded by cfg(target_arch = "aarch64"); NEON is always available on AArch64.
    // All pointer offsets are within bounds: i in 0..128, j in 0..128 step 8,
    // so i+j+7 <= 127+127 = 254 < 256 = RING_DEG.
    unsafe {
        for i in 0..HALF {
            let ai = vdupq_n_u16(a[i]);
            for j in (0..HALF).step_by(8) {
                let bv = vld1q_u16(b.as_ptr().add(j));
                let ov = vld1q_u16(out.as_ptr().add(i + j));
                let res = vmlaq_u16(ov, ai, bv);
                vst1q_u16(out.as_mut_ptr().add(i + j), res);
            }
        }
    }
}

/// NEON-accelerated ring multiply-accumulate using one level of Karatsuba.
///
/// Computes `acc += a * b` in the ring ℤ\[X\]/(X^256 + 1). See the scalar
/// [`super::ring_arith::ring_mul_acc`] for a detailed explanation of the
/// Karatsuba decomposition; this function applies the same algorithm with
/// NEON-vectorised inner loops.
pub(super) fn ring_mul_acc_neon(acc: &mut RingElem, a: &RingElem, b: &RingElem) {
    let a_lo: &[u16; HALF] = a.0[..HALF].try_into().unwrap();
    let a_hi: &[u16; HALF] = a.0[HALF..].try_into().unwrap();
    let b_lo: &[u16; HALF] = b.0[..HALF].try_into().unwrap();
    let b_hi: &[u16; HALF] = b.0[HALF..].try_into().unwrap();

    // Three sub-products
    let mut z0 = [0u16; RING_DEG];
    let mut z2 = [0u16; RING_DEG];
    schoolbook_128_neon(&mut z0, a_lo, b_lo);
    schoolbook_128_neon(&mut z2, a_hi, b_hi);

    // Cross term: (a_lo + a_hi) * (b_lo + b_hi)
    let mut a_sum = [0u16; HALF];
    let mut b_sum = [0u16; HALF];

    // Safety: all accesses within bounds; HALF = 128 is divisible by 8.
    unsafe {
        for i in (0..HALF).step_by(8) {
            let al = vld1q_u16(a_lo.as_ptr().add(i));
            let ah = vld1q_u16(a_hi.as_ptr().add(i));
            vst1q_u16(a_sum.as_mut_ptr().add(i), vaddq_u16(al, ah));

            let bl = vld1q_u16(b_lo.as_ptr().add(i));
            let bh = vld1q_u16(b_hi.as_ptr().add(i));
            vst1q_u16(b_sum.as_mut_ptr().add(i), vaddq_u16(bl, bh));
        }
    }

    let mut z3 = [0u16; RING_DEG];
    schoolbook_128_neon(&mut z3, &a_sum, &b_sum);

    // Accumulate into acc, reducing mod (X^256 + 1).
    //
    // For j in 0..HALF:
    //   acc[j]      += z0[j] - z2[j] - (z3[j+H] - z0[j+H] - z2[j+H])
    //   acc[j+HALF] += z0[j+H] - z2[j+H] + (z3[j] - z0[j] - z2[j])
    //
    // Safety: all accesses within bounds; HALF = 128 is divisible by 8.
    unsafe {
        for j in (0..HALF).step_by(8) {
            // --- low half ---
            let acc_lo = vld1q_u16(acc.0.as_ptr().add(j));
            let z0_lo = vld1q_u16(z0.as_ptr().add(j));
            let z2_lo = vld1q_u16(z2.as_ptr().add(j));
            let z3_hi = vld1q_u16(z3.as_ptr().add(j + HALF));
            let z0_hi = vld1q_u16(z0.as_ptr().add(j + HALF));
            let z2_hi = vld1q_u16(z2.as_ptr().add(j + HALF));

            let z1_wrap = vsubq_u16(vsubq_u16(z3_hi, z0_hi), z2_hi);
            let res_lo = vsubq_u16(vsubq_u16(vaddq_u16(acc_lo, z0_lo), z2_lo), z1_wrap);
            vst1q_u16(acc.0.as_mut_ptr().add(j), res_lo);

            // --- high half ---
            let acc_hi = vld1q_u16(acc.0.as_ptr().add(j + HALF));
            let z3_lo = vld1q_u16(z3.as_ptr().add(j));

            let z1_direct = vsubq_u16(vsubq_u16(z3_lo, z0_lo), z2_lo);
            let res_hi = vaddq_u16(vsubq_u16(vaddq_u16(acc_hi, z0_hi), z2_hi), z1_direct);
            vst1q_u16(acc.0.as_mut_ptr().add(j + HALF), res_hi);
        }
    }
}

// ---------------------------------------------------------------------------
// Element-wise operations
// ---------------------------------------------------------------------------

/// NEON-accelerated right shift of all coefficients.
pub(super) fn shift_right_neon(coeffs: &mut [u16; RING_DEG], shift: usize) {
    // Safety: RING_DEG = 256 is divisible by 8; shift fits in i16 for valid Saber params.
    unsafe {
        let sh = vdupq_n_s16(-(shift as i16));
        for i in (0..RING_DEG).step_by(8) {
            let v = vld1q_u16(coeffs.as_ptr().add(i));
            vst1q_u16(coeffs.as_mut_ptr().add(i), vshlq_u16(v, sh));
        }
    }
}

/// NEON-accelerated left shift of all coefficients.
pub(super) fn shift_left_neon(coeffs: &mut [u16; RING_DEG], shift: usize) {
    // Safety: same as shift_right_neon.
    unsafe {
        let sh = vdupq_n_s16(shift as i16);
        for i in (0..RING_DEG).step_by(8) {
            let v = vld1q_u16(coeffs.as_ptr().add(i));
            vst1q_u16(coeffs.as_mut_ptr().add(i), vshlq_u16(v, sh));
        }
    }
}

/// NEON-accelerated wrapping add of a constant to all coefficients.
pub(super) fn wrapping_add_to_all_neon(coeffs: &mut [u16; RING_DEG], val: u16) {
    // Safety: RING_DEG = 256 is divisible by 8.
    unsafe {
        let vv = vdupq_n_u16(val);
        for i in (0..RING_DEG).step_by(8) {
            let v = vld1q_u16(coeffs.as_ptr().add(i));
            vst1q_u16(coeffs.as_mut_ptr().add(i), vaddq_u16(v, vv));
        }
    }
}

/// NEON-accelerated fused add-then-right-shift of all coefficients.
pub(super) fn wrapping_add_and_shift_right_neon(
    coeffs: &mut [u16; RING_DEG],
    val: u16,
    shift: usize,
) {
    // Safety: RING_DEG = 256 is divisible by 8; shift fits in i16.
    unsafe {
        let vv = vdupq_n_u16(val);
        let sh = vdupq_n_s16(-(shift as i16));
        for i in (0..RING_DEG).step_by(8) {
            let v = vld1q_u16(coeffs.as_ptr().add(i));
            let added = vaddq_u16(v, vv);
            vst1q_u16(coeffs.as_mut_ptr().add(i), vshlq_u16(added, sh));
        }
    }
}

/// NEON-accelerated element-wise wrapping addition: `out = a + b`.
pub(super) fn add_neon(a: &[u16; RING_DEG], b: &[u16; RING_DEG], out: &mut [u16; RING_DEG]) {
    // Safety: RING_DEG = 256 is divisible by 8.
    unsafe {
        for i in (0..RING_DEG).step_by(8) {
            let av = vld1q_u16(a.as_ptr().add(i));
            let bv = vld1q_u16(b.as_ptr().add(i));
            vst1q_u16(out.as_mut_ptr().add(i), vaddq_u16(av, bv));
        }
    }
}

/// NEON-accelerated element-wise wrapping subtraction: `out = a - b`.
pub(super) fn sub_neon(a: &[u16; RING_DEG], b: &[u16; RING_DEG], out: &mut [u16; RING_DEG]) {
    // Safety: RING_DEG = 256 is divisible by 8.
    unsafe {
        for i in (0..RING_DEG).step_by(8) {
            let av = vld1q_u16(a.as_ptr().add(i));
            let bv = vld1q_u16(b.as_ptr().add(i));
            vst1q_u16(out.as_mut_ptr().add(i), vsubq_u16(av, bv));
        }
    }
}

// ---------------------------------------------------------------------------
// Tests — cross-validate NEON against trivial reference implementations
// ---------------------------------------------------------------------------

#[cfg(test)]
mod test {
    use super::*;

    /// Trivial reference schoolbook — used only to validate the NEON version.
    fn reference_schoolbook_128(out: &mut [u16; RING_DEG], a: &[u16; HALF], b: &[u16; HALF]) {
        for i in 0..HALF {
            for j in 0..HALF {
                out[i + j] = out[i + j].wrapping_add(a[i].wrapping_mul(b[j]));
            }
        }
    }

    #[test]
    fn neon_schoolbook_matches_reference() {
        use rand::RngCore;
        let mut rng = rand::rng();

        for _ in 0..50 {
            let mut a = [0u16; HALF];
            let mut b = [0u16; HALF];
            for x in a.iter_mut() {
                *x = (rng.next_u32() & 0x1FFF) as u16; // 13-bit coefficients like Saber
            }
            for x in b.iter_mut() {
                *x = (rng.next_u32() & 0x1FFF) as u16;
            }

            let mut out_ref = [0u16; RING_DEG];
            let mut out_neon = [0u16; RING_DEG];
            reference_schoolbook_128(&mut out_ref, &a, &b);
            schoolbook_128_neon(&mut out_neon, &a, &b);
            assert_eq!(out_ref, out_neon, "schoolbook_128 mismatch");
        }
    }

    #[test]
    fn neon_ring_mul_acc_matches_reference() {
        let mut rng = rand::rng();

        for _ in 0..50 {
            let a = RingElem::rand(&mut rng);
            let b = RingElem::rand(&mut rng);

            // Compute using NEON (Karatsuba + NEON kernels)
            let mut acc_neon = RingElem::default();
            ring_mul_acc_neon(&mut acc_neon, &a, &b);

            // Compute using a completely independent full-degree reference schoolbook,
            // reduced mod X^256 + 1
            let mut product = [0u16; 2 * RING_DEG];
            for i in 0..RING_DEG {
                for j in 0..RING_DEG {
                    product[i + j] = product[i + j].wrapping_add(a.0[i].wrapping_mul(b.0[j]));
                }
            }
            let mut acc_ref = RingElem::default();
            for i in 0..RING_DEG {
                acc_ref.0[i] = product[i].wrapping_sub(product[i + RING_DEG]);
            }

            assert_eq!(acc_neon.0, acc_ref.0, "ring_mul_acc mismatch");
        }
    }

    #[test]
    fn neon_add_sub_roundtrip() {
        let mut rng = rand::rng();
        let a = RingElem::rand(&mut rng);
        let b = RingElem::rand(&mut rng);

        let mut sum = [0u16; RING_DEG];
        add_neon(&a.0, &b.0, &mut sum);

        let mut diff = [0u16; RING_DEG];
        sub_neon(&sum, &b.0, &mut diff);

        assert_eq!(diff, a.0, "add then sub should round-trip");
    }

    #[test]
    fn neon_shift_roundtrip() {
        let mut rng = rand::rng();
        let orig = RingElem::rand(&mut rng);

        // shift right 3 then left 3: the low 3 bits of each coefficient are lost
        let mut coeffs = orig.0;
        shift_right_neon(&mut coeffs, 3);
        shift_left_neon(&mut coeffs, 3);

        for i in 0..RING_DEG {
            assert_eq!(
                coeffs[i],
                orig.0[i] & !(0b111), // low 3 bits lost by right shift
                "shift round-trip mismatch at index {i}"
            );
        }
    }

    #[test]
    fn neon_wrapping_add_and_shift_right_matches_split() {
        use rand::RngCore;
        let mut rng = rand::rng();

        let mut coeffs_fused = [0u16; RING_DEG];
        let mut coeffs_split = [0u16; RING_DEG];
        for x in coeffs_fused.iter_mut() {
            *x = rng.next_u32() as u16;
        }
        coeffs_split.copy_from_slice(&coeffs_fused);

        let val = 0x0100u16;
        let shift = 3;

        // Fused version
        wrapping_add_and_shift_right_neon(&mut coeffs_fused, val, shift);

        // Two-step version
        wrapping_add_to_all_neon(&mut coeffs_split, val);
        shift_right_neon(&mut coeffs_split, shift);

        assert_eq!(coeffs_fused, coeffs_split, "fused vs split mismatch");
    }
}
