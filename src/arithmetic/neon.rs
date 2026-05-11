//! NEON-accelerated polynomial arithmetic for aarch64.
//!
//! This module provides SIMD-accelerated versions of the core arithmetic operations
//! using ARM NEON intrinsics. All functions operate on `u16` coefficient arrays that
//! represent polynomials in the ring ℤ\[X\]/(X^256 + 1).
//!
//! NEON is part of the base AArch64 ISA, so no runtime feature detection is needed.

#![allow(unsafe_code)]

use core::arch::aarch64::*;

use crate::consts::RING_DEG;

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
    use super::super::RingElem;
    use super::*;

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
