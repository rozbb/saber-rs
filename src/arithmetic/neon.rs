//! NEON-accelerated polynomial arithmetic for aarch64.
//!
//! This module provides SIMD-accelerated versions of the core arithmetic operations
//! using ARM NEON intrinsics. All functions operate on `u16` coefficient arrays that
//! represent polynomials in the ring ℤ\[X\]/(X^256 + 1).
//!
//! The ring multiplication uses the Toom-Cook 4-way algorithm (matching the AVX2
//! implementation's evaluation/interpolation scheme), with NEON-vectorised 64×64
//! schoolbook at the leaf level. This reduces total multiply-accumulate operations
//! by ~42% compared to the previous Karatsuba + schoolbook-128 approach.
//!
//! NEON is part of the base AArch64 ISA, so no runtime feature detection is needed.

#![allow(unsafe_code)]

use core::arch::aarch64::*;

use super::RingElem;
use crate::consts::RING_DEG;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// One quarter of the ring degree — TC-4 splits a polynomial into 4 quarters.
const QUARTER: usize = RING_DEG / 4; // 64

/// 3^{-1} mod 2^{16} (3 × 43691 = 131073 ≡ 1 mod 65536).
const INV3: u16 = 43691;
/// 9^{-1} mod 2^{16}.
const INV9: u16 = 36409;
/// 15^{-1} mod 2^{16}.
const INV15: u16 = 61167;
const INT45: u16 = 45;
const INT30: u16 = 30;

// ---------------------------------------------------------------------------
// Toom-Cook 4-way evaluation
// ---------------------------------------------------------------------------

/// Evaluate a polynomial at the 7 Toom-Cook points.
///
/// Splits `poly[0..256]` into 4 quarters q0..q3 (each 64 coefficients) and
/// computes (in the same order as the AVX2 implementation):
///
///   aw\[0\] = P(∞) = q3
///   aw\[1\] = P(2) = q0 + 2q1 + 4q2 + 8q3
///   aw\[2\] = P(1) = q0 + q1 + q2 + q3
///   aw\[3\] = P(-1) = q0 - q1 + q2 - q3
///   aw\[4\] = 8·P(1/2) = 8q0 + 4q1 + 2q2 + q3
///   aw\[5\] = 8·P(-1/2) = 8q0 - 4q1 + 2q2 - q3
///   aw\[6\] = P(0) = q0
#[inline]
unsafe fn tc_eval(poly: &[u16; RING_DEG], aw: &mut [[u16; QUARTER]; 7]) {
    for i in (0..QUARTER).step_by(8) {
        let r0 = vld1q_u16(poly.as_ptr().add(i));
        let r1 = vld1q_u16(poly.as_ptr().add(i + QUARTER));
        let r2 = vld1q_u16(poly.as_ptr().add(i + 2 * QUARTER));
        let r3 = vld1q_u16(poly.as_ptr().add(i + 3 * QUARTER));

        let r4 = vaddq_u16(r0, r2); // q0 + q2
        let r5 = vaddq_u16(r1, r3); // q1 + q3

        // P(1) = q0 + q1 + q2 + q3
        vst1q_u16(aw[2].as_mut_ptr().add(i), vaddq_u16(r4, r5));
        // P(-1) = q0 - q1 + q2 - q3
        vst1q_u16(aw[3].as_mut_ptr().add(i), vsubq_u16(r4, r5));

        // 8·P(1/2): (4q0 + q2)·2 ± (4q1 + q3)
        let r4_half = vshlq_n_u16::<1>(vaddq_u16(vshlq_n_u16::<2>(r0), r2));
        let r5_half = vaddq_u16(vshlq_n_u16::<2>(r1), r3);
        vst1q_u16(aw[4].as_mut_ptr().add(i), vaddq_u16(r4_half, r5_half));
        vst1q_u16(aw[5].as_mut_ptr().add(i), vsubq_u16(r4_half, r5_half));

        // P(2) = q0 + 2q1 + 4q2 + 8q3
        let p2 = vaddq_u16(
            vaddq_u16(r0, vshlq_n_u16::<1>(r1)),
            vaddq_u16(vshlq_n_u16::<2>(r2), vshlq_n_u16::<3>(r3)),
        );
        vst1q_u16(aw[1].as_mut_ptr().add(i), p2);

        // P(0) = q0
        vst1q_u16(aw[6].as_mut_ptr().add(i), r0);
        // P(∞) = q3
        vst1q_u16(aw[0].as_mut_ptr().add(i), r3);
    }
}

// ---------------------------------------------------------------------------
// Schoolbook 64×64 (NEON-vectorised leaf)
// ---------------------------------------------------------------------------

/// NEON-accelerated schoolbook multiplication of two 64-coefficient polynomials.
///
/// Computes `c = a * b` where the product has at most degree 126, stored in
/// `c[0..128]`. Caller must ensure `c` is zeroed on entry.
#[inline]
unsafe fn schoolbook_64(c: &mut [u16; 2 * QUARTER], a: &[u16; QUARTER], b: &[u16; QUARTER]) {
    // Pre-load all 8 chunks of b into registers so they stay in the register
    // file across all 64 outer iterations, saving repeated loads from memory.
    let b0 = vld1q_u16(b.as_ptr());
    let b1 = vld1q_u16(b.as_ptr().add(8));
    let b2 = vld1q_u16(b.as_ptr().add(16));
    let b3 = vld1q_u16(b.as_ptr().add(24));
    let b4 = vld1q_u16(b.as_ptr().add(32));
    let b5 = vld1q_u16(b.as_ptr().add(40));
    let b6 = vld1q_u16(b.as_ptr().add(48));
    let b7 = vld1q_u16(b.as_ptr().add(56));

    for i in 0..QUARTER {
        let ai = vdupq_n_u16(*a.get_unchecked(i));
        let p = c.as_mut_ptr().add(i);
        vst1q_u16(p, vmlaq_u16(vld1q_u16(p as *const _), ai, b0));
        vst1q_u16(p.add(8), vmlaq_u16(vld1q_u16(p.add(8) as *const _), ai, b1));
        vst1q_u16(
            p.add(16),
            vmlaq_u16(vld1q_u16(p.add(16) as *const _), ai, b2),
        );
        vst1q_u16(
            p.add(24),
            vmlaq_u16(vld1q_u16(p.add(24) as *const _), ai, b3),
        );
        vst1q_u16(
            p.add(32),
            vmlaq_u16(vld1q_u16(p.add(32) as *const _), ai, b4),
        );
        vst1q_u16(
            p.add(40),
            vmlaq_u16(vld1q_u16(p.add(40) as *const _), ai, b5),
        );
        vst1q_u16(
            p.add(48),
            vmlaq_u16(vld1q_u16(p.add(48) as *const _), ai, b6),
        );
        vst1q_u16(
            p.add(56),
            vmlaq_u16(vld1q_u16(p.add(56) as *const _), ai, b7),
        );
    }
}

// ---------------------------------------------------------------------------
// Toom-Cook interpolation
// ---------------------------------------------------------------------------

/// Toom-Cook interpolation: recover the full product from 7 evaluation-point
/// products.
///
/// `w[k]` is the 128-coeff product at TC point `k`. Output: the 512-coeff
/// unreduced product in `result`.
///
/// The interpolation formulas match the AVX2 implementation (and the SABER
/// reference C code), expressed with NEON intrinsics on u16 coefficient arrays.
unsafe fn tc_interpol(w: &[[u16; 2 * QUARTER]; 7], result: &mut [u16; 2 * RING_DEG]) {
    let inv3_v = vdupq_n_u16(INV3);
    let inv9_v = vdupq_n_u16(INV9);
    let inv15_v = vdupq_n_u16(INV15);
    let int45_v = vdupq_n_u16(INT45);
    let int30_v = vdupq_n_u16(INT30);
    let zero_v = vdupq_n_u16(0);

    for i in (0..2 * QUARTER).step_by(8) {
        let r0 = vld1q_u16(w[0].as_ptr().add(i)); // product at ∞
        let mut r1 = vld1q_u16(w[1].as_ptr().add(i)); // P(2)·Q(2)
        let mut r2 = vld1q_u16(w[2].as_ptr().add(i)); // P(1)·Q(1)
        let mut r3 = vld1q_u16(w[3].as_ptr().add(i)); // P(-1)·Q(-1)
        let mut r4 = vld1q_u16(w[4].as_ptr().add(i)); // 64·P(½)·Q(½)
        let mut r5 = vld1q_u16(w[5].as_ptr().add(i)); // 64·P(-½)·Q(-½)
        let r6 = vld1q_u16(w[6].as_ptr().add(i)); // P(0)·Q(0)

        // Interpolation — identical arithmetic to the AVX2 path
        r1 = vaddq_u16(r1, r4);
        r5 = vsubq_u16(r5, r4);
        r3 = vsubq_u16(r3, r2);
        r3 = vshrq_n_u16::<1>(r3);
        r4 = vsubq_u16(r4, r0);
        r4 = vsubq_u16(r4, vshlq_n_u16::<6>(r6));
        r4 = vshlq_n_u16::<1>(r4);
        r4 = vaddq_u16(r4, r5);
        r2 = vaddq_u16(r2, r3);
        r1 = vsubq_u16(r1, vshlq_n_u16::<6>(r2));
        r1 = vsubq_u16(r1, r2);
        r2 = vsubq_u16(r2, r6);
        r2 = vsubq_u16(r2, r0);
        r1 = vaddq_u16(r1, vmulq_u16(r2, int45_v));
        r4 = vsubq_u16(r4, vshlq_n_u16::<3>(r2));
        r4 = vmulq_u16(r4, inv3_v);
        r4 = vshrq_n_u16::<3>(r4);
        r5 = vaddq_u16(r5, r1);
        r1 = vaddq_u16(r1, vshlq_n_u16::<4>(r3));
        r1 = vmulq_u16(r1, inv9_v);
        r1 = vshrq_n_u16::<1>(r1);
        r3 = vaddq_u16(r1, r3);
        r3 = vsubq_u16(zero_v, r3); // negate
        let temp_val = vmulq_u16(r1, int30_v);
        let temp2 = vsubq_u16(temp_val, r5);
        let temp3 = vmulq_u16(temp2, inv15_v);
        r5 = vshrq_n_u16::<2>(temp3);
        r2 = vsubq_u16(r2, r4);
        r1 = vsubq_u16(r1, r5);

        // Store into the 512-coeff result array.
        //
        // The 7 interpolated sections are spaced QUARTER apart. For i < QUARTER
        // (first half of each section) we assign directly; for i ≥ QUARTER
        // (second half) we add into the overlap with the previous section,
        // except r0 which starts a new non-overlapping tail and is assigned.
        let p = result.as_mut_ptr();
        if i < QUARTER {
            vst1q_u16(p.add(0 * QUARTER + i), r6);
            vst1q_u16(p.add(1 * QUARTER + i), r5);
            vst1q_u16(p.add(2 * QUARTER + i), r4);
            vst1q_u16(p.add(3 * QUARTER + i), r3);
            vst1q_u16(p.add(4 * QUARTER + i), r2);
            vst1q_u16(p.add(5 * QUARTER + i), r1);
            vst1q_u16(p.add(6 * QUARTER + i), r0);
        } else {
            vst1q_u16(
                p.add(0 * QUARTER + i),
                vaddq_u16(vld1q_u16(p.add(0 * QUARTER + i) as *const _), r6),
            );
            vst1q_u16(
                p.add(1 * QUARTER + i),
                vaddq_u16(vld1q_u16(p.add(1 * QUARTER + i) as *const _), r5),
            );
            vst1q_u16(
                p.add(2 * QUARTER + i),
                vaddq_u16(vld1q_u16(p.add(2 * QUARTER + i) as *const _), r4),
            );
            vst1q_u16(
                p.add(3 * QUARTER + i),
                vaddq_u16(vld1q_u16(p.add(3 * QUARTER + i) as *const _), r3),
            );
            vst1q_u16(
                p.add(4 * QUARTER + i),
                vaddq_u16(vld1q_u16(p.add(4 * QUARTER + i) as *const _), r2),
            );
            vst1q_u16(
                p.add(5 * QUARTER + i),
                vaddq_u16(vld1q_u16(p.add(5 * QUARTER + i) as *const _), r1),
            );
            // r0 (leading coefficient section) is a fresh assignment
            vst1q_u16(p.add(6 * QUARTER + i), r0);
        }
    }
}

// ---------------------------------------------------------------------------
// Core polynomial multiplication
// ---------------------------------------------------------------------------

/// NEON-accelerated ring multiply-accumulate using Toom-Cook 4-way.
///
/// Computes `acc += a * b` in the ring ℤ\[X\]/(X^256 + 1).
///
/// Algorithm:
///  1. TC-4 evaluation of both inputs at 7 points.
///  2. For each point: schoolbook multiply the two 64-coeff evaluations.
///  3. TC-4 interpolation to recover the 512-coeff unreduced product.
///  4. Reduce mod X^256+1 (subtract upper half) and accumulate.
pub(super) fn ring_mul_acc_neon(acc: &mut RingElem, a: &RingElem, b: &RingElem) {
    // Safety: guarded by cfg(target_arch = "aarch64"); NEON is always available.
    unsafe { ring_mul_acc_neon_impl(acc, a, b) }
}

unsafe fn ring_mul_acc_neon_impl(acc: &mut RingElem, a: &RingElem, b: &RingElem) {
    // 1. TC-4 evaluation
    let mut aw = [[0u16; QUARTER]; 7];
    let mut bw = [[0u16; QUARTER]; 7];
    tc_eval(&a.0, &mut aw);
    tc_eval(&b.0, &mut bw);

    // 2. Point-wise multiplication (schoolbook-64 at each evaluation point)
    let mut w = [[0u16; 2 * QUARTER]; 7];
    for k in 0..7 {
        schoolbook_64(&mut w[k], &aw[k], &bw[k]);
    }

    // 3. TC-4 interpolation → 512-coeff unreduced product
    let mut unreduced = [0u16; 2 * RING_DEG];
    tc_interpol(&w, &mut unreduced);

    // 4. Reduce mod (X^256 + 1) and accumulate: acc[i] += lo[i] - hi[i]
    for i in (0..RING_DEG).step_by(8) {
        let lo = vld1q_u16(unreduced.as_ptr().add(i));
        let hi = vld1q_u16(unreduced.as_ptr().add(i + RING_DEG));
        let cur = vld1q_u16(acc.0.as_ptr().add(i));
        vst1q_u16(acc.0.as_mut_ptr().add(i), vaddq_u16(cur, vsubq_u16(lo, hi)));
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

    #[test]
    fn neon_ring_mul_acc_matches_reference() {
        use crate::consts::MODULUS_Q_BITS;
        let mut rng = rand::rng();

        // The Toom-Cook 4-way interpolation uses right-shift divisions (÷2, ÷4, ÷8)
        // that are exact in ℤ but lose the top bits in ℤ/(2^16). The result is
        // correct modulo 2^MODULUS_Q_BITS (= q = 8192), which is all that SABER
        // requires. We compare only the low MODULUS_Q_BITS bits of each coefficient.
        let mask = (1u16 << MODULUS_Q_BITS) - 1; // 0x1FFF

        for _ in 0..50 {
            let a = RingElem::rand(&mut rng);
            let b = RingElem::rand(&mut rng);

            // Compute using NEON (Toom-Cook 4-way)
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

            for i in 0..RING_DEG {
                assert_eq!(
                    acc_neon.0[i] & mask,
                    acc_ref.0[i] & mask,
                    "ring_mul_acc mismatch at index {i}: neon={} ref={}",
                    acc_neon.0[i],
                    acc_ref.0[i],
                );
            }
        }
    }

    #[test]
    fn neon_ring_mul_acc_accumulates() {
        let mut rng = rand::rng();

        let a = RingElem::rand(&mut rng);
        let b = RingElem::rand(&mut rng);
        let c = RingElem::rand(&mut rng);
        let d = RingElem::rand(&mut rng);

        use crate::consts::MODULUS_Q_BITS;
        let mask = (1u16 << MODULUS_Q_BITS) - 1;

        // acc = a*b + c*d  (via two calls to ring_mul_acc)
        let mut acc = RingElem::default();
        ring_mul_acc_neon(&mut acc, &a, &b);
        ring_mul_acc_neon(&mut acc, &c, &d);

        // reference: compute a*b and c*d separately, then add
        let mut ab = RingElem::default();
        ring_mul_acc_neon(&mut ab, &a, &b);
        let mut cd = RingElem::default();
        ring_mul_acc_neon(&mut cd, &c, &d);

        for i in 0..RING_DEG {
            assert_eq!(
                acc.0[i] & mask,
                (ab.0[i].wrapping_add(cd.0[i])) & mask,
                "accumulation mismatch at index {i}"
            );
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
