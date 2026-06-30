//! AVX2-accelerated polynomial arithmetic for x86_64.
//!
//! This module provides SIMD-accelerated versions of the core arithmetic operations
//! using AVX2 intrinsics. All functions operate on `u16` coefficient arrays that
//! represent polynomials in the ring ℤ\[X\]/(X^256 + 1).
//!
//! AVX2 processes 16 × u16 values per instruction (256-bit vectors), double the
//! width of ARM NEON (8 × u16).

#![allow(unsafe_code)]

#[cfg(target_arch = "x86")]
use core::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::*;

use super::RingElem;
use crate::consts::RING_DEG;

/// Half the ring degree — we split 256-coeff polys into two 128-coeff halves for Karatsuba.
const HALF: usize = RING_DEG / 2;

// ---------------------------------------------------------------------------
// Core polynomial multiplication
// ---------------------------------------------------------------------------

/// AVX2-accelerated schoolbook multiplication of two 128-coefficient polynomials.
///
/// Computes `out += a * b` where `a` and `b` are degree-127 polynomials.
/// `out` must be zeroed on entry and have room for 256 coefficients.
///
/// Strategy: broadcast each `a[i]` to all 16 lanes, then sweep across `b` in
/// chunks of 16, using multiply-accumulate to process 16 products per iteration.
///
/// # Safety
///
/// Requires AVX2 support on the executing CPU.
#[target_feature(enable = "avx2")]
#[inline(never)]
unsafe fn schoolbook_128(out: &mut [u16; RING_DEG], a: &[u16; HALF], b: &[u16; HALF]) {
    // i in 0..128, j in {0,16,32,...,112}, so i+j+15 <= 127+127 = 254 < 256 = RING_DEG.
    for i in 0..HALF {
        let ai = _mm256_set1_epi16(a[i] as i16);
        for j in (0..HALF).step_by(16) {
            let bv = _mm256_loadu_si256(b.as_ptr().add(j) as *const __m256i);
            let ov = _mm256_loadu_si256(out.as_ptr().add(i + j) as *const __m256i);
            let res = _mm256_add_epi16(ov, _mm256_mullo_epi16(ai, bv));
            _mm256_storeu_si256(out.as_mut_ptr().add(i + j) as *mut __m256i, res);
        }
    }
}

/// AVX2-accelerated ring multiply-accumulate using one level of Karatsuba.
///
/// Computes `acc += a * b` in the ring ℤ\[X\]/(X^256 + 1). See the scalar
/// [`super::ring_arith::ring_mul_acc`] for a detailed explanation of the
/// Karatsuba decomposition; this function applies the same algorithm with
/// AVX2-vectorised inner loops.
pub(super) fn ring_mul_acc_avx2(acc: &mut RingElem, a: &RingElem, b: &RingElem) {
    // Safety: module gated by cfg(all(feature = "avx2", target_arch = "x86_64")).
    // Callers who enable the avx2 feature assert their hardware supports it.
    unsafe { ring_mul_acc_avx2_impl(acc, a, b) }
}

/// # Safety
///
/// Requires AVX2 support on the executing CPU.
#[target_feature(enable = "avx2")]
unsafe fn ring_mul_acc_avx2_impl(acc: &mut RingElem, a: &RingElem, b: &RingElem) {
    let a_lo: &[u16; HALF] = a.0[..HALF].try_into().unwrap();
    let a_hi: &[u16; HALF] = a.0[HALF..].try_into().unwrap();
    let b_lo: &[u16; HALF] = b.0[..HALF].try_into().unwrap();
    let b_hi: &[u16; HALF] = b.0[HALF..].try_into().unwrap();

    // Three sub-products
    let mut z0 = [0u16; RING_DEG];
    let mut z2 = [0u16; RING_DEG];
    schoolbook_128(&mut z0, a_lo, b_lo);
    schoolbook_128(&mut z2, a_hi, b_hi);

    // Cross term: (a_lo + a_hi) * (b_lo + b_hi)
    let mut a_sum = [0u16; HALF];
    let mut b_sum = [0u16; HALF];

    // HALF = 128 is divisible by 16.
    for i in (0..HALF).step_by(16) {
        let al = _mm256_loadu_si256(a_lo.as_ptr().add(i) as *const __m256i);
        let ah = _mm256_loadu_si256(a_hi.as_ptr().add(i) as *const __m256i);
        _mm256_storeu_si256(
            a_sum.as_mut_ptr().add(i) as *mut __m256i,
            _mm256_add_epi16(al, ah),
        );

        let bl = _mm256_loadu_si256(b_lo.as_ptr().add(i) as *const __m256i);
        let bh = _mm256_loadu_si256(b_hi.as_ptr().add(i) as *const __m256i);
        _mm256_storeu_si256(
            b_sum.as_mut_ptr().add(i) as *mut __m256i,
            _mm256_add_epi16(bl, bh),
        );
    }

    let mut z3 = [0u16; RING_DEG];
    schoolbook_128(&mut z3, &a_sum, &b_sum);

    // Accumulate into acc, reducing mod (X^256 + 1).
    //
    // For j in 0..HALF:
    //   acc[j]      += z0[j] - z2[j] - (z3[j+H] - z0[j+H] - z2[j+H])
    //   acc[j+HALF] += z0[j+H] - z2[j+H] + (z3[j] - z0[j] - z2[j])
    //
    // HALF = 128 is divisible by 16.
    for j in (0..HALF).step_by(16) {
        // --- low half ---
        let acc_lo = _mm256_loadu_si256(acc.0.as_ptr().add(j) as *const __m256i);
        let z0_lo = _mm256_loadu_si256(z0.as_ptr().add(j) as *const __m256i);
        let z2_lo = _mm256_loadu_si256(z2.as_ptr().add(j) as *const __m256i);
        let z3_hi = _mm256_loadu_si256(z3.as_ptr().add(j + HALF) as *const __m256i);
        let z0_hi = _mm256_loadu_si256(z0.as_ptr().add(j + HALF) as *const __m256i);
        let z2_hi = _mm256_loadu_si256(z2.as_ptr().add(j + HALF) as *const __m256i);

        let z1_wrap = _mm256_sub_epi16(_mm256_sub_epi16(z3_hi, z0_hi), z2_hi);
        let res_lo = _mm256_sub_epi16(
            _mm256_sub_epi16(_mm256_add_epi16(acc_lo, z0_lo), z2_lo),
            z1_wrap,
        );
        _mm256_storeu_si256(acc.0.as_mut_ptr().add(j) as *mut __m256i, res_lo);

        // --- high half ---
        let acc_hi = _mm256_loadu_si256(acc.0.as_ptr().add(j + HALF) as *const __m256i);
        let z3_lo = _mm256_loadu_si256(z3.as_ptr().add(j) as *const __m256i);

        let z1_direct = _mm256_sub_epi16(_mm256_sub_epi16(z3_lo, z0_lo), z2_lo);
        let res_hi = _mm256_add_epi16(
            _mm256_sub_epi16(_mm256_add_epi16(acc_hi, z0_hi), z2_hi),
            z1_direct,
        );
        _mm256_storeu_si256(acc.0.as_mut_ptr().add(j + HALF) as *mut __m256i, res_hi);
    }
}

// ---------------------------------------------------------------------------
// Element-wise operations
// ---------------------------------------------------------------------------

/// AVX2-accelerated right shift of all coefficients.
pub(super) fn shift_right_avx2(coeffs: &mut [u16; RING_DEG], shift: usize) {
    // Safety: module gated by cfg(all(feature = "avx2", target_arch = "x86_64")).
    // RING_DEG = 256 is divisible by 16; shift fits in i32 for valid Saber params.
    unsafe {
        shift_right_avx2_impl(coeffs, shift);
    }
}

/// # Safety
///
/// Requires AVX2 support on the executing CPU.
#[target_feature(enable = "avx2")]
unsafe fn shift_right_avx2_impl(coeffs: &mut [u16; RING_DEG], shift: usize) {
    let sh = _mm_cvtsi32_si128(shift as i32);
    for i in (0..RING_DEG).step_by(16) {
        let v = _mm256_loadu_si256(coeffs.as_ptr().add(i) as *const __m256i);
        _mm256_storeu_si256(
            coeffs.as_mut_ptr().add(i) as *mut __m256i,
            _mm256_srl_epi16(v, sh),
        );
    }
}

/// AVX2-accelerated left shift of all coefficients.
pub(super) fn shift_left_avx2(coeffs: &mut [u16; RING_DEG], shift: usize) {
    // Safety: same as shift_right_avx2.
    unsafe {
        shift_left_avx2_impl(coeffs, shift);
    }
}

/// # Safety
///
/// Requires AVX2 support on the executing CPU.
#[target_feature(enable = "avx2")]
unsafe fn shift_left_avx2_impl(coeffs: &mut [u16; RING_DEG], shift: usize) {
    let sh = _mm_cvtsi32_si128(shift as i32);
    for i in (0..RING_DEG).step_by(16) {
        let v = _mm256_loadu_si256(coeffs.as_ptr().add(i) as *const __m256i);
        _mm256_storeu_si256(
            coeffs.as_mut_ptr().add(i) as *mut __m256i,
            _mm256_sll_epi16(v, sh),
        );
    }
}

/// AVX2-accelerated wrapping add of a constant to all coefficients.
pub(super) fn wrapping_add_to_all_avx2(coeffs: &mut [u16; RING_DEG], val: u16) {
    // Safety: RING_DEG = 256 is divisible by 16.
    unsafe {
        wrapping_add_to_all_avx2_impl(coeffs, val);
    }
}

/// # Safety
///
/// Requires AVX2 support on the executing CPU.
#[target_feature(enable = "avx2")]
unsafe fn wrapping_add_to_all_avx2_impl(coeffs: &mut [u16; RING_DEG], val: u16) {
    let vv = _mm256_set1_epi16(val as i16);
    for i in (0..RING_DEG).step_by(16) {
        let v = _mm256_loadu_si256(coeffs.as_ptr().add(i) as *const __m256i);
        _mm256_storeu_si256(
            coeffs.as_mut_ptr().add(i) as *mut __m256i,
            _mm256_add_epi16(v, vv),
        );
    }
}

/// AVX2-accelerated fused add-then-right-shift of all coefficients.
pub(super) fn wrapping_add_and_shift_right_avx2(
    coeffs: &mut [u16; RING_DEG],
    val: u16,
    shift: usize,
) {
    // Safety: RING_DEG = 256 is divisible by 16; shift fits in i32.
    unsafe {
        wrapping_add_and_shift_right_avx2_impl(coeffs, val, shift);
    }
}

/// # Safety
///
/// Requires AVX2 support on the executing CPU.
#[target_feature(enable = "avx2")]
unsafe fn wrapping_add_and_shift_right_avx2_impl(
    coeffs: &mut [u16; RING_DEG],
    val: u16,
    shift: usize,
) {
    let vv = _mm256_set1_epi16(val as i16);
    let sh = _mm_cvtsi32_si128(shift as i32);
    for i in (0..RING_DEG).step_by(16) {
        let v = _mm256_loadu_si256(coeffs.as_ptr().add(i) as *const __m256i);
        let added = _mm256_add_epi16(v, vv);
        _mm256_storeu_si256(
            coeffs.as_mut_ptr().add(i) as *mut __m256i,
            _mm256_srl_epi16(added, sh),
        );
    }
}

/// AVX2-accelerated element-wise wrapping addition: `out = a + b`.
pub(super) fn add_avx2(a: &[u16; RING_DEG], b: &[u16; RING_DEG], out: &mut [u16; RING_DEG]) {
    // Safety: RING_DEG = 256 is divisible by 16.
    unsafe {
        add_avx2_impl(a, b, out);
    }
}

/// # Safety
///
/// Requires AVX2 support on the executing CPU.
#[target_feature(enable = "avx2")]
unsafe fn add_avx2_impl(a: &[u16; RING_DEG], b: &[u16; RING_DEG], out: &mut [u16; RING_DEG]) {
    for i in (0..RING_DEG).step_by(16) {
        let av = _mm256_loadu_si256(a.as_ptr().add(i) as *const __m256i);
        let bv = _mm256_loadu_si256(b.as_ptr().add(i) as *const __m256i);
        _mm256_storeu_si256(
            out.as_mut_ptr().add(i) as *mut __m256i,
            _mm256_add_epi16(av, bv),
        );
    }
}

/// AVX2-accelerated element-wise wrapping subtraction: `out = a - b`.
pub(super) fn sub_avx2(a: &[u16; RING_DEG], b: &[u16; RING_DEG], out: &mut [u16; RING_DEG]) {
    // Safety: RING_DEG = 256 is divisible by 16.
    unsafe {
        sub_avx2_impl(a, b, out);
    }
}

/// # Safety
///
/// Requires AVX2 support on the executing CPU.
#[target_feature(enable = "avx2")]
unsafe fn sub_avx2_impl(a: &[u16; RING_DEG], b: &[u16; RING_DEG], out: &mut [u16; RING_DEG]) {
    for i in (0..RING_DEG).step_by(16) {
        let av = _mm256_loadu_si256(a.as_ptr().add(i) as *const __m256i);
        let bv = _mm256_loadu_si256(b.as_ptr().add(i) as *const __m256i);
        _mm256_storeu_si256(
            out.as_mut_ptr().add(i) as *mut __m256i,
            _mm256_sub_epi16(av, bv),
        );
    }
}

// ---------------------------------------------------------------------------
// Tests — cross-validate AVX2 against trivial reference implementations
// ---------------------------------------------------------------------------

#[cfg(test)]
mod test {
    use super::*;

    /// Trivial reference schoolbook — used only to validate the AVX2 version.
    fn reference_schoolbook_128(out: &mut [u16; RING_DEG], a: &[u16; HALF], b: &[u16; HALF]) {
        for i in 0..HALF {
            for j in 0..HALF {
                out[i + j] = out[i + j].wrapping_add(a[i].wrapping_mul(b[j]));
            }
        }
    }

    #[test]
    fn avx2_schoolbook_matches_reference() {
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
            let mut out_avx2 = [0u16; RING_DEG];
            reference_schoolbook_128(&mut out_ref, &a, &b);
            // Safety: test only runs on x86_64 with AVX2
            unsafe { schoolbook_128(&mut out_avx2, &a, &b) };
            assert_eq!(out_ref, out_avx2, "schoolbook_128 mismatch");
        }
    }

    #[test]
    fn avx2_ring_mul_acc_matches_reference() {
        let mut rng = rand::rng();

        for _ in 0..50 {
            let a = RingElem::rand(&mut rng);
            let b = RingElem::rand(&mut rng);

            // Compute using AVX2 (Karatsuba + AVX2 kernels)
            let mut acc_avx2 = RingElem::default();
            ring_mul_acc_avx2(&mut acc_avx2, &a, &b);

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

            assert_eq!(acc_avx2.0, acc_ref.0, "ring_mul_acc mismatch");
        }
    }

    #[test]
    fn avx2_add_sub_roundtrip() {
        let mut rng = rand::rng();
        let a = RingElem::rand(&mut rng);
        let b = RingElem::rand(&mut rng);

        let mut sum = [0u16; RING_DEG];
        add_avx2(&a.0, &b.0, &mut sum);

        let mut diff = [0u16; RING_DEG];
        sub_avx2(&sum, &b.0, &mut diff);

        assert_eq!(diff, a.0, "add then sub should round-trip");
    }

    #[test]
    fn avx2_shift_roundtrip() {
        let mut rng = rand::rng();
        let orig = RingElem::rand(&mut rng);

        // shift right 3 then left 3: the low 3 bits of each coefficient are lost
        let mut coeffs = orig.0;
        shift_right_avx2(&mut coeffs, 3);
        shift_left_avx2(&mut coeffs, 3);

        for i in 0..RING_DEG {
            assert_eq!(
                coeffs[i],
                orig.0[i] & !(0b111), // low 3 bits lost by right shift
                "shift round-trip mismatch at index {i}"
            );
        }
    }

    #[test]
    fn avx2_wrapping_add_and_shift_right_matches_split() {
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
        wrapping_add_and_shift_right_avx2(&mut coeffs_fused, val, shift);

        // Two-step version
        wrapping_add_to_all_avx2(&mut coeffs_split, val);
        shift_right_avx2(&mut coeffs_split, shift);

        assert_eq!(coeffs_fused, coeffs_split, "fused vs split mismatch");
    }
}
