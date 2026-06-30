//! AVX2-accelerated polynomial arithmetic for x86_64.
//!
//! This module provides SIMD-accelerated versions of the core arithmetic operations
//! using AVX2 intrinsics. All functions operate on `u16` coefficient arrays that
//! represent polynomials in the ring ℤ\[X\]/(X^256 + 1).
//!
//! The ring multiplication uses the Toom-Cook 4-way algorithm from the SABER
//! reference C implementation, with 2-level Karatsuba at the inner layer and
//! 16×16 schoolbook in transposed form for maximum AVX2 throughput.

#![allow(unsafe_code)]

#[cfg(target_arch = "x86")]
use core::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::*;

use super::RingElem;
use crate::consts::RING_DEG;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Number of u16 lanes per __m256i vector.
const CHUNK: usize = 16;
/// Number of AVX vectors per polynomial (256 / 16).
const AVX_N: usize = RING_DEG / CHUNK;
/// Number of AVX vectors per Toom-Cook quarter (16 / 4).
const SMALL_LEN: usize = AVX_N / 4;

/// 3^{-1} mod 2^{16} (3 * 43691 = 131073 ≡ 1 mod 65536).
const INV3: u16 = 43691;
/// 9^{-1} mod 2^{16}.
const INV9: u16 = 36409;
/// 15^{-1} mod 2^{16}.
const INV15: u16 = 61167;
const INT45: u16 = 45;
const INT30: u16 = 30;

// ---------------------------------------------------------------------------
// Inline helpers (used inside #[target_feature(enable = "avx2")] context)
// ---------------------------------------------------------------------------

#[inline(always)]
unsafe fn vadd(a: __m256i, b: __m256i) -> __m256i {
    _mm256_add_epi16(a, b)
}

#[inline(always)]
unsafe fn vsub(a: __m256i, b: __m256i) -> __m256i {
    _mm256_sub_epi16(a, b)
}

#[inline(always)]
unsafe fn vmul(a: __m256i, b: __m256i) -> __m256i {
    _mm256_mullo_epi16(a, b)
}

#[inline(always)]
unsafe fn vzero() -> __m256i {
    _mm256_setzero_si256()
}

#[inline(always)]
unsafe fn vset1(val: u16) -> __m256i {
    _mm256_set1_epi16(val as i16)
}

// ---------------------------------------------------------------------------
// Load / store polynomials as AVX vectors
// ---------------------------------------------------------------------------

/// Load a 256-coefficient polynomial into 16 AVX vectors.
#[target_feature(enable = "avx2")]
unsafe fn load_poly(coeffs: &[u16; RING_DEG]) -> [__m256i; AVX_N] {
    let mut v = [vzero(); AVX_N];
    for i in 0..AVX_N {
        v[i] = _mm256_loadu_si256(coeffs.as_ptr().add(i * CHUNK) as *const __m256i);
    }
    v
}

// ---------------------------------------------------------------------------
// Toom-Cook 4-way evaluation
// ---------------------------------------------------------------------------

/// Evaluate a polynomial at the 7 Toom-Cook points.
///
/// Input: `poly[0..16]` — 16 AVX vectors representing a 256-coeff polynomial,
///        split into 4 quarters q0..q3 of `SMALL_LEN` (4) vectors each.
///
/// Output: `aw[0..28]` — 7 evaluation polynomials × 4 vectors.
#[target_feature(enable = "avx2")]
unsafe fn tc_eval(poly: &[__m256i; AVX_N], aw: &mut [__m256i; 28]) {
    for i in 0..SMALL_LEN {
        let r0 = poly[i]; // q0
        let r1 = poly[i + SMALL_LEN]; // q1
        let r2 = poly[i + 2 * SMALL_LEN]; // q2
        let r3 = poly[i + 3 * SMALL_LEN]; // q3

        let r4 = vadd(r0, r2);
        let r5 = vadd(r1, r3);

        // P(1) = r0 + r1 + r2 + r3
        aw[2 * SMALL_LEN + i] = vadd(r4, r5);
        // P(-1) = r0 - r1 + r2 - r3
        aw[3 * SMALL_LEN + i] = vsub(r4, r5);

        // 8*P(1/2) = 8r0 + 4r1 + 2r2 + r3
        let r4_half = _mm256_slli_epi16(vadd(_mm256_slli_epi16(r0, 2), r2), 1);
        let r5_half = vadd(_mm256_slli_epi16(r1, 2), r3);
        aw[4 * SMALL_LEN + i] = vadd(r4_half, r5_half);
        // 8*P(-1/2) = 8r0 - 4r1 + 2r2 - r3
        aw[5 * SMALL_LEN + i] = vsub(r4_half, r5_half);

        // P(2) = r0 + 2r1 + 4r2 + 8r3
        let r4_two = vadd(
            vadd(r0, _mm256_slli_epi16(r1, 1)),
            vadd(_mm256_slli_epi16(r2, 2), _mm256_slli_epi16(r3, 3)),
        );
        aw[1 * SMALL_LEN + i] = r4_two;

        // P(0) = r0 (constant term)
        aw[6 * SMALL_LEN + i] = r0;
        // P(∞) = r3 (leading coefficient)
        aw[0 * SMALL_LEN + i] = r3;
    }
}

// ---------------------------------------------------------------------------
// 2-level Karatsuba evaluation
// ---------------------------------------------------------------------------

/// Karatsuba evaluation for a single 4-vector polynomial (64 coefficients).
///
/// Splits into 4 chunks of 16 and produces 9 sub-problem terms.
#[target_feature(enable = "avx2")]
unsafe fn kara_eval_single(poly: &[__m256i], bucket: &mut [__m256i], base: usize) {
    let r0 = poly[0];
    let r1 = poly[1];
    let r2 = poly[2];
    let r3 = poly[3];

    bucket[base] = r0;
    bucket[base + 1] = r1;
    bucket[base + 2] = r2;
    bucket[base + 3] = r3;
    bucket[base + 4] = vadd(r0, r1);
    bucket[base + 5] = vadd(r2, r3);
    let r0r2 = vadd(r0, r2);
    let r1r3 = vadd(r1, r3);
    bucket[base + 6] = r0r2;
    bucket[base + 7] = r1r3;
    bucket[base + 8] = vadd(r0r2, r1r3);
}

/// Apply Karatsuba evaluation to all 7 TC evaluation polynomials.
#[target_feature(enable = "avx2")]
unsafe fn kara_eval(aw: &[__m256i; 28], bucket: &mut [__m256i; 64]) {
    for k in 0..7 {
        kara_eval_single(&aw[k * SMALL_LEN..k * SMALL_LEN + SMALL_LEN], bucket, k * 9);
    }
    // Pad entry 63 with zero (7 * 9 = 63 entries used, bucket[63] is padding)
    bucket[63] = vzero();
}

// ---------------------------------------------------------------------------
// 16×16 transpose
// ---------------------------------------------------------------------------

/// Transpose a 16×16 matrix of u16 stored as 16 `__m256i` vectors in-place.
///
/// Before: `m[row]` holds the 16 u16 values of row `row`.
/// After:  `m[col]` holds the 16 u16 values of column `col`.
///
/// All 16 inputs are read into locals before any output is written to avoid
/// aliasing issues (the C reference interleaves reads and writes).
#[target_feature(enable = "avx2")]
unsafe fn transpose_16x16(m: &mut [__m256i]) {
    debug_assert!(m.len() >= 16);

    // Read all 16 inputs first
    let m0 = m[0];
    let m1 = m[1];
    let m2 = m[2];
    let m3 = m[3];
    let m4 = m[4];
    let m5 = m[5];
    let m6 = m[6];
    let m7 = m[7];
    let m8 = m[8];
    let m9 = m[9];
    let m10 = m[10];
    let m11 = m[11];
    let m12 = m[12];
    let m13 = m[13];
    let m14 = m[14];
    let m15 = m[15];

    // Phase 1: unpacklo_epi16 on pairs
    let r0 = _mm256_unpacklo_epi16(m0, m1);
    let r1 = _mm256_unpacklo_epi16(m2, m3);
    let r2 = _mm256_unpacklo_epi16(m4, m5);
    let r3 = _mm256_unpacklo_epi16(m6, m7);
    let r4 = _mm256_unpacklo_epi16(m8, m9);
    let r5 = _mm256_unpacklo_epi16(m10, m11);
    let r6 = _mm256_unpacklo_epi16(m12, m13);
    let r7 = _mm256_unpacklo_epi16(m14, m15);

    // Phase 2: unpacklo/hi_epi32
    let temp = _mm256_unpacklo_epi32(r0, r1);
    let temp0 = _mm256_unpacklo_epi32(r2, r3);
    let temp1 = _mm256_unpacklo_epi32(r4, r5);
    let temp2 = _mm256_unpacklo_epi32(r6, r7);

    let r8 = _mm256_unpackhi_epi32(r0, r1);
    let r9 = _mm256_unpackhi_epi32(r2, r3);
    let r10 = _mm256_unpackhi_epi32(r4, r5);
    let r11 = _mm256_unpackhi_epi32(r6, r7);

    // Phase 3: unpacklo/hi_epi64
    let s0 = _mm256_unpacklo_epi64(temp, temp0);
    let s2 = _mm256_unpackhi_epi64(temp, temp0);
    let s1 = _mm256_unpacklo_epi64(temp1, temp2);
    let s3 = _mm256_unpackhi_epi64(temp1, temp2);

    // Phase 4: high halves
    let h_temp = _mm256_unpackhi_epi16(m0, m1);
    let h_temp0 = _mm256_unpackhi_epi16(m2, m3);
    let h_temp1 = _mm256_unpackhi_epi16(m4, m5);
    let h_temp2 = _mm256_unpackhi_epi16(m6, m7);
    let h_r4 = _mm256_unpackhi_epi16(m8, m9);
    let h_r5 = _mm256_unpackhi_epi16(m10, m11);
    let h_r6 = _mm256_unpackhi_epi16(m12, m13);
    let h_r7 = _mm256_unpackhi_epi16(m14, m15);

    // Write outputs 0, 8, 1, 9 from s0..s3
    m[0] = _mm256_permute2x128_si256(s0, s1, 0x20);
    m[8] = _mm256_permute2x128_si256(s0, s1, 0x31);
    m[1] = _mm256_permute2x128_si256(s2, s3, 0x20);
    m[9] = _mm256_permute2x128_si256(s2, s3, 0x31);

    // Phase 5: epi64 on r8..r11
    let s4 = _mm256_unpacklo_epi64(r8, r9);
    let s5 = _mm256_unpacklo_epi64(r10, r11);
    let s6 = _mm256_unpackhi_epi64(r8, r9);
    let s7 = _mm256_unpackhi_epi64(r10, r11);

    m[3] = _mm256_permute2x128_si256(s6, s7, 0x20);
    m[11] = _mm256_permute2x128_si256(s6, s7, 0x31);
    m[2] = _mm256_permute2x128_si256(s4, s5, 0x20);
    m[10] = _mm256_permute2x128_si256(s4, s5, 0x31);

    // Phase 6: high-half epi32
    let hh0 = _mm256_unpacklo_epi32(h_temp, h_temp0);
    let hh1 = _mm256_unpacklo_epi32(h_temp1, h_temp2);
    let hh2 = _mm256_unpacklo_epi32(h_r4, h_r5);
    let hh3 = _mm256_unpacklo_epi32(h_r6, h_r7);

    let hh8 = _mm256_unpacklo_epi64(hh0, hh1);
    let hh10 = _mm256_unpackhi_epi64(hh0, hh1);
    let hh9 = _mm256_unpacklo_epi64(hh2, hh3);
    let hh11 = _mm256_unpackhi_epi64(hh2, hh3);

    m[4] = _mm256_permute2x128_si256(hh8, hh9, 0x20);
    m[12] = _mm256_permute2x128_si256(hh8, hh9, 0x31);
    m[5] = _mm256_permute2x128_si256(hh10, hh11, 0x20);
    m[13] = _mm256_permute2x128_si256(hh10, hh11, 0x31);

    // Phase 7: high-half epi32 hi
    let hi0 = _mm256_unpackhi_epi32(h_temp, h_temp0);
    let hi1 = _mm256_unpackhi_epi32(h_temp1, h_temp2);
    let hi2 = _mm256_unpackhi_epi32(h_r4, h_r5);
    let hi3 = _mm256_unpackhi_epi32(h_r6, h_r7);

    let hi4 = _mm256_unpacklo_epi64(hi0, hi1);
    let hi6 = _mm256_unpackhi_epi64(hi0, hi1);
    let hi5 = _mm256_unpacklo_epi64(hi2, hi3);
    let hi7 = _mm256_unpackhi_epi64(hi2, hi3);

    m[6] = _mm256_permute2x128_si256(hi4, hi5, 0x20);
    m[14] = _mm256_permute2x128_si256(hi4, hi5, 0x31);
    m[7] = _mm256_permute2x128_si256(hi6, hi7, 0x20);
    m[15] = _mm256_permute2x128_si256(hi6, hi7, 0x31);
}

// ---------------------------------------------------------------------------
// 16×16 schoolbook multiplication (transposed form)
// ---------------------------------------------------------------------------

/// Schoolbook multiplication of two degree-15 polynomials in transposed form.
///
/// Each `__m256i` holds one coefficient across 16 independent problems.
/// Output `c[0..31]` where `c[31] = 0`.
#[target_feature(enable = "avx2")]
unsafe fn schoolbook_16x16(a: &[__m256i], b: &[__m256i], c: &mut [__m256i]) {
    debug_assert!(a.len() >= 16);
    debug_assert!(b.len() >= 16);
    debug_assert!(c.len() >= 32);

    for k in 0..31 {
        c[k] = vzero();
    }
    c[31] = vzero();

    for i in 0..16 {
        for j in 0..16 {
            c[i + j] = vadd(c[i + j], vmul(a[i], b[j]));
        }
    }
    c[31] = vzero();
}

// ---------------------------------------------------------------------------
// Karatsuba interpolation
// ---------------------------------------------------------------------------

/// Karatsuba interpolation for a single TC evaluation point.
///
/// Recovers the 128-coeff product from 9 sub-products stored in `c_bucket`.
/// `base_idx` is `k * 9` for TC point `k`.
#[target_feature(enable = "avx2")]
unsafe fn kara_interpol_single(c_bucket: &[__m256i], base_idx: usize, result: &mut [__m256i; 8]) {
    // Map sub-product index to c_bucket position.
    // After reverse-transpose, sub-product `s` has:
    //   low 16 coeffs at c_bucket[block*32 + offset]
    //   high 16 coeffs at c_bucket[block*32 + 16 + offset]
    // where block = flat/16, offset = flat%16, flat = base_idx + s.
    #[inline(always)]
    fn idx(base_idx: usize, s: usize) -> (usize, usize) {
        let flat = base_idx + s;
        let b = flat / 16;
        let j = flat % 16;
        (b * 32 + j, b * 32 + 16 + j)
    }

    let (p0l, p0h) = idx(base_idx, 0); // r0*s0
    let (p1l, p1h) = idx(base_idx, 1); // r1*s1
    let (p2l, p2h) = idx(base_idx, 2); // r2*s2
    let (p3l, p3h) = idx(base_idx, 3); // r3*s3
    let (p4l, p4h) = idx(base_idx, 4); // (r0+r1)*(s0+s1)
    let (p5l, p5h) = idx(base_idx, 5); // (r2+r3)*(s2+s3)
    let (p6l, p6h) = idx(base_idx, 6); // (r0+r2)*(s0+s2)
    let (p7l, p7h) = idx(base_idx, 7); // (r1+r3)*(s1+s3)
    let (p8l, p8h) = idx(base_idx, 8); // (r0+r1+r2+r3)*(s0+s1+s2+s3)

    let r0 = c_bucket[p0l];
    let mut r1 = c_bucket[p0h];
    let mut r2 = c_bucket[p1l];
    let r3 = c_bucket[p1h];
    let r4 = c_bucket[p2l];
    let mut r5 = c_bucket[p2h];
    let mut r6 = c_bucket[p3l];
    let r7 = c_bucket[p3h];

    let mut c6_l = c_bucket[p6l];
    let mut c6_h = c_bucket[p6h];
    let mut c7_l = c_bucket[p7l];
    let c7_h = c_bucket[p7h];

    // Cross terms
    let cross_z1_l = vsub(vsub(c_bucket[p8l], c6_l), c7_l);
    let cross_z1_h = vsub(vsub(c_bucket[p8h], c6_h), c7_h);
    let cross23_l = vsub(vsub(c_bucket[p5l], r4), r6);
    let cross23_h = vsub(vsub(c_bucket[p5h], r5), r7);
    let cross01_l = vsub(vsub(c_bucket[p4l], r0), r2);
    let cross01_h = vsub(vsub(c_bucket[p4h], r1), r3);

    // Mix: add cross terms to overlapping positions
    r5 = vadd(r5, cross23_l);
    r1 = vadd(r1, cross01_l);
    c6_h = vadd(c6_h, cross_z1_l);
    r6 = vadd(r6, cross23_h);
    r2 = vadd(r2, cross01_h);
    c7_l = vadd(c7_l, cross_z1_h);

    // z1 = z1_full - z0 - z2
    c6_l = vsub(vsub(c6_l, r0), r4);
    c6_h = vsub(vsub(c6_h, r1), r5);
    c7_l = vsub(vsub(c7_l, r2), r6);
    let c7_h_final = vsub(vsub(c7_h, r3), r7);

    // Assemble result = z0 + z1*X^32 + z2*X^64
    result[0] = r0;
    result[1] = r1;
    result[2] = vadd(r2, c6_l);
    result[3] = vadd(r3, c6_h);
    result[4] = vadd(r4, c7_l);
    result[5] = vadd(r5, c7_h_final);
    result[6] = r6;
    result[7] = r7;
}

// ---------------------------------------------------------------------------
// Toom-Cook interpolation
// ---------------------------------------------------------------------------

/// Toom-Cook interpolation: recover the full product from 7 evaluation-point products.
///
/// `w[0..7]` are the 7 products, each 8 vectors (128 coefficients).
/// Output: `result[0..16]` — the final 256-coeff product reduced mod X^256+1.
#[target_feature(enable = "avx2")]
unsafe fn tc_interpol(w: &[[__m256i; 8]; 7], result: &mut [__m256i; AVX_N]) {
    let inv3_avx = vset1(INV3);
    let inv9_avx = vset1(INV9);
    let inv15_avx = vset1(INV15);
    let int45_avx = vset1(INT45);
    let int30_avx = vset1(INT30);

    // res_output has 32 vectors (512 coefficients) for the unreduced product
    let mut res_output = [vzero(); 32];

    for i in 0..8 {
        let r0 = w[0][i]; // product at infinity
        let mut r1 = w[1][i]; // P(2)*Q(2)
        let mut r2 = w[2][i]; // P(1)*Q(1)
        let mut r3 = w[3][i]; // P(-1)*Q(-1)
        let mut r4 = w[4][i]; // 64*P(1/2)*Q(1/2)
        let mut r5 = w[5][i]; // 64*P(-1/2)*Q(-1/2)
        let r6 = w[6][i]; // P(0)*Q(0)

        r1 = vadd(r1, r4);
        r5 = vsub(r5, r4);
        r3 = vsub(r3, r2);
        r3 = _mm256_srli_epi16(r3, 1);
        r4 = vsub(r4, r0);
        r4 = vsub(r4, _mm256_slli_epi16(r6, 6));
        r4 = _mm256_slli_epi16(r4, 1);
        r4 = vadd(r4, r5);
        r2 = vadd(r2, r3);
        r1 = vsub(r1, _mm256_slli_epi16(r2, 6));
        r1 = vsub(r1, r2);
        r2 = vsub(r2, r6);
        r2 = vsub(r2, r0);
        r1 = vadd(r1, vmul(r2, int45_avx));
        r4 = vsub(r4, _mm256_slli_epi16(r2, 3));
        r4 = vmul(r4, inv3_avx);
        r4 = _mm256_srli_epi16(r4, 3);
        r5 = vadd(r5, r1);
        r1 = vadd(r1, _mm256_slli_epi16(r3, 4));
        r1 = vmul(r1, inv9_avx);
        r1 = _mm256_srli_epi16(r1, 1);
        r3 = vadd(r1, r3);
        r3 = vsub(vzero(), r3); // negate
        let temp_val = vmul(r1, int30_avx);
        let temp2 = vsub(temp_val, r5);
        let temp3 = vmul(temp2, inv15_avx);
        r5 = _mm256_srli_epi16(temp3, 2);
        r2 = vsub(r2, r4);
        r1 = vsub(r1, r5);

        // Store into res_output
        if i < SMALL_LEN {
            // First half — assign
            res_output[0 * SMALL_LEN + i] = r6;
            res_output[1 * SMALL_LEN + i] = r5;
            res_output[2 * SMALL_LEN + i] = r4;
            res_output[3 * SMALL_LEN + i] = r3;
            res_output[4 * SMALL_LEN + i] = r2;
            res_output[5 * SMALL_LEN + i] = r1;
            res_output[6 * SMALL_LEN + i] = r0;
        } else {
            // Second half — add to overlap (except r0 which is assign)
            res_output[0 * SMALL_LEN + i] = vadd(res_output[0 * SMALL_LEN + i], r6);
            res_output[1 * SMALL_LEN + i] = vadd(res_output[1 * SMALL_LEN + i], r5);
            res_output[2 * SMALL_LEN + i] = vadd(res_output[2 * SMALL_LEN + i], r4);
            res_output[3 * SMALL_LEN + i] = vadd(res_output[3 * SMALL_LEN + i], r3);
            res_output[4 * SMALL_LEN + i] = vadd(res_output[4 * SMALL_LEN + i], r2);
            res_output[5 * SMALL_LEN + i] = vadd(res_output[5 * SMALL_LEN + i], r1);
            res_output[6 * SMALL_LEN + i] = r0; // assign, not add
        }
    }

    // Reduce mod X^256 + 1: result[i] = res_output[i] - res_output[i + 16]
    for i in 0..AVX_N {
        result[i] = vsub(res_output[i], res_output[i + AVX_N]);
    }
}

// ---------------------------------------------------------------------------
// Core polynomial multiplication
// ---------------------------------------------------------------------------

/// AVX2-accelerated ring multiply-accumulate using Toom-Cook 4-way.
///
/// Computes `acc += a * b` in the ring ℤ\[X\]/(X^256 + 1).
pub(super) fn ring_mul_acc_avx2(acc: &mut RingElem, a: &RingElem, b: &RingElem) {
    // Safety: module gated by cfg(all(feature = "avx2", target_arch = "x86_64")).
    unsafe { ring_mul_acc_avx2_impl(acc, a, b) }
}

/// # Safety
///
/// Requires AVX2 support on the executing CPU.
#[target_feature(enable = "avx2")]
unsafe fn ring_mul_acc_avx2_impl(acc: &mut RingElem, a: &RingElem, b: &RingElem) {
    // 1. Load polynomials into AVX vectors
    let a_avx = load_poly(&a.0);
    let b_avx = load_poly(&b.0);

    // 2. TC evaluation of b + Karatsuba evaluation + transpose
    let mut bw = [vzero(); 28];
    tc_eval(&b_avx, &mut bw);
    let mut b_bucket = [vzero(); 64];
    kara_eval(&bw, &mut b_bucket);
    transpose_16x16(&mut b_bucket[0..16]);
    transpose_16x16(&mut b_bucket[16..32]);
    transpose_16x16(&mut b_bucket[32..48]);
    transpose_16x16(&mut b_bucket[48..64]);

    // 3. TC evaluation of a + Karatsuba evaluation + transpose
    let mut aw = [vzero(); 28];
    tc_eval(&a_avx, &mut aw);
    let mut a_bucket = [vzero(); 64];
    kara_eval(&aw, &mut a_bucket);
    transpose_16x16(&mut a_bucket[0..16]);
    transpose_16x16(&mut a_bucket[16..32]);
    transpose_16x16(&mut a_bucket[32..48]);
    transpose_16x16(&mut a_bucket[48..64]);

    // 4. Schoolbook: 4 blocks of 16×16
    let mut c_bucket = [vzero(); 128];
    schoolbook_16x16(&a_bucket[0..16], &b_bucket[0..16], &mut c_bucket[0..32]);
    schoolbook_16x16(&a_bucket[16..32], &b_bucket[16..32], &mut c_bucket[32..64]);
    schoolbook_16x16(&a_bucket[32..48], &b_bucket[32..48], &mut c_bucket[64..96]);
    schoolbook_16x16(&a_bucket[48..64], &b_bucket[48..64], &mut c_bucket[96..128]);

    // 5. Reverse transpose (8 blocks of 16)
    for blk in 0..8 {
        transpose_16x16(&mut c_bucket[blk * 16..(blk + 1) * 16]);
    }

    // 6. Karatsuba interpolation → 7 products of 8 vectors each
    let mut w = [[vzero(); 8]; 7];
    for k in 0..7 {
        kara_interpol_single(&c_bucket, k * 9, &mut w[k]);
    }

    // 7. TC interpolation → result[16]
    let mut mul_result = [vzero(); AVX_N];
    tc_interpol(&w, &mut mul_result);

    // 8. Accumulate into acc
    for i in 0..AVX_N {
        let cur = _mm256_loadu_si256(acc.0.as_ptr().add(i * CHUNK) as *const __m256i);
        let res = vadd(cur, mul_result[i]);
        _mm256_storeu_si256(acc.0.as_mut_ptr().add(i * CHUNK) as *mut __m256i, res);
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

    #[test]
    fn avx2_ring_mul_acc_matches_reference() {
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

            // Compute using AVX2 (Toom-Cook 4-way)
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

            for i in 0..RING_DEG {
                assert_eq!(
                    acc_avx2.0[i] & mask,
                    acc_ref.0[i] & mask,
                    "ring_mul_acc mismatch at index {i}: avx2={} ref={}",
                    acc_avx2.0[i],
                    acc_ref.0[i],
                );
            }
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
