//! Serialization and deserialization routines for ring elements

// Algorithm 9, BS2POLN
/// Deserializes the given bitstring into a u16 array. Every element of the array has
/// `bits_per_elem` bits (must be ≤ 13), encoded in the lower bits of the word.
pub(crate) fn deserialize<const N: usize>(bytes: &[u8], bits_per_elem: usize) -> [u16; N] {
    // Serialized bitlength must be a multiple of 8, and `bytes` must be the correct length
    debug_assert_eq!((bits_per_elem * N) % 8, 0);
    debug_assert_eq!(bytes.len(), bits_per_elem * N / 8);

    let mut out = [0u16; N];
    let bitmask: u32 = (1 << bits_per_elem) - 1;

    // Sliding window: holds pending bits from the byte stream. We refill from bytes
    // one at a time and extract elements from the bottom.
    let mut window: u32 = 0;
    let mut bits_in_window: usize = 0;
    let mut byte_pos: usize = 0;

    for idx in 0..N {
        // Ensure we have enough bits in the window for one element
        while bits_in_window < bits_per_elem {
            window |= (bytes[byte_pos] as u32) << bits_in_window;
            byte_pos += 1;
            bits_in_window += 8;
        }

        // Extract the lowest bits_per_elem bits as one element
        out[idx] = (window & bitmask) as u16;
        window >>= bits_per_elem;
        bits_in_window -= bits_per_elem;
    }

    out
}

/// Fast specialization of [`deserialize`] for the 13-bit case (matrix/公钥 expansion),
/// which is by far the hottest width. Processes a whole 13-byte group into 8 coefficients
/// with fixed shifts and no per-element branching, so it vectorizes well.
pub(crate) fn deserialize_13(bytes: &[u8; 13 * 256 / 8]) -> [u16; 256] {
    let mut out = [0u16; 256];
    // 256 coeffs = 32 groups of 8, each group packed into 13 bytes.
    for g in 0..32 {
        let b = &bytes[13 * g..13 * g + 13];
        // Widen once so the shifts below can't lose high bits.
        let w = |k: usize| b[k] as u16;
        let o = 8 * g;
        out[o] = w(0) | ((w(1) & 0x1f) << 8);
        out[o + 1] = (w(1) >> 5) | ((w(2)) << 3) | ((w(3) & 0x03) << 11);
        out[o + 2] = (w(3) >> 2) | ((w(4) & 0x7f) << 6);
        out[o + 3] = (w(4) >> 7) | ((w(5)) << 1) | ((w(6) & 0x0f) << 9);
        out[o + 4] = (w(6) >> 4) | ((w(7)) << 4) | ((w(8) & 0x01) << 12);
        out[o + 5] = (w(8) >> 1) | ((w(9) & 0x3f) << 7);
        out[o + 6] = (w(9) >> 6) | ((w(10)) << 2) | ((w(11) & 0x07) << 10);
        out[o + 7] = (w(11) >> 3) | ((w(12)) << 5);
    }
    out
}

/// Fast specialization of [`deserialize`] for the 10-bit case (ciphertext/public-key vector
/// unpacking). Processes a 5-byte group into 4 coefficients with fixed shifts.
pub(crate) fn deserialize_10(bytes: &[u8; 10 * 256 / 8]) -> [u16; 256] {
    let mut out = [0u16; 256];
    // 256 coeffs = 64 groups of 4, each group packed into 5 bytes.
    for g in 0..64 {
        let b = &bytes[5 * g..5 * g + 5];
        let w = |k: usize| b[k] as u16;
        let o = 4 * g;
        out[o] = w(0) | ((w(1) & 0x03) << 8);
        out[o + 1] = (w(1) >> 2) | ((w(2) & 0x0f) << 6);
        out[o + 2] = (w(2) >> 4) | ((w(3) & 0x3f) << 4);
        out[o + 3] = (w(3) >> 6) | (w(4) << 2);
    }
    out
}

#[cfg(test)]
mod ser_test {
    use super::*;

    // The fast 13-bit path must agree with the generic sliding-window deserializer.
    #[test]
    fn deserialize_13_matches_generic() {
        let mut bytes = [0u8; 13 * 256 / 8];
        // Deterministic pseudo-random fill (no rng dep needed here).
        let mut x: u32 = 0x9e3779b9;
        for b in bytes.iter_mut() {
            x = x.wrapping_mul(1664525).wrapping_add(1013904223);
            *b = (x >> 24) as u8;
        }
        let generic: [u16; 256] = deserialize(&bytes, 13);
        let fast = deserialize_13(&bytes);
        assert_eq!(generic, fast);
    }

    // The fast 10-bit path must agree with the generic sliding-window deserializer.
    #[test]
    fn deserialize_10_matches_generic() {
        let mut bytes = [0u8; 10 * 256 / 8];
        let mut x: u32 = 0x12345678;
        for b in bytes.iter_mut() {
            x = x.wrapping_mul(1664525).wrapping_add(1013904223);
            *b = (x >> 24) as u8;
        }
        let generic: [u16; 256] = deserialize(&bytes, 10);
        let fast = deserialize_10(&bytes);
        assert_eq!(generic, fast);
    }
}

// Algorithm 10, POLN2BS
/// Serializes the given u16 array into a bitstring. Every element of the array has `bits_per_elem`
/// bits (must be ≤ 13), encoded in the lower bits of the word.
pub(crate) fn serialize(data: &[u16], out_buf: &mut [u8], bits_per_elem: usize) {
    assert_eq!(out_buf.len(), bits_per_elem * data.len() / 8);

    let bitmask: u32 = (1 << bits_per_elem) - 1;

    // Sliding window: elements are OR'd in at the current position, and complete bytes
    // are flushed out from the bottom.
    let mut window: u32 = 0;
    let mut bits_in_window: usize = 0;
    let mut byte_pos: usize = 0;

    for &elem in data.iter() {
        // Insert this element's bits into the window
        window |= ((elem as u32) & bitmask) << bits_in_window;
        bits_in_window += bits_per_elem;

        // Flush all complete bytes
        while bits_in_window >= 8 {
            out_buf[byte_pos] = window as u8;
            window >>= 8;
            bits_in_window -= 8;
            byte_pos += 1;
        }
    }

    // There should never be any partial bits in the output
    debug_assert_eq!(bits_in_window, 0);
}
