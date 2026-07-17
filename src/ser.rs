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
