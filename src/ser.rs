//! This file implements serialization and deserialization routines for ring elements

// Algorithm 9, BS2POLN
/// Deserializes the given bitstring into a u16 array. Every element of the array has
/// `bits_per_elem` bits (must be ≤ 16), encoded in the lower bits of the word.
pub(crate) fn deserialize<const N: usize>(bytes: &[u8], bits_per_elem: usize) -> [u16; N] {
    assert_eq!(bytes.len(), bits_per_elem * N / 8);

    // We only want the lower bits_per_elem bits to be set in any elem of our output
    let bitmask = (1 << bits_per_elem) - 1;

    // Accumulate all the bits into p
    let mut p = [0u16; N];
    let mut bit_idx = 0;
    while bit_idx < bits_per_elem * N {
        let byte_idx = bit_idx / 8;
        let elem_idx = bit_idx / bits_per_elem;
        let bit_in_byte = bit_idx % 8;
        let bit_in_elem = bit_idx % bits_per_elem;

        // Shift the byte we're reading so the first bit we want is the lowest bit
        let data_to_read = bytes[byte_idx] as u16 >> bit_in_byte;
        // OR the byte into our element, shifitng to align with the first unused bit in the elem
        p[elem_idx] |= data_to_read << bit_in_elem;

        // The above might set some high bits we don't want. Clear them now.
        p[elem_idx] &= bitmask;

        // We have read either: 1) however many bits were remaining in the byte we were
        // reading, or 2) however many unused bits remained in the current element we were
        // writing to. Whichever is smaller.
        let just_read = core::cmp::min(8 - bit_in_byte, bits_per_elem - bit_in_elem);
        bit_idx += just_read;
    }

    p
}

// Algorithm 10, POLN2BS
/// Serializes the given u16 array into a bitstring. Every element of the array has `bits_per_elem`
/// bits (must be ≤ 16), encoded in the lower bits of the word.
pub(crate) fn serialize(data: &[u16], out_buf: &mut [u8], bits_per_elem: usize) {
    assert_eq!(out_buf.len(), bits_per_elem * data.len() / 8);

    // Since we use OR to set the output bits, we must clear the buffer at the beginning
    out_buf.fill(0);

    // We only want to write the lower bits_per_elem bits of any element
    let bitmask = (1 << bits_per_elem) - 1;

    // Write all the bits into the given bytestring
    let mut bit_idx = 0;
    while bit_idx < bits_per_elem * data.len() {
        let byte_idx = bit_idx / 8;
        let elem_idx = bit_idx / bits_per_elem;
        let bit_in_byte = bit_idx % 8;
        let bit_in_elem = bit_idx % bits_per_elem;

        // First clear the unused bits of the element we're going to write
        let elem_to_read = data[elem_idx] & bitmask;
        // Then  shift the element we're writing so the first unwritten bit is the lowest bit
        let elem_to_read = (elem_to_read >> bit_in_elem) as u8;

        // OR the bits into our byte, shifitng to align with the first unused bit in the byte
        out_buf[byte_idx] |= elem_to_read << bit_in_byte;

        // We just wrote either: 1) however many bits remained in the byte we were
        // reading, or 2) however many unused bits remained in the current element we were
        // writing to. Whichever is smaller.
        let just_wrote = core::cmp::min(8 - bit_in_byte, bits_per_elem - bit_in_elem);
        bit_idx += just_wrote
    }
}
