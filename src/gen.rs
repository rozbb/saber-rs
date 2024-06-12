use crate::ring_arith::{deserialize, RingElem, RING_DEG};

use rand_core::CryptoRngCore;

type Secret<const L: usize> = [RingElem; L];

/// The maximum possible μ value is μ=10, set by Lightsaber
const MAX_MU: usize = 10;

/// The maximum possible l value is l=4, set by FireSaber
const MAX_L: usize = 10;

fn gen_secret<const L: usize, const MU: usize>(rng: &mut impl CryptoRngCore) -> Secret<L> {
    let mut backing_buf = [0u8; MAX_L * RING_DEG * MAX_MU / 8];
    let buf = &mut backing_buf[..L * RING_DEG * MU / 8];
    rng.fill_bytes(buf);

    let mut out_elems = [RingElem::default(); L];
    for chunk in buf.chunks(RING_DEG * MU / 8) {
        let data: [u16; 2 * RING_DEG] = deserialize(chunk, MU / 2);
        // Our output element is hamming(data[i]) - hamming(data[i+1]) for every i, skipping by 2
        let mut polyn = RingElem::default();
        for i in 0..RING_DEG {
            let hamming1 = data[2 * i].count_ones() as u16;
            let hamming2 = data[2 * i + 1].count_ones() as u16;
            polyn.0[i] = hamming1.wrapping_sub(hamming2);
        }
    }

    out_elems
}

#[test]
fn test_gen_secret() {
    let mut rng = rand::thread_rng();
    const L: usize = 2;
    const MU: usize = 10;
    gen_secret::<L, MU>(&mut rng);
}
