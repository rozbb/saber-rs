use crate::{
    gen::MAX_L,
    ind_cpa::{self, IndCpaPublicKey, IndCpaSecretKey},
    MODULUS_P_BITS, RING_DEG,
};

use rand_core::CryptoRngCore;
use sha3::{Digest, Sha3_256};
use subtle::{ConditionallySelectable, ConstantTimeEq};

pub type IndCcaPublicKey<const L: usize> = IndCpaPublicKey<L>;
pub type SharedSecret = [u8; 32];

const MAX_MODULUS_T_BITS: usize = 6;

pub struct IndCcaSecretKey<const L: usize> {
    z: [u8; 32],
    hash_pk: [u8; 32],
    cpa_pk: IndCpaPublicKey<L>,
    cpa_sk: IndCpaSecretKey<L>,
}

impl<const L: usize> IndCcaSecretKey<L> {
    pub(crate) fn serialized_len() -> usize {
        32 + 32 + IndCpaPublicKey::<L>::serialized_len() + IndCpaSecretKey::<L>::serialized_len()
    }

    pub(crate) fn to_bytes(&self, out_buf: &mut [u8]) {
        assert_eq!(out_buf.len(), Self::serialized_len());

        let mut idx = 32;
        out_buf[..idx].copy_from_slice(&self.z);

        let delta = 32;
        out_buf[idx..idx + delta].copy_from_slice(&self.hash_pk);
        idx += delta;

        let delta = IndCpaPublicKey::<L>::serialized_len();
        self.cpa_pk.to_bytes(&mut out_buf[idx..idx + delta]);
        idx += delta;

        let delta = IndCpaSecretKey::<L>::serialized_len();
        self.cpa_sk.to_bytes(&mut out_buf[idx..idx + delta]);
        idx += delta;
    }
}

pub fn gen_keypair<const L: usize, const MU: usize>(
    rng: &mut impl CryptoRngCore,
) -> (IndCcaSecretKey<L>, IndCcaPublicKey<L>) {
    let (sk, pk) = ind_cpa::gen_keypair::<L, MU>(rng);

    let mut buf = [0u8; 32 + MAX_L * MODULUS_P_BITS * RING_DEG / 8];
    let pk_bytes = &mut buf[..32 + L * MODULUS_P_BITS * RING_DEG / 8];
    pk.to_bytes(pk_bytes);
    let hash_pk = Sha3_256::digest(pk_bytes);

    let mut z = [0u8; 32];
    rng.fill_bytes(&mut z);

    let sk = IndCcaSecretKey {
        z,
        hash_pk: hash_pk.into(),
        cpa_pk: pk.clone(),
        cpa_sk: sk,
    };

    (sk, pk)
}

pub fn encap<const L: usize, const MU: usize, const MODULUS_T_BITS: usize>(
    rng: &mut impl CryptoRngCore,
    pk: &IndCcaPublicKey<L>,
    out_buf: &mut [u8],
) -> SharedSecret {
    let mut m = [0u8; 32];
    rng.fill_bytes(&mut m);

    // Hash the public key
    let mut buf = [0u8; 32 + MAX_L * MODULUS_P_BITS * RING_DEG / 8];
    let pk_bytes = &mut buf[..32 + L * MODULUS_P_BITS * RING_DEG / 8];
    pk.to_bytes(pk_bytes);
    let hash_pk = Sha3_256::digest(pk_bytes);

    // rk = SHA3-512(hash_pk || m)
    let rk = {
        let mut h = sha3::Sha3_512::new();
        h.update(hash_pk);
        h.update(m);
        h.finalize()
    };

    // r || k = rk
    let (r, k) = rk.split_at(32);
    let r: [u8; 32] = r.try_into().unwrap();
    let k: [u8; 32] = k.try_into().unwrap();

    // ct = Saber.PKE.Enc_pk(m; r)
    ind_cpa::enc_deterministic::<L, MU, MODULUS_T_BITS>(pk, &m, &r, out_buf);

    // r' = SHA3-256(ct)
    let rprime = Sha3_256::digest(out_buf);

    // session key = SHA3-256(r' || k)
    let sess_key = {
        let mut h = Sha3_256::new();
        h.update(rprime);
        h.update(k);
        h.finalize()
    };
    sess_key.into()
}

pub fn decap<const L: usize, const MU: usize, const MODULUS_T_BITS: usize>(
    sk: &IndCcaSecretKey<L>,
    ciphertext: &[u8],
) -> SharedSecret {
    assert_eq!(
        ciphertext.len(),
        MODULUS_T_BITS * RING_DEG / 8 + L * MODULUS_P_BITS * RING_DEG / 8
    );

    let m = ind_cpa::dec::<L, MODULUS_T_BITS>(&sk.cpa_sk, ciphertext);

    // rk = SHA3-512(hash_pk || m)
    let rk = {
        let mut h = sha3::Sha3_512::new();
        h.update(sk.hash_pk);
        h.update(m);
        h.finalize()
    };

    // r || k = rk
    let (r, k) = rk.split_at(32);
    let r: [u8; 32] = r.try_into().unwrap();
    let k: [u8; 32] = k.try_into().unwrap();

    // ct' = Saber.PKE.Enc_pk(m; r)
    let mut buf = [0u8; MAX_MODULUS_T_BITS * RING_DEG / 8 + MAX_L * MODULUS_P_BITS * RING_DEG / 8];
    let reconstructed_ct =
        &mut buf[..MODULUS_T_BITS * RING_DEG / 8 + L * MODULUS_P_BITS * RING_DEG / 8];
    ind_cpa::enc_deterministic::<L, MU, MODULUS_T_BITS>(&sk.cpa_pk, &m, &r, reconstructed_ct);

    // r' = SHA3-256(ct)
    let rprime = Sha3_256::digest(ciphertext);

    // suffix = k if reconstruction matched, else z
    let reconstruction_matched = reconstructed_ct.ct_eq(ciphertext);
    let mut suffix = [0u8; 32];
    for i in 0..32 {
        suffix[i] = u8::conditional_select(&sk.z[i], &k[i], reconstruction_matched);
    }

    // session key = SHA3-256(r' || suffix)
    let sess_key = {
        let mut h = Sha3_256::new();
        h.update(rprime);
        h.update(suffix);
        h.finalize()
    };
    sess_key.into()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn cca_correctness() {
        test_encap_decap::<2, 3, 10>();
        test_encap_decap::<3, 4, 8>();
        test_encap_decap::<4, 6, 6>();
    }

    fn test_encap_decap<const L: usize, const MODULUS_T_BITS: usize, const MU: usize>() {
        let mut rng = rand::thread_rng();

        for _ in 0..100 {
            let (sk, pk) = gen_keypair::<L, MU>(&mut rng);
            let mut ct_buf =
                vec![0u8; MODULUS_T_BITS * RING_DEG / 8 + L * MODULUS_P_BITS * RING_DEG / 8];

            let ss1 = encap::<L, MU, MODULUS_T_BITS>(&mut rng, &pk, &mut ct_buf);
            let ss2 = decap::<L, MU, MODULUS_T_BITS>(&sk, &ct_buf);
            assert_eq!(ss1, ss2);
        }
    }
}
