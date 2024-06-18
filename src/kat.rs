//! This module contains code for running known-answer tests (KATs)

use crate::{
    kem::{decap, encap, IndCcaSecretKey, KemPublicKey},
    pke::ciphertext_len,
};

use std::{string::String, vec::Vec};

use rand_core::{CryptoRng, RngCore};
use serde::{Deserialize, Deserializer, Serialize};

/// The RNG we'll use for KATs. This is an AES-256 CTR DRBG with no personalization strings
struct KatRng(aes_ctr_drbg::DrbgCtx);

impl KatRng {
    fn new(seed: &[u8]) -> Self {
        assert_eq!(seed.len(), 48);
        let mut rng = aes_ctr_drbg::DrbgCtx::new();
        rng.init(&seed, Vec::new());
        KatRng(rng)
    }
}

#[derive(Serialize, Deserialize)]
struct TestVector {
    count: usize,
    #[serde(deserialize_with = "bytes_from_hex")]
    seed: Vec<u8>,
    #[serde(deserialize_with = "bytes_from_hex")]
    pk: Vec<u8>,
    #[serde(deserialize_with = "bytes_from_hex")]
    sk: Vec<u8>,
    #[serde(deserialize_with = "bytes_from_hex")]
    ct: Vec<u8>,
    #[serde(deserialize_with = "bytes_from_hex")]
    ss: Vec<u8>,
}

// Tells serde how to deserialize bytes from the hex representation
fn bytes_from_hex<'de, D>(deserializer: D) -> Result<Vec<u8>, D::Error>
where
    D: Deserializer<'de>,
{
    let mut hex_str = String::deserialize(deserializer)?;
    // Prepend a 0 if it's not even length
    if hex_str.len() % 2 == 1 {
        hex_str.insert(0, '0');
    }
    hex::decode(hex_str).map_err(|e| serde::de::Error::custom(format!("{:?}", e)))
}

// Impl rand_core traits so we can use them with our crate
impl RngCore for KatRng {
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        self.0.get_random(dest)
    }
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.fill_bytes(dest);
        Ok(())
    }

    fn next_u32(&mut self) -> u32 {
        unimplemented!()
    }

    fn next_u64(&mut self) -> u64 {
        unimplemented!()
    }
}

impl CryptoRng for KatRng {}

#[test]
fn drbg() {
    // We have a known vector for seed = [1, 2, 3, ...]
    let mut seed = [0u8; 48];
    for i in 0..48 {
        seed[i] = i as u8;
    }

    let mut rng = KatRng::new(&seed[..]);

    let mut buf = [0u8; 256];
    rng.fill_bytes(&mut buf);

    // Reference string from rusty-saber
    // https://github.com/lkiem/rusty_saber/blob/eac62d05ac4b4b009c82e0bceb02393ccede08e8/src/rng.rs#L152
    let ref1 = [
        0x06u8, 0x15, 0x50, 0x23, 0x4D, 0x15, 0x8C, 0x5E, 0xC9, 0x55, 0x95, 0xFE, 0x04, 0xEF, 0x7A,
        0x25, 0x76, 0x7F, 0x2E, 0x24, 0xCC, 0x2B, 0xC4, 0x79, 0xD0, 0x9D, 0x86, 0xDC, 0x9A, 0xBC,
        0xFD, 0xE7, 0x05, 0x6A, 0x8C, 0x26, 0x6F, 0x9E, 0xF9, 0x7E, 0xD0, 0x85, 0x41, 0xDB, 0xD2,
        0xE1, 0xFF, 0xA1, 0x98, 0x10, 0xF5, 0x39, 0x2D, 0x07, 0x62, 0x76, 0xEF, 0x41, 0x27, 0x7C,
        0x3A, 0xB6, 0xE9, 0x4A, 0x4E, 0x3B, 0x7D, 0xCC, 0x10, 0x4A, 0x05, 0xBB, 0x08, 0x9D, 0x33,
        0x8B, 0xF5, 0x5C, 0x72, 0xCA, 0xB3, 0x75, 0x38, 0x9A, 0x94, 0xBB, 0x92, 0x0B, 0xD5, 0xD6,
        0xDC, 0x9E, 0x7F, 0x2E, 0xC6, 0xFD, 0xE0, 0x28, 0xB6, 0xF5, 0x72, 0x4B, 0xB0, 0x39, 0xF3,
        0x65, 0x2A, 0xD9, 0x8D, 0xF8, 0xCE, 0x6C, 0x97, 0x01, 0x32, 0x10, 0xB8, 0x4B, 0xBE, 0x81,
        0x38, 0x8C, 0x3D, 0x14, 0x1D, 0x61, 0x95, 0x7C, 0x73, 0xBC, 0xDC, 0x5E, 0x5C, 0xD9, 0x25,
        0x25, 0xF4, 0x6A, 0x2B, 0x75, 0x7B, 0x03, 0xCA, 0xB5, 0xC3, 0x37, 0x00, 0x4A, 0x2D, 0xA3,
        0x53, 0x24, 0xA3, 0x25, 0x71, 0x35, 0x64, 0xDA, 0xE2, 0x8F, 0x57, 0xAC, 0xC6, 0xDB, 0xE3,
        0x2A, 0x07, 0x26, 0x19, 0x0B, 0xAA, 0x6B, 0x8A, 0x0A, 0x25, 0x5A, 0xA1, 0xAD, 0x01, 0xE8,
        0xDD, 0x56, 0x9A, 0xA3, 0x6D, 0x09, 0x62, 0x56, 0xC4, 0x20, 0x71, 0x8A, 0x69, 0xD4, 0x6D,
        0x8D, 0xB1, 0xC6, 0xDD, 0x40, 0x60, 0x6A, 0x0B, 0xE3, 0xC2, 0x35, 0xBE, 0xFE, 0x62, 0x3A,
        0x90, 0x59, 0x3F, 0x82, 0xD6, 0xA8, 0xF9, 0xF9, 0x24, 0xE4, 0x4E, 0x36, 0xBE, 0x87, 0xF7,
        0xD2, 0x6B, 0x84, 0x45, 0x96, 0x6F, 0x9E, 0xE3, 0x29, 0xC4, 0x26, 0xC1, 0x25, 0x21, 0xE8,
        0x5F, 0x6F, 0xD4, 0xEC, 0xD5, 0xD5, 0x66, 0xBA, 0x0A, 0x34, 0x87, 0x12, 0x5D, 0x79, 0xCC,
        0x64,
    ];
    assert_eq!(buf, ref1);
}

/// Opens the given RSP file and tests our impl against the given vectors
fn kat_helper<const L: usize, const MODULUS_BITS_T: usize, const MU: usize>(filename: &str) {
    let kat_json_file = std::fs::File::open(filename).unwrap();
    let test_vectors = serde_json::from_reader::<_, Vec<TestVector>>(kat_json_file).unwrap();

    for tv in test_vectors {
        let mut rng = KatRng::new(&tv.seed);
        let (sk, pk) = crate::kem::gen_keypair::<L, MU>(&mut rng);

        let mut sk_buf = vec![0u8; IndCcaSecretKey::<L>::serialized_len()];
        sk.to_bytes(&mut sk_buf);
        assert_eq!(sk_buf, tv.sk, "secret keys do not match");

        let mut pk_buf = vec![0u8; KemPublicKey::<L>::serialized_len()];
        pk.to_bytes(&mut pk_buf);
        assert_eq!(pk_buf, tv.pk, "public keys do not match");

        let mut ct_buf = vec![0u8; ciphertext_len::<L, MODULUS_BITS_T>()];
        let ss1 = encap::<L, MU, MODULUS_BITS_T>(&mut rng, &pk, &mut ct_buf);
        assert_eq!(ct_buf, tv.ct, "ciphertexts do no match");

        let ss2 = decap::<L, MU, MODULUS_BITS_T>(&sk, &ct_buf);
        assert_eq!(ss1, ss2);
        assert_eq!(ss1, tv.ss.as_slice(), "shared secrets do not match");
    }
}

#[test]
fn kat_lightsaber() {
    let filename = "PQCkemKAT_1568.rsp.json";
    const L: usize = 2;
    const MODULUS_BITS_T: usize = 3;
    const MU: usize = 10;

    kat_helper::<L, MODULUS_BITS_T, MU>(filename);
}

#[test]
fn kat_saber() {
    let filename = "PQCkemKAT_2304.rsp.json";
    const L: usize = 3;
    const MODULUS_BITS_T: usize = 4;
    const MU: usize = 8;

    kat_helper::<L, MODULUS_BITS_T, MU>(filename);
}

#[test]
fn kat_firesaber() {
    let filename = "PQCkemKAT_3040.rsp.json";
    const L: usize = 4;
    const MODULUS_BITS_T: usize = 6;
    const MU: usize = 6;

    kat_helper::<L, MODULUS_BITS_T, MU>(filename);
}
