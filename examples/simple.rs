use saber_kem::{
    kem_traits::{Decapsulate, Encapsulate},
    lightsaber::{LightsaberCiphertext, LightsaberPublicKey, LightsaberSecretKey},
};

fn main() {
    let mut rng = rand::thread_rng();

    // Generate a keypair
    let sk = LightsaberSecretKey::generate(&mut rng);
    let pk = sk.public_key();

    // Serialize the secret key, maybe to save on disk
    let mut sk_bytes = [0u8; LightsaberSecretKey::SERIALIZED_LEN];
    sk.to_bytes(&mut sk_bytes);
    let slice_containing_sk = sk_bytes.as_slice();

    // Deserialize the secret key
    // The API only accepts fixed-len slices, so we have to cast it first
    let sk_arr = slice_containing_sk[..LightsaberSecretKey::SERIALIZED_LEN]
        .try_into()
        .unwrap();
    let sk = LightsaberSecretKey::from_bytes(sk_arr);

    // Also serialize and deserialize the public key
    let mut pk_bytes = [0u8; LightsaberPublicKey::SERIALIZED_LEN];
    pk.to_bytes(&mut pk_bytes);
    let slice_containing_pk = pk_bytes.as_slice();
    // The API only accepts fixed-len slices, so we have to cast it first
    let pk_arr = slice_containing_pk[..LightsaberPublicKey::SERIALIZED_LEN]
        .try_into()
        .unwrap();
    let pk = LightsaberPublicKey::from_bytes(pk_arr);

    // Encapsulate a shared secret, ss1, to pk
    let (ct, ss1) = pk.encapsulate(&mut rng).unwrap();
    // The ciphertext is just bytes, so serializing is straightforward
    let ct_bytes = ct.as_ref();

    // Deserializing is also straightforward
    assert_eq!(ct_bytes.len(), LightsaberCiphertext::LEN);
    let mut receiver_ct = LightsaberCiphertext::default();
    receiver_ct.as_mut().copy_from_slice(ct_bytes);

    // Use the secret key to decapsulate the ciphertext
    let ss2 = sk.decapsulate(&receiver_ct).unwrap();

    // Check the shared secrets are equal. NOTE is not a constant-time check (ie not secure). We
    // only do this for testing purposes.
    assert_eq!(ss1.as_bytes(), ss2.as_bytes());

    println!("KEM ran successfully");
}
