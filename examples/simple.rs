use kopis_kem::{
    kopis512::{Kopis512Ciphertext, Kopis512PublicKey, Kopis512SecretKey, KOPIS512_CIPHERTEXT_LEN},
    SharedSecret,
};

fn main() {
    let mut rng = rand::rng();

    // Generate a keypair
    let sk = Kopis512SecretKey::generate(&mut rng);
    let pk = sk.public_key();

    // Serialize the secret key, maybe to save on disk
    let sk_seed = sk.seed();

    // Deserialize the secret key
    // The API only accepts fixed-len slices, so we have to cast it first
    let sk = Kopis512SecretKey::expand_from_seed(&sk_seed);

    // Also serialize and deserialize the public key
    let mut pk_bytes = [0u8; Kopis512PublicKey::SERIALIZED_LEN];
    pk.serialize(&mut pk_bytes);
    let slice_containing_pk = pk_bytes.as_slice();
    // The API only accepts fixed-len slices, so we have to cast it first
    assert_eq!(slice_containing_pk.len(), Kopis512PublicKey::SERIALIZED_LEN);
    let pk_arr = slice_containing_pk.try_into().unwrap();
    let pk = Kopis512PublicKey::from_bytes(pk_arr);

    // Encapsulate a shared secret, ss1, to pk
    let (ct, ss1): (Kopis512Ciphertext, SharedSecret) = pk.encapsulate(&mut rng);
    let slice_containing_ct = ct.as_slice();

    // Deserializing is also straightforward
    assert_eq!(slice_containing_ct.len(), KOPIS512_CIPHERTEXT_LEN);
    let receiver_ct: &Kopis512Ciphertext = slice_containing_ct.try_into().unwrap();

    // Use the secret key to decapsulate the ciphertext
    let ss2 = sk.decapsulate(receiver_ct);

    // Check the shared secrets are equal. NOTE is not a constant-time check (ie not secure). We
    // only do this for testing purposes.
    assert_eq!(ss1.as_bytes(), ss2.as_bytes());

    println!("KEM ran successfully");
}
