use saber_kem::lightsaber::{
    LightsaberCiphertext, LightsaberPublicKey, LightsaberSecretKey, LIGHTSABER_CIPHERTEXT_LEN,
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
    assert_eq!(
        slice_containing_sk.len(),
        LightsaberSecretKey::SERIALIZED_LEN
    );
    let sk_arr = slice_containing_sk.try_into().unwrap();
    let sk = LightsaberSecretKey::from_bytes(sk_arr);

    // Also serialize and deserialize the public key
    let mut pk_bytes = [0u8; LightsaberPublicKey::SERIALIZED_LEN];
    pk.to_bytes(&mut pk_bytes);
    let slice_containing_pk = pk_bytes.as_slice();
    // The API only accepts fixed-len slices, so we have to cast it first
    assert_eq!(
        slice_containing_pk.len(),
        LightsaberPublicKey::SERIALIZED_LEN
    );
    let pk_arr = slice_containing_pk.try_into().unwrap();
    let pk = LightsaberPublicKey::from_bytes(pk_arr);

    // Encapsulate a shared secret, ss1, to pk
    let (_ct, _ss1) = pk.encapsulate(&mut rng);
    // Alternatively, if you have a buffer and want to avoid an extra allocation, encapsulate in
    // place. LightSaberCiphertext is just a byte array, so no conversion necessary:
    let mut ct = [0u8; LIGHTSABER_CIPHERTEXT_LEN];
    let ss1 = pk.encapsulate_in_place(&mut rng, &mut ct);
    let slice_containing_ct = ct.as_slice();

    // Deserializing is also straightforward
    assert_eq!(slice_containing_ct.len(), LIGHTSABER_CIPHERTEXT_LEN);
    let receiver_ct: &LightsaberCiphertext = slice_containing_ct.try_into().unwrap();

    // Use the secret key to decapsulate the ciphertext
    let ss2 = sk.decapsulate(receiver_ct);

    // Check the shared secrets are equal. NOTE is not a constant-time check (ie not secure). We
    // only do this for testing purposes.
    assert_eq!(ss1.as_bytes(), ss2.as_bytes());

    println!("KEM ran successfully");
}
