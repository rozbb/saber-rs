use criterion::{criterion_group, criterion_main, Criterion};
use kopis_kem::kopis768::{Kopis768PublicKey, Kopis768SecretKey, KOPIS768_CIPHERTEXT_LEN};
use rand::{CryptoRng, Rng};
use sha3::{Digest, Sha3_512};
use shake::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
use subtle::ConstantTimeEq;

type Msg1 = [u8; Kopis768PublicKey::SERIALIZED_LEN + 32];
type Msg2 = [u8; KOPIS768_CIPHERTEXT_LEN + 32];

/// Outputs the first message of NoIC and the secret key
fn noic_init_start(rng: &mut impl Rng, sid: &[u8; 32], pw: &[u8; 32]) -> (Msg1, Kopis768SecretKey) {
    let sk: [u8; 32] = rng.random();
    let decap_key = Kopis768SecretKey::expand_from_seed(&sk);
    let encap_key = decap_key.public_key();

    let mut encap_key_bytes = [0u8; Kopis768PublicKey::SERIALIZED_LEN];
    encap_key.serialize(&mut encap_key_bytes);

    let r: [u8; 32] = rng.random();
    // R = G(sid || pw || r)
    let R = {
        let mut h = Shake256::default();
        h.update(sid);
        h.update(pw);
        h.update(&r);

        let mut buf = [0u8; Kopis768PublicKey::SERIALIZED_LEN];
        let mut reader = h.finalize_xof();
        reader.read(&mut buf);
        buf
    };

    // T = R ^ pk
    let mut T = encap_key_bytes;
    for (T_byte, mask_byte) in T.iter_mut().zip(R) {
        *T_byte ^= mask_byte;
    }

    // t = H(sid || pw || T)
    let t = sha3::Sha3_256::new()
        .chain_update(sid)
        .chain_update(pw)
        .chain_update(&T)
        .finalize();

    // s = t ^ r
    let mut s = t;
    for (s_byte, mask_byte) in s.iter_mut().zip(r) {
        *s_byte ^= mask_byte;
    }

    // Output s || T
    let mut out_buf = [0u8; Kopis768PublicKey::SERIALIZED_LEN + 32];
    out_buf[..32].copy_from_slice(&s);
    out_buf[32..].copy_from_slice(&T);

    (out_buf, decap_key)
}

/// Returns the response message and the final key k
fn noic_resp(
    rng: &mut impl CryptoRng,
    sid: &[u8; 32],
    pw: &[u8; 32],
    msg1: &Msg1,
) -> (Msg2, [u8; 32]) {
    let (s, T) = msg1.split_at(32);

    // t = H(sid || pw || T)
    let t = sha3::Sha3_256::new()
        .chain_update(sid)
        .chain_update(pw)
        .chain_update(&T)
        .finalize();

    // r = t ^ s
    let mut r = t;
    for (r_byte, s_byte) in r.iter_mut().zip(s) {
        *r_byte ^= s_byte;
    }

    // R = G(sid || pw || r)
    let R = {
        let mut h = Shake256::default();
        h.update(sid);
        h.update(pw);
        h.update(&r);

        let mut buf = [0u8; Kopis768PublicKey::SERIALIZED_LEN];
        let mut reader = h.finalize_xof();
        reader.read(&mut buf);
        buf
    };

    // pk = R ^ T
    let mut pk_bytes = R;
    for (pk_byte, T_byte) in pk_bytes.iter_mut().zip(T) {
        *pk_byte ^= T_byte;
    }

    let pk = Kopis768PublicKey::from_bytes(&pk_bytes);
    let (ct, k_s) = pk.encapsulate(rng);

    // tag || K = H(K_s,sid,pw,pk,apk,cph)
    let h = Sha3_512::new()
        .chain_update(k_s.as_bytes())
        .chain_update(sid)
        .chain_update(pw)
        .chain_update(pk_bytes)
        .chain_update(msg1)
        .chain_update(&ct)
        .finalize();
    let (tag, k) = h.split_at(32);

    // Output tag || ct
    let mut out = [0u8; KOPIS768_CIPHERTEXT_LEN + 32];
    out[..32].copy_from_slice(tag);
    out[32..].copy_from_slice(&ct);

    (out, k.try_into().unwrap())
}

/// Processes the response message the returns the shared secret or panics on tag error
fn noic_init_end(
    sk: &Kopis768SecretKey,
    sid: &[u8; 32],
    pw: &[u8; 32],
    msg1: &Msg1,
    msg2: &Msg2,
) -> [u8; 32] {
    let (tag, ct) = msg2.split_at(32);
    let k_s = sk.decapsulate(ct.try_into().unwrap());

    let mut pk_bytes = [0u8; Kopis768PublicKey::SERIALIZED_LEN];
    sk.public_key().serialize(&mut pk_bytes);

    // tag || K = H(K_s,sid,pw,pk,apk,cph)
    let h = Sha3_512::new()
        .chain_update(k_s.as_bytes())
        .chain_update(sid)
        .chain_update(pw)
        .chain_update(pk_bytes)
        .chain_update(msg1)
        .chain_update(&ct)
        .finalize();
    let (computed_tag, k) = h.split_at(32);

    if bool::from(tag.ct_eq(computed_tag)) {
        return k.try_into().unwrap();
    } else {
        panic!("invalid tag");
    }
}

fn bench(c: &mut Criterion) {
    use rand::Rng;
    let mut rng = rand::rng();

    let sid: [u8; 32] = rng.random();
    let pw: [u8; 32] = rng.random();

    // Sanity check
    let (msg1, sk) = noic_init_start(&mut rng, &sid, &pw);
    let (msg2, k1) = noic_resp(&mut rng, &sid, &pw, &msg1);
    let k2 = noic_init_end(&sk, &sid, &pw, &msg1, &msg2);
    assert_eq!(k1, k2);

    c.bench_function("noic-kopis-initStart", |b| {
        b.iter(|| noic_init_start(&mut rng, &sid, &pw))
    });

    c.bench_function("noic-kopis-resp", |b| {
        b.iter(|| noic_resp(&mut rng, &sid, &pw, &msg1))
    });

    c.bench_function("noic-kopis-initEnd", |b| {
        b.iter(|| noic_init_end(&sk, &sid, &pw, &msg1, &msg2))
    });
}

criterion_group!(noic, bench);
criterion_main!(noic);
