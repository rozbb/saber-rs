use kopis_kem::{
    kopis1024::Kopis1024SecretKey, kopis512::Kopis512SecretKey, kopis768::Kopis768SecretKey,
};

use criterion::{criterion_group, criterion_main, Criterion};

macro_rules! bench_kopis_variant {
    ($bench_name:ident, $privkey_name:ident) => {
        fn $bench_name(c: &mut Criterion) {
            let randomness = &[0u8; 32];

            let gen_bench_name = format!("{}-gen-keypair-derand", stringify!($bench_name));
            c.bench_function(&gen_bench_name, |b| {
                b.iter(|| $privkey_name::expand_from_seed(randomness))
            });
            let sk = $privkey_name::expand_from_seed(randomness);
            let pk = sk.public_key();

            let encap_bench_name = format!("{}-encap-derand", stringify!($bench_name));
            c.bench_function(&encap_bench_name, |b| {
                b.iter(|| pk.encapsulate_deterministic(randomness))
            });
            let (ct, _) = pk.encapsulate_deterministic(randomness);

            let decap_bench_name = format!("{}-decap", stringify!($bench_name));
            c.bench_function(&decap_bench_name, |b| b.iter(|| sk.decapsulate(&ct)));
        }
    };
}

macro_rules! bench_libcrux_variant {
    ($bench_name:ident, $mod_name:path) => {
        fn $bench_name(c: &mut Criterion) {
            use $mod_name as base_mod;

            use base_mod::portable::unpacked::*;

            let kg_randomness = [0u8; 64];
            let encap_randomness = [0u8; 32];

            let gen_bench_name = format!("{}-gen-keypair-derand", stringify!($bench_name));
            c.bench_function(&gen_bench_name, |b| {
                b.iter(|| generate_key_pair(kg_randomness))
            });
            let kp = generate_key_pair(kg_randomness);
            let pk = kp.public_key();

            let encap_bench_name = format!("{}-encap-derand", stringify!($bench_name));
            c.bench_function(&encap_bench_name, |b| {
                b.iter(|| encapsulate(&pk, encap_randomness))
            });
            let (ct, _) = encapsulate(&pk, encap_randomness);

            let decap_bench_name = format!("{}-decap", stringify!($bench_name));
            c.bench_function(&decap_bench_name, |b| b.iter(|| decapsulate(&kp, &ct)));
        }
    };
}

bench_libcrux_variant!(libcruxmlkem512, libcrux_ml_kem::mlkem512);
bench_libcrux_variant!(libcruxmlkem768, libcrux_ml_kem::mlkem768);
bench_libcrux_variant!(libcruxmlkem1024, libcrux_ml_kem::mlkem1024);

bench_kopis_variant!(kopis512, Kopis512SecretKey);
bench_kopis_variant!(kopis768, Kopis768SecretKey);
bench_kopis_variant!(kopis1024, Kopis1024SecretKey);

fn graviolamlkem768(c: &mut Criterion) {
    use graviola::key_agreement::mlkem768::*;

    let kg_randomness = [0u8; 64];
    let encap_randomness = [0u8; 32];

    c.bench_function("graviolamlkem768-gen-kepair-derand", |b| {
        b.iter(|| DecapKey::keygen_internal(&kg_randomness))
    });

    let sk = DecapKey::generate().unwrap();
    let pk = sk.encapsulation_key();

    c.bench_function("graviolamlkem768-encap-derand", |b| {
        b.iter(|| pk.clone().encaps_internal(Message(encap_randomness)))
    });
    let (_, ct) = pk.encaps().unwrap();

    c.bench_function("graviolamlkem768-decap", |b| {
        b.iter(|| sk.decaps_internal(&ct))
    });
}

criterion_group!(kopis_benches, kopis512, kopis768, kopis1024);
criterion_group!(graviola_benches, graviolamlkem768);
criterion_group!(
    libcrux_benches,
    libcruxmlkem512,
    libcruxmlkem768,
    libcruxmlkem1024
);

criterion_main!(kopis_benches, libcrux_benches, graviola_benches);
