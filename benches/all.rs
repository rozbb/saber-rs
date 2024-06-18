use criterion::{criterion_group, criterion_main, Criterion};
use saber::{ciphertext_len, kem};

pub fn cca(c: &mut Criterion) {
    const L: usize = 2;
    const MODULUS_T_BITS: usize = 3;
    const MU: usize = 10;
    // 3,4,8 and 4,6,6 are the other param sets

    let mut rng = rand::thread_rng();

    c.bench_function("gen-keypair", |b| {
        b.iter(|| kem::gen_keypair::<L, MU>(&mut rng))
    });
    let (sk, pk) = kem::gen_keypair::<L, MU>(&mut rng);

    let mut ct_buf = vec![0u8; ciphertext_len::<L, MODULUS_T_BITS>()];

    c.bench_function("encap", |b| {
        b.iter(|| kem::encap::<L, MU, MODULUS_T_BITS>(&mut rng, &pk, &mut ct_buf))
    });
    kem::encap::<L, MU, MODULUS_T_BITS>(&mut rng, &pk, &mut ct_buf);

    c.bench_function("decap", |b| {
        b.iter(|| kem::decap::<L, MU, MODULUS_T_BITS>(&sk, &ct_buf))
    });
}

criterion_group!(benches, cca);
criterion_main!(benches);
