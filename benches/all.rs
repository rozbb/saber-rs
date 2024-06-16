use criterion::{criterion_group, criterion_main, Criterion};
use rand::thread_rng;
use saber::{ind_cca, RingElem, MODULUS_P_BITS, RING_DEG};

pub fn mul(c: &mut Criterion) {
    let mut rng = thread_rng();
    c.bench_function("schoolbook-mul", |b| {
        b.iter(|| {
            let r1 = RingElem::rand(&mut rng);
            let r2 = RingElem::rand(&mut rng);
            r1.schoolbook_mul(&r2)
        })
    });
    c.bench_function("karatsuba-mul", |b| {
        b.iter(|| {
            let r1 = RingElem::rand(&mut rng);
            let r2 = RingElem::rand(&mut rng);
            r1.karatsuba_mul(&r2)
        })
    });
}

pub fn cca(c: &mut Criterion) {
    const L: usize = 2;
    const MODULUS_T_BITS: usize = 3;
    const MU: usize = 10;
    // 3,4,8 and 4,6,6 are the other param sets

    let mut rng = rand::thread_rng();

    c.bench_function("gen-keypair", |b| {
        b.iter(|| ind_cca::gen_keypair::<L, MU>(&mut rng))
    });
    let (sk, pk) = ind_cca::gen_keypair::<L, MU>(&mut rng);

    let mut ct_buf = vec![0u8; MODULUS_T_BITS * RING_DEG / 8 + L * MODULUS_P_BITS * RING_DEG / 8];

    c.bench_function("encap", |b| {
        b.iter(|| ind_cca::encap::<L, MU, MODULUS_T_BITS>(&mut rng, &pk, &mut ct_buf))
    });
    ind_cca::encap::<L, MU, MODULUS_T_BITS>(&mut rng, &pk, &mut ct_buf);

    c.bench_function("decap", |b| {
        b.iter(|| ind_cca::decap::<L, MU, MODULUS_T_BITS>(&sk, &ct_buf))
    });
}

criterion_group!(benches, cca);
criterion_main!(benches);
