//! ML-KEM-768 round-trip benchmark.
//!
//! Source <https://github.com/FiloSottile/mlkem768/blob/main/mlkem768_test.go> BenchmarkRoundTrip.

use criterion::{
    BenchmarkGroup, Criterion, black_box, criterion_group, criterion_main, measurement::WallTime,
};
use libcrux_ml_kem::mlkem768 as libcrux_mlkem768;
use ml_kem::{Decapsulate, Encapsulate, Kem, MlKem768};
use ml_kem_example::kem::mlkem::EncapsulationKey;
use sha3::{
    Shake128,
    digest::{ExtendableOutput, Update, XofReader},
};

/// Deterministic 32-byte stream so each iteration uses fresh-looking randomness without any
/// syscalls or crypto-rng dependency.
struct Rng(sha3::Shake128Reader);

impl Rng {
    fn new(label: &[u8]) -> Self {
        let mut x = Shake128::default();
        x.update(label);
        Self(x.finalize_xof())
    }

    fn bytes32(&mut self) -> [u8; 32] {
        let mut out = [0u8; 32];
        self.0.read(&mut out);
        out
    }

    fn bytes64(&mut self) -> [u8; 64] {
        let mut out = [0u8; 64];
        self.0.read(&mut out);
        out
    }
}

impl rand_core::TryRng for Rng {
    type Error = core::convert::Infallible;

    fn try_next_u32(&mut self) -> Result<u32, Self::Error> {
        let mut buf = [0u8; 4];
        self.0.read(&mut buf);
        Ok(u32::from_le_bytes(buf))
    }

    fn try_next_u64(&mut self) -> Result<u64, Self::Error> {
        let mut buf = [0u8; 8];
        self.0.read(&mut buf);
        Ok(u64::from_le_bytes(buf))
    }

    fn try_fill_bytes(&mut self, dst: &mut [u8]) -> Result<(), Self::Error> {
        self.0.read(dst);
        Ok(())
    }
}

impl rand_core::TryCryptoRng for Rng {}

fn bench_alice(group: &mut BenchmarkGroup<WallTime>) {
    // --- this crate ---
    let mut setup = Rng::new(b"mlkem768-bench-setup");
    let (ek, dk) = EncapsulationKey::generate(&setup.bytes32(), &setup.bytes32());
    let ct = ek.encaps(&setup.bytes32()).1;

    group.bench_function("this", |b| {
        let mut rng = Rng::new(b"mlkem768-bench-alice");
        b.iter(|| {
            let (ek_s, _) = EncapsulationKey::generate(&rng.bytes32(), &rng.bytes32());
            black_box(&ek_s);
            black_box(dk.decaps(black_box(&ct)));
        });
    });

    // --- ml-kem ---
    let mut setup = Rng::new(b"mlkem768-bench-setup");
    let (ref_dk, ref_ek) = MlKem768::generate_keypair_from_rng(&mut setup);
    let ref_ct = ref_ek.encapsulate_with_rng(&mut setup).0;

    group.bench_function("ml-kem", |b| {
        let mut rng = Rng::new(b"mlkem768-bench-alice");
        b.iter(|| {
            let (ek_s, _) = MlKem768::generate_keypair_from_rng(&mut rng);
            black_box(&ek_s);
            black_box(ref_dk.decapsulate(black_box(&ref_ct)));
        });
    });

    // --- libcrux-ml-kem ---
    let mut setup = Rng::new(b"mlkem768-bench-setup");
    let ref_kp = libcrux_mlkem768::generate_key_pair(setup.bytes64());
    let ref_ct = libcrux_mlkem768::encapsulate(ref_kp.public_key(), setup.bytes32()).0;

    group.bench_function("libcrux-ml-kem", |b| {
        let mut rng = Rng::new(b"mlkem768-bench-alice");
        b.iter(|| {
            let kp = libcrux_mlkem768::generate_key_pair(rng.bytes64());
            black_box(&kp);
            black_box(libcrux_mlkem768::decapsulate(
                ref_kp.private_key(),
                black_box(&ref_ct),
            ));
        });
    });
}

fn bench_bob(group: &mut BenchmarkGroup<WallTime>) {
    // --- this crate ---
    let mut setup = Rng::new(b"mlkem768-bench-setup");
    let (ek, _dk) = EncapsulationKey::generate(&setup.bytes32(), &setup.bytes32());

    group.bench_function("this", |b| {
        let mut rng = Rng::new(b"mlkem768-bench-bob");
        b.iter(|| {
            let (ks, cs) = ek.encaps(black_box(&rng.bytes32()));
            black_box((ks, cs));
        });
    });

    // --- ml-kem ---
    let mut setup = Rng::new(b"mlkem768-bench-setup");
    let (_, ref_ek) = MlKem768::generate_keypair_from_rng(&mut setup);

    group.bench_function("ml-kem", |b| {
        let mut rng = Rng::new(b"mlkem768-bench-bob");
        b.iter(|| {
            let (k, c) = ref_ek.encapsulate_with_rng(&mut rng);
            black_box((k, c));
        });
    });

    // --- libcrux-ml-kem ---
    let mut setup = Rng::new(b"mlkem768-bench-setup");
    let ref_kp = libcrux_mlkem768::generate_key_pair(setup.bytes64());

    group.bench_function("libcrux-ml-kem", |b| {
        let mut rng = Rng::new(b"mlkem768-bench-bob");
        b.iter(|| {
            let (k, c) =
                libcrux_mlkem768::encapsulate(ref_kp.public_key(), black_box(rng.bytes32()));
            black_box((k, c));
        });
    });
}

fn bench_roundtrip(c: &mut Criterion) {
    bench_alice(&mut c.benchmark_group("alice"));
    bench_bob(&mut c.benchmark_group("bob"));
}

criterion_group!(benches, bench_roundtrip);
criterion_main!(benches);
