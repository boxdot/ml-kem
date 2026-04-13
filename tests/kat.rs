//! Accumulated known-answer tests
//!
//! Source <https://github.com/FiloSottile/mlkem768/blob/main/mlkem768_test.go>

use ml_kem_example::kem::mlkem::{Ciphertext, EncapsulationKey};
use sha3::{
    Shake128,
    digest::{ExtendableOutput, Update, XofReader},
};

fn run_accumulated(iterations: usize, expected_hex: &str) {
    let mut rng = Shake128::default().finalize_xof();
    let mut acc = Shake128::default();

    let mut seed = [0u8; 64];
    let mut m = [0u8; 32];
    let mut ct_rand = [0u8; 1088];

    let mut ek_bytes = [0u8; 1184];

    for _ in 0..iterations {
        rng.read(&mut seed);
        rng.read(&mut m);
        rng.read(&mut ct_rand);

        let d: &[u8; 32] = (&seed[0..32]).try_into().unwrap();
        let z: &[u8; 32] = (&seed[32..64]).try_into().unwrap();

        let (ek, dk) = EncapsulationKey::generate(d, z);
        let (k, ct) = ek.encaps(&m);
        let k_reject = dk.decaps(&Ciphertext::from_bytes(&ct_rand));

        ek.encode(&mut ek_bytes);

        acc.update(&ek_bytes);
        acc.update(&ct.to_bytes());
        acc.update(&k);
        acc.update(&k_reject);
    }

    let mut digest = [0u8; 32];
    acc.finalize_xof().read(&mut digest);

    let got = hex_encode(&digest);
    assert_eq!(got, expected_hex, "accumulated KAT hash mismatch");
}

fn hex_encode(bytes: &[u8]) -> String {
    let mut s = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        s.push_str(&format!("{:02x}", b));
    }
    s
}

#[test]
fn mlkem768_accumulated_100() {
    run_accumulated(
        100,
        "1114b1b6699ed191734fa339376afa7e285c9e6acf6ff0177d346696ce564415",
    );
}

#[test]
fn mlkem768_accumulated_10000() {
    run_accumulated(
        10_000,
        "8a518cc63da366322a8e7a818c7a0d63483cb3528d34a4cf42f35d5ad73f22fc",
    );
}

#[ignore = "long test"]
#[test]
fn mlkem768_accumulated_1000000() {
    run_accumulated(
        1_000_000,
        "424bf8f0e8ae99b78d788a6e2e8e9cdaf9773fc0c08a6f433507cb559edfd0f0",
    );
}
