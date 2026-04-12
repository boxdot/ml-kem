use crate::kem::kpke::Ciphertext;

pub(crate) fn shake256(input: &[u8], out: &mut [u8]) {
    use sha3::digest::{ExtendableOutput, Update, XofReader};
    let mut hasher = sha3::Shake256::default();
    hasher.update(input);
    hasher.finalize_xof().read(out);
}

#[allow(non_snake_case)]
pub fn G(input: &[u8]) -> ([u8; 32], [u8; 32]) {
    let hash = sha3_512(input);
    let mut a = [0u8; 32];
    let mut b = [0u8; 32];
    a.copy_from_slice(&hash[0..32]);
    b.copy_from_slice(&hash[32..64]);
    (a, b)
}

fn sha3_512(input: &[u8]) -> [u8; 64] {
    use sha3::Digest;
    let mut hasher = sha3::Sha3_512::new();
    Digest::update(&mut hasher, input);
    hasher.finalize().into()
}

#[allow(non_snake_case)]
pub fn H(input: &[u8]) -> [u8; 32] {
    use sha3::{Digest, Sha3_256};
    let mut hasher = Sha3_256::new();
    hasher.update(input);
    hasher.finalize().into()
}

#[allow(non_snake_case)]
pub fn J(z: &[u8; 32], ct: &Ciphertext) -> [u8; 32] {
    let mut input = [0u8; 32 + 960 + 128];
    input[0..32].copy_from_slice(z);
    input[32..992].copy_from_slice(&ct.c1);
    input[992..1120].copy_from_slice(&ct.c2);
    let mut out = [0u8; 32];
    shake256(&input, &mut out);
    out
}
