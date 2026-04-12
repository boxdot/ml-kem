pub(crate) fn shake256(input: &[u8], out: &mut [u8]) {
    use sha3::digest::{ExtendableOutput, Update, XofReader};
    let mut hasher = sha3::Shake256::default();
    hasher.update(input);
    hasher.finalize_xof().read(out);
}

#[allow(non_snake_case)]
pub(crate) fn G(d: &[u8; 32]) -> ([u8; 32], [u8; 32]) {
    let mut input = [0u8; 33];
    input[0..32].copy_from_slice(d);
    input[32] = 3; // k = 3 for ML-KEM-768

    let hash = sha3_512(&input); // 64 bytes

    let mut rho = [0u8; 32];
    let mut sigma = [0u8; 32];
    rho.copy_from_slice(&hash[0..32]);
    sigma.copy_from_slice(&hash[32..64]);
    (rho, sigma)
}

fn sha3_512(input: &[u8]) -> [u8; 64] {
    use sha3::Digest;
    let mut hasher = sha3::Sha3_512::new();
    Digest::update(&mut hasher, input);
    hasher.finalize().into()
}
