//! Kyber Public Key Encryption

use subtle::{Choice, ConstantTimeEq};

use crate::{
    cyclotomic::ring::Poly,
    hash::G,
    module::{matrix::Mat3x3, vector::Vec3},
};

#[derive(Clone)]
pub struct EncryptionKey {
    t: Vec3,        // Public vector, in NTT domain
    seed: [u8; 32], // Seed to regenerate A on the fly
}

pub struct DecryptionKey {
    s: Vec3, // Secret vector, in NTT domain
}

pub struct Ciphertext {
    pub(crate) c1: [u8; 960], // compressed u with d = 10, 960 bytes
    pub(crate) c2: [u8; 128], // compressed v with d = 4, 128 bytes
}

impl EncryptionKey {
    /// Takes a 32-byte seed, return (ek, dk)
    pub fn generate(seed: &[u8; 32]) -> (EncryptionKey, DecryptionKey) {
        // 1. Expand seed into (rho, sigma) via G = SHA3-512
        let mut g_input = [0u8; 33];
        g_input[0..32].copy_from_slice(seed);
        g_input[32] = 3; // ML-KEM-768 parameter k=3
        let (rho, sigma) = G(&g_input);

        // 2. Generate matrix A from rho, stays in NTT domain
        let a = Mat3x3::sample_ntt(&rho);

        // 3. Sample secret and error from sigma
        let mut s = Vec3::sample_cbd(&sigma, 0); // nonces 0, 1, 2
        let mut e = Vec3::sample_cbd(&sigma, 3); // nonces 3, 4, 5

        // 4. Move secret and error to NTT domain
        s.ntt();
        e.ntt();

        // 5. Compute t = As + e (all in NTT domain)
        let mut t = a.mul_vec_ntt(&s);
        t.add_assign(&e);

        (EncryptionKey { t, seed: rho }, DecryptionKey { s })
    }

    /// Encode to bytes
    pub fn encode(&self, out: &mut [u8; 1184]) {
        self.t.encode((&mut out[0..1152]).try_into().unwrap());
        out[1152..1184].copy_from_slice(&self.seed);
    }

    /// Decode from bytes
    pub fn decode(bytes: &[u8; 1184]) -> Self {
        let mut out = Self {
            t: Vec3::zero(),
            seed: [0u8; 32],
        };
        Vec3::decode(&bytes[0..1152].try_into().unwrap(), &mut out.t);
        out.seed.copy_from_slice(&bytes[1152..1184]);
        out
    }

    /// m is a 32-byte message, r is 32-byte randomness
    pub fn encrypt(&self, m: &[u8; 32], r: &[u8; 32]) -> Ciphertext {
        // 1. Regenerate matrix A from the seed
        let a = Mat3x3::sample_ntt(&self.seed);

        // 2. Sample randomness and errors from r
        let mut r_vec = Vec3::sample_cbd(r, 0); // nonces 0, 1, 2
        let e1 = Vec3::sample_cbd(r, 3); // nonces 3, 4, 5
        let e2 = Poly::sample_cbd(r, 6); // nonces 6

        // 3. NNT transform of r
        r_vec.ntt();

        // 4. u = A^T *r + e1 (u is in coefficient domain)
        let mut u = a.mul_trans_ntt(&r_vec);
        u.intt();
        u.add_assign(&e1);

        // 5. msg polynomial μ: decompress_1(m)
        let mut mu = Poly::zero();
        Poly::decompress::<1, 32>(m, &mut mu);

        // 6. v = t^T * r + e2 + μ (v is in coefficient domain)
        let mut v = self.t.inner_product_ntt(&r_vec);
        v.intt();
        v.add_assign(&e2);
        v.add_assign(&mu);

        // 7. Compress
        let mut c = Ciphertext {
            c1: [0; 960],
            c2: [0; 128],
        };
        u.compress_10(&mut c.c1);
        v.compress_4(&mut c.c2);
        c
    }
}

impl Ciphertext {
    /// Encode to 1088 bytes: c1 (960) || c2 (128)
    pub fn to_bytes(&self) -> [u8; 1088] {
        let mut out = [0u8; 1088];
        out[0..960].copy_from_slice(&self.c1);
        out[960..1088].copy_from_slice(&self.c2);
        out
    }

    /// Decode from 1088 bytes: c1 (960) || c2 (128)
    pub fn from_bytes(bytes: &[u8; 1088]) -> Self {
        let mut c = Ciphertext {
            c1: [0; 960],
            c2: [0; 128],
        };
        c.c1.copy_from_slice(&bytes[0..960]);
        c.c2.copy_from_slice(&bytes[960..1088]);
        c
    }
}

impl DecryptionKey {
    /// Encode the KPKE decryption key (ByteEncode_12 of s in NTT domain).
    pub fn encode(&self, out: &mut [u8; 1152]) {
        self.s.encode(out);
    }

    pub fn decrypt(&self, c: &Ciphertext) -> [u8; 32] {
        // 1. Decompress ciphertext u and v
        let mut u = Vec3::zero();
        Vec3::decompress_10(&c.c1, &mut u);

        let mut v = Poly::zero();
        Poly::decompress_4(&c.c2, &mut v);

        // 2. s^T * u
        u.ntt();
        let mut su = self.s.inner_product_ntt(&u);
        su.intt();

        // 3. w = v - s^T * u
        v.sub_assign(&su);

        // 4. Recover message m from w
        let mut m = [0u8; 32];
        v.compress::<1, 32>(&mut m);
        m
    }
}

impl ConstantTimeEq for Ciphertext {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.c1.ct_eq(&other.c1) & self.c2.ct_eq(&other.c2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kpke_encrypt_decrypt_roundtrip() {
        let d = [0u8; 32];
        let (ek, dk) = EncryptionKey::generate(&d);

        let m = [0u8; 32];
        let r = [0u8; 32];
        let ct = ek.encrypt(&m, &r);
        let recovered = dk.decrypt(&ct);

        assert_eq!(m, recovered);
    }

    #[test]
    fn test_kpke_encrypt_decrypt_nonzero_message() {
        let d = [1u8; 32];
        let (ek, dk) = EncryptionKey::generate(&d);

        let mut m = [0u8; 32];
        m[0] = 0xFF;
        m[15] = 0xAB;
        m[31] = 0x01;

        let r = [2u8; 32];
        let ct = ek.encrypt(&m, &r);
        let recovered = dk.decrypt(&ct);

        assert_eq!(m, recovered);
    }

    #[test]
    fn test_kpke_different_messages_different_ciphertexts() {
        let d = [0u8; 32];
        let (ek, _dk) = EncryptionKey::generate(&d);
        let r = [0u8; 32];

        let m0 = [0u8; 32];
        let mut m1 = [0u8; 32];
        m1[0] = 1;

        let ct0 = ek.encrypt(&m0, &r);
        let ct1 = ek.encrypt(&m1, &r);

        // c1 = A^T*r + e1 does not depend on the message, so only check c2
        assert_eq!(ct0.c1, ct1.c1);
        assert_ne!(ct0.c2, ct1.c2);
    }

    #[test]
    fn test_kpke_different_randomness_different_ciphertexts() {
        let d = [0u8; 32];
        let (ek, _dk) = EncryptionKey::generate(&d);
        let m = [0xABu8; 32];

        let ct0 = ek.encrypt(&m, &[0u8; 32]);
        let ct1 = ek.encrypt(&m, &[1u8; 32]);

        assert_ne!(ct0.c1, ct1.c1);
    }

    #[test]
    fn test_kpke_deterministic() {
        let d = [0u8; 32];
        let (ek, _dk) = EncryptionKey::generate(&d);
        let m = [0xABu8; 32];
        let r = [0u8; 32];

        let ct0 = ek.encrypt(&m, &r);
        let ct1 = ek.encrypt(&m, &r);

        assert_eq!(ct0.c1, ct1.c1);
        assert_eq!(ct0.c2, ct1.c2);
    }
}
