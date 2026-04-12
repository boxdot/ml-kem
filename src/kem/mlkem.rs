use subtle::ConstantTimeEq;

use crate::{
    hash::{G, H, J},
    kem::kpke::{DecryptionKey, EncryptionKey},
};

pub struct EncapsulationKey {
    ek: EncryptionKey,
    ek_hash: [u8; 32], // H(ek)
}

pub struct DecapsulationKey {
    dk: DecryptionKey,
    ek: EncryptionKey,
    ek_hash: [u8; 32], // H(ek)
    z: [u8; 32],       // rejection randomness
}

pub use super::kpke::Ciphertext;

pub type SharedSecret = [u8; 32];

impl EncapsulationKey {
    /// d and z are 32 bytes of fresh randomness
    pub fn generate(d: &[u8; 32], z: &[u8; 32]) -> (EncapsulationKey, DecapsulationKey) {
        let (ek, dk) = EncryptionKey::generate(d);

        let mut ek_bytes = [0u8; 1184];
        ek.encode(&mut ek_bytes);
        let ek_hash = H(&ek_bytes);

        let encaps_key = EncapsulationKey {
            ek: ek.clone(),
            ek_hash,
        };
        let decaps_key = DecapsulationKey {
            dk,
            ek,
            ek_hash,
            z: *z,
        };

        (encaps_key, decaps_key)
    }

    /// m is 32 bytes of fresh randomness
    pub fn encaps(&self, m: &[u8; 32]) -> (SharedSecret, Ciphertext) {
        // 1. Derive (K, r) from m
        // G input: m || H(ek)
        let mut g_input = [0u8; 64];
        g_input[0..32].copy_from_slice(m);
        g_input[32..64].copy_from_slice(&self.ek_hash);
        let (k, r) = G(&g_input);

        // 2. Encrypt m using derived randomness
        let ct = self.ek.encrypt(m, &r);

        (k, ct)
    }
}

impl DecapsulationKey {
    pub fn decaps(&self, ct: &Ciphertext) -> SharedSecret {
        // 1. Decrypt to recover m'
        let m_prime = self.dk.decrypt(ct);

        // 2. Derive (K', r') from m'
        let mut g_input = [0u8; 64];
        g_input[0..32].copy_from_slice(&m_prime);
        g_input[32..64].copy_from_slice(&self.ek_hash);
        let (k_prime, r_prime) = G(&g_input);

        // 3. Re-encrypt and check that the ciphertext matches
        let ct_prime = self.ek.encrypt(&m_prime, &r_prime);

        // 4. Constant-time comparison and rejection
        if ct_prime.ct_eq(ct).into() {
            k_prime
        } else {
            J(&self.z, ct)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mlkem_encaps_decaps_roundtrip() {
        let d = [0u8; 32];
        let z = [1u8; 32];
        let (ek, dk) = EncapsulationKey::generate(&d, &z);

        let m = [2u8; 32];
        let (ss_encaps, ct) = ek.encaps(&m);
        let ss_decaps = dk.decaps(&ct);

        assert_eq!(ss_encaps, ss_decaps);
    }

    #[test]
    fn test_mlkem_tampered_ciphertext_rejected() {
        let d = [0u8; 32];
        let z = [1u8; 32];
        let (ek, dk) = EncapsulationKey::generate(&d, &z);

        let m = [2u8; 32];
        let (ss_encaps, mut ct) = ek.encaps(&m);

        // tamper with ciphertext
        ct.c1[0] ^= 1;

        let ss_decaps = dk.decaps(&ct);

        // shared secrets must differ — implicit rejection fired
        assert_ne!(ss_encaps, ss_decaps);
    }

    #[test]
    fn test_mlkem_rejection_is_deterministic() {
        // same tampered ciphertext must always produce same rejected secret
        let d = [0u8; 32];
        let z = [1u8; 32];
        let (ek, dk) = EncapsulationKey::generate(&d, &z);

        let m = [2u8; 32];
        let (_, mut ct) = ek.encaps(&m);
        ct.c1[0] ^= 1;

        let ss0 = dk.decaps(&ct);
        let ss1 = dk.decaps(&ct);

        assert_eq!(ss0, ss1);
    }

    #[test]
    fn test_mlkem_rejection_depends_on_z() {
        // two keypairs with different z must produce different rejected secrets
        let d = [0u8; 32];
        let (ek, dk0) = EncapsulationKey::generate(&d, &[0u8; 32]);
        let (_, dk1) = EncapsulationKey::generate(&d, &[1u8; 32]);

        let m = [2u8; 32];
        let (_, mut ct) = ek.encaps(&m);
        ct.c1[0] ^= 1;

        let ss0 = dk0.decaps(&ct);
        let ss1 = dk1.decaps(&ct);

        assert_ne!(ss0, ss1);
    }

    #[test]
    fn test_mlkem_different_messages_different_secrets() {
        let d = [0u8; 32];
        let z = [1u8; 32];
        let (ek, dk) = EncapsulationKey::generate(&d, &z);

        let (ss0, ct0) = ek.encaps(&[0u8; 32]);
        let (ss1, ct1) = ek.encaps(&[1u8; 32]);

        assert_ne!(ss0, ss1);
        assert_ne!(ct0.c1, ct1.c1);

        // both must still decapsulate correctly
        assert_eq!(ss0, dk.decaps(&ct0));
        assert_eq!(ss1, dk.decaps(&ct1));
    }

    #[test]
    fn test_mlkem_deterministic() {
        // same inputs must always produce same outputs
        let d = [0u8; 32];
        let z = [1u8; 32];
        let (ek, dk) = EncapsulationKey::generate(&d, &z);
        let m = [2u8; 32];

        let (ss0, ct0) = ek.encaps(&m);
        let (ss1, ct1) = ek.encaps(&m);

        assert_eq!(ss0, ss1);
        assert_eq!(ct0.c1, ct1.c1);
        assert_eq!(ct0.c2, ct1.c2);
        assert_eq!(dk.decaps(&ct0), dk.decaps(&ct1));
    }
}
