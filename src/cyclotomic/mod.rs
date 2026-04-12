//! Implementation of the arithmetic in the cyclotomic ring Rq = Zq[X]/(X^256 + 1),
//! where Zq is the integers modulo Q = 3329.
//!
//! ## Mathematical discussion
//!
//! X^256 + 1 is the irreducible polynomial of degree 256 over the field Zq and is the 512-th
//! cyclotomic polynomial PHI(512). It does not split into linear factors over Zq (because the field
//! does not contain 512th primitive root of unity). However, the field contains a 256th primitive
//! root of unity.
//!
//! The polynomial PHI(512) splits ove Zq in the following irreducible linear factors:
//!
//! (X^2 - zeta_i),
//!
//! where zeta_i are the primitive 256th roots of unity.
//!
//! Since 17^128 = -1 mod Q => 17 is a (primitive) 256th root of unity mod Q. We take it as the
//! generator of the group of 256th roots of unity mod Q.

use crate::{
    cyclotomic::base::BASE_ZETAS,
    hash::shake256,
    integer::{Q, Zq},
};

use base::ZETAS_NTT;

mod base;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Poly {
    /// (2i, 2i+1) -> f_i(X) = a[2i] + a[2i+1] * X (mod X^2 - zeta_i)
    ///
    /// with zeta_i = BASE_ZETAS[i].
    pub(crate) coeffs: [Zq; 256],
}

impl Poly {
    pub fn zero() -> Self {
        Self {
            coeffs: [Zq::zero(); 256],
        }
    }

    /// Sample via CBD(η = 2)
    ///
    /// # Mathematical discussion
    ///
    /// CBD(η) is the centered binomial distribution with parameter eta, that is, is samples each
    /// coefficient as
    ///
    /// x = Σ(i=0..η) a_i - Σ(i=0..η) b_i,
    ///
    /// where a_i, b_i are independent uniform bits.
    ///
    /// For η = 2, this is equivalent to sampling each coefficient as
    ///
    /// x = (a_0 + a_1) - (b_0 + b_1)   ∈ {-2, -1, 0, 1, 2},
    ///
    /// with probabilities
    ///
    /// P(0)  = 6/16 = 0.375
    /// P(±1) = 4/16 = 0.25  each
    /// P(±2) = 1/16 = 0.0625 each.
    pub fn sample_cbd(seed: &[u8; 32], nonce: u8) -> Self {
        // PRF: SHAKE-256(seed || nonce) -> 128 bytes (64η bytes for η=2)
        let mut prf_input = [0u8; 33];
        prf_input[0..32].copy_from_slice(seed);
        prf_input[32] = nonce;

        let mut prf_output = [0u8; 128];
        shake256(&prf_input, &mut prf_output);

        let mut p = Poly::zero();
        // Each byte gives 2 coefficients (4 bits each)
        for (i, byte) in prf_output.into_iter().enumerate() {
            let a0 = byte & 1;
            let a1 = byte >> 1 & 1;
            let b0 = byte >> 2 & 1;
            let b1 = byte >> 3 & 1;
            p.coeffs[2 * i] = Zq::from_int((a0 + a1) as i16 - (b0 + b1) as i16);

            let a0 = byte >> 4 & 1;
            let a1 = byte >> 5 & 1;
            let b0 = byte >> 6 & 1;
            let b1 = byte >> 7 & 1;
            p.coeffs[2 * i + 1] = Zq::from_int((a0 + a1) as i16 - (b0 + b1) as i16);
        }

        p
    }

    pub fn add_assign(&mut self, rhs: &Self) {
        for i in 0..256 {
            self.coeffs[i] += rhs.coeffs[i];
        }
    }

    pub fn sub_assign(&mut self, rhs: &Self) {
        for i in 0..256 {
            self.coeffs[i] = self.coeffs[i] - rhs.coeffs[i];
        }
    }

    pub fn ntt(&mut self) {
        let mut k = 1; // Index into ZETAS_NTT
        let mut len = 128; // Initial stride width
        while len > 1 {
            let mut start = 0;
            while start < 256 {
                let zeta = ZETAS_NTT[k];
                k += 1;

                for j in start..start + len {
                    // t = zeta * a[j + len];
                    let t = self.coeffs[j + len] * zeta;
                    // a[j + len] = a[j] - t;
                    self.coeffs[j + len] = self.coeffs[j] - t;
                    // a[j] = a[j] + t;
                    self.coeffs[j] += t;
                }

                start += len * 2;
            }
            len /= 2;
        }
    }

    /// Inverse number theoretic transform (INTT)
    pub fn intt(&mut self) {
        let mut k = 127;
        let mut len = 2;
        while len <= 128 {
            let mut start = 0;
            while start < 256 {
                let zeta = ZETAS_NTT[k];
                k -= 1;
                for j in start..start + len {
                    let u = self.coeffs[j];
                    let v = self.coeffs[j + len];
                    self.coeffs[j] = u + v;
                    self.coeffs[j + len] = (v - u) * zeta;
                }
                start += len * 2;
            }
            len *= 2;
        }

        // Final scaling by 256^-1 mod Q = 3303
        let f = Zq::from_int(3303);
        for i in 0..256 {
            self.coeffs[i] = self.coeffs[i] * f;
        }
    }

    pub fn pointwise_mul(&self, rhs: &Self) -> Self {
        let mut res = Poly::zero();
        for (i, &zeta) in BASE_ZETAS.iter().enumerate() {
            let a0 = self.coeffs[2 * i];
            let a1 = self.coeffs[2 * i + 1];
            let b0 = rhs.coeffs[2 * i];
            let b1 = rhs.coeffs[2 * i + 1];
            res.coeffs[2 * i] = a0 * b0 + zeta * a1 * b1;
            res.coeffs[2 * i + 1] = a0 * b1 + a1 * b0;
        }
        res
    }

    pub fn mul_add_assign(&mut self, a: &Self, b: &Self) {
        for (i, &zeta) in BASE_ZETAS.iter().enumerate() {
            let a0 = a.coeffs[2 * i];
            let a1 = a.coeffs[2 * i + 1];
            let b0 = b.coeffs[2 * i];
            let b1 = b.coeffs[2 * i + 1];
            self.coeffs[2 * i] += a0 * b0 + zeta * a1 * b1;
            self.coeffs[2 * i + 1] += a0 * b1 + a1 * b0;
        }
    }

    // --- Serialization ---

    pub fn encode(&self, out: &mut [u8; 384]) {
        for i in 0..128 {
            let a = self.coeffs[2 * i].to_int();
            let b = self.coeffs[2 * i + 1].to_int();
            debug_assert!(a < 3329 && b < 3329);
            let a = a as u16;
            let b = b as u16;
            out[3 * i] = (a & 0xFF) as u8; // low 8 bits of a
            out[3 * i + 1] = ((a >> 8) | (b << 4)) as u8; // high 4 bits of a, low 4 bits of b
            out[3 * i + 2] = (b >> 4) as u8; // high 4 bits of b
        }
    }

    pub fn decode(bytes: &[u8; 384], out: &mut Poly) {
        for i in 0..128 {
            let a = bytes[3 * i] as u16;
            let b = bytes[3 * i + 1] as u16;
            let c = bytes[3 * i + 2] as u16;
            let a = a | ((b & 0xF) << 8);
            let b = (b >> 4) | (c << 4);
            out.coeffs[2 * i] = Zq::from_int(a as i16);
            out.coeffs[2 * i + 1] = Zq::from_int(b as i16);
        }
    }

    // --- Compression ---

    pub fn compress<const D: u8, const OUT: usize>(&self, out: &mut [u8; OUT]) {
        debug_assert_eq!(OUT, 256 * D as usize / 8);

        let q = Q as u32;

        let mut buf: u32 = 0; // bits accumulator
        let mut bits: u8 = 0; // number of valid bits in buf
        let mut out_i = 0; // current byte index

        for i in 0..256 {
            // Scale from Zq to Z_{2^D}
            // round(x * 2^D / q) = (x * 2^D + q/2) / q
            let x = self.coeffs[i].to_int() as u32;
            let compressed = (((x << D) + q / 2) / q) & ((1 << D) - 1);

            // Pack D bits into buf
            buf |= compressed << bits;
            bits += D;

            // Flush full bytes
            while bits >= 8 {
                out[out_i] = buf as u8;
                buf >>= 8;
                bits -= 8;
                out_i += 1;
            }
        }

        debug_assert_eq!(bits, 0);
    }

    pub fn decompress<const D: u8, const OUT: usize>(bytes: &[u8; OUT], out: &mut Poly) {
        debug_assert_eq!(OUT, 256 * D as usize / 8);

        let q = Q as u32;

        let mut buf: u32 = 0; // bits accumulator
        let mut bits: u8 = 0; // number of valid bits in buf
        let mut in_i = 0; // current byte index

        let mask: u32 = (1 << D) - 1;

        for i in 0..256 {
            // Extact D bits from buf
            while bits < D {
                buf |= (bytes[in_i] as u32) << bits;
                bits += 8;
                in_i += 1;
            }
            let y = buf & mask;
            buf >>= D;
            bits -= D;

            // Scale from Z_{2^D} to Zq
            // round(y * q / 2^D) = (y * q - 2^(D-1)) / 2^D
            let x = (y * q + (1 << (D - 1))) >> D;
            out.coeffs[i] = Zq::from_int(x as i16);
        }

        debug_assert_eq!(bits, 0);
    }
}

#[cfg(test)]
pub(crate) mod test {
    use quickcheck::{Arbitrary, Gen};
    use quickcheck_macros::quickcheck;

    use super::*;

    impl Arbitrary for Poly {
        fn arbitrary(g: &mut Gen) -> Self {
            let mut p = Poly::zero();
            for i in 0..256 {
                p.coeffs[i] = Zq::from_int(Arbitrary::arbitrary(g));
            }
            p
        }
    }

    #[test]
    fn test_ntt_constant() {
        let mut p = Poly::zero();
        p.coeffs[0] = Zq::from_int(1);

        p.ntt();

        // Constant polynomial 1 should evaluate to 1 in every slot:
        //
        // In degree-2 layout (a0 + a1 * X), a0 = 1, a1 = 0 for all 128 slots.
        for i in 0..128 {
            assert_eq!(p.coeffs[2 * i].to_int(), 1);
            assert_eq!(p.coeffs[2 * i + 1].to_int(), 0);
        }
    }

    #[quickcheck]
    fn test_ntt(a: Poly, b: Poly) -> bool {
        let mut sum = Poly::zero();
        for i in 0..256 {
            sum.coeffs[i] = a.coeffs[i] + b.coeffs[i];
        }

        let mut a = a.clone();
        let mut b = b.clone();

        a.ntt();
        b.ntt();
        sum.ntt();

        (0..256).all(|i| a.coeffs[i] + b.coeffs[i] == sum.coeffs[i])
    }

    #[test]
    fn test_ntt_intt_roundtrip() {
        let mut p = Poly::zero();
        // Fill with a "stress test" pattern
        for i in 0..256 {
            p.coeffs[i] = Zq::from_int(i as i16);
        }

        let original = p.clone();

        p.ntt();
        p.intt();

        for i in 0..256 {
            assert_eq!(
                p.coeffs[i].to_int(),
                original.coeffs[i].to_int(),
                "Mismatch at coefficient {}",
                i
            );
        }
    }

    #[test]
    fn test_full_multiplication_pipeline() {
        // Multiply (1 + 2X) * (3 + 4X) = 3 + 10X + 8X^2

        // 1 + 2X
        let mut a = Poly::zero();
        a.coeffs[0] = Zq::from_int(1);
        a.coeffs[1] = Zq::from_int(2);

        // 3 + 4X
        let mut b = Poly::zero();
        b.coeffs[0] = Zq::from_int(3);
        b.coeffs[1] = Zq::from_int(4);

        // 1. Move to NTT domain
        a.ntt();
        b.ntt();

        // 2. Pointwise multiply
        let mut c = a.pointwise_mul(&b);

        // 3. Move back to coefficient domain
        c.intt();

        assert_eq!(c.coeffs[0].to_int(), 3);
        assert_eq!(c.coeffs[1].to_int(), 10);
        assert_eq!(c.coeffs[2].to_int(), 8);
        for i in 3..256 {
            assert_eq!(c.coeffs[i].to_int(), 0, "Index {} should be 0", i);
        }
    }

    #[test]
    fn test_encode_decode_roundtrip() {
        let mut p = Poly::zero();
        // Use known values that exercise all 12 bits
        p.coeffs[0] = Zq::from_int(0xABC); // 2748
        p.coeffs[1] = Zq::from_int(0x123); // 291
        p.coeffs[2] = Zq::from_int(0);
        p.coeffs[3] = Zq::from_int(3328); // max value

        let mut out = [0u8; 384];
        p.encode(&mut out);

        let mut decoded = Poly::zero();
        Poly::decode(&out, &mut decoded);

        for i in 0..256 {
            assert_eq!(
                p.coeffs[i].to_int(),
                decoded.coeffs[i].to_int(),
                "Mismatch at coefficient {}",
                i
            );
        }
    }

    #[quickcheck]
    fn test_encode_decode_roundtrip_qc(p: Poly) -> bool {
        let mut out = [0u8; 384];
        p.encode(&mut out);

        let mut decoded = Poly::zero();
        Poly::decode(&out, &mut decoded);

        (0..256).all(|i| p.coeffs[i].to_int() == decoded.coeffs[i].to_int())
    }

    // Maximum allowable rounding error per coefficient for each D
    // bound = round(q / 2^(D+1))
    // D=10: round(3329 / 2048) = 2
    // D=4:  round(3329 / 32)   = 104
    pub(crate) fn max_error(d: u8) -> i16 {
        let q = 3329i32;
        ((q + (1 << d)) / (1 << (d + 1))) as i16
    }

    fn check_compress_error<const D: u8, const OUT: usize>(p: &Poly) {
        let mut compressed = [0u8; OUT];
        p.compress::<D, OUT>(&mut compressed);

        let mut decompressed = Poly::zero();
        Poly::decompress::<D, OUT>(&compressed, &mut decompressed);

        let bound = max_error(D);
        for i in 0..256 {
            let orig = p.coeffs[i].to_int();
            let recv = decompressed.coeffs[i].to_int();
            // error is computed mod q on the circle
            let err = (orig - recv).abs().min(3329 - (orig - recv).abs());
            assert!(
                err <= bound,
                "coeff {i}: orig={orig} recv={recv} err={err} bound={bound}"
            );
        }
    }

    #[test]
    fn test_compress_10_zero() {
        check_compress_error::<10, 320>(&Poly::zero());
    }

    #[test]
    fn test_compress_4_zero() {
        check_compress_error::<4, 128>(&Poly::zero());
    }

    #[test]
    fn test_compress_10_known() {
        let mut p = Poly::zero();
        p.coeffs[0] = Zq::from_int(0);
        p.coeffs[1] = Zq::from_int(1664); // q/2, maps to 2^(D-1)
        p.coeffs[2] = Zq::from_int(3328); // q-1, maps back to 0
        check_compress_error::<10, 320>(&p);
    }

    #[test]
    fn test_compress_4_known() {
        let mut p = Poly::zero();
        p.coeffs[0] = Zq::from_int(0);
        p.coeffs[1] = Zq::from_int(1664);
        p.coeffs[2] = Zq::from_int(3328);
        check_compress_error::<4, 128>(&p);
    }

    #[quickcheck]
    fn test_compress_10_qc(p: Poly) -> bool {
        let mut compressed = [0u8; 320];
        p.compress::<10, 320>(&mut compressed);
        let mut decompressed = Poly::zero();
        Poly::decompress::<10, 320>(&compressed, &mut decompressed);
        let bound = max_error(10);
        (0..256).all(|i| {
            let orig = p.coeffs[i].to_int();
            let recv = decompressed.coeffs[i].to_int();
            let err = (orig - recv).abs().min(3329 - (orig - recv).abs());
            err <= bound
        })
    }

    #[quickcheck]
    fn test_compress_4_qc(p: Poly) -> bool {
        let mut compressed = [0u8; 128];
        p.compress::<4, 128>(&mut compressed);
        let mut decompressed = Poly::zero();
        Poly::decompress::<4, 128>(&compressed, &mut decompressed);
        let bound = max_error(4);
        (0..256).all(|i| {
            let orig = p.coeffs[i].to_int();
            let recv = decompressed.coeffs[i].to_int();
            let err = (orig - recv).abs().min(3329 - (orig - recv).abs());
            err <= bound
        })
    }

    #[test]
    fn test_sample_cbd_range() {
        // All coefficients must be in [-2, 2]
        let seed = [0u8; 32];
        let p = Poly::sample_cbd(&seed, 0);
        for i in 0..256 {
            // -2 → 3327
            // -1 → 3328
            //  0 → 0
            //  1 → 1
            //  2 → 2
            let c = p.coeffs[i].to_int();
            assert!(
                c == 0 || c == 1 || c == 2 || c == 3327 || c == 3328,
                "coeff {i} = {c} out of range"
            );
        }
    }

    #[test]
    fn test_sample_cbd_different_nonces() {
        // Different nonces must produce different polynomials
        let seed = [0u8; 32];
        let p0 = Poly::sample_cbd(&seed, 0);
        let p1 = Poly::sample_cbd(&seed, 1);
        assert_ne!(p0, p1);
    }

    #[test]
    fn test_sample_cbd_different_seeds() {
        // Different seeds must produce different polynomials
        let seed0 = [0u8; 32];
        let mut seed1 = [0u8; 32];
        seed1[0] = 1;
        let p0 = Poly::sample_cbd(&seed0, 0);
        let p1 = Poly::sample_cbd(&seed1, 0);
        assert_ne!(p0, p1);
    }

    #[test]
    fn test_sample_cbd_distribution() {
        // Sample many polynomials and check rough frequency of each value.
        // Expected probabilities for η=2:
        // P(0)  = 6/16 = 0.375
        // P(±1) = 4/16 = 0.25
        // P(±2) = 1/16 = 0.0625
        let mut counts = [0i64; 5]; // indices 0..4 for values -2..2
        let n_polys = 100;

        for i in 0..n_polys {
            let seed = [i as u8; 32];
            let p = Poly::sample_cbd(&seed, 0);
            for j in 0..256 {
                let c = p.coeffs[j].to_int();
                let normalized = if c > 3329 / 2 {
                    c as i32 - 3329
                } else {
                    c as i32
                };
                counts[(normalized + 2) as usize] += 1;
            }
        }

        let total = (n_polys * 256) as f64;
        let freq: Vec<f64> = counts.iter().map(|&c| c as f64 / total).collect();

        // Check within 5% of expected
        assert!((freq[2] - 0.375).abs() < 0.05, "P(0)  = {:.3}", freq[2]);
        assert!((freq[1] - 0.25).abs() < 0.05, "P(-1) = {:.3}", freq[1]);
        assert!((freq[3] - 0.25).abs() < 0.05, "P(+1) = {:.3}", freq[3]);
        assert!((freq[0] - 0.0625).abs() < 0.05, "P(-2) = {:.3}", freq[0]);
        assert!((freq[4] - 0.0625).abs() < 0.05, "P(+2) = {:.3}", freq[4]);
    }

    #[test]
    fn test_sample_cbd_deterministic() {
        // Same seed and nonce must always produce the same polynomial
        let seed = [42u8; 32];
        let p0 = Poly::sample_cbd(&seed, 0);
        let p1 = Poly::sample_cbd(&seed, 0);
        assert_eq!(p0, p1);
    }
}
