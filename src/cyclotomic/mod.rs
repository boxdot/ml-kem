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
//!
//!
//! Number theoretic transform (NTT)

use crate::{cyclotomic::base::BASE_ZETAS, integer::Zq};

use base::ZETAS_NTT;

mod base;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Poly {
    /// (2i, 2i+1) -> f_i(X) = a[2i] + a[2i+1] * X (mod X^2 - zeta_i)
    ///
    /// with zeta_i = BASE_ZETAS[i].
    coeffs: [Zq; 256],
}

impl Poly {
    pub fn zero() -> Self {
        Self {
            coeffs: [Zq::zero(); 256],
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
                    self.coeffs[j] = self.coeffs[j] + t;
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

    pub fn pointwise_mul(&mut self, rhs: &Self) -> Self {
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
}

#[cfg(test)]
mod test {
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

    // #[quickcheck]
    // fn test_ntt_intt_roundtrip(mut p: Poly) {
    //     let original = p.clone();
    //
    //     p.ntt();
    //     p.intt();
    //
    //     for i in 0..256 {
    //         assert_eq!(
    //             p.coeffs[i].to_int(),
    //             original.coeffs[i].to_int(),
    //             "Mismatch at coefficient {}",
    //             i
    //         );
    //     }
    // }

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
}
