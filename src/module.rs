use crate::cyclotomic::Poly;

pub struct Vec3([Poly; 3]);

impl Vec3 {
    pub fn zero() -> Self {
        Self::from_polys([Poly::zero(), Poly::zero(), Poly::zero()])
    }

    pub fn from_polys(polys: [Poly; 3]) -> Self {
        Self(polys)
    }

    // --- Arithmetic ---

    pub fn add_assign(&mut self, rhs: &Self) {
        for i in 0..3 {
            self.0[i].add_assign(&rhs.0[i]);
        }
    }

    pub fn sub_assign(&mut self, rhs: &Self) {
        for i in 0..3 {
            self.0[i].sub_assign(&rhs.0[i]);
        }
    }

    pub fn inner_product_ntt(&self, rhs: &Self) -> Poly {
        let mut res = Poly::zero();
        for i in 0..3 {
            res.mul_add_assign(&self.0[i], &rhs.0[i]);
        }
        res.intt();
        res
    }

    // --- NTT ---

    pub fn ntt(&mut self) {
        for i in 0..3 {
            self.0[i].ntt();
        }
    }

    pub fn intt(&mut self) {
        for i in 0..3 {
            self.0[i].intt();
        }
    }

    // --- Serialization ---

    pub fn encode(&self, out: &mut [u8; 1152]) {
        for i in 0..3 {
            let chunk: &mut [u8; 384] = (&mut out[i * 384..(i + 1) * 384]).try_into().unwrap();
            self.0[i].encode(chunk);
        }
    }

    pub fn decode(bytes: &[u8; 1152], out: &mut Vec3) {
        for i in 0..3 {
            let a = bytes[i * 384..(i + 1) * 384].try_into().unwrap();
            Poly::decode(a, &mut out.0[i]);
        }
    }

    // --- Compression ---

    pub fn compress_10(&self, out: &mut [u8; 960]) {
        for i in 0..3 {
            let chunk: &mut [u8; 320] = (&mut out[i * 320..(i + 1) * 320]).try_into().unwrap();
            self.0[i].compress::<10, 320>(chunk);
        }
    }

    pub fn decompress_10(bytes: &[u8; 960], out: &mut Vec3) {
        for i in 0..3 {
            let a = bytes[i * 320..(i + 1) * 320].try_into().unwrap();
            Poly::decompress::<10, 320>(a, &mut out.0[i]);
        }
    }
}

#[cfg(test)]
mod test {
    use quickcheck_macros::quickcheck;

    use crate::{cyclotomic::test::max_error, integer::Zq};

    use super::*;

    #[test]
    fn test_encode_decode_roundtrip() {
        let mut v = Vec3::zero();
        v.0[0].coeffs[0] = Zq::from_int(0xABC); // 2748
        v.0[0].coeffs[1] = Zq::from_int(0x123); // 291
        v.0[1].coeffs[0] = Zq::from_int(3328); // max
        v.0[2].coeffs[255] = Zq::from_int(1); // last coeff of last poly

        let mut out = [0u8; 1152];
        v.encode(&mut out);

        let mut decoded = Vec3::zero();
        Vec3::decode(&out, &mut decoded);

        for p in 0..3 {
            for i in 0..256 {
                assert_eq!(
                    v.0[p].coeffs[i].to_int(),
                    decoded.0[p].coeffs[i].to_int(),
                    "Mismatch at poly {} coeff {}",
                    p,
                    i
                );
            }
        }
    }

    #[test]
    fn test_encode_boundary_between_polys() {
        // Ensure the 384-byte boundary between polys is handled correctly.
        // Last coeff of poly[0] and first coeff of poly[1] should not bleed into each other.
        let mut v = Vec3::zero();
        v.0[0].coeffs[255] = Zq::from_int(3328);
        v.0[1].coeffs[0] = Zq::from_int(1);

        let mut out = [0u8; 1152];
        v.encode(&mut out);

        let mut decoded = Vec3::zero();
        Vec3::decode(&out, &mut decoded);

        assert_eq!(decoded.0[0].coeffs[255].to_int(), 3328);
        assert_eq!(decoded.0[1].coeffs[0].to_int(), 1);
    }

    #[quickcheck]
    fn test_encode_decode_roundtrip_qc(polys: [Poly; 3]) -> bool {
        let v = Vec3::from_polys(polys);

        let mut out = [0u8; 1152];
        v.encode(&mut out);

        let mut decoded = Vec3::zero();
        Vec3::decode(&out, &mut decoded);

        (0..3)
            .all(|p| (0..256).all(|i| v.0[p].coeffs[i].to_int() == decoded.0[p].coeffs[i].to_int()))
    }

    fn check_vec3_compress_error(v: &Vec3) {
        let mut compressed = [0u8; 960];
        v.compress_10(&mut compressed);

        let mut decompressed = Vec3::zero();
        Vec3::decompress_10(&compressed, &mut decompressed);

        let bound = max_error(10);
        for p in 0..3 {
            for i in 0..256 {
                let orig = v.0[p].coeffs[i].to_int();
                let recv = decompressed.0[p].coeffs[i].to_int();
                let err = (orig - recv).abs().min(3329 - (orig - recv).abs());
                assert!(
                    err <= bound,
                    "poly {p} coeff {i}: orig={orig} recv={recv} err={err} bound={bound}"
                );
            }
        }
    }

    #[test]
    fn test_vec3_compress_zero() {
        check_vec3_compress_error(&Vec3::zero());
    }

    #[test]
    fn test_vec3_compress_boundary() {
        // Check that the 320-byte boundary between polys is handled correctly
        let mut v = Vec3::zero();
        v.0[0].coeffs[255] = Zq::from_int(1664); // last coeff of poly[0]
        v.0[1].coeffs[0] = Zq::from_int(1664); // first coeff of poly[1]
        v.0[2].coeffs[255] = Zq::from_int(3328); // last coeff of poly[2]
        check_vec3_compress_error(&v);
    }

    #[quickcheck]
    fn test_vec3_compress_qc(polys: [Poly; 3]) -> bool {
        let v = Vec3::from_polys(polys);
        let mut compressed = [0u8; 960];
        v.compress_10(&mut compressed);

        let mut decompressed = Vec3::zero();
        Vec3::decompress_10(&compressed, &mut decompressed);

        let bound = max_error(10);
        (0..3).all(|p| {
            (0..256).all(|i| {
                let orig = v.0[p].coeffs[i].to_int();
                let recv = decompressed.0[p].coeffs[i].to_int();
                let err = (orig - recv).abs().min(3329 - (orig - recv).abs());
                err <= bound
            })
        })
    }
}
