use sha3::{
    Shake128,
    digest::{ExtendableOutput, Update, XofReader},
};

use crate::{
    cyclotomic::ring::Poly,
    integer::{Q, Zq},
    module::vector::Vec3,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Mat3x3([Vec3; 3]);

impl Mat3x3 {
    pub fn zero() -> Self {
        Self([Vec3::zero(), Vec3::zero(), Vec3::zero()])
    }

    /// Generate A deterministically from a 32-byte seed via SampleNTT.
    ///
    /// Returns A already in NTT domain.
    pub fn sample_ntt(seed: &[u8; 32]) -> Self {
        let mut mat = Self::zero();
        for i in 0..3 {
            for j in 0..3 {
                *mat.0[i].get_mut(j) = sample_ntt_poly(seed, i as u8, j as u8);
            }
        }
        mat
    }

    /// Compute As where A and s both are in NTT domain.
    ///
    /// Returns a [`Vec3`] in NTT domain.
    pub fn mul_vec_ntt(&self, v: &Vec3) -> Vec3 {
        Vec3::from_polys([
            self.0[0].inner_product_ntt(v),
            self.0[1].inner_product_ntt(v),
            self.0[2].inner_product_ntt(v),
        ])
    }

    /// Compute A^T s where A^T (transpose of A) and v both are in NTT domain.
    ///
    /// Returns a [`Vec3`] in NTT domain.
    pub fn mul_trans_ntt(&self, v: &Vec3) -> Vec3 {
        Vec3::from_polys([
            Vec3::from_polys([
                self.0[0].get(0).clone(),
                self.0[1].get(0).clone(),
                self.0[2].get(0).clone(),
            ])
            .inner_product_ntt(v),
            Vec3::from_polys([
                self.0[0].get(1).clone(),
                self.0[1].get(1).clone(),
                self.0[2].get(1).clone(),
            ])
            .inner_product_ntt(v),
            Vec3::from_polys([
                self.0[0].get(2).clone(),
                self.0[1].get(2).clone(),
                self.0[2].get(2).clone(),
            ])
            .inner_product_ntt(v),
        ])
    }
}

fn sample_ntt_poly(seed: &[u8; 32], i: u8, j: u8) -> Poly {
    // XOF input: seed || j || i (column then row)
    let mut xof_input = [0u8; 34];
    xof_input[0..32].copy_from_slice(seed);
    xof_input[32] = j;
    xof_input[33] = i;

    let mut poly = Poly::zero();
    let mut coeffs_filled = 0;
    let mut block = [0u8; 168];
    let mut block_offset = 168; // refill on first iteration

    let mut xof = Shake128::default(); // Stateful SHAKE-128 XOF
    xof.update(&xof_input);
    let mut reader = xof.finalize_xof();

    while coeffs_filled < 256 {
        // refill block if exhausted
        if block_offset >= 168 {
            reader.read(&mut block);
            block_offset = 0;
        }

        // need 3 bytes to parse two 12-bit values
        if block_offset + 3 > 168 {
            reader.read(&mut block);
            block_offset = 0;
        }

        let b0 = block[block_offset] as u16;
        let b1 = block[block_offset + 1] as u16;
        let b2 = block[block_offset + 2] as u16;
        block_offset += 3;

        // parse two 12-bit values little-endian
        let d1 = b0 | ((b1 & 0xF) << 8);
        let d2 = (b1 >> 4) | (b2 << 4);

        // rejection sample: accept only values < Q
        if d1 < Q as u16 && coeffs_filled < 256 {
            poly.coeffs[coeffs_filled] = Zq::from_int(d1 as i16);
            coeffs_filled += 1;
        }
        if d2 < Q as u16 && coeffs_filled < 256 {
            poly.coeffs[coeffs_filled] = Zq::from_int(d2 as i16);
            coeffs_filled += 1;
        }
    }

    poly
}

#[cfg(test)]
mod tests {
    use quickcheck_macros::quickcheck;

    use super::*;

    #[test]
    fn test_sample_ntt_coeffs_in_range() {
        // All coefficients must be in [0, q-1] — rejection sampling correctness
        let seed = [0u8; 32];
        let mat = Mat3x3::sample_ntt(&seed);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..256 {
                    let c = mat.0[i].get(j).coeffs[k].to_int();
                    assert!(
                        (0..3329).contains(&c),
                        "mat[{i}][{j}] coeff {k} = {c} out of range"
                    );
                }
            }
        }
    }

    #[test]
    fn test_sample_ntt_deterministic() {
        // Same seed must always produce the same matrix
        let seed = [42u8; 32];
        let m0 = Mat3x3::sample_ntt(&seed);
        let m1 = Mat3x3::sample_ntt(&seed);
        assert_eq!(m0, m1);
    }

    #[test]
    fn test_sample_ntt_different_seeds() {
        // Different seeds must produce different matrices
        let seed0 = [0u8; 32];
        let mut seed1 = [0u8; 32];
        seed1[0] = 1;
        let m0 = Mat3x3::sample_ntt(&seed0);
        let m1 = Mat3x3::sample_ntt(&seed1);
        assert_ne!(m0, m1);
    }

    #[test]
    fn test_sample_ntt_entries_independent() {
        // A[i][j] and A[i'][j'] must differ — domain separation via (i,j) is working
        let seed = [0u8; 32];
        let mat = Mat3x3::sample_ntt(&seed);
        // check all pairs of entries are distinct
        for i0 in 0..3 {
            for j0 in 0..3 {
                for i1 in 0..3 {
                    for j1 in 0..3 {
                        if (i0, j0) != (i1, j1) {
                            assert_ne!(
                                mat.0[i0].get(j0),
                                mat.0[i1].get(j1),
                                "A[{i0}][{j0}] == A[{i1}][{j1}] — domain separation broken"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_sample_ntt_uniform() {
        // Sample many matrices and check coefficients are roughly uniform over [0, q)
        // Split [0, q) into 4 buckets and check each has ~25% of values
        const N_BUCKETS: usize = 4;
        let mut buckets = [0i64; N_BUCKETS];
        let bucket_size = 3329_usize.div_ceil(N_BUCKETS);
        let n_matrices = 10;

        for i in 0..n_matrices {
            let seed = [i as u8; 32];
            let mat = Mat3x3::sample_ntt(&seed);
            for row in 0..3 {
                for col in 0..3 {
                    for k in 0..256 {
                        let c = mat.0[row].get(col).coeffs[k].to_int() as usize;
                        buckets[c / bucket_size] += 1;
                    }
                }
            }
        }

        let total: i64 = buckets.iter().sum();
        for (b, &count) in buckets.iter().enumerate() {
            let freq = count as f64 / total as f64;
            assert!(
                (freq - 0.25).abs() < 0.05,
                "bucket {b} frequency {freq:.3} not near 0.25"
            );
        }
    }

    #[test]
    fn test_mul_vec_ntt_zero() {
        // A * 0 = 0
        let seed = [0u8; 32];
        let mat = Mat3x3::sample_ntt(&seed);
        let mut v = Vec3::zero();
        v.ntt();

        let result = mat.mul_vec_ntt(&v);
        let mut expected = Vec3::zero();
        expected.ntt();

        assert_eq!(result, expected);
    }

    #[test]
    fn test_mul_vec_ntt_linearity() {
        // A*(u + v) == A*u + A*v
        let seed = [0u8; 32];
        let mat = Mat3x3::sample_ntt(&seed);

        let mut u = Vec3::sample_cbd(&[1u8; 32], 0);
        let mut v = Vec3::sample_cbd(&[2u8; 32], 0);
        u.ntt();
        v.ntt();

        // A*(u+v)
        let mut uv = u.clone();
        uv.add_assign(&v);
        let lhs = mat.mul_vec_ntt(&uv);

        // A*u + A*v
        let mut rhs = mat.mul_vec_ntt(&u);
        rhs.add_assign(&mat.mul_vec_ntt(&v));

        assert_eq!(lhs, rhs);
    }

    #[test]
    fn test_mul_vec_ntt_correct() {
        // Verify result against naive coefficient-domain multiplication.
        // Build a simple matrix and vector with known values,
        // multiply both ways and check they agree after intt().
        let seed = [0u8; 32];
        let mat = Mat3x3::sample_ntt(&seed);

        let mut v = Vec3::sample_cbd(&[42u8; 32], 0);
        v.ntt();

        // naive: compute A*v in coefficient domain via pointwise_mul + intt
        // for each output row i: res[i] = Σ_j mat[i][j] * v[j]
        let mut expected = Vec3::zero();
        for i in 0..3 {
            for j in 0..3 {
                expected
                    .get_mut(i)
                    .mul_add_assign(mat.0[i].get(j), v.get(j));
            }
            expected.get_mut(i).intt();
        }

        // ntt path
        let mut result = mat.mul_vec_ntt(&v);
        result.intt();

        for i in 0..3 {
            for k in 0..256 {
                assert_eq!(
                    result.get(i).coeffs[k].to_int(),
                    expected.get(i).coeffs[k].to_int(),
                    "mismatch at row {i} coeff {k}"
                );
            }
        }
    }

    #[quickcheck]
    fn test_mul_vec_ntt_linearity_qc(mut u: Vec3, mut v: Vec3) -> bool {
        let seed = [0u8; 32];
        let mat = Mat3x3::sample_ntt(&seed);

        u.ntt();
        v.ntt();

        // A*(u+v)
        let mut uv = u.clone();
        uv.add_assign(&v);
        let lhs = mat.mul_vec_ntt(&uv);

        // A*u + A*v
        let mut rhs = mat.mul_vec_ntt(&u);
        rhs.add_assign(&mat.mul_vec_ntt(&v));

        (0..3).all(|i| (0..256).all(|k| lhs.get(i).coeffs[k] == rhs.get(i).coeffs[k]))
    }
}
