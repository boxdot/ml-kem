use crate::integer::Zq;

/// Precomputed primitive 256th roots of unity:
///
/// zeta_i = 17^(2 * brv7(i) + 1) mod Q
///
/// Exponent 2k + 1 ensures that we select only the primitive roots of unity,
/// which are the roots of cyclotomic polynomial PHI(512).
pub const BASE_ZETAS: [Zq; 128] = {
    let mut zetas = [Zq::zero(); 128];
    let mut i = 0;
    let generator = Zq::from_int(17);
    while i < 128 {
        let brv_i = brv7(i as u8);
        let exp = (brv_i as u32) * 2 + 1;
        zetas[i] = generator.pow(exp);
        i += 1;
    }
    zetas
};

/// Twiddle factor for the 7 layers of the NTT (number theoretic transform).
///
/// Note: zeta_ntt_i is not used to keep indexing 1-based.
pub const ZETAS_NTT: [Zq; 128] = {
    let mut zetas = [Zq::zero(); 128];
    let mut i = 0;
    let generator = Zq::from_int(17);
    while i < 128 {
        let brv_i = brv7(i as u8);
        zetas[i] = generator.pow(brv_i as u32);
        i += 1;
    }
    zetas
};

/// 7-bit reversal
const fn brv7(mut n: u8) -> u8 {
    let mut res = 0;
    let mut i = 0;
    while i < 7 {
        res = (res << 1) | (n & 1);
        n >>= 1;
        i += 1;
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zeta_generator() {
        let generator = Zq::from_int(17);
        assert_eq!(BASE_ZETAS[0], generator);
    }
}
