//! Integer arithmetic modulo Q represented in Montgomery form a * R where R = 2^16.

use std::{
    fmt,
    ops::{Add, Mul, Neg, Sub},
};

/// Prime modulus
const Q: i16 = 3329;

/// -Q^-1 mod R
const Q_INV: i16 = -3327;

/// R^2 mod Q
const R2_MOD_Q: i32 = 1353;

/// Internal representation in Montgomery form: x * R mod Q
#[derive(Copy, Clone, PartialEq, Eq)]
pub struct Zq(i16);

impl fmt::Debug for Zq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Zq({}R)", self.0)
    }
}

impl fmt::Display for Zq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_int())
    }
}

impl Zq {
    pub const fn from_int(x: i16) -> Self {
        let a = x as i32 * R2_MOD_Q;
        // a * R^2 * R^-1 = a * R mod Q
        Zq(montgomery_reduce(a))
    }

    pub const fn to_int(self) -> i16 {
        let a = montgomery_reduce(self.0 as i32);
        if a < 0 { a + Q } else { a }
    }

    pub const fn zero() -> Self {
        Zq(0)
    }

    pub const fn one() -> Self {
        Self::from_int(1)
    }

    pub const fn add(self, rhs: Self) -> Self {
        // Normalize inputs from (-Q, Q) to [0, Q)
        let a = self.0 + (Q & (self.0 >> 15));
        let b = rhs.0 + (Q & (rhs.0 >> 15));
        // a + b in [0, 2Q)
        let mut res = a + b;
        res -= Q;
        // Conditional addition if negative
        let mask = res >> 15;
        res += Q & mask;
        Zq(res)
    }

    pub const fn sub(self, rhs: Self) -> Self {
        // Normalize inputs from (-Q, Q) to [0, Q)
        let a = self.0 + (Q & (self.0 >> 15));
        let b = rhs.0 + (Q & (rhs.0 >> 15));
        // a - b in (-(Q-1), Q-1)
        let mut res = a - b;
        // Conditional addition if negative
        let mask = res >> 15;
        res += Q & mask;
        Zq(res)
    }

    pub const fn neg(self) -> Self {
        Zq::zero().sub(self)
    }

    /// Montgomery multiplication: returns a * b * R^-1 mod Q
    ///
    /// <https://en.wikipedia.org/wiki/Montgomery_modular_multiplication>
    pub const fn mul(self, rhs: Self) -> Self {
        let a = self.0 as i32;
        let b = rhs.0 as i32;
        let c = a * b;
        Zq(montgomery_reduce(c))
    }

    pub const fn pow(self, mut exp: u32) -> Self {
        let mut res = Zq::one();
        let mut base = self;
        while exp > 0 {
            if exp % 2 == 1 {
                res = res.mul(base);
            }
            base = base.mul(base);
            exp /= 2;
        }
        res
    }
}

impl Add for Zq {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self.add(rhs)
    }
}

impl Sub for Zq {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self.sub(rhs)
    }
}

impl Mul for Zq {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.mul(rhs)
    }
}

impl Neg for Zq {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Zq::zero() - self
    }
}

/// Montgomery Reduction, that is, division by R mod Q
///
/// Maps a 32-bit value x to x * R^-1 mod Q
pub const fn montgomery_reduce(a: i32) -> i16 {
    let u = (a as i16).wrapping_mul(Q_INV) as i32;
    let t = u * Q as i32;
    let res = (a - t) >> 16;
    res as i16
}

#[cfg(test)]
mod tests {
    use quickcheck::{Arbitrary, Gen};

    use super::*;

    impl Arbitrary for Zq {
        fn arbitrary(g: &mut Gen) -> Self {
            Zq::from_int(Arbitrary::arbitrary(g))
        }
    }

    #[test]
    fn test_from_to_int() {
        // Test round-trip for all elements in the field
        for i in 0..Q {
            let z = Zq::from_int(i);
            assert_eq!(z.to_int(), i, "Round-trip failed for {}", i);
        }
    }

    #[test]
    fn test_additive_identity() {
        let a = Zq::from_int(1234);
        let zero = Zq::zero();
        assert_eq!(a.add(zero).to_int(), 1234);
        assert_eq!(a.sub(a).to_int(), 0);
    }

    #[test]
    fn test_multiplicative_identity() {
        let a = Zq::from_int(1234);
        let one = Zq::from_int(1);
        assert_eq!(a.mul(one).to_int(), 1234);
    }

    #[test]
    fn test_arithmetic_properties() {
        let a = Zq::from_int(1000);
        let b = Zq::from_int(2500);
        let c = Zq::from_int(3);

        // Modular addition
        assert_eq!(a.add(b).to_int(), 171);

        // Distributivity: (a + b) * c = a*c + b*c
        let lhs = a.add(b).mul(c);
        let rhs = a.mul(c).add(b.mul(c));
        assert_eq!(lhs.to_int(), rhs.to_int());
    }

    #[test]
    fn test_montgomery_reduction_bug() {
        let x = Zq::from_int(Q - 1);
        let y = Zq::from_int(Q - 1);

        // (Q-1)*(Q-1) mod Q should be (-1)*(-1) mod Q = 1
        assert_eq!(x.mul(y).to_int(), 1);
    }

    #[test]
    fn test_addition_overflow() {
        // 3328 + 1 = 3329 ≡ 0 mod Q
        let a = Zq::from_int(Q - 1);
        let b = Zq::from_int(1);
        assert_eq!(a.add(b).to_int(), 0);

        // Test a large overflow: (Q-1) + (Q-1)
        assert_eq!(a.add(a).to_int(), Q - 2);
    }

    #[test]
    fn test_subtraction_underflow() {
        let a = Zq::zero();
        let b = Zq::from_int(1);
        assert_eq!(a.sub(b).to_int(), Q - 1);

        // Test a larger underflow: 1 - 100
        let c = Zq::from_int(1);
        let d = Zq::from_int(100);
        assert_eq!(c.sub(d).to_int(), Q - 99);
    }

    #[test]
    fn test_internal_consistency() {
        // Verify that the internal representation (Montgomery form) also obeys the bounds [0, Q)
        let a = Zq::from_int(Q - 1);
        let b = Zq::from_int(1);
        let sum = a.add(b);

        // Even the raw Montgomery bits should be 0, not Q
        assert!(sum.0 >= 0 && sum.0 < Q);
    }

    #[test]
    fn test_montgomery_pow_consistency() {
        let base = Zq::from_int(3);
        // 3^4 = 81
        let res = base.pow(4);
        assert_eq!(res.to_int(), 81);

        // 17^128 ≡ 3328 ≡ -1 (mod 3329)
        let generator = Zq::from_int(17);
        let root = generator.pow(128);
        assert_eq!(root.to_int(), 3328);
    }
}
