//! Integer arithmetic modulo Q represented in Montgomery form a * R where R = 2^16.

use std::ops::{Add, Mul, Sub};

/// Prime modulus
const Q: i16 = 3329;

/// -Q^-1 mod R
const Q_INV: i16 = -3327;

/// R^2 mod Q
const R2_MOD_Q: i32 = 1353;

/// Internal representation in Montgomery form: x * R mod Q
#[derive(Copy, Clone, PartialEq, Eq)]
pub struct Zq(i16);

impl Zq {
    pub fn from_int(x: i16) -> Self {
        let a = x as i32 * R2_MOD_Q;
        // a * R^2 * R^-1 = a * R mod Q
        Zq(montgomery_reduce(a))
    }

    pub fn to_int(self) -> i16 {
        let a = montgomery_reduce(self.0 as i32);
        if a < 0 { a + Q } else { a }
    }

    pub fn zero() -> Self {
        Zq(0)
    }
}

impl Add for Zq {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut res = self.0 + rhs.0;
        res -= Q;
        // Conditional addition if negative
        let mask = res >> 15;
        res += Q & mask;
        Zq(res)
    }
}

impl Sub for Zq {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut res = self.0 - rhs.0;
        // Conditional addition if negative
        let mask = res >> 15;
        res += Q & mask;
        Zq(res)
    }
}

impl Mul for Zq {
    type Output = Self;

    /// Montgomery multiplication: returns a * b * R^-1 mod Q
    ///
    /// <https://en.wikipedia.org/wiki/Montgomery_modular_multiplication>
    fn mul(self, rhs: Self) -> Self::Output {
        let product = self.0 as i32 * rhs.0 as i32;
        Zq(montgomery_reduce(product))
    }
}

/// Montgomery Reduction, that is, division by R mod Q
///
/// Maps a 32-bit value x to x * R^-1 mod Q
fn montgomery_reduce(a: i32) -> i16 {
    let u = (a as i16).wrapping_mul(Q_INV) as i32;
    let t = u * Q as i32;
    let res = (a - t) >> 16;
    res as i16
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
