use core::fmt;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Goldilocks field: p = 2^64 - 2^32 + 1.
///
/// We deliberately use a simple u128 reduction here so every basis uses the
/// same field implementation. Once the benchmark logic is stable we can swap
/// this for a production Goldilocks/M31 implementation without changing the
/// basis algorithms.
#[derive(Copy, Clone, Default, PartialEq, Eq)]
#[repr(transparent)]
pub struct F(pub u64);

impl F {
    pub const MODULUS: u64 = 0xffff_ffff_0000_0001;
    pub const ZERO: Self = Self(0);
    pub const ONE: Self = Self(1);
    pub const GENERATOR: Self = Self(7);

    #[inline(always)]
    pub const fn from_canonical_u64(x: u64) -> Self {
        Self(x)
    }

    #[inline(always)]
    pub fn new(x: u128) -> Self {
        Self(Self::reduce_u128(x))
    }

    #[inline(always)]
    fn reduce_u128(x: u128) -> u64 {
        // Goldilocks:
        // p = 2^64 - 2^32 + 1
        // therefore 2^64 = 2^32 - 1 (mod p).
        const C: u128 = (1u128 << 32) - 1;
        const MASK: u128 = (1u128 << 64) - 1;

        let lo = x & MASK;
        let hi = x >> 64;

        let y = lo + hi * C;

        let lo2 = y & MASK;
        let hi2 = y >> 64;

        let mut z = lo2 + hi2 * C;

        let p = Self::MODULUS as u128;
        if z >= p {
            z -= p;
        }

        z as u64
    }

    #[inline(always)]
    pub fn square(self) -> Self {
        self * self
    }

    #[inline]
    pub fn pow(mut self, mut e: u64) -> Self {
        let mut acc = Self::ONE;
        while e != 0 {
            if e & 1 == 1 {
                acc *= self;
            }
            self = self.square();
            e >>= 1;
        }
        acc
    }

    #[inline]
    pub fn inverse(self) -> Self {
        assert!(self != Self::ZERO);
        self.pow(Self::MODULUS - 2)
    }

    /// Primitive N-th root of unity for N = 2^k, k <= 32.
    pub fn primitive_root_of_unity(n: usize) -> Self {
        assert!(n.is_power_of_two());
        assert!(n <= (1usize << 32));
        let exp = (Self::MODULUS - 1) / (n as u64);
        Self::GENERATOR.pow(exp)
    }
}

impl fmt::Debug for F {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Add for F {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        let sum = self.0 as u128 + rhs.0 as u128;
        let p = Self::MODULUS as u128;

        if sum >= p {
            Self((sum - p) as u64)
        } else {
            Self(sum as u64)
        }
    }
}
impl AddAssign for F {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) { *self = *self + rhs; }
}
impl Sub for F {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        if self.0 >= rhs.0 {
            Self(self.0 - rhs.0)
        } else {
            Self(Self::MODULUS - (rhs.0 - self.0))
        }
    }
}
impl SubAssign for F {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Self) { *self = *self - rhs; }
}
impl Neg for F {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self { if self == Self::ZERO { self } else { Self(Self::MODULUS - self.0) } }
}
impl Mul for F {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: Self) -> Self {
        Self(Self::reduce_u128(self.0 as u128 * rhs.0 as u128))
    }
}
impl MulAssign for F {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Self) { *self = *self * rhs; }
}
