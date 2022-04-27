#![feature(generic_const_exprs)]

use std::ops::{Add, Mul, Sub};

type Mask = u32;

/// Calculate factorial for binominal calculation.
const fn factorial(n: usize) -> usize {
    let mut acc = 1;
    let mut index = 1;
    while index < n {
        index += 1;
        acc *= index;
    }
    acc
}

/// Calculate binominal to estimate k-blade basis count in space N.
pub const fn binominal(n: usize, k: usize) -> usize {
    factorial(n) / factorial(n - k) / factorial(k)
}

/// Generate set of basis combinations.
const fn combinations<const N: usize, const K: usize>() -> [[usize; K]; binominal(N, K)] {
    let mut result = [[0; K]; binominal(N, K)];
    let mut current = [0; K];
    let mut index = 0;
    while index < K {
        current[index] = index;
        index += 1;
    }

    let mut index = 0;
    while index < result.len() {
        result[index] = current;

        let mut pointer = K - 1;

        while pointer != 0 && current[pointer] == N - K + pointer {
            pointer -= 1;
        }
        current[pointer] += 1;
        let mut counter = pointer + 1;
        while counter < K {
            current[counter] = current[counter - 1] + 1;
            counter += 1;
        }

        index += 1;
    }
    result
}

/// Calculate basis mask.
const fn combinations_mask<const N: usize, const K: usize>() -> [Mask; binominal(N, K)] {
    let generated = combinations::<N, K>();
    let mut result = [0; binominal(N, K)];
    let mut index = 0;
    while index < result.len() {
        let mut pointer = 0;
        while pointer < K {
            result[index] ^= 1 << generated[index][pointer];
            pointer += 1;
        }
        index += 1;
    }
    result
}

/// Permutation sign.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Sign {
    Negative,
    Zero,
    Positive,
}

/// Estimate basis mask and permutation sign after basis wedge.
fn basis_mask_and_sign(left_mask: Mask, right_mask: Mask) -> (Mask, Sign) {
    let mask = left_mask | right_mask;
    if left_mask & right_mask != 0 {
        (0, Sign::Zero)
    } else {
        let mut counter = 0;

        let mut right_index = 0;
        while right_index < Mask::BITS - 1 {
            if right_mask & 1 << right_index != 0 {
                let test_mask = !((1 << right_index) - 1);
                counter += (left_mask & test_mask).count_ones();
            }

            right_index += 1;
        }
        (
            mask,
            if counter % 2 == 0 {
                Sign::Positive
            } else {
                Sign::Negative
            },
        )
    }
}

/// K-blade in N-space.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Blade<const K: usize, T, const N: usize>
where
    [(); binominal(N, K)]:,
{
    value: [T; binominal(N, K)],
}

impl<const K: usize, T, const N: usize> Blade<K, T, N>
where
    [(); binominal(N, K)]:,
{
    /// Create blade from the array.
    pub fn new(value: [T; binominal(N, K)]) -> Self {
        Self { value }
    }

    /// Get basis amount required to describe this blade.
    pub fn size(&self) -> usize {
        binominal(N, K)
    }

    /// Get blade grade.
    pub fn grade(&self) -> usize {
        K
    }

    /// Get blade space dimension.
    pub fn dimension(&self) -> usize {
        N
    }

    /// Get individual value.
    pub fn value(&self, index: usize) -> Option<&T> {
        self.value.get(index)
    }

    /// Get individual value mutably.
    pub fn value_mut(&mut self, index: usize) -> Option<&mut T> {
        self.value.get_mut(index)
    }
}

impl<const K: usize, T, const N: usize> Blade<K, T, N>
where
    [(); binominal(N, K)]:,
    T: Default + Copy + Mul<Output = T> + Add<Output = T> + Sub<Output = T>,
{
    /// Perform the wedge operation.
    pub fn wedge<const L: usize>(&self, other: Blade<L, T, N>) -> Blade<{ K + L }, T, N>
    where
        [(); binominal(N, L)]:,
        [(); binominal(N, K + L)]:,
    {
        let mut data = [T::default(); binominal(N, K + L)];

        let left_mask = combinations_mask::<N, K>();
        let right_mask = combinations_mask::<N, L>();

        let result_mask = combinations_mask::<N, { K + L }>();

        for (k_index, k_value) in left_mask.iter().enumerate() {
            for (l_index, l_value) in right_mask.iter().enumerate() {
                let (mask, sign) = basis_mask_and_sign(*k_value, *l_value);
                if sign != Sign::Zero {
                    // This one is not optimal, but it was easy to implement.
                    for (index, value) in data.iter_mut().enumerate() {
                        if result_mask[index] == mask {
                            if sign == Sign::Positive {
                                *value = *value + self.value[k_index] * other.value[l_index]
                            } else {
                                *value = *value - self.value[k_index] * other.value[l_index];
                            }
                        }
                    }
                }
            }
        }

        Blade::new(data)
    }
}

impl<const K: usize, T, const N: usize> Default for Blade<K, T, N>
where
    [(); binominal(N, K)]:,
    T: Copy + Default,
{
    fn default() -> Self {
        let values = [T::default(); binominal(N, K)];
        Self::new(values)
    }
}

pub type Scalar<T, const N: usize> = Blade<0, T, N>;
pub type PseudoScalar<T, const N: usize> = Blade<N, T, N>;

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn blade_len_is_sane() {
        let scalar = Blade::<0, _, 3>::new([1.0]);
        assert_eq!(scalar.size(), 1);

        let vector = Blade::<1, _, 3>::new([1.0, 2.0, 3.0]);
        assert_eq!(vector.size(), 3);

        let bivector = Blade::<2, _, 3>::new([1.0, 2.0, 3.0]);
        assert_eq!(bivector.size(), 3);

        let trivector = Blade::<3, _, 3>::new([1.0]);
        assert_eq!(trivector.size(), 1);
    }

    #[test]
    fn wedge_is_sane() {
        let e1 = Blade::<1, _, 3>::new([1, 0, 0]);
        let e2 = Blade::<1, _, 3>::new([0, 1, 0]);
        let e3 = Blade::<1, _, 3>::new([0, 0, 1]);
        let e12 = e1.wedge(e2);
        assert_eq!(e12, Blade::new([1, 0, 0]));
        let e123 = e12.wedge(e3);
        assert_eq!(e123, Blade::new([1]));

        let zero = e12.wedge(e1);
        assert_eq!(zero, Blade::new([0]));

        let a = Blade::<1, _, 2>::new([3, 2]);
        let b = Blade::<1, _, 2>::new([2, 3]);
        let i = a.wedge(b);
        assert_eq!(i, Blade::new([5]));

        let a = Blade::<1, _, 2>::new([3, 2]);
        let b = Blade::<1, _, 2>::new([-1, 3]);
        let i = a.wedge(b);
        assert_eq!(i, Blade::new([11]));
    }

    #[test]
    fn combinations_is_sane() {
        let expected: [[usize; 3]; 10] = [
            [0, 1, 2],
            [0, 1, 3],
            [0, 1, 4],
            [0, 2, 3],
            [0, 2, 4],
            [0, 3, 4],
            [1, 2, 3],
            [1, 2, 4],
            [1, 3, 4],
            [2, 3, 4],
        ];
        assert_eq!(combinations::<5, 3>(), expected);
    }

    #[test]
    fn combinations_mask_is_sane() {
        let expected = [
            // 8765_43210
            0b_0000_0111, //
            0b_0000_1011, //
            0b_0001_0011, //
            0b_0000_1101, //
            0b_0001_0101, //
            0b_0001_1001, //
            0b_0000_1110, //
            0b_0001_0110, //
            0b_0001_1010, //
            0b_0001_1100, //
        ];
        assert_eq!(combinations_mask::<5, 3>(), expected);
    }

    #[test]
    fn basis_mask_sign_is_sane() {
        let e0 = 0b_0000_0001;
        let e01 = 0b_0000_0011;
        let e0123 = 0b_0000_1111;
        let e2356 = 0b_0110_1100;
        let e235 = 0b_0010_1100;

        assert_eq!(basis_mask_and_sign(e0123, e0).1, Sign::Zero);
        assert_eq!(basis_mask_and_sign(e0, e0123).1, Sign::Zero);

        assert_eq!(basis_mask_and_sign(e0, e2356).1, Sign::Positive);
        assert_eq!(basis_mask_and_sign(e2356, e0).1, Sign::Positive);

        assert_eq!(basis_mask_and_sign(e0, e235).1, Sign::Positive);
        assert_eq!(basis_mask_and_sign(e235, e0).1, Sign::Negative);

        assert_eq!(basis_mask_and_sign(e01, e235).1, Sign::Positive);
        assert_eq!(basis_mask_and_sign(e235, e01).1, Sign::Positive);
    }
}
