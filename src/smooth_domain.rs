/// Returns `true` if `n` is a smooth-{2,3} number, i.e. `n = 2^a * 3^b`.
///
/// `0` is NOT smooth.
#[inline]
pub const fn is_smooth(n: usize) -> bool {
    if n == 0 {
        return false;
    }
    let mut v = n;
    while v % 2 == 0 {
        v /= 2;
    }
    while v % 3 == 0 {
        v /= 3;
    }
    v == 1
}

/// Returns the odd part of `n`: after removing all factors of 2,
/// the remaining value `3^b`.
#[inline]
pub const fn odd_part(n: usize) -> usize {
    assert!(n > 0);
    let mut v = n;
    while v % 2 == 0 {
        v /= 2;
    }
    v
}

/// Number of extra sumcheck rounds needed for the odd residual:
/// `ceil(log2(odd_part))`.
#[inline]
pub const fn extra_rounds(n: usize) -> usize {
    let odd = odd_part(n);
    if odd == 1 {
        0
    } else {
        (usize::BITS - (odd - 1).leading_zeros()) as usize
    }
}
