# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
# Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
"""Small exact Z[phi] kernel for the inverse-Euclidean tail checks.

A pair (a,b) means a+b*phi, where phi**2=phi+1. No floating-point
arithmetic is used. A normalized ratio is stored as its numerator and
positive denominator in this ring; division is unnecessary.
"""

from functools import lru_cache
from math import gcd

ZERO = (0, 0)
ONE = (1, 0)
PHI = (0, 1)
Q = (-2, 1)


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def neg(x):
    return -x[0], -x[1]


def sub(x, y):
    return add(x, neg(y))


def scale(x, n):
    return n*x[0], n*x[1]


def mul(x, y):
    a, b = x
    c, d = y
    return a*c+b*d, a*d+b*c+b*d


def power(x, n):
    if n < 0:
        raise ValueError("only nonnegative powers are needed")
    result = ONE
    while n:
        if n & 1:
            result = mul(result, x)
        x = mul(x, x)
        n //= 2
    return result


def sign(x):
    """Exact real sign, including zero; sqrt(5) is never approximated."""
    a, b = x
    c = 2*a+b
    if b == 0:
        return (c > 0)-(c < 0)
    if c == 0:
        return (b > 0)-(b < 0)
    if (c > 0) == (b > 0):
        return (c > 0)-(c < 0)
    difference = c*c-5*b*b
    assert difference != 0  # Irrationality of sqrt(5), with b nonzero.
    return ((c > 0)-(c < 0))*((difference > 0)-(difference < 0))


@lru_cache(maxsize=4096)
def q_power(n):
    return power(Q, n)


def source_limit(row):
    """C_row=phi² [3/(sqrt(5)phi)]**(row-1), exactly."""
    numerator = mul((1, 1), power((9, -3), row-1))
    denominator = 5**(row-1)
    common = gcd(gcd(abs(numerator[0]), abs(numerator[1])), denominator)
    return (numerator[0]//common, numerator[1]//common), denominator//common


SOURCE_CEILINGS = {
    8: ((707, 0), 1000),
    9: ((2929, 0), 5000),
    10: ((4857, 0), 10000),
    11: ((20111, 0), 50000),
    12: ((333511, 0), 1000000),
    13: ((2765404267, 0), 10000000000),
}


@lru_cache(maxsize=8192)
def ratio(gaps):
    """Exact matched-corner ratio; None denotes a gap at infinity."""
    k = len(gaps)-1
    assert k >= 1 and all(g is None or g >= 0 for g in gaps)
    x, y = (1, 1), (3, -2)
    for index, gap in enumerate(gaps):
        y = ZERO if gap is None else mul(y, q_power(gap))
        if index < k:
            x, y = add(scale(x, 2), y), mul(Q, sub(x, scale(y, 2)))
    numerator = add(mul(x, x), mul(y, y))
    denominator = scale((2, 3), 5**k)
    if all(g is not None for g in gaps):
        denominator = mul(denominator, sub(ONE, q_power(2*k+5+2*sum(gaps))))
    assert sign(numerator) > 0 and sign(denominator) > 0
    return numerator, denominator


def compare(value, threshold):
    numerator, denominator = value
    target, divisor = threshold
    assert divisor > 0
    return sign(sub(scale(numerator, divisor), mul(denominator, target)))


def fibonacci(n):
    if n < 0:
        raise ValueError("nonnegative Fibonacci index required")
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a+b
    return a


def endpoint(gaps):
    """Independent physical integer evaluation of a finite gap tuple."""
    u, v = 1, 2
    for index, gap in enumerate(gaps):
        assert gap is not None and gap >= 0
        if gap:
            fg, fnext = fibonacci(gap), fibonacci(gap+1)
            u, v = (fnext-fg)*u+fg*v, fg*u+fnext*v
        if index+1 < len(gaps):
            u, v = u, u+v
    return u, v


def fixed_checks():
    """Prescribed schema controls only; not a label or word enumeration."""
    assert mul(PHI, PHI) == add(PHI, ONE)
    assert sign(Q) < 0 and sign(add(Q, ONE)) > 0
    assert sign((0, 1)) > 0 and sign((2, -1)) > 0
    assert sign((1, -1)) < 0 and sign(ZERO) == 0
    prescribed = ((0, 0), (1, 0), (0, 1), (2, 1, 0),
                  (1, 0, 3, 2), (0, 2, 1, 0, 4, 3),
                  (1, 120, 2, 3, 16, 7, 0))
    for gaps in prescribed:
        a, b = endpoint(gaps)
        c = sum(gaps)+len(gaps)
        numerator, denominator = ratio(gaps)
        assert scale(numerator, fibonacci(2*c+3)) == scale(denominator, a*a+b*b)
    # Radical coefficients independently stated in the accepted row results.
    expected = {
        10: ((1850202, -1141614), 6250),
        11: ((4015332, -2480058), 6250),
        12: ((21789081, -13463172), 15625),
        13: ((94596498, -58458510), 31250),
    }
    for row, (value, divisor) in expected.items():
        actual, den = source_limit(row)
        assert scale(actual, divisor) == scale(value, den)
        assert compare((value, (divisor, 0)), SOURCE_CEILINGS[row]) < 0
    for row, ceiling in SOURCE_CEILINGS.items():
        value, denominator = source_limit(row)
        assert compare((value, (denominator, 0)), ceiling) < 0
    return len(prescribed)
