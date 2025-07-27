import math


def ncr(n: int, r: int) -> int:
    f = math.factorial
    return f(n) // f(r) // f(n - r)
