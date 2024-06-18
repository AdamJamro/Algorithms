from functools import cache


@cache
def pow_c(num, exp):
    return num ** exp


def polynomial_roll_hash(s, p=53, m=1e9 + 9):
    n = len(s)
    h = 0
    for i in range(n):
        x = pow(p, i) % m
        h += (x * (ord(s[i]) - ord('a')))
        h = h % m
    return h


def paired_hash(s):
    key = (polynomial_roll_hash(s, 31, 1e9 + 9), polynomial_roll_hash(s, 53, 1e9 + 9))
    return key


for word in ["hello", "world", "hello world"]:
    print(polynomial_roll_hash(word))
    print(paired_hash(word))
    print("---")
