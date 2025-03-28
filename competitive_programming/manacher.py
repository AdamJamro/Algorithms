def manacher(s):
    # reduce the problem to odd length palindromes
    s = '#' + '#'.join(s) + '#'
    res = manacher_odd(s)
    return res


def manacher_odd(s):
    size = len(s)
    s = "$" + s + "^"
    p = [0] * (size + 2)
    left, right = 0, 1
    for i in range(size):
        p[i] = min(right - i, p[left + right - i]) if right > i else 0
        while s[i + p[i] + 1] == s[i - p[i] - 1]:
            p[i] += 1
        if i + p[i] > right:
            left = i - p[i]
            right = i + p[i]
    return p[1:-1]  # remove the $ and ^
    # alternatively, return max(p)



