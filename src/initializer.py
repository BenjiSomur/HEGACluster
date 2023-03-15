import math


# def get_theta(nodes):
#     n = len(nodes)
#     if n >= 10000:
#         return 0.9
#     else:
#         return 1 - (1 / math.sqrt(n))

def get_theta(nodes):
    n = len(nodes)
    if n >= 10000:
        return 0.9
    else:
        return 1 - ((1 / math.log2(n))**2)


def get_pop_size(nodes):
    n = len(nodes)
    if n <= 300:
        return 2*(n+1)
    elif n <= 3000:
        return n
    elif n > 3000:
        return n // 4
    else:
        return None


def get_no_gen(nodes):
    n = len(nodes)
    if n <= 300:
        return 50 * n
    elif n <= 3000:
        return 20 * n
    elif n <= 10000:
        return 5 * n
    elif n > 10000:
        return n
    else:
        return None


def get_cp(nodes):
    n = len(nodes)
    if n <= 100:
        return 0.83
    elif n < 1000:
        x = 0.2 * (n - 100)
        return 0.8 + (x / 899)
    elif n >= 10000:
        return 1
    else:
        return None


# def get_mp(nodes):
#     n = len(nodes)
#     x = math.log2(n)
#     return 15 / (17 * x)

def get_mp(nodes):
    n = len(nodes)
    x = math.sqrt(n)
    return 16 / (10 * x)

# def get_mp(nodes):
#     n = len(nodes)
#     x = math.log2(n)
#     return 16 / (0.004 * x)
