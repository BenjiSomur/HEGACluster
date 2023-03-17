from decoder import decode


def memoize_fitness(func):
    cache = {}

    def memoized_func(chrom, ref):
        key = (tuple(map(tuple, chrom)), tuple(map(tuple, ref)))
        if key in cache:
            return cache[key]
        else:
            result = func(chrom, ref)
            cache[key] = result
            return result

    return memoized_func


def memoize(func):
    cache = {}

    def memoized_func(*args):
        _ref, _chrom = args
        # convert _chrom to tuple of tuples
        key = (_ref, tuple(map(tuple, _chrom)))
        if key in cache:
            return cache[key]
        else:
            result = func(*args)
            cache[key] = result
            return result
    return memoized_func


@memoize
def get_clus(_ref, _chrom):
    return next((idx for idx, clust in enumerate(_chrom) if _ref in clust), None)


@memoize_fitness
def turbomq(chrom, ref):
    alpha = [0] * len(chrom)
    beta = [0] * len(chrom)
    chrom_sets = [set(clust) for clust in chrom]
    for idx, clust in enumerate(chrom):
        for idi in clust:
            for idj, val in enumerate(ref[idi-1]):
                if val == 0:
                    continue
                idj_aux = idj + 1
                if idj_aux in chrom_sets[idx]:
                    alpha[idx] += val
                else:
                    id_clusj = get_clus(idj_aux, chrom)
                    if id_clusj is None:
                        continue
                    beta[id_clusj] += val
                    beta[idx] += val
    return sum(alpha[i]/(alpha[i] + beta[i]/2) if alpha[i] + beta[i]/2 > 0 else 0 for i in range(len(chrom)))


def get_clusmax(_chrom):
    size = 0
    for clus in _chrom:
        _sz = len(clus)
        if _sz >= size:
            size = _sz
    return size


def get_clusmin(_chrom):
    size = 10000000
    for clus in _chrom:
        _sz = len(clus)
        if _sz < size:
            size = _sz
    return size


def fitness(_indiv, _ref, _theta):
    _chrom = decode(_indiv, len(_ref))
    tmq = turbomq(_chrom, _ref)
    nc = len(_chrom) / len(_ref)
    _deltaclus = (max(_indiv[1]) - min(_indiv[1])) / len(_ref)
    return [(tmq * _theta) + (((1-_theta)/2) * (nc + (1-_deltaclus))), tmq]


def fitness5(_chrom, _ref, _theta):
    tmq = turbomq2(_chrom, _ref)
    nc = len(_chrom) / len(_ref)
    _deltaclus = (get_clusmax(_chrom) - get_clusmin(_chrom)) / len(_ref)
    return [(tmq * _theta) + (((1-_theta)/2) * (nc + (1-_deltaclus))), tmq]


if __name__ == '__main__':
    _streval = 'from parser import create_table; import numpy as np'
    exec(_streval)
    with open('../mdgs/compiler.mdg', 'r') as f:
        data = f.readlines()
    ref = create_table(data)
    chrom = [[1, 3, 5], [2, 7, 12, 13], [4, 10], [6, 8, 11], [9]]
    # chrom = [[3, 5, 1, 2], [7, 12, 13], [4, 10], [6, 8, 11], [9]]
    # chrom = [[9,11,12,13],[1,3],[8,10],[2,5,7,8,10],[6,14,16],[4,15,17,18,19,20]]
    # chrom = [[4, 9], [2, 10], [5, 6, 1, 14, 3], [
    #    13, 7, 12, 11, 15], [8, 16, 17, 18, 20, 19]]
    # chrom = [[10], [3], [5, 11, 6], [8, 7, 13], [12]]

    # chrom = [[10, 6, 3, 14], [1, 2], [4, 5], [7, 11], [
    #    9, 8, 13, 15], [12, 20, 17, 18, 16, 19], [21, 23, 22, 24]]
    # print(np.array(ref))
    print(turbomq2(chrom, ref))
    print(fitness5(chrom, ref, 0.65))
