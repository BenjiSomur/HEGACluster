from decoder import decode
from turbomq import TurboMQ


class Fitness:
    def __init__(self, ref, theta, cachesize):
        self.cachesize = cachesize
        self.turbo_mq = TurboMQ(ref, cachesize)
        self.ref = tuple(map(tuple, ref))
        self.lenref = len(ref)
        self.theta = theta

    def __call__(self, _indiv):
        if self.turbo_mq.cache.cache_info().currsize >= self.cachesize:
            self.turbo_mq.cache.cache_clear()
        _chrom = decode(_indiv, self.lenref)
        chrom_tuple = tuple(map(tuple, _chrom))
        tmq = self.turbo_mq(chrom_tuple)
        nc = len(_chrom) / self.lenref
        _deltaclus = (max(_indiv[1]) - min(_indiv[1])) / self.lenref
        return [(tmq * self.theta) + (((1-self.theta)/2) * (nc + (1-_deltaclus))), tmq]


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
