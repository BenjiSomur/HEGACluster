from functools import lru_cache


class TurboMQ:
    def __init__(self, ref, cache_size):
        self.ref = tuple(map(tuple, ref))
        self.cache = lru_cache(maxsize=cache_size)(self.__call__)

    @staticmethod
    def get_clus(refid, chrom):
        for idx, clust in enumerate(chrom):
            if refid in clust:
                return idx
        return None

    def __call__(self, chrom):
        alpha = [0] * len(chrom)
        beta = [0] * len(chrom)
        chrom_sets = [set(clust) for clust in chrom]
        for idx, clust in enumerate(chrom):
            for idi in clust:
                for idj, val in enumerate(self.ref[idi-1]):
                    if val == 0:
                        continue
                    idj_aux = idj + 1
                    if idj_aux in chrom_sets[idx]:
                        alpha[idx] += val
                    else:
                        id_clusj = self.get_clus(idj_aux, chrom)
                        if id_clusj is None:
                            continue
                        beta[id_clusj] += val
                        beta[idx] += val
        return sum(alpha[i]/(alpha[i] + beta[i]/2) if alpha[i] + beta[i]/2 > 0 else 0 for i in range(len(chrom)))
