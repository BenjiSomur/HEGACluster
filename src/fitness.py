from decoder import decode


def calc_intraedges(_clust, _ref):
    _sum = 0
    for _idc in _clust:
        for _aux in _clust:
            if _idc == _aux:
                continue
            _sum += _ref[_idc - 1][_aux - 1]
    return _sum


def get_edge_sum(_clus1, _clus2, _ref):
    _eps = 0
    for _idi in _clus1:
        for _idj in _clus2:
            _eps += _ref[_idi - 1][_idj - 1]
            _eps += _ref[_idj - 1][_idi - 1]
    return _eps


def calc_interedges(_clust, _chrom, _ref):
    _sum = 0
    for _aux in _chrom:
        if _aux[0] in _clust:
            continue
        _sum += get_edge_sum(_clust, _aux, _ref)
    return _sum


def get_clus(_ref, _chrom):
    aux = list(range(len(_chrom)))
    for x in aux:
        if _ref in _chrom[x]:
            return x


def turbomq2(_chrom, _ref):
    alpha = [0]*len(_chrom)
    beta = [0]*len(_chrom)
    for _idx, _clust in enumerate(_chrom):
        for _idi in _clust:
            for _idj, val in enumerate(_ref[_idi-1]):
                if val == 0:
                    continue
                _idj_aux = _idj+1
                if _idj_aux in _clust:
                    alpha[_idx] += val
                else:
                    id_clusj = get_clus(_idj_aux, _chrom)
                    beta[_idx] += val
                    beta[id_clusj] += val
    _sum = 0
    for x in range(len(_chrom)):
        _sum += alpha[x]/(alpha[x] + (beta[x]/2))
        #_sum += (2*alpha[x])/((2*alpha[x]) + beta[x])
    return _sum


def turbomq(_chrom, _ref):
    _sum = 0
    for _clust in _chrom:
        _mu = calc_intraedges(_clust, _ref)
        if _mu == 0:
            continue
        _eps = calc_interedges(_clust, _chrom, _ref)
        _mf = _mu / (_mu + (_eps / 2))
        _sum += _mf
    return _sum


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


def fitness3(_indiv, _ref, _theta):
    _chrom = decode(_indiv, len(_ref))
    tmq = turbomq(_chrom, _ref)
    nc = len(_chrom) / len(_ref)
    _deltaclus = (max(_indiv[1]) - min(_indiv[1])) / len(_ref)
    return [(tmq * _theta) + (((1-_theta)/2) * (nc - _deltaclus)), tmq]


def fitness(_indiv, _ref, _theta):
    _chrom = decode(_indiv, len(_ref))
    tmq = turbomq2(_chrom, _ref)
    nc = len(_chrom) / len(_ref)
    _deltaclus = (max(_indiv[1]) - min(_indiv[1])) / len(_ref)
    return [(tmq * _theta) + (((1-_theta)/2) * (nc + (1-_deltaclus))), tmq]


def fitness4(_indiv, _ref, _theta):
    _chrom = decode(_indiv, len(_ref))
    tmq = turbomq(_chrom, _ref)
    nc = len(_chrom) / len(_ref)
    _deltaclus = (max(_indiv[1]) - min(_indiv[1])) / len(_ref)
    return [(tmq * _theta) + (nc * ((1 - _theta)/2)) + ((1 - _deltaclus)*((1-_theta)/2)), tmq]


def fitness2(_chrom, _ref, _theta):
    tmq = turbomq(_chrom, _ref)
    nc = len(_chrom) / len(_ref)
    _deltaclus = (get_clusmax(_chrom) - get_clusmin(_chrom)) / len(_ref)
    return [(tmq * _theta) + (nc * ((1 - _theta)/2)) + ((1 - _deltaclus)*((1-_theta)/2)), tmq]


def fitness5(_chrom, _ref, _theta):
    tmq = turbomq2(_chrom, _ref)
    nc = len(_chrom) / len(_ref)
    _deltaclus = (get_clusmax(_chrom) - get_clusmin(_chrom)) / len(_ref)
    return [(tmq * _theta) + (((1-_theta)/2) * (nc + (1-_deltaclus))), tmq]


if __name__ == '__main__':
    _streval = 'from parser import create_table'
    exec(_streval)
    with open('../mdgs/rcs.mdg', 'r') as f:
        data = f.readlines()
    ref = create_table(data)
    chrom = [[14, 29], [21, 26, 27], [5, 20, 25], [9, 13], [23, 24, 28], [
        1, 2, 3, 4, 6, 7, 8, 10, 11, 12, 15, 16, 17, 18, 19, 22]]
    print(fitness5(chrom, ref, 0.65))
