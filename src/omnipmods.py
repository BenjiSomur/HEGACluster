from decoder import decode
from math import ceil
from fitness import get_clus, turbomq, fitness
from itertools import chain
from copy import deepcopy


def related_to(ref, chrom, modn):
    rels = []
    for no_cl, clus in enumerate(chrom):
        if modn in clus:
            continue
        for modaux in clus:
            if no_cl in rels:
                break
            if ref[modn][modaux-1] > 0 or ref[modaux-1][modn] > 0:
                rels.append(no_cl)
    return rels


def detomns(gbls, chrom, ref):
    for idx in glbls:
        refs = related_to(ref, chrom, idx - 1)
        if len(refs) >= 2:
            yield idx


def from_clus(omnid, ref, chrom):
    depends = [0]*len(chrom)
    nodeps = 0
    for idx, clus in enumerate(chrom):
        for modaux in clus:
            if ref[omnid][modaux-1] > 0 or ref[modaux-1][omnid] > 0:
                depends[idx] += ref[omnid][modaux-1]
                depends[idx] += ref[modaux - 1][omnid]
                nodeps += 1
    auxsums = [*map(lambda x: x / nodeps, depends)]
    return [(idx, score) for idx, score in enumerate(auxsums)]


def extract_omnis(omnis, chrom):
    _aux = chrom[:]
    origs = []
    for nmod in omnis:
        for idx in range(len(_aux)):
            if nmod in _aux[idx]:
                origs.append((nmod, idx))
                _aux[idx].pop(_aux[idx].index(nmod))
                break
    aux2 = remempt(_aux)
    return (aux2, origs)


def remempt(chrom):
    aux = list()
    for clus in chrom:
        if len(clus) == 0:
            continue
        aux.append(list(clus))
    return aux


def get_highest(dscores, gamma):
    aux = list(dscores)
    aux.sort(key=lambda x: x[1])
    return aux[len(aux) - gamma:]


def encode(chrom, nodes):
    max_bits = len("{0:b}".format(nodes))
    ratio = "0{}b".format(max_bits)
    indiv = list()
    intpart = [len(x) for x in chrom]
    binflat = list(chain(*chrom))
    aux = [format(x, ratio) for x in binflat]
    indiv.append(''.join(aux))
    indiv.append(intpart)
    return indiv


def rep_solut(sol, gbls, ref, theta):
    auxsol = list(sol)
    chrom = decode(auxsol[0], len(ref))
    gamma = ceil(0.2 * len(chrom))
    omns = [x for x in detomns(gbls, chrom, ref)]
    chrprim, orgs = extract_omnis(omns, chrom)
    mqprim = turbomq(chrprim, ref)
    for idx, omn in enumerate(omns):
        dscores = from_clus(omn, ref, chrprim)
        gamma = ceil(0.2 * len(chrom))
        highscores = get_highest(dscores, gamma)
        deltamq = [0] * len(highscores)
        for idy, (clusno, _) in enumerate(highscores):
            auxchr = deepcopy(chrprim)
            auxchr[clusno].append(omn)
            auxmq = turbomq(auxchr, ref)
            deltamq[idy] += (auxmq - mqprim)
        maxid = deltamq.index(max(deltamq))
        if deltamq[maxid] > 0:
            chrprim[highscores[maxid][0]].append(omn)
            mqprim = deltamq[maxid]
        else:
            orig = orgs[idx][1]
            chrprim[orig].append(omn)
    ressol = encode(chrprim, len(ref))
    fitprim = fitness(ressol, ref, theta)
    if fitprim[1] > sol[1][1]:
        return [ressol, fitprim]
    else:
        return auxsol


if __name__ == '__main__':
    _streval = 'from parser import create_table, get_globals'
    exec(_streval)
    with open('../mdgs/compiler.mdg', 'r') as f:
        data = f.readlines()
    ref = create_table(data)
    glbls = [*get_globals(ref)]
    chrom = [[1, 3, 5], [2, 7, 12, 13], [4, 10], [6, 8, 11], [9]]
    sol = encode(chrom, len(ref))
    testfit = fitness(sol, ref, 0.65)
    print([*sol, testfit])

    sol2 = rep_solut([sol, testfit], glbls, ref, 0.65)
    print(sol2)
