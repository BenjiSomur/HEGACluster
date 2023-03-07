from decoder import decode
from math import ceil
from fitness import get_clus, turbomq, fitness
from itertools import chain
from copy import deepcopy
from parser import transpose_mat
from random import shuffle
import numpy as np


def related_to(ref, chrom, modn):
    rels = list()
    for no_cl, clus in enumerate(chrom):
        if modn in clus:
            continue
        for modaux in clus:
            if no_cl in rels:
                break
            if ref[modn - 1][modaux-1] > 0 or ref[modaux-1][modn - 1] > 0:
                rels.append(no_cl)
    return rels


def detomns(glbls, chrom, ref):
    for idx in glbls:
        refs = related_to(ref, chrom, idx)
        if len(refs) >= 2:
            yield idx


def from_clus(omnid, ref, chrom):
    depends = [0]*len(chrom)
    auxref = transpose_mat(ref)
    nodeps = sum(ref[omnid - 1]) + sum(auxref[omnid - 1])
    for idx, clus in enumerate(chrom):
        for modaux in clus:
            if ref[omnid - 1][modaux-1] > 0 or ref[modaux-1][omnid - 1] > 0:
                depends[idx] += ref[omnid - 1][modaux-1]
                depends[idx] += ref[modaux - 1][omnid - 1]
    try:
        auxsums = [*map(lambda x: x / nodeps, depends)]
    except ZeroDivisionError:
        print(np.array(ref))
        print(*depends)
        print(omnid)
        print(chrom)
        raise
    return [(idx, score) for idx, score in enumerate(auxsums)]


def extract_omnis(omnis, chrom):
    _aux = deepcopy(chrom)
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
        aux.append(deepcopy(clus))
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


def get_omnilocals(ref, chrom):
    for clus in chrom:
        for nomod in clus:
            omnirefs = related_to(ref, chrom, nomod)
            if len(omnirefs) > 2:
                yield nomod


def extendomns(omns, lcls):
    aux = deepcopy(omns)
    for auxid in lcls:
        if auxid in aux:
            continue
        aux.append(auxid)
    return aux


def rep_solut(sol, gbls, ref, theta):
    auxsol = list(sol)
    chrom = decode(auxsol[0], len(ref))
    gamma = ceil(0.2 * len(chrom)) + 1
    omns = [x for x in detomns(gbls, chrom, ref)]
    omnilocals = [y for y in get_omnilocals(ref, chrom)]
    omns = extendomns(omns, omnilocals[:gamma])
    chrprim, orgs = extract_omnis(omns, chrom)
    for idx, omn in enumerate(omns):
        mqprim = turbomq(chrprim, ref)
        dscores = from_clus(omn, ref, chrprim)
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
        else:
            try:
                orig = orgs[idx][1]
                chrprim[orig].append(omn)
            except IndexError:
                chrprim.append([omn])
    ressol = encode(chrprim, len(ref))
    fitprim = fitness(ressol, ref, theta)
    if fitprim[1] >= sol[1][1]:
        return [*[ressol], fitprim]
    else:
        return deepcopy(sol)


if __name__ == '__main__':
    _streval = 'from parser import create_table, get_globals'
    exec(_streval)
    with open('../mdgs/mini-tunis.mdg', 'r') as f:
        # with open('../mdgs/ispell.mdg', 'r') as f:
        data = f.readlines()
    ref = create_table(data)
    glbls = [*get_globals(ref)]
    # chrom = [[10, 6, 3, 14], [1, 2], [4, 5], [7, 11], [
    #     9, 8, 13, 15], [12, 20, 17, 18, 16, 19], [21, 23, 22, 24]]
    # chrom = [[18, 15, 8, 5], [2, 10, 16, 3], [1, 9, 20], [
    #     11, 7], [13, 14, 4, 17, 6, 19, 12], [21, 22, 23, 24]]

    # chrom = [[9, 11, 12, 13], [1, 3], [8, 10], [
    #     2, 5, 7, 8, 10], [6, 14, 16], [4, 15, 17, 18, 19, 20]]
    chrom = [[9, 7, 6], [1, 11, 3], [2, 8, 4, 5, 10],
             [13, 12], [14, 15, 16], [18, 17, 19, 20]]

    # chrom = [[4, 10], [9, 2, 3], [1, 5, 11, 6], [8, 7, 13], [12]]
    # chrom = [[9, 6, 4], [8, 11], [13, 7, 10, 2, 12], [1, 3, 5]]
    sol = encode(chrom, len(ref))
    testfit = fitness(sol, ref, 0.65)
    print([*[sol], testfit])
    sol2 = rep_solut([sol, testfit], glbls, ref, 0.65)
    print(sol2)
