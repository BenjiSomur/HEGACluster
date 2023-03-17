from decoder import decode, encode
from math import ceil
from fitness import get_clus, turbomq, fitness
from copy import deepcopy
from parser import transpose_mat
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


def get_omnilocals(ref, chrom):
    for clus in chrom:
        for nomod in clus:
            omnirefs = related_to(ref, chrom, nomod)
            if len(omnirefs) > 2:
                yield nomod


def detomnlcls(ref, omnlocals):
    aux_ref = transpose_mat(ref)
    locls = list()
    for idx in omnlocals:
        locls.append((idx, sum(ref[idx - 1]) + sum(aux_ref[idx - 1])))
    aux_locls = sorted(locls, key=lambda x: x[1])
    return aux_locls


def extendomns(omns, lcls):
    aux = set(omns)
    for (auxid, _) in lcls:
        aux.add(auxid)
    return list(aux)


def rep_solut(sol, gbls, ref, theta):
    auxsol = list(sol)
    chrom = decode(auxsol[0], len(ref))
    gamma = ceil(0.2 * len(chrom)) + 1
    omns = [x for x in detomns(gbls, chrom, ref)]
    omnilocals = [y for y in get_omnilocals(ref, chrom)]
    highestomnilocals = detomnlcls(ref, omnilocals)
    omns = extendomns(omns, highestomnilocals[len(
        highestomnilocals) - (gamma + 1):])
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


def mutate_indiv(sol, gbls, ref, theta):
    auxfit = fitness(sol, ref, theta)
    chrom = decode(sol, len(ref))
    gamma = ceil(0.2 * len(chrom)) + 1
    omns = [x for x in detomns(gbls, chrom, ref)]
    omnilocals = [y for y in get_omnilocals(ref, chrom)]
    highestomnilocals = detomnlcls(ref, omnilocals)
    omns = extendomns(omns, highestomnilocals[len(
        highestomnilocals) - (gamma + 1):])
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
    return ressol


if __name__ == '__main__':
    _streval = 'from parser import create_table, get_globals'
    exec(_streval)
    with open('../mdgs/compiler.mdg', 'r') as f:
        # with open('../mdgs/icecast.mdg', 'r') as f:
        data = f.readlines()
    ref = create_table(data)
    glbls = [*get_globals(ref)]
    # chrom = [[10, 6, 3, 14], [1, 2], [4, 5], [7, 11], [
    #     9, 8, 13, 15], [12, 20, 17, 18, 16, 19], [21, 23, 22, 24]]
    # chrom = [[18, 15, 8, 5], [2, 10, 16, 3], [1, 9, 20], [
    #     11, 7], [13, 14, 4, 17, 6, 19, 12], [21, 22, 23, 24]]

    # chrom = [[9, 11, 12, 13], [1, 3], [8, 10], [
    #     2, 5, 7, 8, 10], [6, 14, 16], [4, 15, 17, 18, 19, 20]]
    # chrom = [[9, 7, 6], [1, 11, 3], [2, 8, 4, 5, 10],
    #          [13, 12], [14, 15, 16], [18, 17, 19, 20]]

    chrom = [[4, 10], [9, 2, 3], [1, 5, 11, 6], [8, 7, 13], [12]]
    # chrom = [[9, 6, 4], [8, 11], [13, 7, 10, 2, 12], [1, 3, 5]]
    # chrom = [[1, 40, 3, 35, 25, 2, 14, 27, 7], [4, 5, 6, 10, 8, 13, 21], [9, 15, 28, 12, 16, 23, 18, 24, 44], [20, 11, 19, 26, 31, 46, 47, 51],
    #          [22, 29, 30, 32, 37, 45, 54], [33, 38, 36, 39, 34, 50], [42, 17], [41, 43, 52], [57, 49], [53, 48, 55, 56, 58, 59], [60]]
    # chrom = [[1, 7, 40, 3], [35, 4, 11, 5, 6, # 10, 8, 13, 2, 14, 27], [21, 9, 31, 15, 28, 12, 16, 17, 23], [18, 51, 19, 20, 22],
    # #          [24, 26, 25, 29, 30, 32, 37, 33, 38, 36, 39, 34, 42, 47, 41, 46, 43, 45, 44, 52, 57, 49, 53, 50, 48, 55, 54], [56, 58, 59, 60]]
    # chrom = [[2, 14, 1, 27, 7, 40, 3, 35, 4], [11], [5, 6], [10, 8, 13, 21, 9, 31], [15, 28], [12], [16, 17, 23, 18, 51, 19, 20, 22],
    #          [24, 26, 25, 29, 30, 32, 37, 33, 38, 36, 39, 34, 42, 47, 41, 46, 43, 45, 44, 52, 57, 49, 53, 50, 48, 55, 54, 56, 58, 59, 60]]
    sol = encode(chrom, len(ref))

    # sol = ['000001101000000011100011011001000010001110011011000111000100000101000110001010001000001101010101001001001111011100001100010000010111010010011000101100010100001011010011011010011111101110101111110011010110011101011110100000100101101101110110100001100110100100100111100010110010101010010001101001101011110100111001110001110101110000110111111000111010111011111100', [9, 7, 9, 8, 7, 6, 2, 3, 2, 6, 1]]
    testfit = fitness(sol, ref, 0.65)
    print([*[sol], testfit])
    sol2 = rep_solut([sol, testfit], glbls, ref, 0.65)
    # sol2 = mutate_indiv(sol, glbls, ref, 0.65)
    # fit2 = fitness(sol2, ref, 0.65)
    # print([*[sol2], fit2])
    print(sol2)
