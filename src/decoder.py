from itertools import chain, islice
from struct import unpack
import timeit


# def encode(chrom, nodes):
#     max_bits = len("{0:b}".format(nodes))
#     ratio = "0{}b".format(max_bits)
#     indiv = list()
#     intpart = [len(x) for x in chrom]
#     binflat = list(chain(*chrom))
#     aux = [format(x, ratio) for x in binflat]
#     indiv.append(''.join(aux))
#     indiv.append(intpart)
#     return indiv

def encode(chrom, nodes):
    max_bits = len("{0:b}".format(nodes))
    ratio = "0{}b".format(max_bits)
    intpart = [len(x) for x in chrom]
    binflat = [bit for sublist in chrom for bit in sublist]
    aux = [format(x, ratio) for x in binflat]
    indiv = [''.join(aux), intpart]
    return indiv


# def decode(chrom, nodes):
#     max_size = len(chrom[0])
#     max_bits = len("{0:b}".format(nodes))
#     div_aux = [chrom[0][i:i+max_bits] for i in range(0, max_size, max_bits)]
#     # int_aux = {int(x, 2) for x in div_aux}
#     int_aux = list(map(lambda x: int(x, 2), div_aux))
#     decoded = []
#     for size in chrom[1]:
#         aux = [int_aux.pop(0) for _ in range(size)]
#         decoded.append(aux)
#     return decoded


def decode(chrom, nodes):
    max_size = len(chrom[0])
    max_bits = len("{0:b}".format(nodes))
    div_aux = [chrom[0][i:i+max_bits] for i in range(0, max_size, max_bits)]
    int_aux = [int(x, 2) for x in div_aux]
    decoded = []
    for size in chrom[1]:
        aux = [int_aux.pop(0) for _ in range(size)]
        decoded.append(aux)
    return decoded


if __name__ == '__main__':
    _streval = 'from parser import create_table, get_globals'
    exec(_streval)
    with open('../mdgs/compiler.mdg', 'r') as f:
        data = f.readlines()
    ref = create_table(data)
    sol = [[4, 10], [9, 2, 3], [1, 5, 11, 6], [8, 7, 13], [12]]
    chrom = encode(sol, len(ref))
    nodes = len(ref)
    print(chrom)

    print(timeit.timeit(lambda: decode(chrom, nodes),
          globals=globals(), number=1000000))
    print(timeit.timeit(lambda: decode2(chrom, nodes), number=1000000))
    # print(decode(chrom, nodes))
    # print(decode2(chrom, nodes, max_size, max_bits))
