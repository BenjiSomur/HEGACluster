from operator import itemgetter
from parser import get_nodes, create_table
from initializer import *
from file_manager import write_headers, write_gen, write_csv
from operators import init_population, tournament, crossover, mutation, dominance_sort
from fitness import fitness
from decoder import decode
from graph_maker import create_graph
import sys
import os
import timeit


def main():
    data = None
    it_no = sys.argv[1]
    filename = sys.argv[2]
    _cop = sys.argv[3]
    _intx = sys.argv[4]
    with open('../mdgs/{}.mdg'.format(filename), 'r') as f:
        data = f.readlines()
    nodes = get_nodes(data)
    ref = create_table(data)
    pop_size = get_pop_size(nodes)
    cp = get_cp(nodes)
    mp = get_mp(nodes)
    _theta = get_theta(nodes)
    gens = get_no_gen(nodes)
    _args = {'filename': filename,
             'pop_size': pop_size,
             'cp': cp,
             'mp': mp,
             'it_no': it_no,
             'theta': _theta,
             'type': _cop,
             'intx': _intx
             }
    write_headers(_args)
    raw_pop = init_population(pop_size, len(nodes))
    _pop = []
    for chrom in raw_pop:
        ft = fitness(chrom, ref, _theta)
        _pop.append([chrom, ft])
    del raw_pop
    # _pop.sort(key=lambda x: (x[1][1], -x[1][0]))
    # best = _pop[-1]
    _pop = dominance_sort(_pop)
    best = _pop[0]
    print(best)
    no_gen = 0
    _kwargs = {'filename': filename,
               'it_no': it_no,
               'type': _cop,
               'best': best[0],
               'fitness': best[1],
               'nodes': len(ref),
               'gen': no_gen,
               'intx': _intx
               }

    write_gen(_kwargs)
    write_csv(_kwargs)
    while no_gen < gens:
        parents = list(tournament(_pop, pop_size))
        offsp = list()
        for i in range(0, len(parents), 2):
            gen_chld = crossover(
                [_pop[parents[i]][0], _pop[parents[i+1]][0]], len(nodes), cp, _cop, _intx)
            chld1 = mutation(gen_chld[0], mp)
            chld2 = mutation(gen_chld[1], mp)
            offsp += [chld1, chld2]
        aux = []
        for chrom in offsp:
            ft = fitness(chrom, ref, _theta)
            aux.append([chrom, ft])
        aux.append(best)
        del _pop
        _pop = aux
        del aux
        # _pop.sort(key=lambda x: (x[1][1], -x[1][0]))
        # while len(_pop) < pop_size:
        #     _pop.pop(0)
        # best = _pop[-1]
        _pop = dominance_sort(_pop)
        while len(_pop) < pop_size:
            _pop.pop(-1)
        best = _pop[0]
        print("{}: ".format(no_gen + 1), best)
        no_gen += 1
        _kwargs = {'filename': filename,
                   'it_no': it_no,
                   'type': _cop,
                   'best': best[0],
                   'fitness': best[1],
                   'gen': no_gen,
                   'nodes': len(ref),
                   'intx': _intx}
        write_gen(_kwargs)
        write_csv(_kwargs)

    _kwargs = {'filename': filename,
               'type': _cop,
               'it_no': it_no,
               'chrom': best[0],
               'nodes': nodes,
               'raw_data': data,
               'intx': _intx
               }
    create_graph(_kwargs)

    _path = f'./{filename}/{_cop}/{_intx}'
    if not os.path.exists('{}/bests.csv'.format(_path)):
        with open('{}/bests.csv'.format(_path), 'w') as f:
            f.write('It No., Best MQ, No clusters, max clus, min clus\n')
    with open('{}/bests.csv'.format(_path), 'a') as f:
        f.write('{},{},{},{},{}\n'.format(it_no,
                                          best[1][1],
                                          len(best[0][1]),
                                          max(best[0][1]),
                                          min(best[0][1])))


def extract_globals(ref):
    _globals = []
    for idx in range(len(ref)):
        uncalled = 0
        for line in ref:
            if line[idx] != 0:
                continue
            uncalled += 1
        if uncalled == len(ref):
            _globals.append(idx)
    return _globals


if __name__ == '__main__':
    # main()

    filename = input('Enter the filename: ')
    # _cop = input('Enter crossover method (cpx, onepx): ')
    # _intx = input('Enter the int crossover method (onepx, exchange): ')
    with open('../mdgs/{}.mdg'.format(filename), 'r') as f:
        data = f.readlines()
    nodes = get_nodes(data)
    ref = create_table(data)
    print(extract_globals(ref))

    # pop_size = get_pop_size(nodes)
    # cp = get_cp(nodes)
    # mp = get_mp(nodes)
    # _theta = get_theta(nodes)
    # gens = get_no_gen(nodes)
    # raw_pop = init_population(pop_size, len(nodes))
    # _pop = []
    # for chrom in raw_pop:
    #     ft = fitness(chrom, ref, _theta)
    #     _pop.append([chrom, ft])
    # _pop = dominance_sort(_pop)
    # print(_pop)

    # test_str = """
# raw_pop = init_population(pop_size, len(nodes))
# _pop = []
# for chrom in raw_pop:
#     ft = fitness(chrom, ref, _theta)
#     _pop.append([chrom, ft])
# _pop.sort(key=lambda x: x[1][0])
# best = _pop[-1]
# parents = [x for x in tournament(_pop, pop_size)]
# offsp = list()
# for i in range(0, len(parents), 2):
#     gen_chld = crossover([_pop[parents[i]][0], _pop[parents[i+1]][0]], len(nodes), cp, _cop, _intx)
#     chld1 = mutation(gen_chld[0],mp)
#     chld2 = mutation(gen_chld[1],mp)
#     offsp += [chld1,chld2]
# aux = []
# for chrom in offsp:
#     ft = fitness(chrom, ref, _theta)
#     aux.append([chrom, ft])
# aux.append(best)
# del _pop
# _pop = aux
# _pop.sort(key=lambda x: x[1][0])
# while len(_pop) < pop_size:
#     _pop.pop(0)
# """
#     _setup = '''from operators import crossover, mutation, init_population, tournament; from fitness import fitness'''
#     print(timeit.timeit(test_str, setup=_setup, number=1, globals={'pop_size': pop_size,
#                                                                    'nodes': nodes,
#                                                                    'cp': 1,
#                                                                    '_cop': _cop,
#                                                                    '_intx': _intx,
#                                                                    'mp': mp,
#                                                                    'ref': ref,
#                                                                    '_theta': _theta}))
