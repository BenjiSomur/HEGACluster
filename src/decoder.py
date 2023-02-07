def decode(chrom, nodes):
    max_size = len(chrom[0])
    max_bits = len("{0:b}".format(nodes))
    div_aux = [chrom[0][i:i+max_bits] for i in range(0, max_size, max_bits)]
    # int_aux = {int(x, 2) for x in div_aux}
    int_aux = list(map(lambda x: int(x, 2), div_aux))
    decoded = []
    for size in chrom[1]:
        aux = [int_aux.pop(0) for _ in range(size)]
        decoded.append(aux)
    return decoded
