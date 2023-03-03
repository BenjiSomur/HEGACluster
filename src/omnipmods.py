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
        refs = related_to(ref, chrom, idx)
        if len(refs) >= 2:
            yield idx
