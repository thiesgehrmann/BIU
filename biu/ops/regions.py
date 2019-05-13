"""
A library of functions to deal with genomic regions
"""

from . import lst

def merge(regions):
    """
    Merge genomic regions based
    parameters:
    -----------
    regions: List[3-tuples]
        A list of regions in the form of (seqid, start, end)
    
    output:
    -------
    Returns a list of merged genomic regions
    """
    regions = lst.group(regions, key=lambda x: x[0]) # Group by seqid
    regions = { s : sorted(regions[s], key=lambda x: x[1]) for s in regions }
    merged = {}
    for seqid in regions:
        regs = regions[seqid]
        m = [ regs[0] ]
        for seqid, start, end in regs[1:]:
            if start <= m[-1][2]:
                m[-1] = (seqid, m[-1][1], end)
            else:
                m.append((seqid, start, end))
            #fi
        #efor
        merged[seqid] = m
    #efor
    return lst.flatten(merged.values())
#edef