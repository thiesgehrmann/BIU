import pandas as pd

from .. import utils

###############################################################################

def query(VCFstruct, chrom, start, end):
  """ Take as input a PYVCF Reader object"""

  try:
    return VCFstruct.fetch(str(chrom), int(start), int(end))
  except Exception as e:
    print(e)
    return []
  #etry
#edef

###############################################################################

def makeIdentifier(record, altPos=None):
  if altPos is not None:
    alts = [ record.ALT[altPos].sequence if hasattr(record.ALT[altPos], "sequence") else '-' ]
  else:
    alts = [ a.sequence if hasattr(a, "sequence") else "-" for a in record.ALT ]
  #fi
  return "%s-%d-%s-%s" % (record.CHROM, record.POS, record.REF, '/'.join(alts))
#edef

###############################################################################

def filter(arr, filters=None, gtFilters=None):

  # Filter out variants and samples
  # To filter out samples, we need to modify the _Call.samples and _sample_indexes variables.
  # In my version of pyvcf, the is_filtered function is not defined for some reason, so I add it here.

  def is_filtered(call):
    """ Return True for filtered calls.
        Taken from
    """
    try: # no FT annotation present for this variant
          filt = call.data.FT
    except AttributeError:
      return False
    #etry
    if filt is None or len(filt) == 0: # FT is not set or set to PASS
      return False
    else:
      return True
    #fi
  #edef

  if filters is not None:
    filters = set(filters)
  #fi

  if gtFilters is None:
    gtFilters = filters
  else:
    gtFilters = set(gtFilters)
  #fi

  res = arr
  if (filters is not None) or (gtFilters is not None):
    filtRes = []
    for r in res:
      if (filters is not None) and (r.FILTER is not None):
        if len(set(r.FILTER) & filters) > 0:
          continue # We have matched a filter, do not add this one.
        #fi
      #fi
      if gtFilters is not None:
        samples = r.samples
        filtSamples = []
        for s in samples:
          if is_filtered(s):
            sfilt = set(s.data.FT.split(';'))
            if len(sfilt & gtFilters) == 0:
              filtSamples.append(s)
            #fi
          else:
            filtSamples.append(s)
          #fi
        #efor
        r.samples = filtSamples
        r._sample_indexes = { s.sample : i for i,s in enumerate(r.samples) }
      #fi
      filtRes.append(r)
    #efor
    res = filtRes
  #fi
  return res
#edef

###############################################################################

def filterType(arr, types=None):
  res = arr
  if types is not None:
    types = set(types) & set(['snp', 'indel', 'sv'])
    # types must be set of snp|indel|sv
    typeRes = []
    for r in res:
      if r.var_type in types:
        typeRes.append(r)
      #fi
    #efor
    res = typeRes
  #fi
  return res
#edef

###############################################################################

def filterSubTypes(arr, subTypes=None):
  res = arr
  # types must be set of types defined in https://pyvcf.readthedocs.io/en/latest/_modules/vcf/model.html#_Record var_subtypes
  if subTypes is not None:
    subTypes = set(subTypes)
    subTypeRes = []
    for r in res:
      if r.var_subtype in subTypes:
        subTypeRes.append(r)
      #fi
    #efor
    res = subTypeRes
  #fi
  return arr
#edef

###############################################################################

def filterSamples(arr, sampleFilters=None):
  res = arr
  if sampleFilters is not None:
    sampleFilters = set(sampleFilters)
    sampleFiltRes = []
    for r in res:
      filtSamples = [ r.genotype(s) for s in sampleFilters if s in r._sample_indexes ]
      r.samples = filtSamples
      r._sample_indexes = { s.sample : i for i,s in enumerate(r.samples) }
      sampleFiltRes.append(r)
    #efor
    res = sampleFiltRes
  #efor
  return res
#edef

###############################################################################

def extract(arr, extract=None):
  extractions = {
    "raw" : lambda x: x,
    "?" : lambda x: utils.dbm("Unknown extraction option '%s'." % str(extract)),
    "summary" : lambda x: summary(x) }

  dExtract = "raw" if (extract is None) else (extract.lower() if isinstance(extract, str) else "?")

  return extractions[dExtract](arr)
#edef

###############################################################################

def summary(arr, altPos=None, refPos=None):
  # Basically copied from Erik's getGenoCountsFromMat function,
  # except that I WANT TO ALLOW FOR SUMMARIES FROM NON 0/1 genotype calls, but also for 1/2 0/2 etc...
  # NEED TO ADAPT THIS!

  def singleSummary(record, refp, altp):
    gtypes = [ s.data.GT if hasattr(s.data, "GT") else '-' for s in  record.samples ]
    #S = ( makeIdentifier(record, ),
    #      len([ x for x in gtypes if x in [ "0/0","0|0" ] ]),
    #      len([ x for x in gtypes if x in [ "0" ] ]),
    #      len([ x for x in gtypes if x in [ "0/1","1/0","1|0","0|1" ] ]),
    #      len([ x for x in gtypes if x in [ "1" ] ]),
    #      len([ x for x in gtypes if x in [ "1/1","1|1" ] ]),
    #      len([ x for x in gtypes if x not in [ "0","0/0","0|0", "1","0/1","0|1","1/0","1|0","1/1","1|1" ] ]) )
    S = ( makeIdentifier(record, altp-1),
          len([ x for x in gtypes if x in [ "%d/%d" % (refp, refp),"%d|%d" % (refp,refp) ] ]),
          len([ x for x in gtypes if x in [ "%d" % refp ] ]),
          len([ x for x in gtypes if x in [ "%d/%d" % (refp, altp),"%d/%d" % (altp, refp), "%d|%d" % (altp, refp),"%d|%d" % (refp, altp) ] ]),
          len([ x for x in gtypes if x in [ "%d" % (altp) ] ]),
          len([ x for x in gtypes if x in [ "%d/%d" % (altp, altp),"%d|%d" % (altp, altp) ] ]),
          len([ x for x in gtypes if x not in [ "%d" % refp, "%d/%d" % (refp, refp) ,"%d|%d" % (refp, refp), "%d" % altp, "%d/%d" % (refp, altp),"%d|%d" % (refp, altp),"%d/%d" % (altp, refp),"%d|%d" % (altp, refp), "%d/%d" % (altp, altp),"%d|%d" % (altp, altp) ] ]) )
    return S
  #edef

  return pd.DataFrame( [ singleSummary(v, 0 if refPos is None else refPos[i], 1 if altPos is None else altPos[i] ) for i,v in enumerate(arr)], 
                       columns=[ "id", "RR", "R", "RA", "A", "AA", "U"])
#edef

###############################################################################

def genotypeInfoFieldIndexes(altAlleleID, ref=0):
    """ For multi-allelic sites, genotype INFO fields have a particular order.
        If we want to extract information from these fields for a specific variant allele, we need to be able to reconstruct this order for an arbitrary number of alleles.
        The order is defined in the VCF specification, and implemented here.
        Search for 'GL : genotype likelihoods' in VCF spec.
        F(j/k) = (k*(k+1)/2)+j. I """
    
    F = lambda k,j : int(k*(k+1)/2)+j
    return [ F(ref, ref), F(ref, altAlleleID), F(altAlleleID,altAlleleID) ]
#edef

###############################################################################

#class Vcf(object):
#
#  def __init__(self, data *kwargs):
#
#  #edef
#
#  def filter(self, filters):
#
#  #edef
#
##eclass
