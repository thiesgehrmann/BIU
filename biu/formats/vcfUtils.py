import pandas as pd

import vcf
import intervaltree

from .. import utils

###############################################################################

class VCF(object):

  __reader = None
  __tabix = None
  __fileName = None
  __vcfArgs = None
  __readerIndex = None

  def __init__(self, data, tabix=False, **kwargs):
    self.__tabix = tabix
    self.__vcfArgs = kwargs

    if not(isinstance(data, str)):
      utils.dbm("VCF Input source is list of Records.")
      self.__reader = data
      self.__makeIndex(data)
    elif self.__tabix:
      utils.dbm("VCF Input source is tabixed file.")
      self.__fileName = data
      self.__reader = vcf.Reader(filename=self.__fileName, compressed=True, **kwargs)
    else:
      utils.dbm("VCF Input source is unindexed file.")
      self.__fileName =  data
      self.__reader = [ v for v in vcf.Reader(open(self.__fileName, 'r'), **self.__vcfArgs) ]
      self.__makeIndex(self.__reader)
    #fi
  #edef

  @property
  def samples(self):
    if self.__reader is None:
      if len(self.__readerIndex) == 0:
        return []
      #fi
      for item in self.__readerIndex[list(self.__readerIndex.keys())[0]]:
        return list(item.data._sample_indexes.keys())
      #efor
    else:
      return self.__reader.samples
    #fi
      return []
  #edef

  def __iter__(self):
    return self.__reader.__iter__()
  #edef

  def __next__(self):
    return self.__reader.__next__()
  #edef

  def __str__(self):
    dstr  = "VCF object\n"
    dstr += " Where: %s\n" % (self.__fileName if self.__fileName is not None else hex(id(self)))
    if self.__tabix is False:
      dstr += "Entries: %d\n" % sum( [len(self.__readerIndex[k]) for k in self.__readerIndex] )
    #fi
    dstr += "Number of genotypes: %d\n" % len(self.samples)
    dstr += "Tabix: %s\n" % ("Yes" if self.__tabix else "No")
    return dstr
  #edef

  def __makeIndex(self, entries):
    utils.dbm("Building VCF Index. May take a while.")

    self.__readerIndex = {}

    for record in entries:
      if record.CHROM not in self.__readerIndex:
        self.__readerIndex[record.CHROM] = intervaltree.IntervalTree()
      #fi
      self.__readerIndex[record.CHROM].addi(record.POS-1, record.POS, record)
    #efor
  #edef

  def __singleQuery(self, seqid, start, end):
    if self.__tabix:
      res = self.__tabixQuery(seqid, start, end)
    else:
      res = self.__nonTabixQuery(seqid, start, end)
    #fi
    return res
  #edef

  def query(self, seqid, start, end, **kwargs):
    return self.queryRegions( [ (seqid, start, end) ], **kwargs)
  #edef

  def queryRegions(self, regions,
                   filters=None,
                   gtFilters=None,
                   sampleFilters=None,
                   types=None,
                   subTypes=None,
                   nAlleles=None,
                   extract=None):
    """ This function is hacked together from internals of the pyVCF _Record and _Call classes. I hope that their definitions don't change anytime soon.
    """

    res = []
    for (seqid, start, end) in regions:
      for r in self.__singleQuery(seqid, start, end):
        res.append(r)
    #efor

    res = self.filterType(res, types=types)
    res = self.filterSubTypes(res, subTypes=subTypes)
    res = self.filter(res, filters=filters, gtFilters=gtFilters)
    res = self.filterSamples(res, sampleFilters=sampleFilters)
    res = self.filterAlleles(res, nAlleles=nAlleles)

    res = self.extract(res, extract=extract)

    return res
  #edef

  def __nonTabixQuery(self, chrom, start, end):
    chrom = str(chrom)
    if chrom not in self.__readerIndex:
      utils.dbm("Provided chromosome '%s' is not present in VCF." % chrom)
      return []
    else:
      return [ res.data for res in self.__readerIndex[chrom].search(start, end) ]
    #fi
  #edef
  
  def __tabixQuery(self, chrom, start, end):
    """ Take as input a PYVCF Reader object"""
  
    try:
      return self.__reader.fetch(str(chrom), int(start), int(end))
    except Exception as e:
      print(e)
      return []
    #etry
  #edef
  
  ###############################################################################

  @staticmethod
  def makeIdentifier(record, altPos=None):
    if altPos is not None:
      alts = [ record.ALT[altPos].sequence if hasattr(record.ALT[altPos], "sequence") else '-' ]
    else:
      alts = [ a.sequence if hasattr(a, "sequence") else "-" for a in record.ALT ]
    #fi
    return "%s-%d-%s-%s" % (record.CHROM, record.POS, record.REF, '/'.join(alts))
  #edef
  
  ###############################################################################

  @staticmethod  
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

  @staticmethod
  def filterAlleles(arr, nAlleles=None):
    res =  arr
    if nAlleles is not None:
      res = [ v for v in arr if len(v.ALT) <= nAlleles ]
    #fi
    return res
  #edef

  @staticmethod
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

  @staticmethod
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

  @staticmethod
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

  @staticmethod
  def extract(arr, extract=None):
    extractions = {
      "raw" : lambda x: x,
      "?" : lambda x: utils.dbm("Unknown extraction option '%s'." % str(extract)),
      "summary" : lambda x: VCF.summary(x) }
  
    dExtract = "raw" if (extract is None) else (extract.lower() if isinstance(extract, str) else "?")
  
    return extractions[dExtract](arr)
  #edef
  
  ###############################################################################

  @staticmethod
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
      S = ( VCF.makeIdentifier(record, altp-1),
            len([ x for x in gtypes if x in [ "%d/%d" % (refp, refp),"%d|%d" % (refp,refp) ] ]),
            len([ x for x in gtypes if x in [ "%d" % refp ] ]),
            len([ x for x in gtypes if x in [ "%d/%d" % (refp, altp),"%d/%d" % (altp, refp), "%d|%d" % (altp, refp),"%d|%d" % (refp, altp) ] ]),
            len([ x for x in gtypes if x in [ "%d" % (altp) ] ]),
            len([ x for x in gtypes if x in [ "%d/%d" % (altp, altp),"%d|%d" % (altp, altp) ] ]),
            len([ x for x in gtypes if x not in [ "%d" % refp, "%d/%d" % (refp, refp) ,"%d|%d" % (refp, refp), "%d" % altp, "%d/%d" % (refp, altp),"%d|%d" % (refp, altp),"%d/%d" % (altp, refp),"%d|%d" % (altp, refp), "%d/%d" % (altp, altp),"%d|%d" % (altp, altp) ] ]) )
      return S
    #edef
  
    return pd.DataFrame( [ singleSummary(v, 0 if refPos is None else refPos[i], 1 if altPos is None else altPos[i] ) for i,v in enumerate(arr)], 
                         columns=[ "id", "RR", "R", "RA", "A", "AA", "O"])
  #edef
  
  ###############################################################################

  @staticmethod
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

#eclass
