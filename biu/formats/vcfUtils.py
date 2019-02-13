from .. import utils

pd           = utils.py.loadExternalModule("pandas")
np           = utils.py.loadExternalModule("numpy")
vcf          = utils.py.loadExternalModule("vcf")
intervaltree = utils.py.loadExternalModule("intervaltree")

###############################################################################

class VCF(object):

  __slots__ = [ '__reader', '__tabix', '__fileName', '__template', '__vcfArgs', '__readerIndex', '__records' ]

  #__reader = None
  #__tabix = None
  #__fileName = None
  #__template = None
  #__vcfArgs = None
  #__readerIndex = None

  def __init__(self, data, template=None, tabix=False, **kwargs):
    self.__tabix = tabix
    self.__vcfArgs = kwargs
    self.__readerIndex = None
    self.__tabix = tabix
    self.__fileName = None
    self.__template = None
    self.__records  = None

    if not(isinstance(data, str)):
      utils.dbm("VCF Input source is list of Records.")
      self.__reader = list(data)
      self.__template = template
      self.__fileName = None
    elif self.__tabix:
      utils.dbm("VCF Input source is tabixed file.")
      self.__fileName = data
      self.__reader = vcf.Reader(filename=self.__fileName, compressed=True, **kwargs)
      self.__template = self.__fileName
    else:
      utils.dbm("VCF Input source is unindexed file.")
      self.__fileName =  data
      self.__reader = list([ v for v in vcf.Reader(open(self.__fileName, 'r'), **self.__vcfArgs) ])
      self.__template = self.__fileName
    #fi
  #edef

  def write(self, outFile, template=None):
    """
    write: Write the VCF structure to a VCF file
    Inputs: outFile: The location to which it should it be written
            template: VCF template to write with (needed for meta data). Can be specified here if the VCF file doesn't already have one
    """

    # NOTE: This may not work if you have used gtFilters or sampleFilters
    template = self.__template if template is None else template
    if template is None:
      utils.error("There is no template defined for this VCF file.")
      return False
    #fi

    with open(outFile, 'wb') as ofd:
      writer = vcf.Writer(ofd, template)
      for record in self.__reader:
        writer.write_record(record)
      #efor
    #ewith
    return True
  #edef

  @property
  def filename(self):
    """
    filename: The filename at which the VCF object exists (if at all)
    """
    return self.__fileName

  @property
  def samples(self):
    """
    samples: A list of samples for which genotypes are defined
    """
    if isinstance(self.__reader, list):
      if len(self.__getReaderIndex()) == 0:
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

  def __getitem__(self, k):
    if isinstance(self.__reader, list):
      return self.__reader[k]
    else:
      utils.warning("Tabix Indexed VCF files cannot be indexed numerically.")
      return []
    #fi
  #edef

  def __iter__(self):
    return self.__reader.__iter__()
  #edef

  def __next__(self):
    return self.__reader.__next__()
  #edef

  def __len__(self):
    return len(self.__reader)
  #edef

  @property
  def records(self):
    """
    records: Return a list of records
    This may take some time if the file is large...
    """
    if self.__records is None:
      self.__records = [ r for r in self.__reader ]
    #edef
    return self.__records
  #edef

  @property
  def template(self):
    """
    template: The VCF template
    """
    return self.__template
  #edef

  def __str__(self):
    dstr  = "VCF object\n"
    dstr += " Where: %s\n" % (self.__fileName if self.__fileName is not None else hex(id(self)))
    dstr += " Template: %s\n" % (self.__template if self.__template is not None else "No Template!")
    if self.__tabix is False:
      dstr += " Entries: %d\n" % sum( [len(self.__getReaderIndex(k)) for k in self.__getReaderIndex()] )
    #fi
    dstr += " Number of Samples: %d\n" % len(self.samples)
    dstr += " Tabix: %s\n" % ("Yes" if self.__tabix else "No")
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

  def __getReaderIndex(self, k=None):
    if self.__readerIndex is None:
      self.__makeIndex(self.records)
    #fi

    if k is None:
      return self.__readerIndex
    elif k not in self.__readerIndex:
      return None
    else:
      return self.__readerIndex[k]
    #fi
  #edef

  def __singleQuery(self, seqid, start, end):
    if self.__tabix:
      res = self.__tabixQuery(seqid, start, end)
    else:
      res = self.__nonTabixQuery(seqid, start, end)
    #fi
    return res
  #edef

  def whoHas(self, chromosome, pos, alt, ref=0):
    """
    whoHas: Return a list of samples in which the variant is present
    Inputs: chromosome: The Chromosome ID
            pos: The nucleotide position
            alt: The alternative allele (A/C/G/T)
            ref: The reverence allele (default 0)
    Output: List of sample IDs
    """
    var = self.getVar(chromosome, pos, alt)
    if var is None:
      return []
    #fi
    var, altPos = var
    altGTs = [ '%d/%d' % (altPos, ref), '%d/%d' % (ref, altPos), '%d|%d' % (altPos, ref), '%d|%d' % (ref, altPos) ]
    theyHave = []
    for sample in var.samples:
      if not(hasattr(sample.data, 'GT')):
        continue
      #fi
      if sample.data.GT in altGTs:
        theyHave.append(sample.sample)
      #fi
    #efor
    return theyHave
  #edef

  def getVar(self, chromosome, pos, alt, **kwargs):
    """
    getVar: Get variant at a specific location
    Inputs: chromosome: The Chromosome ID
            pos: The nucleotide position
            alt: The alternative allele (A/C/G/T)
    Outputs: Return (VCF._Record, altPos) if the alternative variant exists, otherwise None
             AltPos is the index of the alternative allele (in the case of a multi-allelic site)
    """
    V = self.query(chromosome, int(pos)-1, pos, extract='raw', **kwargs)
    for v in V:
      try:
        altPos = [ a.sequence for a in v.ALT ].index(alt)
        return (v, altPos+1)
      except ValueError:
        return None
      #etry
    #efor
    return None
  #edef

  def query(self, seqid, start, end, **kwargs):
    """
    query: Query the VCF file in a given region
           If the file is tabix-indexed, then the tabix index is used. Otherwise, an index is constructed using intervaltrees.
    Inputs: seqid: The chromosome ID
            start: The start of the region of the query
            end:   The end of the region of the query
            **kwargs: See queryRegions for these options (filters, gtFilters, sampleFilters, types, subTypes, nAlleles, extract)
    Output: depends on value of extract. See queryRegions

    This function is a wrapper for queryRegions, where query regions is called with parameters ( [ (seqid, start, end) ], **kwargs)
    """
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
    """
    queryRegions: Query multiple regions at once
    Inputs: regions:       A list of regions of the form [ (seqid, start, end) ... ]
            filters:       [list of str ]A list of filters to apply on the FILTER column of the VCF record (Removes all records that match the filter)
            gtFilters:     [list of str]A list of filters to apply to genotype calls
            sampleFilters: [list of str ]Only keep these samples
            types:         [list of str] Only keep these types of variants (see pyvcf documentation for types)
            subtypes:      [list of str] Only keep these subtypes of variants (see pyvcf docs)
            nAlleles:      [int] Maximum number of alleles
            extract:       None -> Return VCF structure
                           "raw" -> Return a list of _Record objects
                           "summary" -> Return a summary
    Outputs: See extract option

    This function is hacked together from internals of the pyVCF _Record and _Call classes. I hope that their definitions don't change anytime soon.
    """

    res = []
    for (seqid, start, end) in regions:
      for r in self.__singleQuery(seqid, start, end):
        res.append(r)
    #efor

    res = self.filterType(res, types=types)
    res = self.filterSubTypes(res, subTypes=subTypes)
    res = self.filterFilter(res, filters=filters, gtFilters=gtFilters)
    res = self.filterAlleles(res, nAlleles=nAlleles)
    res = self.filterSamples(res, sampleFilters=sampleFilters)

    res = self.extract(res, extract=extract)

    if extract is None:
      res = VCF(data=res, template=self.__template)
    #fi

    return res
  #edef

  def filter(self, filters=None,
                   gtFilters=None,
                   sampleFilters=None,
                   types=None,
                   subTypes=None,
                   nAlleles=None):
    """
    filter: See options to queryRegions
    """

    if not(isinstance(self.__reader, list)):
      utils.error("Cannot filter tabix indexed VCF files. First perform a query on it, and then you can filter it, or load without the tabix file")
      return None
    #fi

    res = self.__reader
    res = self.filterType(res, types=types)
    res = self.filterSubTypes(res, subTypes=subTypes)
    res = self.filterSamples(res, sampleFilters=sampleFilters)
    res = self.filterFilter(res, filters=filters, gtFilters=gtFilters)
    res = self.filterAlleles(res, nAlleles=nAlleles)

    return VCF(data=res, template=self.__template)
  #edef

  #edef

  def __nonTabixQuery(self, chrom, start, end):
    chrom = str(chrom)

    # Build the index if we haven't already - Delay building it if we don't need to!
    idx = self.__getReaderIndex(chrom)
    if idx is None:
      utils.dbm("Provided chromosome '%s' is not present in VCF." % chrom)
      return []
    else:
      return [ res.data for res in idx.search(start, end) ]
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
    """
    makeIdentifier: Make a unique identifier for variant.
    Inputs: record : [VCF._Record] object (e.g. from getVar, or query)
            altPos: [int] In the case of a multiallelic locus, select the a specific other alternative (altPos > 1)
    """
    if altPos is not None:
      alts = [ record.ALT[altPos].sequence if hasattr(record.ALT[altPos-1], "sequence") else '-' ]
    else:
      alts = [ a.sequence if hasattr(a, "sequence") else "-" for a in record.ALT ]
    #fi
    return "%s-%d-%s-%s" % (record.CHROM, record.POS, record.REF, '/'.join(alts))
  #edef
  
  ###############################################################################

  @staticmethod  
  def filterFilter(arr, filters=None, gtFilters=None):
  
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

          # If there are no VARIANT genotype calls that pass this filter, then we should remove the variant.
          if (len([ s for s in filtSamples if s.is_variant ]) == 0) and (len(samples) > 0):
            continue
          #fi
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
          F(j/k) = (k*(k+1)/2)+j. I
      """
      
      F = lambda k,j : int(k*(k+1)/2)+j
      return [ F(ref, ref), F(ref, altAlleleID), F(altAlleleID,altAlleleID) ]
  #edef

  ###############################################################################

  def genotype_matrix(self):
    """
    Create a genotype matrix from the records in the VCF

    Returns:
        G: a DataFrame with columns as samples, and rows as variants
        Cells are encoded as 0/1/2 for homozygous ref/heterozygous/homozygous variant.

    NOTE: ASSUMES REFERENCE is reference.
    If you wish to use a different variant as reference, this function will not work!

    NOTE: catches some weird errors in the VCF structure, whereby a person might be missing.
    The internal gm object makes sure that this doesnt interfere with the results.
    """
    
    class gm(object):
        def __init__(self, samples, n_variants):
            self.values = None
            self.n_variants = n_variants            
            self.sample_idx = { s : i for (i,s) in enumerate(samples) }
            
            self.gt = np.zeros((self.n_samples, self.n_variants))
        #edef
        
        @ property
        def n_samples(self):
            return len(self.sample_idx)
        #edef
        
        @property
        def shape(self):
            return self.gt.shape
        #edef
        
        @property
        def samples(self):
            return list(self.sample_idx.keys())
        
        def __setitem__(self, key, value):
            sample_id, variant_id = key
            if sample_id not in self.sample_idx:
                self.sample_idx[sample_id] = len(self.sample_idx)
                self.gt.resize((self.n_samples, self.n_variants))
            #fi
            
            if variant_id >  self.n_variants:
                self.gt.resize((self.n_samples, variant_id+1))
                self.n_variants = variant_id
            #fi
            
            sample_key = self.sample_idx[sample_id]
            self.gt[sample_key, variant_id] = value
        #edef
    #eclass

    G = gm(self.samples, len(self.records))
    for i, record in enumerate(self.records):
        for sample in record.samples:
            gt = sample.data.GT.replace('|', ' ').replace('/',' ').split(' ')
            if len(gt) == 1:
                gt = 1 if gt[0] == '1' else 0
            elif gt[0] != gt[1]:
                gt = 1
            elif (gt[0] == '0') or (gt[1] == '0'):
                gt = 0
            else:
                gt = 2
            #fi
            G[sample.sample,i] = gt
        #efor
    #efor

    return pd.DataFrame(G.gt, index=G.samples, columns=[self.makeIdentifier(r) for r in self.records])
  #edef

###############################################################################

#eclass
