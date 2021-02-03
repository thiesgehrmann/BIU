import os

from .. import utils

cyvcf2 = utils.py.loadExternalModule('cyvcf2')
np     = utils.py.loadExternalModule('numpy')
pd     = utils.py.loadExternalModule('pandas')

class MyOwnVariant(object):
    def __init__(self, data):
        
        if not isinstance(data, str):
            self._str, self._CHROM, self._POS, self._ID, self._REF, self._ALT, self._QUAL, self._FILTER, self._INFO, self._FORMAT, self._samples = data
        else:
        
            fields = data.strip().split('\t')
            self._str    = data
            self._CHROM  = fields[0]
            self._POS    = int(fields[1])
            self._ID     = fields[2]
            self._REF    = fields[3]
            self._ALT    = fields[4].split(',')
            self._QUAL   = fields[5]
            self._FILTER = fields[6].split(',')

            self._FILTER = None if self._FILTER[0].lower() in ['.', 'pass'] else self._FILTER

            self._INFO   = dict([ x.split('=') for x in fields[7].split(';') ])
            if len(fields) > 8:
                self._FORMAT  = fields[8].split(':')
                self._samples = dict(zip(self._FORMAT, zip(*[p.split(':') for p in fields[9:] ])))
            else:
                self._FORMAT  = []
                self._samples = {}
            #fi
        #fi
    #edef
    
    def copy(self):
        return MyOwnVariant([ x.copy() for x in [ self._str, self._CHROM, self._POS, self._ID, self._REF, self._ALT, self._QUAL, self._FILTER, self._INFO, self._FORMAT, self._samples] ])
    #edef

    @property
    def INFO(self):
        return self._INFO
    #edef
    
    @property
    def FILTER(self):
        return self._FILTER
    #edef
    
    @property
    def FORMAT(self):
        return self._FORMAT
    #edef
    
    @property
    def CHROM(self):
        return self._CHROM
    #edef
    
    @property
    def POS(self):
        return self._POS
    
    @property
    def REF(self):
        return self._REF
    #edef
    
    @property
    def ALT(self):
        return self._ALT
    #edef
    
    @property
    def ID(self):
        return self._ID
    #edef
    
    @property
    def QUAL(self):
        return self._QUAL
    #edef
    
    @property
    def genotypes(self):
        gt = self.format('GT')
        ret = []
        
        def fmt(h):
            if h == '.':
                return -1
            else:
                return int(h)
            #fi
        #edef
        
        for g in gt:
            g = g[0]
            phased = '|' in g
            ret.append([ fmt(h) for h in g.split('|' if phased else '/')] + [phased])
        #efor
        return ret
    #edef
    
    @property
    def ploidy(self):
        gt = self.format('GT')
        ret = []
        
        for g in gt:
            g = g[0]
            phased = '|' in g
            ret.append(len(g.split('|' if phased else '/')))
        #efor
        return ret
    #edef
    
    @property
    def gt_phases(self):
        return [ '|' in gt[0] for gt in self.format('GT') ]
    #edef
    
    def format(self, field, vtype=None):
        if field in self._FORMAT:
            return ( v.split(';') for v in self._samples[field] )
        else:
            raise Exception("FUCK")
        #fi
    #edef
    
    def set_format(self, field, value):
        """
        PARAMETERS:
        field: String
        value: list[String|Int|object]
        """
        nsamples = len(self._samples[self._FORMAT[0]])
        if len(value) != nsamples:
            raise Exception("NOT THE RIGHT # of SAMPLES")
        else:
            value = [
                [v] if isinstance(v,str) or (not hasattr(v, '__iter__')) else v for v in value
            ]
            self._samples[field] = [ ';'.join([str(v) for v in val]) for val in value ]
        #fi
    #edef
    
    @property
    def is_snp(self):
        "boolean indicating if the variant is a SNP."
        if len(self.REF) > 1: return False
        for alt in self.ALT:
            if alt not in ["A", "C", "G", "T"]:
                return False
            #fi
        #efor
        return len(self.ALT) >= 1
    #edef

    @property
    def is_indel(self):
        "boolean indicating if the variant is an indel."
        is_sv = self.is_sv

        if len(self.REF) > 1 and not is_sv: return True

        for alt in self.ALT:
            if alt == ".":
                return True
            #fi
            if len(alt) != len(self.REF):
                if not is_sv:
                    return True
                #fi
            #fi
        return False

    @property
    def is_transition(self):
        "boolean indicating if the variant is a transition."
        if len(self.ALT) != 1:
            return False
        #fi

        if not self.is_snp:
            return False
        #fi

        ref = self.REF
            # just one alt allele
        alt_allele = self.ALT[0]
        if ((ref == 'A' and alt_allele == 'G') or
            (ref == 'G' and alt_allele == 'A') or
            (ref == 'C' and alt_allele == 'T') or
            (ref == 'T' and alt_allele == 'C')):
                return True
        #fi
        return False
    #edef

    @property
    def is_deletion(self):
        "boolean indicating if the variant is a deletion."
        if len(self.ALT) > 1:
            return False
        #fi

        if not self.is_indel:
            return False
        #fi
        if len(self.ALT) == 0:
            return True
        alt = self.ALT[0]
        if alt is None or alt == ".":
            return True
        #fi
        
        if len(self.REF) > len(alt):
            return True
        #fi
        return False
        #fi
    #edef

    @property
    def is_sv(self):
        "boolean indicating if the variant is an SV."
        return self.INFO.get('SVTYPE',None) is not None
    #edef


    @property
    def var_type(self):
        "type of variant (snp/indel/sv)"
        if self.is_snp:
            return "snp"
        elif self.is_indel:
            return "indel"
        elif self.is_sv:
            return "sv"
        else:
            return "unknown"
        #fi
    #edef

    @property
    def var_subtype(self):
        if self.is_snp:
            if self.is_transition:
                return "ts"
            #fi
            if len(self.ALT) == 1:
                return "tv"
            #fi
            return "unknown"

        elif self.is_indel:
            if self.is_deletion:
                return "del"
            if len(self.ALT) == 1:
                return "ins"
            else:
                return "unknown"

        svt = self.INFO.get("SVTYPE")
        if svt is None:
            return "unknown"
        if svt == "BND":
            return "complex"
        if self.INFO.get('IMPRECISE') is None:
            return svt
        return self.ALT[0].strip('<>')

    @property
    def start(self):
        "0-based start of the variant."
        return self.POS
    #edef

    @property
    def end(self):
        "end of the variant. the INFO field is parsed for SVs."
        return self.POS + len(self.ALT[0])
    #edef
    
    def __str__(self):
        info = ';'.join('%s=%s' % (k,v) for (k,v) in self._INFO.items())
        people = '\t'.join([ ':'.join(x) for x in zip(*[self._samples[k] for k in self._FORMAT])])
        return '\t'.join([self.CHROM, str(self.POS), self._ID, self._REF, ','.join(self.ALT),
                         self._QUAL, '.' if self._FILTER is None else ','.join(self._FILTER),
                         info, ':'.join(self._FORMAT), people])
    #edef


    #aaf
    #
    #    alternate allele frequency across samples in this VCF.
    #
    #call_rate
    #
    #    proprtion of samples that were not UKNOWN.
    #
    #
    #gt_alt_depths
    #
    #    get the count of alternate reads as a numpy array.
    #
    #gt_alt_freqs
    #
    #    get the freq of alternate reads as a numpy array.
    #
    #gt_bases
    #
    #    numpy array indicating the alleles in each sample.
    #
    #gt_depths
    #
    #    get the read-depth for each sample as a numpy array.
    #
    #gt_phases
    #
    #    get a boolean indicating wether each sample is phased as a numpy array.
    #
    #gt_phred_ll_het
    #
    #    get the PL of het for each sample as a numpy array.
    #
    #gt_phred_ll_homalt
    #
    #    get the PL of hom_alt for each sample as a numpy array.
    #
    #gt_phred_ll_homref
    #
    #    get the PL of Hom ref for each sample as a numpy array.
    #
    #gt_quals
    #
    #    get the GQ for each sample as a numpy array.
    #
    #gt_ref_depths
    #
    #    get the count of reference reads as a numpy array.
    #
    #gt_types
    #
    #    gt_types returns a numpy array indicating the type of each sample.
    #
    #    HOM_REF=0, HET=1. For gts012=True HOM_ALT=2, UKNOWN=3
    #
    #
    #num_called
    #
    #    number of samples that were not UKNOWN.
    #
    #num_het
    #
    #    number heterozygous samples at this variant.
    #
    #num_hom_alt
    #
    #    number homozygous alternate samples at this variant.
    #
    #num_hom_ref
    #
    #    number homozygous reference samples at this variant.
    #
    #num_unknown
    #
    #    number unknown samples at this variant.
    #
    #ploidy
    #
    #    get the ploidy of each sample for the given record.
    #
    #
    #set_pos(self, int pos0)
    #
    #    set the POS to the given 0-based position

#eclass

###################################################################

class BiuVariant(object):
    
    _var = None
    
    def __init__(self, data):
        if isinstance(data, cyvcf2.cyvcf2.Variant) or isinstance(data, MyOwnVariant):
            self._var = data
        elif isinstance(data, BiuVariant):
            self._var = data._var
        elif isinstance(data, str):
            self._var = MyOwnVariant(data)
        else:
            print('Unknown variant type: %s' % type(data))
        #fi
    #edef
    
    def __getattr__(self, name):
        if name == 'set_format':
            self.switch()
        #fi
        
        return getattr(self._var, name)
    #edefp
    
    def switch(self):
        """Switch the internal representation to the slower, more flexible one"""
        if not isinstance(self._var, MyOwnVariant):
            self._var = MyOwnVariant(str(self._var).strip())
        #fi
        return self
    #edef
    
    def copy(self):
        return self.__class__(self._var)
    #edef
    
    def __dir__(self):
        """
        To allow tab-completion for the registered objects.
        """
        return object.__dir__(self) + object.__dir__(self._var)
    #edef
    
    def __repr__(self):
        return "%s(%s:%d %s/%s)" % (self.__class__.__name__, self._var.CHROM,
                                    self._var.POS, self._var.REF, ",".join(self._var.ALT))
    #edef
    
    def __str__(self):
        return str(self._var)
    #edef
    
    def make_identifier(self, alt_pos=0):
        """
        Make an identifier for a given variant.
        
        parameters:
        self: a BiuVariant object.
        alt_pos: Integer. Which alternative allele to use in the identifier
        
        returns:
        A string identifier for the variant
        """
        return '%s-%d-%s-%s' % (str(self.CHROM), self.POS, self.REF, self.ALT[alt_pos])
    #edef
    
    def summary(self, altpos=None):
        """
        Make a summary of the variants at this variant
        ASSUMES A MONOPLOID/DIPLOID ORGANISM!
        parameters:
        self: BiuVariant Variant object
        altpos: Integer, List[Integer]
            Which alternative allele indexes to consider
            
        Returns: DataFrame
        """
        
        counts  = []
        idxs = []
        
        if altpos is None:
            altpos = list(range(len(self.ALT)))
        elif isinstance(altpos, int):
            altpos = [altpos]
        #fi
        
        for i in altpos:
            alt = self.ALT[i]
            ident = self.make_identifier(i)
            pos = i+1

            #                   RR, RA, AA, R, A, O
            v_count = np.array([ 0,  0,  0, 0, 0, 0])
            n_genotypes = len(self.genotypes)

            for gt in self.genotypes:
                if len(gt) == 2: # Monoploid!
                    v_count[3+(gt[0]==pos)] += 1
                elif len(gt) == 3: # Diploid
                    v_count[(gt[0]==pos)+(gt[1]==pos)] += 1
                else: #Polyploid! Don't count this sample!
                    pass
                #fi
            #efor

            # Set the Other count
            v_count[-1] = n_genotypes - sum(v_count)

            counts.append(v_count)
            idxs.append(ident)
        #efor

        return pd.DataFrame(counts, index=idxs, columns=['RR','RA','AA','R','A', 'O'])
    #edef
#eclass

###################################################################

class VCF_filter(object):
    """
    A structure to describe a VCF filtering step.
    Used internally.
    """
    def __init__(self, name, params):
        self._name   = name
        self._params = params
    #edef
    
    @property
    def name(self):
        return self._name
    #edef
    
    @property
    def params(self):
        return self._params
    #edef
    
    def __str__(self):
        """
        String representation of the filter
        """
        return '%s(%s)' % (self._name, ', '.join([str(p) for p in self._params]))
    #edef
#eclass

###################################################################

class VCF2(object):
    """
    A cyvcf2 wrapper for general querying of tabixed, and non-tabixed files.
    
    Internally, the data is represented in one of three ways:
     * VCF_records: Essentially a wrapper around a list of VCF records.
     * VCF_cyvcf2:  A wrapper around the usual cyvcf2 interface
     * VCF_tabix_cyvcf2: A wrapper which provides a tabix query-able interface to the cyvcf2 object.
     
    There are five different filters, and the internal representation changes depending on the filter.
    These transitions are explained below.
     * filter_samples   (Only retrieve records for a specific set of samples)
        VCF_cyvcf2       -> VCF_cyvcf2
        VCF_records      -> Raises Exception (You can't select samples at this stage)

     * filter_region    (Retrieve records for a specific region)
        VCF_cyvcf2       -> VCF_records
        VCF_records      -> VCF_records
        
     * filter_filter    (Retrieve records that do not pass certain filters)
        VCF_cyvcf2       -> VCF_records
        VCF_records      -> VCF_records
        
     * filter_n_alleles (Retrieve only records that have at most a certain number of alternative alleles)
        VCF_cyvcf2       -> VCF_records
        VCF_records      -> VCF_records

     * filter_vartypes  (Retrieve only records of a specific type)
        VCF_cyvcf2       -> VCF_records
        VCF_records      -> VCF_records
     
     Typical usage:
     
     x = VCF2('my.vcf.tgz')               # Internal representation: VCF_cyvcf2
     x = x.filter(12, 50012321, 50042321) # Internal representation: VCF_records
     x = x.filter(samples=[a,b,c,d])      # Internal representation: VCF_records
     
     x = VCF2('my.vcf.tgz')               # Internal representation: VCF_cyvcf2
     x = x.filter(samples=[a,b,c,d])      # Internal representation: VCF_cyvcf2
     x = x.filter(12, 50012321, 50042321) # Internal representation: VCF_records
     x = x.filter(12, 50012321, 50032321) # Internal representation: VCF_records
    
    """
    
    HOM_REF = 0
    HET     = 1
    HOM_ALT = 2
    UNKNOWN = 3
    
    def __init__(self, file=None, tabix=True, samples=None,
                 filter_stack=None, vcf_object=None,
                 internal_variant_representation=BiuVariant):
        """
        Initialize a VCF2 object.
        
        Typical usage:
        
        x = VCF2("my.vcf")
        
        By default, we will try to use a tabix index if it exists.
        
        If you have tabix-indexed the file, but do not wish to make use of this, then you can do the following:
        
        x = VCF2("my.vcf", tabix=False)
        
        parameters:
        ------------
        file: String. Path to the file you want to use
        tabix: Boolean. Use a tabix index, if one exists?
        samples: list of strings. Which samples to use?

        filter_stack: internal use
        vcf_object:   internal use representation of the VCF data.
        """
        
        if isinstance(vcf_object, (VCF2_records, VCF2_cyvcf2) ):
            self._vcf = vcf_object
        elif file is not None:
            #if tabix & os.path.exists(file + '.tbi'):
            #    self._vcf = VCF2_cyvcf2_tabix(file, samples=samples)
            #else:
            self._vcf = VCF2_cyvcf2(file, samples=samples, internal_variant_representation=internal_variant_representation)
            #fi
        else:
            raise ValueError("Incorrectly specified initialization.")
        #fi
        
        self._internal_variant_representation = internal_variant_representation
        self._filter_stack = [] if filter_stack is None else filter_stack
    #edef

    @property
    def records(self):
        """
        A list of Variant records that are described in this VCF object.
        """
        self.to_records()
        
        return list(self._vcf.records)
    #edef
    
    @property
    def samples(self):
        """
        A list of the samples that are described in this VCF object
        """
        return self._vcf.samples
    #edef
    
    @property
    def filter_stack(self):
        """
        A list of filters that have been applied to result in this VCF2 object.
        """
        return self._filter_stack
    #edef
        
    
    def filter(self, chrom=None, start=None, end=None, samples=None, filters=None,
               vartypes=None, n_alleles=None, samples_format=None):
        """
        Perform a filtering operation.
        You can filter on:
        
        region:            vcf.filter(chrom, start, end)
        samples:           vcf.filter(samples=[a,b,c,d])
        filters:           vcf.filter(filters=[a,b,c])
        variant types:     vcf.filter(vartypes=[a,b,c])
            vartypes must be in ['snp', 'indel', 'cv' ]
        # variant alleles: vcf.filter(n_alleles=2)
        
        Filter samples based on a format filtering operation
        vcf.filter(samples_format=(['FT', lambda x: len(set(x.split('|')) & set(['VQLOW'])) == 0]))
        
        You can chain filtering steps as you wish:
            vcf.filter(chrom, start, end).filter(n_alleles=2)
          
        Or you can perform them all in one go (preferred):
            vcf.filter(chrom, start, end, n_alleles=2)
        
        Returns VCF2 type object.
        """
        filter_list = []
        
        if samples is not None:
            filt = VCF_filter("samples", [samples])
            filter_list.append(filt)
        
        if (chrom is not None) and (start is not None) and (end is not None):
            filt = VCF_filter("region", (chrom, start, end))
            filter_list.append(filt)
        #fi
        
        if filters is not None:
            filt = VCF_filter("filters", [filters])
            filter_list.append(filt)
        #fi
        
        if vartypes is not None:
            filt = VCF_filter("vartypes", [vartypes])
            filter_list.append(filt)
        #fi
        
        if n_alleles is not None:
            filt = VCF_filter("n_alleles", [n_alleles])
            filter_list.append(filt)
        #fi
        
        if samples_format is not None:
            for (field, test) in samples_format:
                filt = VCF_filter('samples_format', [field, test])
                filter_list.append(filt)
            #efor
        #fi
        
        filter_functions = {
            "samples" :   lambda obj: obj.filter_samples,
            "region" :    lambda obj: obj.filter_region,
            "filters" :   lambda obj: obj.filter_filter,
            "vartypes" :  lambda obj: obj.filter_vartype,
            "n_alleles" : lambda obj: obj.filter_n_alleles,
            "samples_format" : lambda obj: obj.filter_samples_format
        }
        
        vcf_object = self._vcf
        for filt in filter_list:
            vcf_object = filter_functions[filt.name](vcf_object)(*filt.params)
        #efor
        
        return self.__class__(vcf_object=vcf_object, filter_stack=self._filter_stack + filter_list)
    #edef
    
    def select_records(self, record_indexes, filter_name="select_records", filter_params=[]):
        """
        Given a list of record indexes (integers), return a new VCF object which contains only those records
        """
        
        records_obj = VCF2_records([ self._vcf.records[i] for i in record_indexes ], self._vcf.samples)
        filt = VCF_filter(filter_name, filter_params)
        return self.__class__(vcf_object=records_obj, filter_stack=self._filter_stack + [ filt ])
    #edef
    
    def summary(altPos=None, refPos=None):
        """
        Make a table, per variant, of the frequency, etc.
        """
        self.to_records()
        return self._vcf.summary(altPos, refPos)
    #edef
    
    def to_records(self):
        """
        Convert the internal representation to a VCF_records representation
        """
        if not isinstance(self._vcf, VCF2_records):
            #records = [ self._internal_variant_representation(v) for v in self._vcf.records]
            #print(records)
            records_obj = VCF2_records(self._vcf.records, self._vcf.samples)
            self._vcf = records_obj
        #fi
    #edef

    def who_has(self, chrom=None, pos=None, ref=None, alt=None, var=None):
        """
        Determine who has a specific variant

        parameters:
        -----------
        chrom: str|int. Which chromosome the variant is on
        pos:   int. What position the variant is on.
        ref:   str. What is the reference allele?
        alt:   str. What is the alternative allele?

        returns:
        --------
        List of sample IDs for who has the variant.
        """
        
        if (chrom is not None) and (pos is not None) and (ref is not None) and (alt is not None) and (var is None):
            var = self.get_var(chrom, pos, ref, alt)
            if var is None:
                return []
            #fi
        elif (alt is not None) and (var is not None) and isinstance(var, cyvcf2.Variant):
            ref = var.REF
        else:
            raise ValueError("Incorrect parameter specification.")
        #fi

        ref = ref.upper()
        alt = alt.upper()

        alleles = [ var.REF.upper() ] + [ a.upper() for a in var.ALT ]

        alt_pos = alleles.index(alt)
        ref_pos = alleles.index(ref)

        rel_pos = alt_pos

        if ref_pos != 0:
            # In the case of an allele switch, we need to swap the reference and allele positions
            rel_pos = ref_pos
        #fi

        person_has = [ (a1 == rel_pos) or (a2 == rel_pos) for (a1,a2, phased) in var.genotypes ]

        return np.array(self.samples)[person_has]
    #edef
    
    def get_var(self, chrom, pos, ref, alt):
        """
        Get the variant record for a specific variant
        
        parameters:
        -----------
        chrom: str|int. Which chromosome the variant is on
        pos:   int. What position the variant is on.
        ref:   str. What is the reference allele?
        alt:   str. What is the alternative allele?
        
        returns:
        --------
        A cyvcf2.VCF.Variant object if the variant exists. Otherwise None
        """
        
        pos = int(pos)
        
        V = self._vcf.filter_region(chrom, pos, pos)
        
        ref = ref.lower()
        alt = alt.lower()
        for v in V:
            if (str(v.CHROM) == str(chrom)) and (v.POS == pos):
                vref = v.REF.lower()
                valt = [ a.lower() for a in v.ALT ]
                if ((vref == ref) and (alt in valt)):
                    return v
                elif ((vref == alt) and (ref in valt)):
                    utils.msg.dbm("There has been an allele switch at %s" % self.make_identifier(v))
                    return v
                #fi
            #fi
        #efor
        return None
    #edef
    
    def genotype_matrix(self):
        """
        Create a genotype matrix for each variant in the object
        
        NOTE:
        For a polyploid organism, the gt counts will be (e.g. for )
        
        """
        self.to_records()
        
        gts  = []
        idxs = []

        for var in self.records:
            for i, alt in enumerate(var.ALT):
                ident = self.make_identifier(var, i)
                pos = i+1
                var_gts = [ sum([a == pos for a in gt[:-1] ]) for gt in var.genotypes ]
                gts.append(var_gts)
                idxs.append(ident)
            #efor
        #efor

        return pd.DataFrame(gts, index=idxs, columns=self.samples).transpose()
    #edef
    
    @utils.decorators.class_or_instance_method
    def summary(obj, variants=None):
        """
        Create a summary matrix for each variant in the object
        
        ASSUMES a MONO/DIPLOID ORGANISM!
        """
        
        if variants is None:
            if obj.is_instance:
                obj.self.to_records()
                variants = obj.self.records
            else:
                raise ValueError("To use this function as a classmethod, you must specify a list of variants you wish to summarize.")
            #fi
        #fi
        
        S = [ v.summary() for v in variants ]
        return pd.concat(S)
    #edef
    

    @classmethod
    def make_identifier(cls, variant, alt_pos=0):
        """
        Make an identifier for a given variant.
        
        parameters:
        variant: a cyvcf2.Variant object.
        alt_pos: Integer. Which alternative allele to use in the identifier
        
        returns:
        A string identifier for the variant
        """
        return '%s-%d-%s-%s' % (str(variant.CHROM), variant.POS, variant.REF, variant.ALT[alt_pos])
    #edef

    @staticmethod
    def genotype_info_field_indexes(altAlleleID, ref=0):
        """
        For multi-allelic sites, genotype INFO fields have a particular order.
        If we want to extract information from these fields for a specific variant allele,
          we need to be able to reconstruct this order for an arbitrary number of alleles.
        The order is defined in the VCF specification, and implemented here.
        Search for 'GL : genotype likelihoods' in VCF spec.
        F(j/k) = (k*(k+1)/2)+j. I
        """
        
        F = lambda k,j : int(k*(k+1)/2)+j
        return [ F(ref, ref), F(ref, altAlleleID), F(altAlleleID,altAlleleID) ]
    #edef

    def __str__(self):
        """
        String representation of the object
        """
        dstr = "%s Object\n" % self.__class__.__name__
        dstr += "Internal representation: %s\n\n" % (self._vcf.__class__.__name__)
        
        dstr += 'Samples: %d\n\n' % len(self._vcf.samples)
        
        if len(self.filter_stack) > 0:
            dstr += 'Applied filters:\n'
            for filt in self.filter_stack:
                dstr += ' * %s\n' % str(filt)
            #efor
        #fi
        return dstr
    #edef
    
    def __repr__(self):
        """
        String representation of the object.
        """
        return str(self)
    #edef
    
    def __iter__(self):
        """
        Iterate over the records in this VCF object
        """
        return self._vcf.__iter__()
    #edef
    
    def __len__(self):
        """
        Return the number of records in this VCF object.
        """
        return len(self.records)
    #edef
    
    
    def __add__(self, other):
        """
        Add two VCF objects.
        Esentially concatenates the VCF records, on the basis that the samples described in both are identical.
        
        parameters:
        ------------
        other: A VCF2 object.
        
        Returns:
        a VCF2 object with VCF_records internal representation containing all VCF records.
        """
        if not isinstance(other, VCF2):
            raise ValueError("Operation not permitted for other object of type %s." % type(other))
        #fi
        
        if not all([ p[0] == p[1] for p in zip(self.samples, other.samples) ]):
            raise ValueError("The two VCF objects do not describe the same samples!!!")
        #fi
        
        records_obj = VCF2_records( list(self.records) + list(other.records), self.samples)
        filt = VCF_filter("ADD", [[ str(f) for f in self.filter_stack], [ str(f) for f in other.filter_stack ]])
        return self.__class__(vcf_object=records_obj, filter_stack=[filt])
    #edef

    
    @classmethod
    def merge(cls, to_merge):
        """
        Merge a set of VCF objects
        In general this is faster than adding multiple times
        """
        samples    = to_merge[0].samples
        inst_class = to_merge[0].__class__
        
        for o in to_merge[1:]:
            if not isinstance(o, inst_class):
                raise ValueError("Operation not permitted for other object of type %s." % type(other))
            #fi
            
            if not all([ p[0] == p[1] for p in zip(samples, o.samples) ]):
                raise ValueError("The VCF objects do not describe the same samples!!!")
            #fi
        #efor

        records = []
        for o in to_merge:
            records.extend(o.records)
        #efor
        
        records_obj = VCF2_records(records, samples)
        filt = VCF_filter("MERGE", ["%d objects" % len(to_merge)])
        return inst_class(vcf_object=records_obj, filter_stack=[filt])
#eclass

##############################################################################

class VCF2_master_type(object):
    def __init__(self, records=None, samples=None, internal_variant_representation=BiuVariant):
        self._vcf          = None
        self._records      = records
        self._samples      = samples
        self._internal_variant_representation = internal_variant_representation
    #edef
    
    @property
    def records(self):
        return self._records
    #edef
    
    @property
    def samples(self):
        return self._samples
    #edef
    
    def __iter__(self):
        return self.records.__iter__()
    #edef
    
    def __len__(self):
        return len(self.records)
    #edef
    
    def __str__(self):
        dstr = "%s object\n" % self.__class__.__name__
        return dstr
    #edef
    
    def __repr__(self):
        dstr = "%s object\n" % self.__class__.__name__
        return dstr
    #edef
#eclass

##############################################################################

class VCF2_records(VCF2_master_type):
    def __init__(self, *pargs, **kwargs):
        super(VCF2_records, self).__init__(*pargs, **kwargs)
        
        self._records = [ self._internal_variant_representation(r) if isinstance(r, cyvcf2.cyvcf2.Variant) else r
                             for r in self.records ]
        
    #edef
    
    def filter_vartype(self, vartypes):
        return VCF2_records([ v for v in self.records if v.var_type in vartypes ], self.samples,
                            internal_variant_representation=self._internal_variant_representation)
    #edef

    def filter_filter(self, filters):
        def test(var):
            if var.FILTER is None:
                return True
            #fi
            
            return not any([f in filters for f in var.FILTER])
        #edef
        
        return VCF2_records([ v for v in self.records if test(v) ], self.samples,
                            internal_variant_representation=self._internal_variant_representation)
    #edef

    def filter_samples(self, samples):
        raise Exception("Cannot filter samples in the VCF2_records representation. Please filter by samples first! See documentation.")
    #edef

    def filter_region(self, chrom, start, end):
        # Quick and dirty first
        test = lambda v: (str(v.CHROM) == str(chrom)) and (v.POS >= start) and (v.POS <= end)
        return VCF2_records([ v for v in self.records if test(v) ], self.samples,
                            internal_variant_representation=self._internal_variant_representation)
    #edef

    def filter_n_alleles(self, n_alleles):
        test = lambda v: len(v.ALT) <= n_alleles
        return VCF2_records([ v for v in self.records if test(v) ], self.samples,
                            internal_variant_representation=self._internal_variant_representation)
    #edef
    
    def filter_samples_format(self, field, test):
        """
        Filter samples on a test on a format field.
        All samples that do not pass this filter are set to the missing genotype
        Variants in which no individual has the variant are removed from the list also.
        
        NOTE: The value of the field in question will ALWAYS BE A LIST (split on ';')
              THIS IS REGARDLESS OF HOW MANY VALUES THERE SHOULD BE IN THE FORMAT
        
        parameters:
        self: VCF2_records object
        field: string
        test: function
        
        """
        
        new_records = []
        
        for r in self.records:
            # It is IMPERATIVE THAT THE VARIANT IS IN "MyOwnVariant" format here!!!
            new = r.switch().copy()
            
            GT = r.format('GT')
            PL = r.ploidy
            F  = r.format(field)
            
            new_gt = []
            for gt, pl, f in zip(GT, PL, F):
                if test(f):
                    new_gt.append(gt)
                else:
                    new_gt.append('/'.join(['.']*pl))
                #fi
            #efor
            
            new.set_format('GT', new_gt)
            
            for gt in new.genotypes:
                if any([x > 0 for x in gt[:-1]]):
                    new_records.append(new)
                    break
                #fi
            #efor
        #efor
        
        return VCF2_records(new_records, self.samples,
                            internal_variant_representation=self._internal_variant_representation)
    #edef

#eclass

##############################################################################

class VCF2_cyvcf2(VCF2_master_type):
    def __init__(self, file, *pargs, **kwargs):
        super(VCF2_cyvcf2, self).__init__(*pargs, **kwargs)
        
        self._vcf = cyvcf2.VCF(file, lazy=True, gts012=True, samples=self._samples)
        self._file = file
    #edef
    
    @property
    def records(self):
        if self._records is None:
            self._records = [ self._internal_variant_representation(v) for v in self._vcf ]
        #fi
        
        return self._records
    #edef
    
    @property
    def samples(self):
        return self._vcf.samples
    #edef
    
    def filter_vartype(self, vartypes):
        return VCF2_records([ r for r in self.records], self.samples,
                            internal_variant_representation=self._internal_variant_representation).filter_vartype(vartypes)
    #edef

    def filter_filter(self, filters):
        return VCF2_records([ r for r in self.records], self.samples,
                            internal_variant_representation=self._internal_variant_representation).filter_filter(filters)
    #edef

    def filter_samples(self, samples):
        return VCF2_cyvcf2(self._file, samples,
                           internal_variant_representation=self._internal_variant_representation)
    #edef

    def filter_region(self, chrom, start, end):
        region = '%s:%s-%s' % (str(chrom), str(int(start)), str(int(end)))
        return VCF2_records([ r for r in self._vcf(region) ], self.samples,
                            internal_variant_representation=self._internal_variant_representation)
    #edef

    def filter_n_alleles(self, n_alleles):
        return VCF2_records([ r for r in self.records], self.samples,
                            internal_variant_representation=self._internal_variant_representation).filter_n_alleles(n_alleles)
    #edef
    
    def filter_samples_format(self, field, test):
        """
        Filter samples on a test on a format field.
        All samples that do not pass this filter are set to the missing genotype
        
        parameters:
        self: VCF2_records object
        field: string
        test: function
        
        """
        return VCF2_records([ r for r in self.records], self.samples,
                            internal_variant_representation=self._internal_variant_representation).filter_samples_format(field, test)
    #edef
#eclass

###################################################################