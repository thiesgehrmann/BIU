import pandas as pd

from ..structures import Dataset2
from .. import utils
from .. import formats

pd = utils.py.loadExternalModule('pandas')

class HGNC(Dataset2):
    """
    The HUGO Gene Name Consortium Annotations.
    
    This is an interface to the HGNC Annotations.
    All annotations are provided here.
    Given a symbol, for example, you can look up the other identifiers it has:
      hgnc = biu.db.HGNC()
      hgnc.symbol['MTOR'].ensembl_gene_id
      
    There are many annotations/labels, and they can all be accessed the same way:
      hgnc.xxx
    Where `xxx` is the name of the annotation.
    
    Each attribute `hgnc.xxx` is encoded as a MappingIndex Object.
    See the documentation for this object for more details on how to access it.
    
    You can look up the values that are indexed in each annotation with, e.g.:
      hgnc.symbol.keys()
      
    Take care in how you use these annotations.
    Note that if several exist, they are separated by a '|'
      
    These are the annotations that can be accessed with `hgnc.xxx`
      hgnc_id: The HGNC Identified
      symbol:  The symbol it is attached to
      name:    The name that describes it
      locus_group: The group of locus it is in. Note that
      locus_type: 
      status:  Is it approved/waiting
      location:  What is the location of this locus
      location_sortable: 
      alias_symbol: What alias symbols exist?
      alias_name:   What alias names exist.
      prev_symbol:  What previous symbols existed?
      prev_name:    What previous names existed?
      gene_family:  
      gene_family_id: 
      date_approved_reserved: 
      date_symbol_changed: 
      date_name_changed: 
      date_modified: 
      entrez_id: 
      ensembl_gene_id: 
      vega_id: 
      ucsc_id: 
      ena: 
      refseq_accession: 
      ccds_id: 
      uniprot_ids: 
      pubmed_id: 
      mgd_id: 
      rgd_id: 
      lsdb: 
      cosmic: 
      omim_id: 
      mirbase: 
      homeodb: 
      snornabase: 
      bioparadigms_slc: 
      orphanet: 
      pseudogene_org: 
      horde_id: 
      merops: 
      imgt: 
      iuphar: 
      kznf_gene_catalog: 
      mamit_trnadb: 
      cd: 
      lncrnadb: 
      enzyme_id: 
      intermediate_filament_db: 
      rna_central_ids: 
      lncipedia: 
      gtrnadb: 
    
    """
    def __init__(self, *pargs, **kwargs):
        super(HGNC, self).__init__("HGNC", *pargs, **kwargs)
    
        def process_table(fin, fout):
    
            T = pd.read_csv(fin, sep='\t', dtype=object)

            def fillzip(*pargs):
                mlen = max([len(x) for x in pargs])
                farg = [ list(x) + ([None] * (mlen-len(x))) for x in pargs ]
                return zip(*farg)
            #edef

            mult_cols = T.columns[T.applymap(lambda x: '|' in x if hasattr(x, '__contains__') else False).any()]
            sel = A.applymap(lambda x: '|' in x if hasattr(x, '__contains__') else False).any(axis=1)

            T_no_mult = T[~sel]
            T_mult    = T[sel]
            T_unmult = []
            for i, row in T_mult.iterrows():
                if i % 100 == 0:
                    print('\r%d/%d' % (i, T.shape[0]), end='')
                #fi
                mult = [ row[col].split('|') if not(pd.isna(row[col])) else [] for col in mult_cols ]
                for m in fillzip(*mult):
                    row[mult_cols] = m
                    T_unmult.append(row.copy())
                #efor
            #efor
            T = pd.concat([T_no_mult, pd.DataFrame(T_unmult, columns=T.columns.values)])
            #efor

            T.to_pickle(fout)
            return utils.Acquire2.STATUS_SUCCESS
        #edef
    
        f = biu.utils.Acquire2().curl('ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt')
        f = f.func(lambda fin, fout: process_table(fin[0], fout))
        self._obj.add_file("hgnc.pkl", f)
        
        self._obj.register('_table', ["hgnc.pkl"], lambda f: pd.read_pickle(f["hgnc.pkl"]))
        
        ##################
        # Here there is a little trick to make sure that items like this are properly registered:
        # lambda f, const1=const1, const2=const2: func(f, const1, const2)
        # Internally in the Dataset, it is called as func(f), and the consts are constant.
        ##################
        attributes = [(0, 'hgnc_id'),
                     (1, 'symbol'),
                     (3, 'locus_group'),
                     (4, 'locus_type'),
                     (5, 'status'),
                     (6, 'location'),
                     (7, 'location_sortable'),
                     (8, 'alias_symbol'),
                     (10, 'prev_symbol'),
                     (12, 'gene_family'),
                     (13, 'gene_family_id'),
                     (14, 'date_approved_reserved'),
                     (15, 'date_symbol_changed'),
                     (16, 'date_name_changed'),
                     (17, 'date_modified'),
                     (18, 'entrez_id'),
                     (19, 'ensembl_gene_id'),
                     (20, 'vega_id'),
                     (21, 'ucsc_id'),
                     (22, 'ena'),
                     (23, 'refseq_accession'),
                     (24, 'ccds_id'),
                     (25, 'uniprot_ids'),
                     (26, 'pubmed_id'),
                     (27, 'mgd_id'),
                     (28, 'rgd_id'),
                     (29, 'lsdb'),
                     (30, 'cosmic'),
                     (31, 'omim_id'),
                     (32, 'mirbase'),
                     (33, 'homeodb'),
                     (34, 'snornabase'),
                     (35, 'bioparadigms_slc'),
                     (36, 'orphanet'),
                     (37, 'pseudogene_org'),
                     (38, 'horde_id'),
                     (39, 'merops'),
                     (40, 'imgt'),
                     (41, 'iuphar'),
                     (42, 'kznf_gene_catalog'),
                     (43, 'mamit_trnadb'),
                     (44, 'cd'),
                     (45, 'lncrnadb'),
                     (46, 'enzyme_id'),
                     (47, 'intermediate_filament_db'),
                     (48, 'rna_central_ids'),
                     (49, 'lncipedia'),
                     (50, 'gtrnadb')]
        for idx, attr in attributes:
            self._obj.register(attr, [], lambda f, idx=idx: formats.MappingIndex(self._table, idx))
        #efor
    #edef
#eclass