from .. structures import Dataset2
from .. import formats
from .. import utils
from .. import stats

pd = utils.py.loadExternalModule("pandas")
np = utils.py.loadExternalModule('numpy')

import os

###############################################################################

# https://www.biostars.org/p/71737/


class TFCAT(Dataset2):
    """
    The TFCat dataset contains a list of putative and verified Transcription Factors.
    
    You can access the Dataset with an entrezID
    """

    def __init__(self, min_strength=1, *pargs, **kwargs):
        """
          Load the TFCat dataset

          Inputs:
           - min_strength: Integer[0,3] The minimum evidence strength to allow a TF in the list
           - *pargs, **kwargs : Dataset2 structure arguments
          Output:
           KEGG data structure
        """
        super(TFCAT, self).__init__("TFCAT", *pargs, **kwargs)
        
        url = 'http://cisreg.ca/tfcat/index.php'
        args = [
            ('-H', 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8'),
            ('-H', 'Accept-Language: en-GB,en;q=0.5'),
            ('--compressed',),
            ('-H', 'Referer: http://cisreg.caindex.php'),
            ('-H', 'Content-Type: application/x-www-form-urlencoded'),
            ('-H', 'Connection: keep-alive'),
            ('-H', 'Upgrade-Insecure-Requests: 1'),
            ('-H', 'DNT: 1'),
            ('--data', 'gene_taxonomy_classification=ALL%2C&gene_judgement_classification=ALL%2C&gene_DBD_protein_group_classification=ALL%2C&gene_DBD_protein_family_classification=ALL%2C&judgement_list=ALL&taxonomy_list=ALL&protein_group_list=ALL&protein_family_list=ALL&cluster_id_sel=&gene_id_sel=&action=get_TFCat_data') ]
        f = utils.Acquire2().curl(url, args=args)
        f = f.cmd("""curl $(cat {} | grep -e "javascript:newWindow=window.open('tmp/" | cut -d\\' -f2 | uniq | sed -e 's;.*;http://cisreg.ca/tfcat/&;')""",
                 placeholder='{}')
        f = f.cmd("sed -e '/^[^0-9a-zA-Z]/d'")

        self._obj.add_file('db_file', f)
        
        def load_table(f):
            T = pd.read_csv(f["db_file"], sep='\t', encoding='ISO-8859-1', dtype=object,
                           names=['entrez', 'description', 'group', 'reviewer',
                                  'judgement', 'taxa', 'reviewer_commments', 'pubmed_id', 'function',
                                  'species', 'evidence_strength', 'pubmed_comments'], skiprows=1)

            smap = { 'weak' : 1,
                     'medium' : 2,
                     'strong' : 3,
                      'nan': 0,
                      'unknown' : 0,
                      'not selected' : 0}
            T['score'] = T.evidence_strength.apply(lambda s: smap[str(s).lower()])
            
            return T[T.score >= min_strength]
        
        self._obj.register('_table', ["db_file"], lambda f: load_table(f))
        
        for idx, attr in [(0, 'entrez'),
                          (2, 'group'),
                          (7, 'pubmed_id'),
                          (8, 'function'),
                          (9, 'species')]:
            self._obj.register(attr, [], lambda f, idx=idx: formats.MappingIndex(self._table, key=idx) )
        #efor
        
        self._add_str_func(lambda x: "Minimum strength: %d" % min_strength)
    #edef
#eclass