from ..structures import Dataset2
from .. import formats
from .. import utils

import xml.etree.ElementTree as ET

###############################################################################

class GO2(Dataset2):
    """
    Interface to the Gene Ontology annotation dataset
    """

    versions = {
      "human" : "http://geneontology.org/gene-associations/goa_human.gaf.gz",
      "mouse" : "http://geneontology.org/gene-associations/gene_association.mgi.gz",
      "drosophilia" : "http://geneontology.org/gene-associations/gene_association.fb.gz" }
    
    _version = None

    def __init__(self, version=list(versions.keys())[0], *pargs, **kwargs):
        """
        Initialize the object.
        
        parameters:
        -----------
        version: [human|mouse|drosophilia] The annotation set you want to use\
        *pargs, **kwargs: Additional arguments to biu.structures.Dataset2
        """
        super(GO2, self).__init__("GeneOntology/%s" % version, *pargs, **kwargs)
        self._version = version
        
        def parseGOXML(infile, outfile):
            xml = ET.parse(infile[0])
            root = xml.getroot()
            terms    = [ child for child in root if child.tag == 'term' ]
            termInfo = [ { tc.tag : (tc.text if (tc.tag in ['id','name','namespace']) else tc.find('defstr')) for tc in term if (tc.tag in ['id','name','namespace', 'def'])  } for term in terms]
            termInfo = [ ( t['id'], t['namespace'], t['name'], t['def'].text if hasattr(t['def'], 'text') else '') for t in termInfo ]
   
            with open(outfile, 'w') as ofd:
                ofd.write('\t'.join(['id', 'namespace', 'name', 'desc']) + '\n')
                for ti in termInfo:
                    ofd.write('\t'.join(ti) + '\n')
                #efor
              #ewith
            return 0
        #edef
        
        self._obj.add_file("terminfo.tsv", utils.Acquire2().curl('http://archive.geneontology.org/latest-termdb/go_daily-termdb.obo-xml.gz').gunzip().func(parseGOXML))
        
        self._obj.add_file("annots.gaf", utils.Acquire2().curl(self.versions[self._version]).gunzip())
        
        self._obj.register("_gaf", [ "annots.gaf" ],  lambda f: formats.GAF(f['annots.gaf'], comment='!', delimiter='\t'))
        self._obj.register("terminfo", [ "terminfo.tsv" ], lambda f: formats.Map(f["terminfo.tsv"], header=True, 
                                                                                 names=['id', 'namespace', 'name', 'desc'],
                                                                                 delimiter='\t') )
        

        self._add_str_func(lambda s: "Version: %s" % self._version)
      #edef


    ###############################################################################

    @property
    def annotations(self):
      """
        Get a list of all possible annotations (GO terms)
      """
      return self._gaf.annotations
    #edef

    @property
    def objects(self):
        """
          Get a list of all objects that are annotated (usually uniprot protein IDs)
        """
        return self._gaf.objects
    #edef

    def get_annots(self, objectID):
        """
          get_annots : Get a list of annotations for a specific object (uniprot protein ID)
          Input:
           - ObjectID : (Uniprot protein ID)
          Output:
           - List of GO annotations for given objectID
        """
        return self._gaf.get_annots(objectID)
    #edef

    def get_annotated(self, annotID):
        """
        get_annotated : Get all objects annotated with a speficic GO term
        Input:
          - annotID : GO term
        Output:
          - List of objectIDs (uniprot protein IDs)
        """
        return self._gaf.get_annotated(annotID)
    #e  def

    def enrich(self, your_set, pathway=None, abcd_values=False,  method=None, **kwargs):
        """
        Enrich: Check enrichment of GO pathways in a given set
  
        Inputs:
          - yourSet: List of Uniprot protein IDs to test
          - pathway: List of GO terms (or single term) to test (Defaults to all terms that your objectIDs are present in)
          - abcd_values: Boolean. If True, it will return the actual element values in the contingency table, rather than just counts
          - method: Type of multiple testing correction procedure to use
          - **kwargs: Additional arguments for multple testing procedure
  
        Outputs:
         - df : Pandas Data Frame of test results
        """
        return self._gaf.enrich(your_set, pathway, abcd_values, method, **kwargs)
    #edef
                           
    def summary(self, your_sets):
        """
        For a set of proteins, count the occurrence of each pathway.
  
        parameters:
        -----------
        your_sets: List[String] | List[List[String]] | dict[string->List[String]]
            A list of uniprot protein identifiers, or a list/dict of lists of uniprot identifiers
  
        Returns:
        Pandas Dataframe of pathway counts
          
        example usage:
          
        """
  
        S = self._gaf.summary(your_sets)
        S['description'] = [ self.terminfo.id[p][0].name for p in S.index ]
        return S
    #edef
      
    def annotate(self, your_sets):
        """
        For a set of proteins, Annotate them
  
        parameters:
        -----------
        your_sets: List[String] | List[List[String]] | dict[string->List[String]]
            A list of uniprot protein identifiers, or a list/dict of lists of uniprot identifiers
  
        Returns:
        Pandas Dataframe of pathway counts
          
        example usage:
          
        """

        A = self._gaf.annotate(your_sets)
        A['description'] = [ self.terminfo.id[p][0].name for p in A.pathway ]
        return A
    #edef

#eclass