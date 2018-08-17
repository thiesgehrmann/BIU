from ..structures import Dataset
from .. import formats
from .. import utils
from ..config import settings

import xml.etree.ElementTree as ET

###############################################################################

class GO(Dataset):

  versions = {
    "human" : "http://geneontology.org/gene-associations/goa_human.gaf.gz",
    "mouse" : "http://geneontology.org/gene-associations/gene_association.mgi.gz",
    "drosophilia" : "http://geneontology.org/gene-associations/gene_association.fb.gz" }

  def __init__(self, version=list(versions.keys())[0], where=None, **kwargs):
    fileIndex = self.__genFileIndex(version, where=where)
    Dataset.__init__(self, fileIndex, **kwargs)
    self.version = version

    self._registerObject('_gaf', formats.GAF, ['gaf'], fileIndex["gaf"].path, skiprows=1, delimiter='\t')
    self._registerObject('terminfo', formats.Map, ['terminfo'], fileIndex["terminfo"].path, header=True, names=['id', 'namespace', 'name', 'desc'], delimiter='\t')

    self._addStrFunction(lambda s: "Version: %s" % self.version)
  #edef

  def __genFileIndex(self, version, where=None):
     finalPath = '%s/geneOntology/%s' % ( (settings.getDataDir() if where is None else where), version)
     url = self.versions[version]

     def parseGOXML(infile, outfile):
       xml = ET.parse(infile)
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

     files = {}
     files['gaf'] = utils.Acquire(where=where).curl(url).gunzip().finalize('%s/annots.gaf' % finalPath)
     files['terminfo'] = utils.Acquire().curl('http://archive.geneontology.org/latest-termdb/go_daily-termdb.obo-xml.gz').gunzip().func(parseGOXML).finalize('%s/terminfo.tsv' % finalPath)
     return files
  #edef


  ###############################################################################

  @property
  def gaf(self):
    return self._getObject("gaf")
  #edef

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

  def getAnnots(self, objectID):
    """
      getAnnots : Get a list of annotations for a specific object (uniprot protein ID)
      Input:
       - ObjectID : (Uniprot protein ID)
      Output:
       - List of GO annotations for given objectID
    """
    return self._gaf.getAnnots(objectID)
  #edef

  def getAnnotated(self, annotID):
    """
    getAnnotated : Get all objects annotated with a speficic GO term
    Input:
      - annotID : GO term
    Output:
      - List of objectIDs (uniprot protein IDs)
    """
    return self._gaf.getAnnotated(annotID)
  #edef

  def enrich(self, *pargs, **kwargs):
    """
    Enrich: Check enrichment of GO terms on a given set.
      See help for biu.formats.GAF.enrich
    """
    return self._gaf.enrich(*pargs, **kwargs)
  #edef

#eclass
