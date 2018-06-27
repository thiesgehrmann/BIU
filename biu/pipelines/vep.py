from . import Pipeline
from .. import formats

import pandas as pd

import inspect, os
snakemakeFile  = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + '/vep/Snakefile'

###############################################################################

class VEP(Pipeline):

  _defaultConfig = {
    "install_species" : "homo_sapiens",
    "install_assembly" : "GRCh38",
    "vep_options" : "",
    "output_file_name" : "vep_annotations.tsv"
  }

  _output = None

  def __init__(self, data, config=None, **kwargs):

    Pipeline.__init__(self, snakemakeFile, **kwargs)

    smConfig = self._defaultConfig
    if isinstance(config, dict):
      smConfig.update(config)
    #fi

    if isinstance(data, formats.VCF):
      data = [ record for record in data ]
    #fi

    smConfig["vcf_file"] = self.__writeTemporaryFile(data)

    self.setConfig(smConfig)
    self.run(["output"])
  #edef

  def __writeTemporaryFile(self, vcfArray, hashName=True):
    fileName, exists = self._generateInputFileName(vcfArray)

    if not(exists):
      with open(fileName, "w") as ofd:
        for record in vcfArray:
           alts = '/'.join([ a.sequence if hasattr(a, 'sequence') else '.' for a in record.ALT])
           ofd.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (record.CHROM, record.POS, record.POS, '%s/%s' % ( record.REF, '/'.join(alts)), '+', formats.VCF.makeIdentifier(record))) 
        #efor
      #ewith
    #fi

    return fileName
  #edef

  def getAnnotations(self):
    if not(self.success):
      return None
    #fi

    if self._output is None:
      self._output = pd.read_csv(self.getAnnotationFileName(), 
                                 delimiter='\t', comment='#', header=None, 
                                 names=["variant_id","location","allele","gene","feature","feature_type","consequence",
                                        "cdna_position","cds_position","protein_position","amino_acids","codons",
                                        "existing_variation","extra"])
    #fi
    return self._output
  #edef

  def getAnnotationFileName(self):
    return '%s/%s' % (self.config["outdir"], self.config["output_file_name"])
  #edef

#eclass
