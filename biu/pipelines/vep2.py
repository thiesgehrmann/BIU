from ..structures import Pipeline
from .. import formats
from .. import utils

pd = utils.py.loadExternalModule("pandas")

import inspect, os
snakemakeFile  = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + '/vep/Snakefile'

###############################################################################

class VEP2(Pipeline):

    __defaultConfig = {
      "install_species" : "homo_sapiens",
      "install_assembly" : "GRCh38",
      "vep_options" : "",
      "output_file_name" : "vep_annotations.tsv"
    }

    _output = None

    def __init__(self, data, config=None, **kwargs):
        """
        Initialize the VEP pipeline, designed to work with the VCF2 formatted VCF records
        
        parameters:
        -----------
        data: List[VCF2 records]
        config: Additional options for the VEP pipeline
        **kwargs: 
        """

        Pipeline.__init__(self, snakemakeFile, {**self.__defaultConfig, **config}, **kwargs)

        self.setConfig(vcf_file=self.__writeTemporaryFile(data))
        if self.autorun:
            self.run(["output"])
        #fi
    #edef

    def __writeTemporaryFile(self, data, hashName=True):
        fileName, exists = self._generateInputFileName(data)

        if not(exists):
            with open(fileName, "w") as ofd:
                for record in data:
                     alts = '/'.join([ a if a else '.' for a in record.ALT ])
                     ofd.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (record.CHROM, record.POS, record.POS, '%s/%s' % ( record.REF, alts), '+', formats.VCF2.make_identifier(record))) 
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