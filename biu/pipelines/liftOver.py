from . import Pipeline
from .. import utils

import vcf
import pandas as pd

import inspect, os
snakemakeFile  = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + '/liftOver/Snakefile'

chainFiles = {
  ("hg19", "hg38") : "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
  ("hg38", "hg19") : "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
}

def listVersions():
  print("Available liftOver chain Files:")
  for (fromAssembly, toAssembly) in chainFiles:
    print(" * %s -> %s" % (fromAssembly, toAssembly))
  #efor

#edef


class LiftOver(Pipeline):

  _defaultConfig = {
    "output_file_name" : "liftedOutput.bed"
  }

  _output = None

  def __init__(self, data, fromAssembly="hg19", toAssembly="hg38", config=None, where=None, **kwargs):

    Pipeline.__init__(self, snakemakeFile, **kwargs)

    smConfig = self._defaultConfig
    if isinstance(config, dict):
      smConfig.update(config)
    #fi

    smConfig["fromAssembly"] = fromAssembly.lower()
    smConfig["toAssembly"]   = toAssembly.lower()
    smConfig["chainURL"]     = chainFiles[(smConfig["fromAssembly"], smConfig["toAssembly"])]
    if isinstance(data, str):
      smConfig["input_file"] = data
    else:
      smConfig["input_file"] = self._writeTemporaryFile(data)
    #fi

    self.setConfig(smConfig)
    self.run(["output"])
  #edef

  def _writeTemporaryFile(self, data):
    fileName, exists = self._generateInputFileName(data)

    def chrFormatter(s):
      if s[:3] == 'chr':
        return s
      else:
        return 'chr%s' % s
      #fi
    #edef

    if not(exists):
      with open(fileName, "w") as ofd:
        for i, region in enumerate(data):
          if len(region) == 2:
            ofd.write("%s\t%s\t%s\t%d\n" % (chrFormatter(region[0]), str(region[1]), str(region[1]), i))
          elif len(region) == 3:
            ofd.write("%s\t%s\t%s\t%d\n" % (chrFormatter(region[0]), str(region[1]), str(region[2]), i))
          else:
            utils.error("The provided region %s is not valid. Must have tuple format (chr,start[,end]).")
          #fi
        #efor
      #ewith
    #fi

    return fileName
  #edef

  def getLiftOver(self):
    if not(self.success):
      utils.error("The pipeline did not complete successfully, you cannot retrieve results yet.")
      return None
    #fi
    if self._output is None:
      output = pd.read_csv(self.getAnnotationFileName(),
                           delimiter='\t', comment='#', header=None)
      output = output.rename(columns={0:"chr", 1:"start", 2:"end", 3:"id"})
      output["chr"] = output["chr"].apply(lambda x: x[3:])

      if "id" in output.columns:
        indexedOutput = []
        for i, row in enumerate(output.values):
          if i == row[3]:
            indexedOutput.append( (row[0], row[1], row[2]) )
          else:
            indexedOutput.append( None )
          #fi
        #efor
        output = indexedOutput
      else:
        output = output.values
      #fi
      self._output = output
    #fi
    return self._output
  #edef

  def getAnnotationFileName(self):
    return '%s/%s' % (self.config["outdir"], self.config["output_file_name"])
  #edef

#eclass
###############################################################################
