
from . import genomes as genomes
from . import db as db
from . import maps as maps
from . import formats as formats
from . import structures as structures
from . import utils as utils
from . import stats as stats

from . import processing as processing

from .config import config
from .config import settings

from . import pipelines as pipelines

def __version__():
  print("BIU (Bio Utilities) python module")
  print(config.settings.dumps())
  print(" Current config hash: %s" % processing.lst.hash(config.settings.dumps()))
#edef

__missingDependencies = settings.missingDependencies()
if len(__missingDependencies[0]) > 0:
  utils.msg.warning("The following dependencies of BIU are missing. Functionality of BIU will be affected.\n  %s" % ', '.join(__missingDependencies[0]))
#fi

if len(__missingDependencies[1]) > 0:
  utils.msg.warning("Some optional dependencies of BIU are missing. Functionality of BIU may be affected.\n  %s" % ', '.join(__missingDependencies[1]))
#fi
