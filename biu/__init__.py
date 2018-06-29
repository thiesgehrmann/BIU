
from . import genomes as genomes
from . import db as db
from . import maps as maps
from . import formats as formats
from . import structures as structures
from . import utils as utils
from . import stats as stats

from . import processing as processing

from .config import config

from . import pipelines as pipelines

def __version__():
  print("BIU (Bio Utilities) python module")
  print(config.settings.dumps())
  print(" Current config hash: %s" % processing.lst.hash(config.settings.dumps()))
#edef
