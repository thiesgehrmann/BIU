import sys

from ..config import settings as settings

###############################################################################

def dbm(message):
  if settings.getDebugState():
    for line in str(message).split('\n'):
      if settings.getDebugStream() == 'stdout':
        sys.stdout.write('D: %s\n' % line)
      else:
        sys.stderr.write('D: %s\n' % line)
      #fi
    #efor
  #fi
#edef

###############################################################################

def error(message):
  if settings.getErrorState():
    for line in str(message).split('\n'):
      sys.stderr.write('E: %s\n' % line)
    #efor
  #fi
#edef

###############################################################################

def warning(message):
  if settings.getWarningState():
    for line in str(message).split('\n'):
      sys.stderr.write('W: %s\n' % line)
    #efor
  #fi
#edef

###############################################################################
