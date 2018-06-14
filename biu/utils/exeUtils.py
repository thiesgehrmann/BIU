import os
import shlex, subprocess, os

from . import msgUtils as msg

###############################################################################

def runCommand(cmd, bg=False, stdin=None, stdout=None, stderr=None, shell=False, verbose=False):
  if verbose:
    msg.dbm(cmd)
  #fi

  if not shell:
    cmd = shlex.split(cmd)
  #fi

  p = subprocess.Popen(cmd, stdin=stdin, stdout=stdout, stderr=stderr, shell=shell);
  if bg:
    return p;
  else:
    (pid, r) = os.waitpid(p.pid, 0);
    return r;
  #fi
#edef

###############################################################################

def getCommandOutput(cmd, stdin=None, stderr=False, shell=False, verbose=False):
  if verbose:
    msg.dbm(cmd)
  #fi

  if not shell:
    cmd = shlex.split(cmd)
  #fi

  #p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell);
  p = subprocess.run(cmd, stdin=stdin, stdout=subprocess.PIPE, stderr=(subprocess.PIPE if stderr else None), shell=shell)

  if stderr:
    return (p.stdout, p.stderr)
  #fi
  return p.stdout
#edef
