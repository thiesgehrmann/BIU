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

  p = subprocess.Popen(cmd, stdin=stdin, stdout=stdout, stderr=stderr, shell=shell, env=os.environ);
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
  p = subprocess.run(cmd, stdin=stdin, stdout=subprocess.PIPE, stderr=(subprocess.PIPE if stderr else None), shell=shell, env=os.environ)

  if stderr:
    return (p.stdout, p.stderr)
  #fi
  return p.stdout
#edef

###############################################################################

def which(program):
    def is_exe(fpath):
      return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
      if is_exe(program):
        return program
      #fi
    else:
      for path in os.environ["PATH"].split(os.pathsep):
        exe_file = os.path.join(path, program)
        if is_exe(exe_file):
          return exe_file
        #fi
      #efor
    #fi

    return None
#edef

###############################################################################

def exists(cmd):
  if which(cmd) is None:
    return False
  #fi
  return True
#edef
  
