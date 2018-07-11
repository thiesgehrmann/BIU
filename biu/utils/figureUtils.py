from . import pyUtils as py

plt = py.loadExternalModule('matplotlib.pylab')
pltvenn = py.loadExternalModule('matplotlib_venn')

###############################################################################

def subplots(dpi=300, **kwargs):
  fig, axes = plt.subplots(**kwargs)
  if hasattr(axes, '__len__'):
    axes = axes.flatten()
  else:
    axes = [axes]
  return fig, axes
#edef

###############################################################################

def venn(*sets, ax=None, names=None):
  import  matplotlib_venn as pltvenn
  if len(sets) == 2:
    return pltvenn.venn2(sets, ax=ax, set_labels=names)
  elif len(sets) == 3:
    return pltvenn.venn3(sets, ax=ax, set_labels=names)
  else:
    return
  #fi
#edef
