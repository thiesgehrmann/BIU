import matplotlib.pylab as plt

def subplots(**kwargs):
  fig, axes = ply.subplots(**kwargs)
  if hasattr(axes, '__len__'):
    axes = axes.flatten()
  #fi
  return fig, axes
#edef
  
