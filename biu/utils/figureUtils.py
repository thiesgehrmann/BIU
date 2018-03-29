import matplotlib.pylab as plt

def subplots(dpi=300, **kwargs):
  fig, axes = plt.subplots(**kwargs)
  if hasattr(axes, '__len__'):
    axes = axes.flatten()
  else:
    axes = [axes]
  return fig, axes
#edef
  
