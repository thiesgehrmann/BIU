from .. import utils

np = utils.py.loadExternalModule("numpy")

def linearInterpolation(curveX, curveY, p, isSorted=False):
    """
    linearInterpolation: Linearly Interpolate a value on a curve
    
    Inputs:
      curveX : The X axis coordinates for the curve
      curveY : The Y axis coordinates for the curve
      p: The location to interpolate to (Value or list of values)
      isSorted : Is it sorted? (Otherwise I sort it)
      
    Output:
      - Interpolated point (or points if c is list)
    """

    if not isSorted:
        curveX, curveY = zip(*sorted(zip(curveX, curveY), key=lambda x: x[0]))
    #fi
    
    if hasattr(p, '__len__'):
        return np.array([ linearInterpolation(curveX, curveY, v) for v in p ])
    #fi

    pos = np.searchsorted(curveX, p, side='right', sorter=None)
    
    if pos == 0:
        return curveY[0]
    elif pos >= len(curveX)-1:
        return curveY[-1]
    #fi
    
    x1, y1 = curveX[pos], curveY[pos]
    x2, y2 = curveX[pos+1], curveY[pos+1]
    return y1 + (y2 - y1) * (p - x1) / (x2 - x1)
#edef
