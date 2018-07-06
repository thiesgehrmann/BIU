import hashlib

def freq(lst):
  """
    Count the number of occurances in a list.
    Input: an iterable list
    Output: Dictionary with count of each item
  """
  F = {}
  for item in lst:
    F[item] = F.get(item, 0) + 1
  #efor
  return F
#edef

def group(lst, key=lambda x: x[0]):
  """
    Group items based on a certain key.
    Input: lst: an iterable list of tuple (or indexable values)
           key: A function to determine grouping
    Output: Dictionary with items grouped by key.
  """
  F = {}
  for item in lst:
    k = key(item)
    F[k] = F.get(k, []) + [item]
  #efor
  return F
#edef

def flatten(lst):
  return [ item for group in lst for item in group ]
#edef

def hash(arr, strategy="tmb", f=hashlib.md5):
  h = f()

  if strategy == 'tmb': # Top Middle Bottom
    middle = int(len(arr) / 2)
    for o in arr[:10] + arr[middle:middle+10] + arr[:-10]:
      h.update(str(o).encode())
    #efor
  elif strategy == 'all':
    for o in arr:
      h.update(str(o).encode())
    #efor
  #fi

  return h.hexdigest()
#edef
