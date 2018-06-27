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
    F[k] = F.get(k, []) + [append(item)]
  #efor
  return F
#edef
