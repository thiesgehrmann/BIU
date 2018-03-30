# Adapted from https://codereview.stackexchange.com/q/10481

class LazyProxy(object):
  def __init__(self, cls, *params, **kwargs):
    self.__dict__["_cls"]=cls
    self.__dict__["_params"]=params
    self.__dict__["_kwargs"]=kwargs

    self.__dict__["_obj"]=None
  #edef

  def __init_obj(self):
    self.__dict__["_obj"]=object.__new__(self.__dict__["_cls"])
    self.__dict__["_obj"].__init__(*self.__dict__["_params"],
                                   **self.__dict__["_kwargs"])
  #edef

  def __getattr__(self, name):
    if self.__dict__["_obj"] is None:
      self.__init_obj()
    #fi
    return getattr(self.__dict__["_obj"], name)
  #edef

  def __setattr__(self, name, value):
    if self.__dict__["_obj"] is None:
      self.__init_obj()
    #fi
    setattr(self.__dict__["_obj"], name, value)
  #edef

  def __str__(self):
    if self.__dict__["_obj"] is None:
      self.__init_obj()
    #fi
    return self.__dict__["_obj"].__str__()
  #edef

  def __call__(self):
    if self.__dict__["_obj"] is None:
      self.__init_obj()
    #fi
    return self.__dict__["_obj"].__call__()
  #edef

  def __len__(self):
    if self.__dict__["_obj"] is None:
      self.__init_obj()
    #fi
    return self.__dict__["_obj"].__len__()
  #edef

  def __eq__(self):
    if self.__dict__["_obj"] is None:
      self.__init_obj()
    #fi
    return self.__dict__["_obj"].__eq__()
  #edef

  def __ne__(self):
    if self.__dict__["_obj"] is None:
      self.__init_obj()
    #fi
    return self.__dict__["_obj"].__ne__()
  #edef

  def __add__(self):
    if self.__dict__["_obj"] is None:
      self.__init_obj()
    #fi
    return self.__dict__["_obj"].__add__()
  #edef

  def __radd__(self):
    if self.__dict__["_obj"] is None:
      self.__init_obj()
    #fi
    return self.__dict__["_obj"].__radd__()
  #edef

  def __getitem__(self, key):
    if self.__dict__["_obj"] is None:
      self.__init_obj()
    #fi
    return self.__dict__["_obj"].__getitem__(key)
  #edef

  def __setitem__(self, key, value):
    if self.__dict__["_obj"] is None:
      self.__init_obj()
    #fi
    return self.__dict__["_obj"].__setitem__(key, value)
  #edef

  def __iter__(self):
    if self.__dict__["_obj"] is None:
      self.__init_obj()
    #fi
    return self.__dict__["_obj"].__iter__()
  #edef

  def __next__(self):
    if self.__dict__["_obj"] is None:
      self.__init_obj()
    #fi
    return self.__dict__["_obj"].__next__()
  #edef

  @property
  def lazyInitialized(self):
    return (self.__dict__["_obj"] is not None)
  #edef

#eclass

class LazyObject(object):
  def __new__(cls, *params, **kwargs):
    return LazyProxy(cls, *params, **kwargs)
  #edef
#eclass

#can be used as such:
#class A(LazyObject): # classes meant to be lazy loaded are derived from LazyInit
#    def __init__(self, x, y):
#        print("Init A")
#        self.x=14+x
#        print("Init B")
#        self.y=15+y
#    
#
#a=A(1,2)
#print("Go")
#print("x=", a.x)
#print("x=", a.y)
#print("cls=", a.__dict__["_cls"])
