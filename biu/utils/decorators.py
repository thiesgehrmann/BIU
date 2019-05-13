import functools
import inspect
import warnings

##########################################################################

string_types = (type(b''), type(u''))

##########################################################################

def deprecated(reason):
    """
    This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used.
    
    Taken from: https://stackoverflow.com/a/40301488
    """

    if isinstance(reason, string_types):

        # The @deprecated is used with a 'reason'.
        #
        # .. code-block:: python
        #
        #    @deprecated("please, use another function")
        #    def old_function(x, y):
        #      pass

        def decorator(func1):

            if inspect.isclass(func1):
                fmt1 = "Call to deprecated class {name} ({reason})."
            else:
                fmt1 = "Call to deprecated function {name} ({reason})."
            #fi

            @functools.wraps(func1)
            def new_func1(*args, **kwargs):
                warnings.simplefilter('always', DeprecationWarning)
                warnings.warn(
                    fmt1.format(name=func1.__name__, reason=reason),
                    category=DeprecationWarning,
                    stacklevel=2
                )
                warnings.simplefilter('default', DeprecationWarning)
                return func1(*args, **kwargs)
            #edef

            return new_func1
        #edef

        return decorator

    elif inspect.isclass(reason) or inspect.isfunction(reason):

        # The @deprecated is used without any 'reason'.
        #
        # .. code-block:: python
        #
        #    @deprecated
        #    def old_function(x, y):
        #      pass

        func2 = reason

        if inspect.isclass(func2):
            fmt2 = "Call to deprecated class {name}."
        else:
            fmt2 = "Call to deprecated function {name}."
        #fi

        @functools.wraps(func2)
        def new_func2(*args, **kwargs):
            warnings.simplefilter('always', DeprecationWarning)
            warnings.warn(
                fmt2.format(name=func2.__name__),
                category=DeprecationWarning,
                stacklevel=2
            )
            warnings.simplefilter('default', DeprecationWarning)
            return func2(*args, **kwargs)
        #edef

        return new_func2
    else:
        raise TypeError(repr(type(reason)))
    #fi
#edef

##################################################################

from collections import namedtuple

class class_or_instance_method(object):
    """
    Defines a decorator that provides a class instance object if a function is called from a class instance,
    or a class object if it is called as a classmethod.
    
    e.g.
    
    class test(object):
        def __init__(self):
            self.a = 10
        @biu.utils.decorators.class_or_instance_method
        def coi(obj, a=0):
            if obj.is_instance:
                print(obj.self.a)
            else:
                print(a)
                
    NOTE: the format of the coi() function:
        def coi(obj,        # A named tuple (see below)
                *pargs,      # Additional arguments to the function
                **kwargs     # Additional, named arguments to the function)
    
    """
    def __init__(self, f):
        self.f = f
        self._nt = namedtuple("class_or_instance", ["self", "is_instance"])
    #edef

    def __get__(self, instance, owner):
        obj = self._nt(instance, True) if instance is not None else self._nt(owner, False)

        def newfunc(*args, **kwargs):
            return self.f(obj, *args, **kwargs)
        #edef
        return newfunc
    #edef
#eclass

##################################################################