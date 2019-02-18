from .. import utils
from ..structures import Dataset2

pd = utils.py.loadExternalModule("pandas")

class Iris(Dataset2):
    """
    The IRIS Dataset
    """
    def __init__(self, *pargs, **kwargs):
        super(Iris, self).__init__("Iris", *pargs, **kwargs)

        file = utils.Acquire2().curl('https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data')

        self._obj.add_file("iris.tsv", file)
        self._obj.register("iris", ["iris.tsv"],
                           lambda x: pd.read_csv(x["iris.tsv"], index_col=False, names=['a','b','c','d','class']),
                           docstring="An Pandas DataFrame of the IRIS data")
        
        # An example of how to make a more complicated loading scheme
        def iris_mod_func(d):
            ir = pd.read_csv(d["iris.tsv"], index_col=False, names=['a','b','c','d','class'])
            ir['e'] = ir[['a','b','c','d']].sum(axis=1)
            return ir
        #edef
        self._obj.register("iris_mod", [], iris_mod_func,
                           docstring="An Pandas DataFrame with modified the IRIS data")
        
        # An example of how to prevent loading a file unnecessarily...
        def iris_mod2_func(d):
            ir = self.iris.copy()
            ir['f'] = ir[['a','b','c','d']].prod(axis=1)
            return ir
        #edef
        self._obj.register("iris_mod2", [], iris_mod2_func,
                           docstring="An Pandas DataFrame with modified the IRIS data")
        
        self._add_str_func(lambda x: "The IRIS flower petal dataset.")
    #edef

#eclass
