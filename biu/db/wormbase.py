from ..structures import Dataset2
from .. import formats
from .. import utils

requests = utils.py.loadExternalModule("requests")

class WormBase(Dataset2):
    def __init__(self, *pargs, **kwargs):
        """
        Initialize the WormBaseRest data structure.
        
        version: Only one version exists at the moment.
        *pargs, **kwargs. See arguments for biu.structures.Dataset2
        
        """
        super(WormBase, self).__init__("WormBaseRest", *pargs, **kwargs)
        
        self._obj.add_file('queries.sqlite', utils.Acquire2().touch('queries.sqlite'))
        self._obj.register("queries", ["queries.sqlite"], lambda f: formats.SQLDict(f["queries.sqlite"]))

    #edef
    
    def query(self, Type, Class, obj_id, widget_or_field):
        """
        See: https://wormbase.org/about/userguide/for_developers#3--10
        and: http://rest.wormbase.org/index.html
        """
        url = "http://rest.wormbase.org/rest/%s/%s/%s/%s" % (Type, Class, obj_id, widget_or_field)
        if url not in self.queries:
            res = requests.get(url)
            if res.status_code == 200:
                self.queries[url] = res.json()
            #fi
        #fi
        return self.queries[url]
    #edef
    
    def rnai_gene(self, rnai):
        try:
            return wbr.query("widget", "rnai", rnai, "overview")['fields']['targets']['data'][0]['gene']['id']
        except Exception as e:
            return None
        #edef
    #edef
    
#eclass