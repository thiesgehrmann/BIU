"""
Tools to aid in the exploration of hierarchical clustering

Includes:
    * Improved tree structure
    * Modularity-based hierarchical clustering
    * Dendrogram drawing
    * Cluster coloring

"""

from collections import namedtuple

from .. import utils
from .. import ops

np  = utils.py.loadExternalModule("numpy")
plt = utils.py.loadExternalModule('matplotlib.pylab')

####################################################################

def prepare_tree(T):
    """
    The output of the linkage is confusing, and doesn't contain the leaf nodes.
    We add them here, and also do some pre-processing so that it is easy to keep
    track of the leaf nodes that descend from each internal node (an optimization).
    
    parameters:
    -----------
    T: The output of scipy.cluster.hierarchy.linkage
    
    returns:
    -------
    A tree with leaves included.
    Each node is represented as a namedtuple with attributes:
     * children: List[Integer] (Empty if leaf)
     * parents: Integer (None if root)
     * height: Float
     * leaves: List[Integer]
     * info: Dict of additional tags you can add
         Contains also:
         is_leaf: Boolean
    """
      
    nt = namedtuple("TreeNode", ["children", "parent", "height", "leaves", "info"])

    # The prepared tree
    # Nodes have values (child_left_id, child_right_id, height, leaves, info)
    P = [ nt([], None, 0, [i], dict(is_leaf=True)) for i in range(len(T) + 1) ];
    
    
    for t in T:
        c_left, c_right, height, leaf_count = t;
        c_left  = int(c_left) 
        c_right = int(c_right)
        
        node_leaves = P[c_right].leaves + P[c_left].leaves;
        
        P.append( nt([c_left, c_right], None, height, node_leaves, dict(is_leaf=False)) );
        
        P[c_left]  = nt(P[c_left].children, len(P)-1, P[c_left].height, P[c_left].leaves, P[c_left].info)
        P[c_right] = nt(P[c_right].children, len(P)-1, P[c_right].height, P[c_right].leaves, P[c_right].info)
        
    #efor
    
    return P;
#edef 

####################################################################

def deprepare_tree(T, height=lambda x: x.height):
    """
    Translate a tree back to the original scipy format
    
    parameters:
    -----------
    
    T: The tree from a prepare_tree call
    height: callable.
        A function to define the height of a node (given the node ID)
        default is the original height.
        
    returns:
    --------
    Tree compatible with the scipy.hierarchy tree structure
    """
    P = []
    for n in T:
        if n.info["is_leaf"]:
            continue
        #fi
        if len(n.children) != 2:
            return []
        #fi
        P.append([n.children[0], n.children[1], height(n), len(n.leaves)])
    #efor
    return P
#edef

####################################################################

def modularity_ayroles_module(C, M):
    """
    Calculate the per-module modularity
    
    parameters:
    ----------
    C: numpy.ndarray
        The correlation network between nodes in the tree
    sigma: Float Default:1
        How much to de-exaggerate the correlation scores with (rij-1)/(sigma^2)
    M: List[Integer]
        For which nodes should the modularity be calculated?
    
    returns:
    --------
    Modularity score for this module
    
    """
    mask    = np.zeros(C.shape[0],dtype=bool)
    mask[M] = True
    X = np.sum(np.triu(C[mask,:][:,mask]))
    Y = np.sum(np.triu(C))
    Z = np.sum(C[mask,:][:,])    

    modularity = X/Y - (Z/(2*Y))**2
    
    return modularity
#edef

####################################################################

def modularity_haq_module(C, sigma, M):
    """
    Calculate Modularity with weighted modularity from
    Haq,N.F. et al. (2019) Community structure detection from networks with weighted modularity. Pattern Recognit. Lett., 122, 14–22.
    
    parameters:
    -----------
    C: numpy.ndarray
        The correlation network between nodes in the tree
    sigma: Float Default:1
        How much to de-exaggerate the correlation scores with (rij-1)/(sigma^2)
    M: List[Integer]
        For which nodes should the modularity be calculated?
        
    returns:
    --------
    Modularity score for this module

    """
    n = len(M)
    li = 1 + (2 * np.sum(np.triu(C[M,:][:,M]))) / ( np.exp(1/sigma**2) * n * (n-1)) if n > 1 else 1 
    qi = modularity_ayroles_module(C, M)

    return li * qi
#edef

####################################################################

def modularity_cut(C, T, sigma=1, minsize=20, percentile=95, method='ayroles'):
    """
    Identify cutting points in a tree that satisfy certain modularity conditions
    
    parameters:
    -----------
    C: numpy.ndarray
        The correlation network between nodes in the tree
    T: tree from prepare_tree
    sigma: Float Default:1
        How much to de-exaggerate the correlation scores with (rij-1)/(sigma^2)
    minsize: Integer Default:20
        The minimum size of cluster to return
    percentile: float [0-100] Default:95
        Return only clusters with a modularity score in the top x% percentile
    method: String in { 'ayroles', 'haq' }
        The type of modularity to use:
            ayroles: Modularity
            haq: Weighted modularity
        
    returns:
    --------
    List[Integer]
        The nodes which satisfy the conditions
    """
    
    C = np.exp((np.abs(C) - 1) / (sigma**2))
    
    for node_id in range(len(T)):
        if method.lower() == 'haq':
            T[node_id].info["modularity"] = modularity_haq_module(C, sigma, T[node_id].leaves)
        else:
            T[node_id].info["modularity"] = modularity_ayroles_module(C, T[node_id].leaves)
        #fi
    #efor
    
    modules = [ -1 ]
    
    threshold = np.percentile([ n.info["modularity"] for n in T if not n.info["is_leaf"] ], percentile)
    
    while True:
        next_modules = []
        for node in modules:
            if T[node].info.get("visited",False):
                next_modules.append(node)
                continue
            #fi
            valid_children = [ c for c in T[node].children if
                              (T[c].info["modularity"] >= threshold) & (len(T[c].leaves) > minsize)]
            
            if len(valid_children) > 0:
                next_modules.extend(valid_children)
            else:
                next_modules.append(node)
            #fi
            #print(node, valid_children, next_modules)
            T[node].info["visited"] = True
        #efor
        next_modules = set(next_modules)
        if next_modules == modules:
            break
        else:
            modules = next_modules
        #fi
    #ewhile
    return modules
#edef

####################################################################

def cutpoints_to_partitions(T, CP, order=None):
    """
    Given a set of nodes, derives a cluster assignment for each node
    
    parameters:
    -----------
    T: the tree, from process_tree
    CP: List[Integer]
        The nodes defining the clusters
    order: A dendrogram ordering to inform the re-numbering of clusters
    
    returns:
    --------
    P: List[Integer]
        A list of cluster IDs per leaf node in the tree. 0 indicates no cluster assignment
    """
    P = np.zeros(len([n for n in T if n.info["is_leaf"]]))
    for n in CP:
        P[T[n].leaves] = n
    #efor
    
    if order is not None:
        rename = { n: i+1 for i,n in enumerate(sorted(CP, key=lambda x: np.min(np.where(P[order] == x)))) }
        P = [ rename.get(c,0) for c in P ]
    #fi
    
    return np.array(P)
#edef

####################################################################

def cutpoints_to_colors(T, CP, cmap=None):
    """
    Determine colors for clusters in the dendrogram.
    
    parameters:
    -----------
    T: the tree, from process_tree
    CP: List[Integer]
        The nodes defining the clusters
    cmap: matplotlib color map
        Which colors to use for the clusters
        If None -> plt.get_cmap('Dark2')
    
    returns:
    --------
    C: List[matplotlib color] per node in tree (includes leaves)
    """
    C = [ 'b' for n in T ]
    nleaves = len([n for n in T if n.info["is_leaf"]])
    
    cmap = plt.get_cmap('Dark2') if cmap is None else cmap
    
    clusters = CP
    for i, cluster in enumerate(clusters):
        color = cmap(i % 20)
        nodes = [ cluster ]
        while len(nodes) > 0:
            node = nodes.pop()
            C[node] = cmap(i)
            for c in T[node].children:
                nodes.append(c)
            #efor
        #ewhile
    #efor
    
    return C
#edef

####################################################################

def draw_dend(dend, ax, draw_leaves=False, colors=None):
    """
    Draw a dendrogram, with additional options
    
    parameters:
    -----------
    dend: The dendrogram
        output from scipy dendrogram
    ax: matplotlib axis
        The axis to plot to
    draw_leaves: boolean
        Draw the leaves, or do not?
    colors: List[matplotlib color code]
        The colors to use for the nodes. If None, then use colors in dendrogram. If False, do not draw colors
        
    returns:
        matplotlib axis.
    """
    colors = dend['color_list'] if colors is None else ['k']*len(dend['icoord']) if colors == False else colors
    for X, Y, c in zip(dend['icoord'], dend['dcoord'], colors):
        XY = [ (x,y) for (x,y) in zip(X,Y) if draw_leaves or (y != 0)]
        X, Y = zip(*XY)
        ax.plot(X, Y, c=c)
    #efor
    ax.set_xlim([min(ops.lst.flatten(dend['icoord'])), max(ops.lst.flatten(dend['icoord']))])
    return ax
#edef