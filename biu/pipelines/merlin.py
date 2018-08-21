from ..structures import Pipeline
from .. import formats
from .. import processing
from .. import utils
from .. import db

import inspect, os
import datetime
import csv

pd = utils.py.loadExternalModule("pandas")
plt = utils.py.loadExternalModule("matplotlib.pylab")
patches = utils.py.loadExternalModule("matplotlib.patches")
sns = utils.py.loadExternalModule("seaborn")
np = utils.py.loadExternalModule("numpy")

####################################################################

snakemakeFile =  os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + '/merlin/Snakefile'

class Merlin(Pipeline):

  __defaultConfig = {}

  def __init__(self, pedObject=None, mapFile=None, pedFile=None, datFile=None, ibd=True, npl=True, pairs=False, perFamily=True, heritability=True, merlinOptions="",
               config={}, dbsnp=db.DBSNP(), lodThresh=3.0, lodDrop=1.0, **kwargs):

    Pipeline.__init__(self, snakemakeFile, {**self.__defaultConfig, **config}, **kwargs)

    pedFileName = None
    datFileName = None
    if (pedFile is not None) and (datFile is not None):
      pedFileName = os.path.abspath(pedFile)
      datFileName = os.path.abspath(datFile)
    elif pedObject is not None:
      pedFileName, exists = self._generateInputFileName([ str(id(pedObject)), datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"), 'ped' ])
      datFileName, exists = self._generateInputFileName([ str(id(pedObject)), datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"), 'dat' ])
      pedObject.write(pedFileName, datFileName)
    else:
      utils.msg.error("You must specify a PED object, or PED and DAT files")
      return self
    #fi

    if mapFile is None:
      utils.msg.error("You must specify a MAP file.")
      return self
    #fi
      
    if (pairs and npl) or not(pairs or npl):
      utils.msg.error("You must choose either --npl or --pairs")
      return self
    #fi

    self.setConfig(pedfile=pedFileName, datfile=datFileName, mapfile=mapFile,
                   ibd=(1 if ibd else 0),
                   npl=(1 if npl else 0),
                   pairs=(1 if pairs else 0),
                   perfam=(1 if perFamily else 0),
                   heritability=(1 if heritability else 0),
                   merlinOptions=merlinOptions)

    if self.autorun:
      self.run(["output"])
    #fi
    self.__output = {}

    self.__dbsnp     = dbsnp
    self.__lodThresh = lodThresh
    self.__lodDrop   = lodDrop
  #edef

  #############################################################################

  def _outputFileNames(self):
    files = {}
    files["linkage"] = '%s/linkage-nonparametric.tbl' % self.config["outdir"]
    files["ibd"] = '%s/linkage.ibd' % self.config["outdir"]
    files["perfamlod"] = '%s/linkage.lod' % self.config["outdir"] 
    return files
  #edef

  @property
  def ibd(self):
    if not self.success:
      return None
    #fi

    if 'ibd' not in self.__output:
      self.__output['ibd'] = self.ibdReader(self._outputFileNames()['ibd'], self.__dbsnp)
    #fi

    return self.__output['ibd']
  #edef

  @property
  def linkage(self):
    if not self.success:
      return None
    #fi

    if 'linkage' not in self.__output:
      res = pd.read_csv(self._outputFileNames()['linkage'], sep='\t')
       # Remove rows with NA fields
      res = res[(res == "na").apply(lambda row: not any(row), axis=1)]
      res["POS"] = res["POS"].astype("float")
      res["nt_position"] = res.LABEL.apply(lambda rs: self.__dbsnp[rs][1])
      res['cluster'] = self.annotateLDRegions(res, self.__lodThresh, self.__lodDrop)
      self.__output['linkage'] = res
    #fi

    return self.__output['linkage']
  #edef

  @property
  def perfamlod(self):
    if not self.success:
      return None
    #fi

    if 'perfamlod' not in self.__output:
      res = pd.read_csv(self._outputFileNames()['perfamlod'], delim_whitespace=True,
                        skiprows=1, header=None,
                        names=("family", "trait", "analysis", "location", "zscore", "plod", "delta", "lod"))
      res["chr"] = res.location.apply(lambda rs: self.__dbsnp[rs][0])
      res["nt_position"] = res.location.apply(lambda rs: self.__dbsnp[rs][1])
      res = res.sort_values(["chr", "nt_position"])
      self.__output['perfamlod'] = res
    #fi

    return self.__output['perfamlod']
  #edef

  #############################################################################

  def peaks(self, pos="order"):
    """
    peaks: Get peak regions in either the order of markers (default), or in nucleotide positions (pos='nt_position')
    Inputs: pos: [order|POS|nt_position]
                  order: Order of markers
                  POS: POS of marker
                  nt_position : nucleotide position
    Outputs: Dictionary of peak_ID : (peak_chromosome, peak_start, peak_end)
    
    """
    linkageRes = self.linkage

    if pos not in [ 'order', 'POS', 'nt_position' ]:
        raise ValueError
    #fi
    
    if pos == 'order':
        order = linkageRes.groupby("CHR").POS.transform(lambda v: [ p[0] for p in sorted(enumerate([ float(i) for i in v]), key=lambda x: x[1]) ])
        linkageRes["order"] = order
        print
    #fi
        
    groups = linkageRes.groupby("cluster").agg(["min", 'max'])[["CHR", pos]].iterrows()
    linkagePeaks = { row[0] : ( row[1].CHR.min(), row[1][pos].min(), row[1][pos].max() ) for row in groups if row[0] > 0 }
        
    return linkagePeaks
  #edef

  #############################################################################


  def plotLOD(self, chromosomeID=None, axis=None, pos="POS"):
      """
      input:
       * res: Linkage result
       * chromosomeID: Name of chromosome to plot (If none is provided, all are plotted)
       * axis: Which axis to plot on
       * pos: x-axis [POS|nt_position|order]
              POS : marker position
              nt_position: nucleotide position
              order: order of marker in output (integer positions)
      
      output:
       * plt.show()
      """
      
      if pos not in [ 'order', 'POS', 'nt_position' ]:
        raise ValueError
      #fi
        
      res = self.linkage
      
      if chromosomeID is None:
          chromosomes = res.CHR.unique()
          nrows = int(np.ceil(len(chromosomes)/3))
          fig, axes = utils.figure.subplots(figsize=(9*3,3*nrows), ncols=3, nrows=nrows, dpi=500)
          for i, chromosome in enumerate(chromosomes):
              self.plotLOD(chromosome, axes[i], pos=pos)
          #efor
          return plt.show()
      #fi
      
      if axis is None:
          fig, axes = utils.figure.subplots(figsize=(9,3), ncols=1, nrows=1, dpi=500)
          axis = axes[0]
      #fi
  
      
      relRes = res[(res.CHR == chromosomeID)]
      relRes["range"] = range(len(relRes.LOD))
      minLOD, maxLOD = (min(res.LOD), max(res.LOD))
      if pos == "POS":
          positions = relRes.POS
          relRes["range"] = relRes.POS
      elif pos == 'nt_position':
          relRes["range"] = relRes.nt_position
      elif pos == 'order':
          [ p[0] for p in sorted(enumerate(relRes.POS.values), key=lambda x: x[1]) ]
          relRes["range"] = [ p[0] for p in sorted(enumerate(relRes.POS.values), key=lambda x: x[1]) ]
      #fi
      
      # Plot the clusters
      for (chromosome, start, end) in self.peaks(pos=pos).values():
          if chromosome != chromosomeID:
              continue
          #fi
          rect = patches.Rectangle((start, minLOD),end-start, maxLOD - minLOD, linewidth=1,facecolor='r',alpha=.25)
          axis.add_patch(rect)
      #efor
      axis.plot(relRes.range, relRes.LOD)
      axis.plot([0, max(relRes.range)], [3, 3], '-', c='orange')
      axis.plot([0, max(relRes.range)], [2, 2], '-', c='orange')
      axis.plot([0, max(relRes.range)], [0, 0], c='k')
      axis.set_ylim([minLOD, maxLOD])
      axis.set_xlim([0, max(relRes.range)])
      axis.set_title("Chromosome %s" % chromosomeID)
  #edef

  #############################################################################

  def plotFamLOD(self, chromosomeID):
    """
    plotFamLOD: Plot the family LOD contributions for a given chromosome.
    Inputs: ChromosomeID: The ID of the chromosome to plot
    Outputs: Figure
    """
    linkageRes = self.linkage
    famLOD     = self.perfamlod
    
    famLODGroup = famLOD[famLOD.chr == chromosomeID].groupby("family")[["family", "lod", "chr", "nt_position"]]

    D = []
    for i, family in enumerate(famLOD.family.unique()):
        fg = famLODGroup.get_group(family)
        D.append(fg.lod)
    #efor

    # Reorder the elements in the matrix
    D = np.matrix(D)
    order = processing.matrix.order(D, 'correlation', 'complete')
    newD = D[order,:]
    
    # Make the plot
    fig = plt.figure(figsize=(20, 20), dpi=500)

    axes = [ plt.subplot2grid((20,1), (0,0), colspan=1, rowspan=2),
             plt.subplot2grid((20,1), (2,0), colspan=1, rowspan=18) ]

    # Plot the chromosome LOD score, and clusters.
    # Plot using the order, so that we can match with the heatmap
    merlin.plotLOD(chromosomeID, axes[0], pos="order")
    
    # Plot the heatmap
    sns.heatmap(newD, ax=axes[1], cbar=False, cmap="Blues")

    plt.show()
  #edef

  #############################################################################

  def familiesPerPeak(self, lodThreshold=0.3):
    """
    familiesPerPeak: Identify families that have a LOD contribution > than a certain value in each linkage peak
    Inputs: lodThreshold: The LOD threshold to use (default 0.3)
    Outputs: Dictionary of clusterID: [Family IDs contributing to cluster]
    """
    linkagePeaks = self.peaks(pos='nt_position')
    famLOD       = self.perfamlod
    
    clusterFams = {}
    for clusterID, (chromosome, start, end) in sorted(linkagePeaks.items(), key=lambda x: (x[1][1], x[1][2])):
        families = famLOD[(famLOD.chr == chromosome) & (famLOD.lod > lodThreshold) & 
                          (famLOD.nt_position > start) & (famLOD.nt_position < end)].family.values
        clusterFams[clusterID] = set(families)
    #efor
    return clusterFams
  #edef
    
  #############################################################################

  @staticmethod
  def ibdReader(ibdFile, dbsnp=db.DBSNP()):
    return IBDStatus(ibdFile, dbsnp)
  #edef

  @staticmethod
  def annotateLDRegions(linkageRes, lodThresh=3.0, lodDrop=1.0):
      dropLodThresh = lodThresh - lodDrop
      clusterIndex = 0
      inCluster = False
      clusterChr = None
      
      clusterIDs = [ -1 for i in linkageRes.LOD ]
      maxLodPerCluster = {}
      for i, lodScore in enumerate(linkageRes.LOD):
          if clusterChr != clusterChr:
              inCluster = False
          #fi
          
          if inCluster:
              if lodScore < dropLodThresh:
                  inCluster = False
              #fi
              clusterIDs[i] = clusterIndex
          elif lodScore >= dropLodThresh:
              inCluster = True
              clusterIndex = clusterIndex + 1
              maxLodPerCluster[clusterIndex] = lodScore
              clusterIDs[i-1] = clusterIndex
              clusterIDs[i] = clusterIndex
          #fi
          
          if inCluster:
              if maxLodPerCluster[clusterIndex] < lodScore:
                  maxLodPerCluster[clusterIndex] = lodScore
      #efor
      
      invalidClusters = set([ ci for ci in maxLodPerCluster if maxLodPerCluster[ci] < lodThresh])
      clusterIDs = [ -1 if ci in invalidClusters else ci for ci in clusterIDs]
      
      renamedClusterIDs = { cid: (-1 if cid == -1 else i) for (i, cid) in enumerate(sorted(set(clusterIDs))) }
      
      return [ renamedClusterIDs[cid] for cid in clusterIDs ]
  #edef

#eclass

###############################################################################

class IBDStatus(object):

  families = None
  markers  = None
  markerLookup = None

  def __init__(self, fileName, dbsnp):
    self.families = {}
    self.markers = set([])
    self.markerLookup = {}
    with open(fileName, 'r') as ifd:
      reader = csv.reader(ifd, delimiter=' ')
      header = next(reader)
      for row in reader:
        family, person1, person2, marker, _, p0, p1, p2 = row
        self.markers.add(marker)
        if family not in self.families:
          self.families[family] = { '__familymembers' : set([]) }
        #fi
        self.families[family]['__familymembers'].update([person1, person2])
        k = self.__lookupkey(person1, person2)
        if k not in self.families[family]:
          self.families[family][k] = {}
        #fi
        self.families[family][k][marker] = max([(0, p0), (1, p1), (2,p2)], key=lambda x: x[1])[0]
      #efor
    #ewith
    
    for marker in self.markers:
      pos = dbsnp[marker]
      if pos is None:
        continue
      #fi
      chrID, location = pos
      if (chrID is None) or (location is None):
        continue
      #fi
      if chrID not in self.markerLookup:
        self.markerLookup[chrID] = []
      #fi
      self.markerLookup[chrID].append((location, marker))
    #efor
    
    self.markerLookup = { chrID : sorted(self.markerLookup[chrID], key=lambda x: x[0]) for chrID in self.markerLookup}
      
  #edef
  
  def __lookupkey(self, person1, person2):
    return (person1, person2) if person1 < person2 else (person2, person1)
  #edef
          
          
  def lookup(self, famID, person, marker):
    sibs = [ s for s in self.families[famID]['__familymembers'] if s != person ]
    return { sib : self.families[famID][self.__lookupkey(person, sib)].get(marker, 0) for sib in sibs}
  #edef
 
  def neighboringMarkers(self, chromosome, pos):
    chromosome = str(chromosome)
    beforeMarker = None
    afterMarker  = None
    for i, (markerPos, marker) in enumerate(self.markerLookup[chromosome]):
      if markerPos > pos:
        beforeMarker = self.markerLookup[chromosome][i-1][1]
        afterMarker  = marker
        break
      #fi
    #efor
    return beforeMarker, afterMarker
  #edef

  def personIBDAtPos(self, famID, person, chromosome, pos):
    beforeMarker, afterMarker = self.neighboringMarkers(chromosome, pos)
    return self.personIBDAtMarker(famID, person, beforeMarker, afterMarker)
  #edef

  def personIBDAtMarker(self, famID, person, marker1, marker2):
    ibdBefore = self.lookup(famID, person, marker1)
    ibdAfter  = self.lookup(famID, person, marker2)

    return{ sibling : (ibdBefore[sibling], ibdAfter[sibling]) for sibling in ibdBefore }
  #edef

  def __str__(self):
    dstr = "Merlin IBD status object\n"
    dstr += " Families: %d\n" % len(self.families)
    dstr += " Individuals: %d\n" % sum([ len(self.families[f]['__familymembers']) for f in self.families])
    dstr += " Markers: %d\n" % len(self.markers)
    return dstr

#eclass
