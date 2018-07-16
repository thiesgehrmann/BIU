from ..structures import Pipeline
from .. import formats
from .. import utils
from .. import db

import inspect, os
import datetime
import csv

pd = utils.py.loadExternalModule("pandas")

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
      self.__output['perfamlod'] = res
    #fi

    return self.__output['perfamlod']
  #edef

  #############################################################################

  def peaks(self, pos="order"):
    """ Get peaks in either the order of markers (default), or in nucleotide positions (pos='nt_position')"""
    linkageRes = self.linkage

    if pos == 'order':
        order = linkageRes.groupby("CHR").POS.transform(lambda v: [ p[0] for p in sorted(enumerate([ float(i) for i in v]), key=lambda x: x[1]) ])
        linkageRes["order"] = order
        print
    #fi
        
    groups = linkageRes.groupby("cluster").agg(["min", 'max'])[["CHR", pos]].iterrows()
    linkagePeaks = [ (row[0], (row[1].CHR.min(),  row[1][pos].min(), row[1][pos].max()))
                for row in groups
                 if row[0] > 0]
        
    return dict(linkagePeaks)
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
