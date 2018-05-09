from .. import utils

import errno
import re
import os
import csv

###############################################################################
# TO MATCH THE MERLIN documentation

class Individual(object):
  __slots__ = [ '__famID', '__indivID', '__fatherID', '__motherID', '__gender', '__features' ]
  def __init__(self, famID, indivID, fatherID, motherID, gender, features):
    self.__famID    = str(famID)
    self.__indivID  = str(indivID)
    self.__fatherID = fatherID
    self.__motherID = motherID
    self.__gender   = gender
    self.__features = features
  #edef

  def __str__(self):
    dstr = "Pedigree Individual\n"
    dstr += " Family ID: %s\n" % self.__famID
    dstr += " Individual ID: %s\n" % self.__indivID
    dstr += " Mother/Father ID: %s/%s\n" % (self.__motherID, self.__fatherID)
    dstr += " Gender: %s\n" % self.__gender
    dstr += " Affection status:\n"
    for affection in self.__features.affections:
      dstr += "  %s : %s" % (self.__features.getFeatureName(affection), self.__features[affection])
    #efor

    return dstr
  #edef

  def toRow(self):
    row = [self.__famID, self.__indivID, 
           self.__fatherID if self.__fatherID is not None else '0',
           self.__motherID if self.__motherID is not None else '0',
           self.__gender] + self.__features.raw
    return row
  #edef

  @staticmethod
  def fromRow(row, datFormat):
    row = [ e.strip() for e in row ]
    famID   = str(row[0])
    indivID = str(row[1])
    father  = None if row[2] == '0' else row[2]
    mother  = None if row[3] == '0' else row[3]
    gender  = 'm' if int(row[4]) == 1 else 'f'

    features = PhenoGenotype(row[5:], datFormat)
    return Individual(famID, indivID, father, mother, gender, features)
  #edef

  @property
  def famID(self):
    return self.__famID
  #edef

  @property
  def ID(self):
    return self.__indivID
  #edef

  def setID(self, newID):
    self.__indivID = newID
  #edef

  @property
  def fatherID(self):
    return self.__fatherID
  #edef

  def setFather(self, newID):
    self.__fatherID = newID
  #edef

  @property
  def motherID(self):
    return self.__motherID
  #edef

  def setMother(self, newID):
    self.__motherID = newID
  #edef

  @property
  def gender(self):
    return self.__gender
  #edef

  @property
  def features(self):
    return self.__features
  #edef

  @property
  def isFounder(self):
    return (self.fatherID is None) and (self.motherID is None)
  #edef

  def getFeature(self, feature):
    return self.__features[feature]
  #edef

  def setFeature(self, feature, value):
    self.__features[feature] = value
  #edef

  def copy(self, datFormat):
    features = self.__features.copy(datFormat)
    return Individual(self.famID, self.__indivID, self.__fatherID, self.motherID, self.gender, features)
  #edef

#eclass

###############################################################################

class PhenoGenotype(object):
  __slots__ = [ '__row', '__datFormat' ]
  def __init__(self, row, datFormat):

    self.__row = row
    self.__datFormat = datFormat

    nFeaturesRow = len(row)
    nFeaturesDat = len(datFormat)

    fixedRow = { datFormat.getFeatureName(i) : datFormat.emptyValueOfField(i) for i in range(nFeaturesDat) }
    rowIndex = 0
    datIndex = 0
    while datIndex < nFeaturesDat:
      if rowIndex >= nFeaturesRow:
        break
      #fi
      featureName = datFormat.getFeatureName(datIndex)
      featureType = self.__datFormat.getFeatureType(datIndex)
      if (featureType == 'm') and ('/' not in row[rowIndex]):
        fixedRow[featureName] = '%s/%s' % (row[rowIndex], row[rowIndex+1])
        rowIndex += 2
        datIndex += 1
      else:
        fixedRow[featureName] = row[rowIndex]
        rowIndex += 1
        datIndex += 1
      #fi
    #ewhile

    self.__row = fixedRow

    if rowIndex < nFeaturesRow:
      utils.error("More provided features than in DAT file. IGNORED!")
    elif datIndex < nFeaturesDat:
      utils.error("More DAT features than provided. Filling with unknown values. Verify your PED/DAT file.")
    #fi
  #edef

  @property
  def raw(self):
    row = [ str(self.__row[featureName]) for (featureType, featureName) in self.__datFormat ]
    return row
  #edef 

  def __getitem__(self, identifier):
    if isinstance(identifier, int):
      fieldName = self.__datFormat.getFeatureName(identifier)
      if fieldName is None:
        return None
      #fi
      return self.__row[fieldName]
    elif isinstance(identifier, str):
      return self.__row[identifier]
    #fi
    return None
  #edef

  def __len__(self):
    return len(self.__row) if self.__row is not None else 0
  #edef

  def __setitem__(self, key, value):
    if isinstance(key, int):
      fieldName = self.__datFormat.getFeatureName(key)
      if fieldName is None:
        return None
      #fi
      self.__row[fieldName] = value
    elif isinstance(key, str):
      self.__row[key] = value
    #fi
    return None
  #edef

  def getFeatureName(self, identifier):
    return self.__datFormat[identifier][2]
  #edef

  @property
  def affections(self):
    return self.__datFormat.affections
  #edef

  @property
  def covariates(self):
    return self.__datFormat.covariates
  #edef

  @property
  def traits(self):
    return self.__datFormat.traits
  #edef

  @property
  def markers(self):
    return self.__datFormat.markers
  #edef

  def emptyValueOfField(self, fieldID):
    return self.__datFormat.emptyValueOfField(self, fieldID)
  #edef

  def __str__(self):
    dstr = "Genotype and Phenotype Object\n"
    for fieldName in self.__row:
      dstr += " %s: %s\n" % (fieldName, self.__row[fieldName])
    #efor
    return dstr
  #edef

  def copy(self, datFormat):
    return PhenoGenotype(self.raw, datFormat)
  #edef

#eclass

###############################################################################

class Family(object):
  slots = [ '__famID', '__members', '__datFormat' ]
  def __init__(self, famID, datFormat, members=[]):
    self.__famID   = str(famID)
    self.__datFormat = datFormat
    if isinstance(members, dict):
      self.__members = members
    else:
      self.__members = { member.ID : member for member in members }
    #fi
  #edef

  def __str__(self):
    dstr = ""
    dstr += "Pedigree Family\n"
    dstr += " Members: %d\n" % len(self.__members)
    dstr += " Founders: %d\n" % self.nFounders
    return dstr
  #edef

  def __contains__(self, memberID):
    return str(memberID) in self.__members
  #fi

  @property
  def famID(self):
    return self.__famID
  #edef

  def add(self, individual):
    if isinstance(individual, Individual):
      if individual.ID in self.__members:
        utils.dbm("Overwriting individual '%s' in family '%s'" % (individual.ID, self.famID))
      #fi
      self.__members[individual.ID] = individual
    else:
      utils.error("Cannot add this object to family")
    #fi
  #edef

  def newMember(self, indivID, fatherID, motherID, gender):
    if isinstance(fatherID, Individual):
      if fatherID.famID != self.famID:
        utils.error("Family ID of father (%s) is not this family: '%s'" %(father.famID, self.famID))
        return None
      #fi
      fatherID = fatherID.ID
    elif (fatherID is not None) and (fatherID not in self.__members):
      utils.warning("FatherID '%s' is not present in this family." % fatherID)
    #fi
    if isinstance(motherID, Individual):
      if motherID.famID != self.famID:
        utils.error("Family ID of mother (%s) is not this family: '%s'" %(mother.famID, self.famID))
        return None
      #fi
      motherID = motherID.ID
    elif (motherID is not None) and (motherID not in self.__members):
      utils.warning("MotherID '%s' is not present in this family." % motherID)
    #fi

    features = PhenoGenotype([], self.__datFormat)
    newPerson = Individual(self.famID, indivID, fatherID, motherID, gender, features)
    self.add(newPerson)
    return newPerson
  #edef

  def delMember(self, memberID):
    memberID = str(memberID)
    if memberID in self.__members:
      del self.__members[memberID]
      for member in self.__members:
        if self.__members[member].motherID == memberID:
          self.__members[member].setMother(None)
        elif self.__members[member].fatherID == memberID:
          self.__members[member].setFather(None)
        #fi
      #efor
    #fi
  #edef

  def __len__(self):
    return len(self.__members)
  #edef

  def __iadd__(self, individual):
    self.add(individual)
  #edef

  @property
  def nFounders(self):
    return len([m for m in self.__members.values() if m.isFounder])
  #edef

  @property
  def members(self):
    return self.__members
  #edef

  def __iter__(self):
    return self.__members.__iter__()
  #edef

  def __getitem__(self, key):
    return self.members[key]
  #edef

  def changeMemberID(self, currentID, newID):
    if currentID not in self.__members:
      utils.warning("Cannot change name of '%s', no such member." % currentID)
      return
    #fi
    self.__members[newID] = self.__members[currentID]
    self.__members[newID].setID(newID)
    del self.__members[currentID]

    for memberID in self.__members:
      if memberID != newID:
        if self.__members[memberID].motherID == currentID:
          self.__members[memberID].setMother(newID)
        elif self.__members[memberID].fatherID == currentID:
          self.__members[memberID].setFather(newID)
        #fi
      #fi
    #efor
  #edef

  def copy(self, datFormat):
    members = { memberID: self.__members[memberID].copy(datFormat) for memberID in self.__members }
    return Family(self.famID, datFormat, members)
  #edef

#eclass
    
###############################################################################    

class DAT(object):

  __slots__ = [ '__fields', '__mask', '__types', '__names', '__fileName' ]

  def __init__(self, data, delimiter=' ', quotechar='#', mask=None):

    self.__fields = None
    self.__types  = None
    self.__names  = None
    self.__fileName = None
    self.__mask   = None

    fields = []
    if isinstance(data, str):
      with open(data, 'r') as ifd:
        reader = csv.reader(ifd, delimiter=delimiter, quotechar=quotechar)
        for row in reader:
          row = [ col for col in row if col != "" ]
          if len(row) != 2:
            continue
          #fi
          fieldType  = row[0].lower()
          fieldValue = row[1]
          if fieldType == 'e':
            break
          #fi
          if fieldType not in [ 'a', 'c', 't', 'm', 's' ]:
            utils.warning("Unknown fieldType '%s'" % fieldType)
            fieldType = 's'
          #fi
          fields.append((fieldType, fieldValue))
        #efor
      #ewith
      self.__fields = fields
      self.__fileName = data
    else:
      self.__fields = data
    #fi

    if (mask is None) or len(mask) != len(self.__fields):
      self.__mask = [ False for field in self.__fields ]
    else:
      self.__mask = mask
    #fi
    self.__indexFields()
  #edef

  def __indexFields(self):
    self.__types = {}
    self.__names = {}
    for i, (fieldType, fieldName) in enumerate(self.__fields):
      fieldType = fieldType.lower()
      if fieldType not in self.__types:
        self.__types[fieldType] = []
      #fi
      self.__types[fieldType].append(i)
      self.__names[fieldName] = i
    #efor

  def add(self, fieldType, fieldName):
    if fieldName in self.__names:
      utils.error("Field '%s' already exists." % fieldName)
      return
    #fi
    fieldType = fieldType.lower()
    if fieldType not in 'actms':
      utils.error("'%s' not a valid field type." % fieldType)
    else:
      self.__fields.append((fieldType, fieldName))
      self.__indexFields()
      self.__mask.append(False)
    #fi
  #edef

  @property
  def affections(self):
    return self.__types['a'] if 'a' in self.__types else []
  #edef

  def emptyFeatures(self):
    return [ self.emptyValueOfField(fieldType) for (fieldType, fieldName) in self.__fields ]
  #edef

  @property
  def covariates(self):
    return self.__types['c'] if 'c' in self.__types else []
  #edef

  @property
  def traits(self):
    return self.__types['t'] if 't' in self.__types else []
  #edef

  @property
  def markers(self):
    return self.__types['m'] if 'm' in self.__types else []
  #edef

  def __getitem__(self, identifier):
    if isinstance(identifier, int):
      fieldType, fieldName = self.__fields[identifier]
      return (identifier, fieldType, fieldName)
    elif isinstance(identifier, str):
      fieldID = self.__names[identifier]
      fieldType, fieldName = self.__fields[fieldID]
      return (fieldID, fieldType, fieldName)
    else:
      return None
    #fi
  #edef

  def __len__(self):
    return len(self.__fields)

  def __iter__(self):
    return self.__fields.__iter__()
  #edef

  def __contains__(self, field):
    return field in self.__names
  #edef

  def getFeatureType(self, field):
    fieldID, fieldType, fieldName = self[field]
    return fieldType
  #eidef

  def setFeatureType(self, field, newFieldType):
    fieldID, fieldType, fieldName = self[field]
    newFieldType = newFieldType.lower()
    if newFieldType in [ 'a', 'c', 't', 'm', 's' ]:
      self.__fields[fieldID] = (newFieldType, fieldName)
    else:
      utils.error("Invalid feature type '%s'." % newFieldType)
  #edef

  def getFeatureName(self, field):
    fieldID, fieldType, fieldName = self[field]
    return fieldName
  #eidef

  def getFeatureID(self, fieldName):
    fieldID, fieldType, fieldName = self[fieldName]
    return fieldID
  #edef

  def maskFeature(self, featureName):
    featureID, featureType, featureName = self[featureName]
    self.__mask[featureID] = True
  #edef

  def unmaskFeature(self, featureName, newFeatureType=None):
    featureID, featureType, featureName = self[featureName]
    self.__mask[featureID] = False
    if newFeatureType is not None:
      self.setFeatureType(featureName, newFeatureType)
    self.__mask[featureID] = False
  #edef

  def detectType(self, value):
    if ('/' in value) or ( ' ' in value):
      return 'm'
    elif value.upper() == 'X':
      return 't'
    elif value in [ '0', '1', '2' ]:
      return 'a'
    else:
      return 's'
    #fi
  #edef

  def emptyValueOfField(self, fieldID):
    fieldID, fieldType, fieldName = self[fieldID]

    if fieldType == 'a':
      return 'X'
    elif fieldType == 'c':
      return 'X'
    elif fieldType == 't':
      return 'X'
    elif fieldType == 'm':
      return '0/0'
    elif fieldType == 's':
      return 'X'
    else:
      return '0'
    #fi
  #edef

  def write(self, fileName):
    with open(fileName, "w") as ofd:
      for i, (fieldType, fieldName) in enumerate(self.__fields):
        fieldType = 'S' if self.__mask[i] else fieldType.upper()
        ofd.write("%s\t%s\n" % (fieldType, fieldName))
      #efor
      ofd.write("E\tEND-OF-DATA\n")
      self.__fileName = fileName
    #ewith
  #edef

  def copy(self):
    return DAT([ (fieldType, fieldName) for (fieldType, fieldName) in self.__fields], mask=self.__mask)
  #edef

#eclass

###############################################################################

class PED(object):

  __slots__ = [ 'families', '__fileName', '__datFile', '__datFormat' ]

  def __init__(self, data, datFile=None, datFormat=None, **kwargs):

    self.families  = {}
    self.__fileName  = None
    self.__datFile   = datFile
    self.__datFormat = None

    if datFormat is not None:
      self.__datFormat = datFormat
    elif isinstance(datFile, str):
      self.__datFormat = DAT(datFile)
      self.__datFile = datFile
    else:
      utils.warning("You must provide a DAT file to match. I will try to guess them!")
      self.__datFormat = DAT([])
      self.__datFile = None
    #fi

    if isinstance(data, str):
      self.families = PED.fromFile(data, self.__datFormat)
      self.__fileName = data
    elif isinstance(data, dict):
      self.families = data
    elif isinstance(data, Family):
      self.families = { data.famID : data }
    else:
      self.families = { f.famID for f in data }
    #fi
  #edef

  def __contains__(self, famID):
    return str(famID) in self.families
  #edef

  def __getitem__(self, famID):
    famID = str(famID)
    if famID in self.families:
      return self.families[famID]
    #fi
    return None
  #edef

  def __delitem__(self, famID):
    famID = str(famID)
    if famID in self.families:
      del self.families[famID]
    #fi
  #edef

  def subset(self, famIDs):
    datFormatCopy = self.__datFormat.copy()
    if isinstance(famIDs, str) or isinstance(famIDs, int):
      famIDs = [ famIDs ]
    #fi
    famIDs = [ str(famID) for famID in famIDs ]
    retFams = { famID : self.families[famID].copy(datFormatCopy) for famID in famIDs if (famID in self.families)}
    return PED(data = retFams, datFormat=datFormatCopy)
  #edef

  def emptyFeatures(self):
    self.__datFormat.emptyFeatures()
  #edef

  #def delFeature(self, featureName)
  #def renameFeature(self, featureName)
  def maskFeature(self, feature):
    self.__datFormat.maskFeature(feature)
  #edef
    
  def unmaskFeature(self, feature, newFeatureType=None):
    self.__datFormat.unmaskFeature(feature, newFeatureType)
  #edef


  def addFeature(self, featureType, featureName, defaultValue=None):
    self.__datFormat.add(featureType, featureName)
    emptyValue = self.__datFormat.emptyValueOfField(featureName) if defaultValue is None else defaultValue

    for famID in self:
      for memberID in self[famID]:
        self[famID][memberID].setFeature(featureName, emptyValue)
      #efor
    #efor
  #edef

  def getFeature(self, featureName):
    values = {}
    if featureName not in self.__datFormat:
      utils.error("Feature '%s' doesn't exist." % featureName)
      return values
    #fi

    for famID in self.families:
      family = self[famID]
      for memberID in family:
        values[(famID, memberID)] = family[memberID].getFeature(featureName)
      #efor
    #efor
    return values
  #edef

  def __iter__(self):
    return self.families.__iter__()
  #edef

  def newFamily(self, famID):
    famID = str(famID)
    if famID in self.families:
      utils.warning("Overwriting family '%s'." % famID)
    #fi
    self.families[famID] = Family(famID, self.__datFormat, [])
    return self.families[famID]
  #edef

  @staticmethod
  def fromFile(fileName, datFormat, delimiter='\t', quotechar='#', nrows=None, **kwargs):
    irow = 0
    families = {}
    with open(fileName, 'r') as ifd:
      rePattern = re.compile(r"[\s]+")
      
      for line in ifd:
        row = [ col for col in rePattern.split(line) if col != "" ]
        irow = irow + 1
        if (nrows is not None) and (irow > nrows):
          utils.dbm("I will; break here")
          break
        #fi
        if len(row) < 5:
          break
        #fi
        
        indiv = Individual.fromRow(row, datFormat)
        if indiv.famID not in families:
          families[indiv.famID] = Family(indiv.famID, datFormat, [indiv])
        else:
          families[indiv.famID].add(indiv)
        #fi
      #efor
    #ewith
    return families
  #edef

  @property
  def nFeatures(self):
    return len(self.__datFormat)
  #edef

  def __str__(self):
    dstr  = "PED object\n"
    dstr += " Where: %s\n" % (self.__fileName if self.__fileName is not None else hex(id(self)))
    dstr += " DAT file: %s\n" % (self.__datFile if self.__datFile is not None else hex(id(self.__datFormat)))
    dstr += " Families: %d\n" % len(self.families)
    dstr += "  Founders: %d\n" % sum([ self.families[famID].nFounders for famID in self.families ])
    dstr += "  Total: %d\n" % ( sum([ len(self.families[famID]) for famID in self.families ]) )
    dstr += " Features: %d\n" % (self.nFeatures)
    dstr += "  Affections: %d\n" % len(self.__datFormat.affections)
    dstr += "  Covariates: %d\n" % len(self.__datFormat.covariates)
    dstr += "  Traits: %d\n" % len(self.__datFormat.traits)
    dstr += "  Markers: %d\n" % len(self.__datFormat.markers)

    return dstr
  #edef

  def write(self, fileName, datFileName=None):

    if datFileName is None:
      if fileName[-3:].lower() == 'ped':
        datFileName = fileName[:-3] + 'dat'
      else:
        datFileName = fileName + '.dat'
      #fi
    #fi

    with open(fileName, 'w') as ofd:
      for famID in self.families:
        for memberID in self.families[famID].members:
          member = self.families[famID][memberID]
          row = member.toRow()
          ofd.write('\t'.join(row) + '\n')
        #efor
      #efor
      self.__fileName = fileName
    #ewith

    self.__datFormat.write(datFileName)
  #edef

  standardValues = {
    'affected_true' : 2,
    'affected_false' : 1,
    'affected_unknown' : 0,
  }

#eclass
