from ..structures import Dataset
from ..config import settings as settings
from .. import formats
from .. import utils
from .. import ops
from .. import medical

pd = utils.py.loadExternalModule("pandas")
np = utils.py.loadExternalModule("numpy")

import re

class GOTO(Dataset):

    versions = { "current":
        { "chrs" : [ "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"] }
    }
    version = None
    vcf = None
    
    __samplesPerIndividual = None
    
    def __init__(self, version=list(versions.keys())[0], where="/exports/molepi/GOTO_HRC1.1", localCopy={}):
        fileIndex = self.__genFileIndex(version, where)

        Dataset.__init__(self, fileIndex, localCopy=localCopy)
        self.version = version

        for chrID in self.versions[self.version]["chrs"]:
            self._registerObject('vcf_%s' % chrID, formats.VCF, [ "vcf_%s" % chrID, "vcf_%s_tbi" % chrID ], fileIndex["vcf_%s" % chrID].path, tabix=True)
        #efor
        
        self._registerObject('expr', pd.read_csv, ['expr'], fileIndex['expr'].path, index_col=0)
        self._registerObject('blood', pd.read_csv, ['blood_cov'], fileIndex['blood_cov'].path)
        self._registerObject('muscle', pd.read_csv, ['muscle_cov'], fileIndex['muscle_cov'].path)

        
        self._addStrFunction(lambda s: "Version: %s" % self.version)
    #edef

    def __genFileIndex(self, version, where=None):
        files = {}
        for chrID in self.versions[version]["chrs"]:
            files['vcf_%s' % chrID] = utils.Acquire("%s/chr%s.dose.vcf.gz" % (where, chrID), where=where)
            files['vcf_%s_tbi' % chrID] = utils.Acquire("%s/chr%s.dose.vcf.gz.tbi" % (where, chrID), where=where)
        #efor
        files['expr'] = utils.Acquire('/exports/molepi/tgehrmann/data/GOTO_data/goto_expr.csv')
        files['blood_cov'] = utils.Acquire('/exports/molepi/tgehrmann/data/GOTO_data/datasheet_RNAseq_blood_V2.csv')
        files['muscle_cov'] = utils.Acquire('/exports/molepi/tgehrmann/data/GOTO_data/muscle_QC_covariates_filesv3.csv')

        #files['phen'] = utils.Acquire("/dev/null")

        return files
    #edef
    
    @property
    def blood(self):
        """
        Blood sample covariates.
        Categorical variables are strings, provided as a categorical dtype 
        Continuous variables are provided as floats
        """
        if not self._objectLoaded('blood'):
            self._getObject('blood')['Sample'] = self._getObject('blood').Sample.apply(lambda s: s.replace('-', '_'))
            #self._getObject('blood')['flowcell'] = self._getObject('blood').Sample.apply(lambda s: s.split('_')[1])
            self._getObject('blood')['flowcell_lane'] = self._getObject('blood').Sample.apply(lambda s: s.split('_')[2])

            def castStr(value):
                return str(value) if str(value).strip() != "" else None
            #edef

            def castFloat(value):
                try:
                    if not str(value).strip():
                        return np.nan
                    #fi
                    return float(value)
                except Exception as e:
                    return np.nan
                #etry
            #edef


            FRS = []
            for i, row in self._getObject('blood').iterrows():
                frs = medical.health.indicators.framingham_risk_score(castFloat(row.Sex) == 1, castFloat(row.Age_Baseline),
                                                                      castFloat(row.Cholesterol), castFloat(row.Cholesterol_HDL),
                                                                      castFloat(row.Systolic_Blood_Pressure),
                                                                      castFloat(row.Antihypertensive_medication) == 1,
                                                                      castFloat(row.Current_Smoker) == 1,
                                                                      False)
                FRS.append(frs)
            #edef
            self._getObject('blood')['FRS'] = np.array(FRS)

            discreteCovariates = [ 'IOP2_ID', 'timepoint', 'Sex', 'Status', 'Lipid_lowering_medication',
                                   'Current_Smoker', 'Antihypertensive_medication', 'Date_of_sample_collection',
                                   'intervention', 'nutridrink', 'sampID', 'Labnr', 'Sample', 'studID',
                                   'flowcell', 'Lane', 'Index', 'blood.or.muscle', 'plate',
                                   'RNA_blood_Isolatieseries', 'RNA_Isolation_series', 'flowcell_lane' ]
            continuousCovariates = [ 'Age_Baseline', 'Glucose', 'Cortisol', 'IGF_BP3', 'LN_Insulin', 
                                     'IGF_1', 'IGF1_IGFBP3', 'LN_ALAT', 'Albumin', 'ALP', 'LN_ASAT', 'Cholesterol',
                                     'Cholesterol_HDL', 'Creatinine', 'GFR', 'LN_Triglycerides', 'LN_TSH',
                                     'Uric_acid', 'VitaminD', 'LDL_Cholesterol', 'HDL_Cholesterol', 'DHEA_S', 'fT3', 'LN_Leptin',
                                     'LN_Adiponectin', 'LN_Leptin_Adiponectin', 'LN_FFA', 'VitaminE', 'LN_TNFa', 'LN_IL6',
                                     'ASAT_ALAT', 'BMI', 'Systolic_Blood_Pressure', 'Hb', 'Ht', 'Erythrocytes', 'MCV', 'MCH',
                                     'MCHC', 'Leukocytes', 'Eosinophiles', 'Basophiles', 'Neutrophiles',
                                     'Lymfocytes', 'Monocytes', 'Eosinophil_perc', 'Basophil_perc', 'fT4', 'LN_gGT', 'LN_CRP',
                                     'Neutrophil_perc', 'Lymfocyte_perc', 'Monocyte_perc', 'Thrombocytes',
                                     'Total_reads', 'Mapped_reads', 'mappedreads_perc',
                                     'PFALIGNEDBASES', 'MEDIAN5PRIMEBIAS', 'MEDIAN3PRIMEBIAS',
                                     'MEDIAN5PRIMETO3PRIMEBIAS', 'meaninsertsize', 'standarddeviation',
                                     'medianinsertsize', 'Nanodrop260280', 'passedQC_perc',
                                     'LabonchipRIN', 'Labonchip28S18S', 'totalyieldug', 'FRS' ]
            
            for c in discreteCovariates:
                self._getObject('blood')[c] = self._getObject('blood')[c].map(castStr).astype('category')
            #efor
            
            for c in continuousCovariates:
                self._getObject('blood')[c] = self._getObject('blood')[c].map(castFloat).astype(float)
            #ef r

        #fi
        return self._getObject('blood').set_index('Sample')
    #edef
    
    @property
    def muscle(self):
        if not self._objectLoaded('muscle'):
            self._getObject('muscle')['Sample'] = self._getObject('muscle').Sample.apply(lambda s: s.replace('-', '_'))
        #fi
        return self._getObject('muscle').set_index('Sample')
    #edef
    
    @property
    def individuals(self):
        return list(map(str, ops.lst.uniq(list(self.blood.IOP2_ID.values) + list(self.muscle.IOP2_ID.values))))
    #edef

    @property
    def samplesPerIndividual(self):
        if self.__samplesPerIndividual is None:
            muscle = ops.lst.group(self.muscle[["IOP2_ID", "Sample", "Visitnr"]].values)
            blood  = ops.lst.group(self.blood[["IOP2_ID", "Sample", "timepoint"]].values)
            self.__samplesPerIndividual = { str(p) : { 'blood' : { str(s[2]) : s[1] for s in blood.get(str(p),[]) }, 
                                                       'muscle' : { str(s[2]) : s[1] for s in muscle.get(int(p),[]) } }
                                           for p in self.individuals }
        #fi
        return self.__samplesPerIndividual
    #edef
    
    def getGeneVariants(self, hg, ensemblGeneID, upstream=1000000, downstream=1000000, mafThreshold=None, **queryArgs):
        gffEntry = hg.gff['gene:%s' % ensemblGeneID]
        variants = self.query(gffEntry.seqid, gffEntry.start - upstream, gffEntry.end + downstream, extract='raw', **queryArgs)
        nvars    = len(variants)
        
        if mafThreshold is not None:
            def satisfiesMAF(value):
                maf = min(value, 1-value)
                return maf >= mafThreshold
            #edef
            variants = [ v for v in variants if satisfiesMAF(float(v.INFO.get('AF', 0))) ]
        #fi
        
        vars = formats.VCF(variants)
        
        print("Found %d variants for %s. Reduced to: %d" % (nvars, ensemblGeneID, len(variants)))
        
        return vars
    #edef
        
    
    def getGeneVariantGenotypes(self, *pargs, **kwargs):

        vars = self.getGeneVariants(*pargs, **kwargs)

        matrix = np.zeros((len(vars.samples), len(vars)))
        for i, var in enumerate(vars):
            matrix[:,i] = np.array([ sum([ int(a) for a in re.split(r'[\|\/]', s.data.GT) ]) for s in var.samples ])
        #edef

        genotypes = pd.DataFrame(matrix,
                                 index=[ s.split('_')[0] for s in vars.samples ],
                                 columns=[ 'var_' + formats.VCF.makeIdentifier(v) for v in vars ])

        return genotypes
    #edef
    
    #def getGeneExprInSample(self, ensemblGeneID):
        

    def query(self, chrID, start, end, **kwargs):
        chrID = str(chrID)
        oname = "vcf_%s" % chrID
        if self._objectExists(oname):
            return self._getObject(oname).query(chrID, start, end, **kwargs)
        else:
            utils.error("Could not find chromosome '%s'" % chrID)
            return iter(())
        #fi
    #edef
    
    def queryRegions(self, regions, extract=None, **kwargs):
        R = []
        for (c,s,e) in regions:
            R.extend(list(self.query(c,s,e, **kwargs)))
        #efor
        return formats.VCF.extract(R, extract=extract)
    #edef

    def getVar(self, chromosome, *pargs, **kwargs):
        chromosome = str(chromosome)
        oname = "vcf_%s" % chromosome
        if not self._objectExists(oname):
            utils.error("Could not find chromosome '%s'" % chromosome)
            return None
        #fi
        return self._getObject(oname).getVar(chromosome, *pargs, **kwargs)
    #edef

    def whoHas(self, chromosome, *pargs, **kwargs):
        chromosome = str(chromosome)
        oname = "vcf_%s" % chromosome
        if not self._objectExists(oname):
            utils.error("Could not find chromosome '%s'" % chromosome)
            return None
        #fi
        return self._getObject(oname).whoHas(chromosome, *pargs, **kwargs)
    #edef
    
#eclass
