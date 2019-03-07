from ..structures import Dataset2
from .. import formats
from .. import utils
from .. import ops
from .. import medical

pd = utils.py.loadExternalModule("pandas")
np = utils.py.loadExternalModule("numpy")

class GOTO2(Dataset2):
    """
    An interface to the GOTO dataset.
    
    Typical usage:
    --------------
    goto = biu.db.GOTO2()
    
    # RNA-Seq expression data
    goto.expr   # Expression counts
    goto.blood  # Covariates for blood samples
    goto.muscle # Covariates for muscle samples
    
    # Genetic variation
    goto.filter(15, 50001000, 50005000)
    goto.filter_regions([(15, 50001000, 50005000),(13, 50001000, 50005000)])
    
    # Metabolomics data
    goto.metabolomics # Metabolyte measurements
    goto.metabolomics_cov # Covariates for metabolyte data
    """

    versions = { "current":
        { "chrs" : [ "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"] }
    }
    
    version = None
    
    __samplesPerIndividual = None
    
    def __init__(self, version=list(versions.keys())[0], *pargs, **kwargs):

        super(GOTO2, self).__init__("GOTO/%s" % version , *pargs, **kwargs)
        self.version = version
        
        for chrom in self.versions[version]["chrs"]:
            self._obj.add_file("vcf_%s.vcf.tgz" % chrom, utils.Acquire2("/exports/molepi/GOTO_HRC1.1/chr%s.dose.vcf.gz" % (chrom)), finalize=False)
            self._obj.add_file("vcf_%s.vcf.tgz.tbi" % chrom, utils.Acquire2("/exports/molepi/GOTO_HRC1.1/chr%s.dose.vcf.gz.tbi" % (chrom)), finalize=False)
            self._obj.register('vcf_%s' % chrom, ["vcf_%s.vcf.tgz" % chrom, "vcf_%s.vcf.tgz.tbi" % chrom],
                              lambda f, k=chrom: formats.VCF2(f["vcf_%s.vcf.tgz" % k], tabix=True))
        #efor
        
        self._obj.add_file("expr.csv", utils.Acquire2('/exports/molepi/tgehrmann/data/GOTO_data/goto_expr.csv'), finalize=False)
        self._obj.register("expr", ['expr.csv'], lambda f: pd.read_csv(f['expr.csv'], index_col=0) )
        
        def load_blood_cov(file_name):
            D = pd.read_csv(file_name)
            D['Sample'] = D.Sample.apply(lambda s: s.replace('-', '_'))
            D['flowcell_lane'] = D.Sample.apply(lambda s: s.split('_')[2])

            discrete_covariates = [ 'IOP2_ID', 'timepoint', 'Sex', 'Status', 'Lipid_lowering_medication',
                                   'Current_Smoker', 'Antihypertensive_medication', 'Date_of_sample_collection',
                                   'intervention', 'nutridrink', 'sampID', 'Labnr', 'Sample', 'studID',
                                   'flowcell', 'Lane', 'Index', 'blood.or.muscle', 'plate',
                                   'RNA_blood_Isolatieseries', 'RNA_Isolation_series', 'flowcell_lane' ]
            continuous_covariates = [ 'Age_Baseline', 'Glucose', 'Cortisol', 'IGF_BP3', 'LN_Insulin', 
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
                                     'LabonchipRIN', 'Labonchip28S18S', 'totalyieldug' ]
            
            for c in discrete_covariates:
                D[c] = ops.series.cast_str(D[c]).astype('category')
            #efor
            
            for c in continuous_covariates:
                D[c] = ops.series.cast_float(D[c])
            #efor
            
            FRS = [ medical.health.indicators.framingham_risk_score(row.Sex == 1,
                                                                    row.Age_Baseline,
                                                                    row.Cholesterol, 
                                                                    row.Cholesterol_HDL,
                                                                    row.Systolic_Blood_Pressure,
                                                                    row.Antihypertensive_medication == 1,
                                                                    row.Current_Smoker == 1,
                                                                    False)
                   for i, row in D.iterrows() ]
            D['FRS'] = np.array(FRS)

            return D.set_index('Sample')
        #edef
        
        self._obj.add_file("blood_cov.csv", utils.Acquire2('/exports/molepi/tgehrmann/data/GOTO_data/datasheet_RNAseq_blood_V2.csv'), finalize=False)
        self._obj.register("blood", ["blood_cov.csv"], lambda f: load_blood_cov(f['blood_cov.csv']))
        
        def load_muscle_cov(file_name):
            D = pd.read_csv(file_name)
            D['Sample'] = D.Sample.apply(lambda s: s.replace('-', '_'))
            
            categorical = [ 'sampID', 'IOP2_ID', 'Sex', 'Lipid_lowering_medication', 'Antihypertensive_medication', 
                            'Status', 'studID', 'flowcell', 'Lane', 'Index', 'Visitnr', 'Biopsynumberused',
                            'Isolationsseries', 'plate', 'blood.or.muscle' ]
            numeric     = [ 'Age_Baseline','Nanodrop260280', 'LabonchipRIN', 'Labonchip28S18S', 'totalyieldug',
                            'passedQC_perc', 'Total_reads', 'Mapped_reads', 'mappedreads_perc', 'PFALIGNEDBASES',
                            'MEDIAN5PRIMEBIAS', 'MEDIAN3PRIMEBIAS', 'MEDIAN5PRIMETO3PRIMEBIAS', 'meaninsertsize',
                            'standarddeviation', 'medianinsertsize', 'Glucose', 'Creatinine', 'Cholesterol', 'LN_Insulin',
                            'WholeBody_percentage_fat', 'WholeBody_percentage_lean', 'total_legs_fat_percentage',
                            'total_legs_lean_percentage', 'BMI', 'ln_Collagen', 'Pax7_per_nucleus', 'ln_Pax7_per_fiber',
                            'Ki67_per_nucleus', 'ln_Ki67_per_fiber', 'MyHC1', 'MyHC1_2A', 'MyHC2A', 'MyHC1_2X', 'MyHC2A_2X',
                            'MyHC2X', 'ln_CSA_um2', 'ln_lipiddropletarea', 'HDL_cholesterol', 'LN_Triglycerides',
                            'LDL_cholesterol', 'LN_Leptin_Adiponectin', 'XXLVLDLL_LNscaled', 'XLVLDLL_LNscaled',
                            'LVLDLL_LNscaled', 'MVLDLL_LNscaled', 'SVLDLL_LNscaled', 'XSVLDLL_LNscaled', 'IDLL_LNscaled',
                            'LLDLL_LNscaled', 'SLDLL_LNscaled', 'XLHDLL_LNscaled', 'MHDLL_LNscaled', 'SHDLL_LNscaled',
                            'VLDLD_LNscaled', 'LDLD_LNscaled', 'HDLD_LNscaled', 'SerumC_LNscaled', 'VLDLC_LNscaled',
                            'RemNAtC_LNscaled', 'LDLC_LNscaled', 'HDLC_LNscaled', 'EstC_LNscaled', 'FreeC_LNscaled',
                            'TotPG_LNscaled', 'PC_LNscaled', 'SM_LNscaled', 'TotCho_LNscaled', 'ApoA1_LNscaled',
                            'ApoB_LNscaled', 'TotFA_LNscaled', 'FALen_LNscaled', 'UnsatDeg_LNscaled', 'LA_LNscaled',
                            'CLA_LNscaled', 'FAw3_LNscaled', 'FAw6_LNscaled', 'PUFA_LNscaled', 'MUFA_LNscaled',
                            'SFA_LNscaled', 'Glc_LNscaled', 'Lac_LNscaled', 'Pyr_LNscaled', 'Cit_LNscaled', 'Glol_LNscaled',
                            'Ala_LNscaled', 'Gln_LNscaled', 'Gly_LNscaled', 'His_LNscaled', 'Ile_LNscaled', 'Leu_LNscaled',
                            'Val_LNscaled', 'Phe_LNscaled', 'Tyr_LNscaled', 'Ace_LNscaled', 'AcAce_LNscaled',
                            'bOHBut_LNscaled', 'Crea_LNscaled', 'Alb_LNscaled', 'Gp_LNscaled', 'Gripstrength' ]
            
            for c in categorical:
                D[c] = ops.series.cast_str(D[c]).astype('category')
            #efor
            
            for c in numeric:
                D[c] = ops.series.cast_float(D[c])
            #efor

            return D.set_index('Sample')
        #edef

        self._obj.add_file("muscle_cov.csv", utils.Acquire2('/exports/molepi/tgehrmann/data/GOTO_data/muscle_QC_covariates_filesv2.csv'), finalize=False)
        self._obj.register("muscle", ["muscle_cov.csv"], lambda f: load_muscle_cov(f['muscle_cov.csv']))

        self._obj.add_file('metab_biomat.csv', utils.Acquire2('~/repos/GOTO_analysis/data/biomaterial.csv'), finalize=False)
        self._obj.add_file('metab_visit.csv', utils.Acquire2('~/repos/GOTO_analysis/data/visit.csv'), finalize=False)
        self._obj.add_file('metab.csv', utils.Acquire2('~/repos/GOTO_analysis/data/nightingale_metabolomics.csv'), finalize=False)

        self._obj.register("metabolomics", ["metab.csv"], lambda f: pd.read_csv(f['metab.csv']).drop(columns=['labnr']).set_index('biomaterial_id'))
        
        def load_metab_cov(f):
            BM = pd.read_csv(f["metab_biomat.csv"])
            VS = pd.read_csv(f["metab_visit.csv"])
            
            D = BM.set_index(['person_id','intervention']).join(VS.set_index(['person_id','intervention'])).reset_index()
            D['time_point'] = 1 + 2 * D.intervention + (1 - D.fasted)
            
            categorical = [ 'person_id', 'intervention', 'fasted', 'labnr', 'type', 'biopsy_nr',
                            'time_fatbiopsy', 'time_musclebiopsy', 'time_nutridrink', 'LLnr', 'biomaterial_id',
                            'visit_date', 'lipid_lowering_medication', 'antihypertensive_medication', 'time_point' ]

            numeric = [ 'length', 'weight', 'waist_circumference', 'hip_circumference', 'armspan', 'length_leg',
                        'AF_score', 'systolic_bp', 'diastolic_bp', 'SPPB', 'walking', 'chair_stand', 'remaining' ]

            for col in categorical:
                D[col] = ops.series.cast_str(D[col]).astype('category')
            #efor

            for col in numeric:
                D[col] = ops.series.cast_float(D[col])
            #efor
            
            return D.set_index('biomaterial_id')
        #edef
        self._obj.register("metabolomics_cov", ["metab_biomat.csv", "metab_visit.csv"], lambda f: load_metab_cov(f))
        
        # Add all the Database tables, for now.
        tgz = '/exports/molepi/tgehrmann/data/GOTO_data/GOTO-database-711cf1c7e3ff39b8f44a9015759c702752af3525.tar.gz'
        tgz = utils.Acquire2(tgz)
        tgz = tgz.gunzip().untar('GOTO-database-711cf1c7e3ff39b8f44a9015759c702752af3525/tables')
        tables = [ 'accelerometry', 'aliquots', 'BIA_mobile', 'BIA_standing', 'biomaterial', 'cell_counts', 
                   'ckcl', 'dexa_hip', 'dexa_spine', 'dexa_whole_body', 'glycans', 'grip_strength', 'histology',
                   'MRI_spectroscopy', 'nightingale_metabolomics', 'person', 'respiratory', 'rna_seq', 'visit' ]

        [self._obj.add_file("%s.csv" % table, tgz.select('%s.csv' % table)) for table in tables]
        [self._obj.register("db_%s" % t, ["%s.csv" % t], lambda f, t=t: pd.read_csv(f["%s.csv" % t])) for t in tables]
        
        self._add_str_func(lambda s: "Version: %s" % self.version)
    #edef

    @property
    def individuals(self):
        """
        Return a list of all IOP2_IDs in the rna-seq and metabolomics datasets
        """
        return list(map(str, ops.lst.uniq(list(self.blood.IOP2_ID.values) + list(self.muscle.IOP2_ID.values) + list(self.metabolomics_cov.person_id.values))))
    #edef

    @property
    def samples_per_individual(self):
        if self.__samplesPerIndividual is None:
            muscle = ops.lst.group(self.muscle.reset_index()[["IOP2_ID", "Sample", "Visitnr"]].values)
            blood  = ops.lst.group(self.blood.reset_index()[["IOP2_ID", "Sample", "timepoint"]].values)
            self.__samplesPerIndividual = { str(p) : { 'blood' : { str(s[2]) : s[1] for s in blood.get(str(p),[]) }, 
                                                       'muscle' : { str(s[2]) : s[1] for s in muscle.get(str(p),[]) } }
                                           for p in self.individuals }
        #fi
        return self.__samplesPerIndividual
    #edef

    def filter(self, chrom, start, end, *pargs, **kwargs):
        """
        Perform a filter for several regions
        
        parameters:
        -----------
        chrom: String|int. Chromosome of interest
        start: int. Start of region of interest
        end: int. End of region of interest
        *pargs, **kwargs: See additional arguments for VCF2.filter
        
        Returns: VCF2 object
        """
        oname = "vcf_%s" % str(chrom)

        if oname not in self._obj:
            raise AttributeError("Could not find chromosome '%s'" % chrom)
        #fi
        
        return self._obj[oname].filter(chrom, start, end, *pargs, **kwargs)
    #edef
    
    def filter_regions(self, regions, chrom=None, start=None, end=None, *pargs, **kwargs):
        """
        Perform a filter for several regions
        
        parameters:
        -----------
        regions: A list of 3-tuples (chrom, start, end) for each region of interest
        chrom, start, end: Ignored
        *pargs, **kwargs: See additional arguments for VCF2.filter
        
        Returns: VCF2 object
        """
        rets = [ self.filter(c, s, e, *pargs, **kwargs) for (c,s,e) in regions ]
        return rets[0].merge(rets)
    #edef


    def get_var(self, chrom, *pargs, **kwargs):
        """
        Get the variant record for a specific variant
        
        parameters:
        -----------
        chrom: str|int. Which chromosome the variant is on
        pos:   int. What position the variant is on.
        ref:   str. What is the reference allele?
        alt:   str. What is the alternative allele?
        
        returns:
        --------
        A cyvcf2.VCF.Variant object if the variant exists. Otherwise None
        """
        oname = "vcf_%s" % str(chrom)

        if oname not in self._obj:
            raise AttributeError("Could not find chromosome '%s'" % chrom)
        #fi

        return self._obj[oname].get_var(chrom, *pargs, **kwargs)
    #edef

    def who_has(self, chrom, *pargs, **kwargs):
        """
        Determine who has a specific variant

        parameters:
        -----------
        chrom: str|int. Which chromosome the variant is on
        pos:   int. What position the variant is on.
        ref:   str. What is the reference allele?
        alt:   str. What is the alternative allele?

        returns:
        --------
        List of sample IDs for who has the variant.
        """
        oname = "vcf_%s" % str(chrom)

        if oname not in self._obj:
            raise AttributeError("Could not find chromosome '%s'" % chrom)
        #fi

        return self._obj[oname].who_has(chrom, *pargs, **kwargs)
    #edef

    
#eclass