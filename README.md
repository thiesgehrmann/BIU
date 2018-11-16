# BIU: BioInformatics Utilities.

## Functionality

BIU is a library of utilities and tools that provides:

  * Access to several standard (and not so standard) Bioinformatics formats (`biu.formats`)

  * Automated retrieval and interface to data from several databases (`biu.db`):
    * Genomes (Ensembl, Flybase, Wormbase)
    * Database/Datasets (ClinVar, DBSnp, GnomAD, Kegg, GO, UniProt, CADD, etc.)

  * Mapping between genome Identifiers (`biu.maps`):
    * Within a genome, where possible (gene -> Transcript -> Protein)
    * Between databases (with Ensembl BioMart)
    * Between species (DIOPT)

  * Bioinformatic pipelines (`biu.pipelines`):
    * BLAST/Diamond (`biu.pipelines.Blast`)
    * Variant Effect Predictor (biu.pipelines.VEP)
    * Write your own! (`biu.structures.Pipeline`)

  * Analysis of bioinformatic data (`biu.analysis`)
    * RNASeq data: (`biu.analysis.rnaseq`)
    * Covariates: (`biu.analysis.covariates`)

  * Medical tools (`biu.medical`):
    * health parameters: (`biu.medical.health`) (Framingham Risk Score, etc)

  * Data Structure processing (`biu.processing`):
    * Python lists
    * numpy Matrices
    * pandas DataFrames
   
  * Statistical analysis (`biu.stats`)
    * Regression (`biu.stats.regression`):
      * LOWESS regression
    * Enrichment (`biu.stats.enrichment`)
    * Multiple Testing Correction (FWER, FDR) (`biu.stats.correction`)
    * Statistical genetics (`biu.stats.genetics`)
      * Hardy Weinberg Equilibrium test)

  * Mathematical functionality (`biu.math`)
    * Interpolation (`biu.math.interpolation`)

## Examples

### Map a human gene to the worm, starting only with HGNC symbol
   ```python
   import biu
   
   # Load a genome from Ensembl
   hg = biu.genomes.Ensembl(grch37=True)
   wg = biu.genomes.Wormbase('Caenorhabditis_elegans')
   
   # Load maps
   biomart = biu.maps.BioMart(grch37=True)
   diopt = biu.maps.DIOPT()
   
   h_symbol = 'MTOR'
   h_id     = biomart.hgnc_symbol[h_symbol][0].entrezgene
   
   # Using DIOPT, Identify high quality mappings
   mtor_homologs = diopt.batch([mtor_geneid])
   w_symbol, w_id = mtor_homologs.c_elegans[0]
   
   print('%s(%s) -> %s(%s)' % (h_symbol, h_id, w_symbol, w_id))
   # MTOR(2475) -> let-363(WBGene00002583)
   ```

### Get all genes in the KEGG pathways that MTOR is in
   ```python
   import biu
   
   biomart = biu.maps.BioMart(grch37=True)
   kegg    = biu.db.KEGG()
   
   mtor_entrez    = biomart.hgnc_symbol["MTOR"][0].entrezgene
   mtor_pathways  = kegg.getGenePathways(mtor_entrez)
   mtor_neighbors = set([ g for p in mtor_pathways for g in kegg.getPathwayGeneIDs(p) ])
   
   print(len(mtor_neighbors))
   # 1934
   ```

### Retrieve information about a DBSNP variant
  ```python
  import biu
  where = '/exports/molepi/tgehrmann/biu/'
  biu.config.settings.setWhere(where)
  dbsnp = biu.db.DBSNP("human_9606_b150_GRCh37p13")
  
  chromosome, position = dbsnp["rs4894"]
  record = dbsnp.query(chromosome, position, position+1).records[0]
  
  print("Lookup rs4894: Chr(%s) Pos(%d) %s>%s" % (chromosome, position, record.REF, record.ALT[0] ))
  print("Reverse lookup: %s" % dbsnp(chromosome, position+1, alt))
  
  #Lookup rs4894: Chr(22) Pos(39917514) A>C
  #Reverse lookup: rs4894
  ```

### Lookup all the variants in gene region
   ```python
   import biu
   
   hg = biu.genomes.Ensembl(grch37=True)
   biomart = biu.maps.BioMart(grch37=True)
   dbsnp   = biu.db.DBSNP("human_9606_b150_GRCh37p13")
   
   symbol = 'MTOR'
   ensembl_gene_id = biomart.hgnc_symbol[symbol][0].ensembl_gene_id
   gff_entry = hg.gff.getID('gene:%s' % ensembl_gene_id)
   mtor_vars = dbsnp.query(gff_entry.seqid, gff_entry.start, gff_entry.end)
   print(len(mtor_vars))
   # 18577
   ```

## Installation    

Currently, there is no proper pythonic way to install the package.
Just clone the github repository and add BIU to your path.
  ```python
     os.path.append('path/to/BIU')
  ```

### Dependencies

BIU depends on a lot of modules.
However, much functionality in BIU does not depend on many other packages.
Additionally, BIU is implemented in such a way that it will load, even if none of those dependencies are actually installed.
This means that you can load BIU without any dependencies installed, and check if the functionality you are interested in works without them.
A recommended baseline for most functionality would be:
  * numpy
  * pandas
  * scipy
  * tabix (needed for some datasets, and VCF indexing)
  * vcf (needed for some datasets and VCF loading)
  * intervaltree (needed for indexing formats such as VCF and GFF files)

All external dependencies currently used are:
  * fastcluster
  * intervaltree
  * matplotlib
  * matplotlib_venn
  * numpy
  * openpyxl
  * pandas
  * scipy
  * seaborn
  * sklearn
  * snakemake
  * sqlite3
  * statsmodels
  * tabix
  * vcf
  * xlrd

## Documentation

See [Documentation](docs#biu-documentation) for documentation.

## Example usage

See [example.ipynb](docs/example.ipynb) for an example usage.
=======
