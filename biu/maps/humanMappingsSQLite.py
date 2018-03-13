from ..structures import fileManager as fm
from ..structures import resourceManager as rm

from .. import utils

import os

###############################################################################

def urlFileIndex():
  files = {}
  files["gene2ensembl"] = ("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz", "gene2ensembl.tsv", {})
  files["gene2refseq"]  = ("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz", "gene2refseq.tsv", {})
  files["geneInfo"]     = ("ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/Homo_sapiens.gene_info.gz", "geneinfo.tsv", {})
  files["uniprotmap"]   = ("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz", "uniprotmap.tsv", {})
  files["sqlite_db"]    = (None, 'db.sqlite', {})
  return { k : (u, 'humanMappings_%s' % (l), o) for (k, (u, l, o)) in files.items() }
#edef

###############################################################################

class HumanMapping(fm.FileManager):

  gene2ensembl = None
  #gene2refseq  = None
  geneInfo     = None
  uniprotmap   = None
  sqlite       = None

  def __init__(self, **kwargs):
    fm.FileManager.__init__(self, urlFileIndex(), objects=["sqlite"], **kwargs)

    dbFile = self.getFileName("sqlite_db")

    self.touchFile("sqlite_db")
    self.sqlite = rm.SQLiteResourceManager(self, "sqlite_db")

    tableDetails = {
      "geneIDInfo": {
        "what" : "geneInfo",
        "names" : [ "tax_id", "geneid", "symbol", "locustag", "synonyms", "dbxrefs", "chromosome", "map_location", "description", "type_of_gene", "symbol_from_nomenclature_authority", "full_name_from_nomenclature_authority", "nomenclature_status", "other_designations", "modification_date", "feature_type" ],
        "index" : [ "geneid", "symbol" ],
        "delimiter" : '\t' },

      "gene2Ensembl" : {
        "what" : "gene2ensembl",
        "names" : [ "tax_id", "geneid", "ensembl_gene_identifier", "rna_nucleotide_accession.version", "ensembl_rna_identifier", "protein_accession.version", "ensembl_protein_identifier" ],
        "index" : [ "geneid", "ensembl_gene_identifier" ],
        "delimiter" : '\t' },

      "gene2refseq" : {
        "what" : "gene2refseq",
        "names" : [ "tax_id", "geneid", "status", "rna_nucleotide_accession.version", "rna_nucleotide_gi", "protein_accession.version", "protein_gi", "genomic_nucleotide_accession.version", "genomic_nucleotide_gi", "start_position_on_the_genomic_accession", "end_position_on_the_genomic_accession", "orientation", "assembly", "mature_peptide_accession.version", "mature_peptide_gi", "symbol" ],
        "index" : ["geneid" ],
        "delimiter" : '\t' }
    }

    existingTables = self.sqlite.getTableNames()
    for table in tableDetails:
      if table not in existingTables:
        utils.dbm("Adding tables into mapping database. This might take a while...")
        self.sqlite.dropTable("mapping")
        what = tableDetails[table]["what"]
        if not(self.satisfyRequiredFiles([what])):
          utils.error("Unable to download '%s'. Cannot make Table." % what)
        else:
          utils.dbm("Adding '%s' into SQLite db as table '%s'" % (what, table))
          self.sqlite.createTableFromFile(self.getFileName(what), table, names=tableDetails[table]["names"], delimiter=tableDetails[table]["delimiter"])
        #fi
      #fi
    #efor

    if "mapping" not in existingTables:
      self.createMappingTable()
    #fi 
  #edef

  ############################################################################### 

  def createMappingTable(self):
    self.sqlite.dropTable("mapping")
    self.sqlite.execute("""
       CREATE TABLE mapping AS
       SELECT DISTINCT
         geneIDInfo.geneid,
         geneIDInfo.symbol,
         gene2Ensembl.ensembl_gene_identifier AS ensemblid
       FROM geneIDInfo 
       LEFT JOIN gene2Ensembl ON gene2Ensembl.geneid = geneIDInfo.geneid
       UNION ALL
       SELECT DISTINCT
         geneIDInfo.geneid,
         geneIDInfo.symbol,
         gene2Ensembl.ensembl_gene_identifier AS ensemblid
       FROM gene2Ensembl
       LEFT JOIN geneIDInfo ON gene2Ensembl.geneid = geneIDInfo.geneid
       where geneIDInfo.symbol IS NULL;""")
   #edef

  def getEnsemblGeneID(self, ensemblID):
    query = """
      SELECT DISTINCT mapping.geneid
      FROM mapping
      WHERE mapping.ensemblid = ?"""
    return [ r[0] for r in self.sqlite.execute(query, [ ensemblID ]) ]

  def getEnsemblSymbol(self, ensemblID):
    query = """
      SELECT DISTINCT mapping.symbol
      FROM mapping
      WHERE mapping.ensemblid = ?"""
    return [ r[0] for r in self.sqlite.execute(query, [ ensemblID ]) ]
  #edef

  def getSymbolEnsembl(self, geneSymbol):
    query = """
      SELECT DISTINCT mapping.ensemblid
      FROM mapping
      WHERE mapping.symbol = ?"""
    return [ r[0] for r in self.sqlite.execute(query, [ geneSymbol ]) ]
  #edef

  def getSymbolGeneID(self, geneSymbol):
    query = """
      SELECT DISTINCT mapping.geneid
      FROM mapping
      WHERE mapping.symbol = ?"""
    return [ r[0] for r in self.sqlite.execute(query, [ geneSymbol ]) ]

  def getGeneIDEnsembl(self, geneID):
    query = """
      SELECT DISTINCT mapping.ensemblid
      FROM mapping
      WHERE mapping.geneid = ?"""
    return [ r[0] for r in self.sqlite.execute(query, [ geneID ]) ]
  #edef

  def getGeneIDSymbol(self, geneID):
    query = """
      SELECT DISTINCT mapping.symbol
      FROM mapping
      WHERE mapping.geneid = ?"""
    return [ r[0] for r in self.sqlite.execute(query, [ geneID ]) ]
  #edef

#eclass
