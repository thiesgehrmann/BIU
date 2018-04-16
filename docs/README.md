# BIU documentation

The BIU toolkit contains many many utilities that help my analysis. Here I provide documentation for them. I choose to do this via Jupyer Notebooks because it allows me to document not only the functionality but also the usage.

## Configuration

 * Changing global settings: biu.config.settings

## Accessing/Handling different formats
 * Standard Bioinformatics formats
  * Fasta Files : [biu.formats.Fasta](biu.formats.Fasta.ipynb)
  * GFF3 Files :  [biu.formats.GFF3](biu.formats.GFF3.ipynb)
  * VCF Files :   [biu.formats.VCF](biu.formats.VCF.ipynb)
  * XLSX Files :  [biu.formats.XLSX](biu.formats.XLSX.ipynb)
  * XLS Files : (Works like XLSX files)
  * GAF Filies :  [biu.formats.GAF](biu.formats.GAF.ipynb)

 * Additional formats
  * SQLite Databases :  [biu.formats.SQLite](biu.formats.SQLite.ipynb)
  * SQLite Dictionary : [biu.formats.SQLDict](biu.formats.SQLDict.ipynb)

## Databases

### Internal datasets

 Internal datasets contain data that cannot be accessed from the internet.
 They primarily can only be accessed from the SHARK cluster of the LUMC.

 * BBMRI : [biu.db.BBMRI](biu.db.BBMRI.ipynb)
 * LLS   : [biu.db.LLS](biu.db.LLS.ipynb)

### External Databases

 External datasets can be downloaded from the internet.
 Usually they can be retrieved by BIU, but for some, such as GTeX, you must provide the files yourself.

 * CADD : biu.db.CADD 
 * Clinvar : biu.db.ClinVar
 * Cosmic : biu.db.Cosmic
 * Genomes : [biu.db.Genome](biu.db.Genome.ipynb)
 * GnomAD : [biu.db.Gnomad](biu.db.Gnomad.ipynb)
 * Gene Ontology : [biu.db.GO](biu.db.GO.ipynb)
 * GTeX : biu.db.GTeX
 * HAGR : biu.db.HAGR
 * KEGG : biu.db.KEGG
 * MiRmine : biu.db.MiRmine
 * Neo4jDB : biu.db.Neo4jDB
 * Reactome : biu.db.Reactome
 * Uniprot : biu.db.UniProt

 * Human Gene Essentiality metrics
  * Residual Variation Intolerance Score (RVIS): biu.db.RVIS
  * Gene Damage Index (GDI): biu.db.GDI
  * De Novo Excess (DNE): biu.db.DNE

## Mapping utilities

 * Human mappings
  * Dictionary maps: biu.maps.Human
  * SQLite maps: biu.maps.HumanS

## Pipelines

 * Simple pipelines
  * VEP: biu.pipelines.VEP
  * LiftOver: biu.pipelines.LiftOver

## Internals

 * Utils : biu.utils
 * FileManager : biu.structures.FileManager
 * ResourceManager : biu.structures.ResourceManager
 * LazyObject : biu.structures.LazyObject
 * Pipeline: biu.pipelines.Pipeline
