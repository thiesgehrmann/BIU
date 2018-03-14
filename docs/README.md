# BIU documentation

The BIU toolkit contains many many utilities that help my analysis. Here I provide documentation for them. I choose to do this via Jupyer Notebooks because it allows me to document not only the functionality but also the usage.

## Configuration

 * Changing global settings: biu.config.settings

## Accessing/Handling different formats
 * Standard Bioinformatics formats
  * Fasta Files : biu.formats.Fasta
  * GFF3 Files : biu.formats.GFF3
  * VCF Files : biu.formats.VCF
  * XLSX Files : biu.formats.XLSX
  * GAF Filies : biu.formats.GAF

 * Additional formats
  * SQLite Databases : biu.formats.SQLite
  * SQLite Dictionary : biu.formats.SQLDict

## Databases

 * Internal LUMC Datasets
  * BBMRI
  * LLS

 * External Databases
  * CADD : biu.db.CADD 
  * Clinvar : biu.db.ClinVar
  * Cosmic : biu.db.Cosmic
  * Genomes : biu.db.Genome
  * GnomAD : biu.db.Gnomad
  * Gene Ontology : biu.db.GO
  * GTeX : biu.db.GTeX
  * HAGR : biu.db.HAGR
  * KEGG : biu.db.KEGG
  * MiRmine : biu.db.MiRmine
  * Neo4jDB : biu.db.Neo4jDB
  * Reactome : biu.db.Reactome
  * Uniprot : biu.db.UniProt

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
