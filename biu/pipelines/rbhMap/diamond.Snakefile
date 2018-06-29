
import inspect, os
__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

__OUTDIR__ = config["outdir"]

###############################################################################

localrules: mappingBlastDB, mapping, output

  # Generate blast databases for each genome
rule mappingBlastDB:
  input:
    fa = lambda wildcards: config["genomes"][wildcards.genome]["fasta"]
  output:
    db = "%s/blastdb.{genome}.db"% __OUTDIR__
  conda: "conda.yaml"
  shell: """
    diamond makedb --in "{input.fa}" -d "{output.db}"
    touch "{output.db}"
  """

  # Blast genome_1 vs genome_2 and vice-versa
rule mappingBlastQuery:
  input:
    db    = lambda wildcards: expand("%s/blastdb.{genome}.db" % (__OUTDIR__), genome=[wildcards.genome_1]),
    trans = lambda wildcards: config["genomes"][wildcards.genome_2]["fasta"]
  output:
    res = "%s/result.{genome_1}.{genome_2}.tsv"% __OUTDIR__
  conda: "conda.yaml"
  params:
    blast_fields = config["blast_fields"],
    evalue       = config["e_threshold"],
    max_target_seqs = config["max_target_seqs"]
  threads: 5
  shell: """
    diamond blastp -p "{threads}" -f "6" {params.blast_fields} --max-target-seqs "{params.max_target_seqs}" -e "{params.evalue}" -d "{input.db}" -q "{input.trans}" -o "{output.res}"
  """

  # Perform a reciprocal best blast hit to identify karyollele pairs
rule mapping:
  input:
    g1v2 = expand("%s/result.{genome_2}.{genome_1}.tsv" % __OUTDIR__, genome_1=["A"], genome_2=["B"]),
    g2v1 = expand("%s/result.{genome_1}.{genome_2}.tsv" % __OUTDIR__, genome_1=["A"], genome_2=["B"]),
    g1fa = config["genomes"]["A"]["fasta"],
    g2fa = config["genomes"]["B"]["fasta"]
  output:
    mapping  = '%s/mappingFile.tsv' % __OUTDIR__
  params:
    blast_fields = config["blast_fields"]
  run:
    import utils as utils

      # Load the blast hits
    g12Hits = utils.indexListBy(utils.readBlastFile(input.g1v2[0], params.blast_fields), lambda x: x.qseqid)
    g21Hits = utils.indexListBy(utils.readBlastFile(input.g2v1[0], params.blast_fields), lambda x: x.qseqid)

      # For each gene, select the best hit in the other genome
    bestg12Hits = dict([ (qseqid, max(L, key=lambda x: x.bitscore)) for (qseqid,L) in g12Hits.items() ])
    bestg21Hits = dict([ (qseqid, max(L, key=lambda x: x.bitscore)) for (qseqid,L) in g21Hits.items() ])

      # Match the best hits to identify reciprocal best hits, and write the mappings to file
    matches = []
    with open(output.mapping, "w") as ofd:
      for key in bestg12Hits.keys():
        besthit = bestg12Hits[key].sseqid
        if (besthit in bestg21Hits) and (key == bestg21Hits[besthit].sseqid):
          ofd.write("%s\t%s\t%f\t%f\n" % (key, besthit, bestg21Hits[besthit].evalue, bestg21Hits[besthit].bitscore))
        #fi
      #efor
    #ewith

rule output:
  input:
    mapping = rules.mapping.output.mapping
  output:
    res = "%s/%s" % (config["outdir"], config["output_file_name"])
  shell: """
    ln -s "{input.mapping}" "{output.res}"
  """
