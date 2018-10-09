
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

  # Blast genome_db vs genome_query and vice-versa
rule mappingBlastQuery:
  input:
    db    = lambda wildcards: "%s/blastdb.%s.db" % (__OUTDIR__, wildcards.genome_db),
    trans = lambda wildcards: config["genomes"][wildcards.genome_query]["fasta"]
  output:
    res = "%s/result.{genome_db}.{genome_query}.tsv"% __OUTDIR__
  conda: "conda.yaml"
  params:
    blast_fields = config["blast_fields"],
    evalue       = config["e_threshold"],
    options      = config["options"]
  threads: 5
  shell: """
    diamond blastp -p "{threads}" -f "6" {params.blast_fields} {params.options} -e "{params.evalue}" -d "{input.db}" -q "{input.trans}" -o "{output.res}"
  """

  # Perform a reciprocal best blast hit to identify karyollele pairs
rule output:
  input:
    g1v2 = expand("%s/result.{genome_db}.{genome_query}.tsv" % __OUTDIR__, genome_db=["A"], genome_query=["B"]),
    g1fa = config["genomes"]["A"]["fasta"],
    g2fa = config["genomes"]["B"]["fasta"]
  output:
    res = "%s/%s" % (config["outdir"], config["output_file_name"])
  shell: """
    ln -s "{input.g1v2}" "{output.res}"
  """
