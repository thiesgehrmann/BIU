import inspect, os
__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

__OUTDIR__ = config["outdir"]

###############################################################################

localrules: output

  # Blast genome_1 vs genome_2 and vice-versa
rule mappingBlastQuery:
  input:
    db    = config["genomes"]['A']["fasta"],
    query = config["genomes"]['B']["fasta"]
  output:
    res = "%s/result.{genome_1}.{genome_2}.tsv"% __OUTDIR__
  conda: "conda.yaml"
  params:
    blast_fields = config["blast_fields"],
    evalue       = config["e_threshold"],
    options   = config["options"],
    dbType    = 'prot' if config["genomes"]['A']["is_prot"] else 'dna',
    queryType = 'prot' if config["genomes"]['B']["is_prot"] else 'dna'
  threads: 5
  shell: """
    blat "{input.db}" "{input.query}" -t={params.dbType} -q={params.queryType} {params.options} -out=psl "{output.res}"
  """

  # Perform a reciprocal best blast hit to identify karyollele pairs
rule output:
  input:
    g1v2 = "%s/result.B.A.tsv" % __OUTDIR__
  output:
    res = "%s/%s" % (config["outdir"], config["output_file_name"])
  shell: """
    ln -s "{input.g1v2}" "{output.res}"
  """
