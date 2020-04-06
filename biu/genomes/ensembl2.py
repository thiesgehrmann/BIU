from ..structures import Dataset2
from ..utils import Acquire2
from .. import formats

import ftplib

class Ensembl2(Dataset2):
    def __init__(self, release=92, organism='homo_sapiens', grch37=False, *pargs, **kwargs):
        version = '%s.%s%s' % (organism, str(release), '.grch37' if grch37 else '')
        super(Ensembl2, self).__init__("ensembl/%s" % version, *pargs, **kwargs)
        
        self._obj.add_file('ensembl_index.sqlite', Acquire2().touch('ensembl_index'))
        self._obj.register('ensembl_index', ['ensembl_index.sqlite'],
                           lambda f: formats.SQLDict(f['ensembl_index.sqlite']))
        
        files = self._get_from_ensembl_index(grch37, str(release), organism)
        
        if files['gff'] is not None:
            self._obj.add_file('gff3.gff', Acquire2().curl(files['gff']).gunzip())
            self._obj.register('gff', ['gff3.gff'], lambda f: formats.GFF3(f['gff3.gff']))
        #fi
        
        if files['dna'] is not None:
            indivs = [ Acquire2().curl(f) for f in files['dna'] ]
            self._obj.add_file('genome.fa', Acquire2().merge(indivs, method='zcat'))
            self._obj.register('genome', ['genome.fa'], lambda f: formats.Fasta(f['genome.fa'], 
                                                                                seqType=formats.Sequence.DNATYPE))
        #fi
        
        if files['cds'] is not None:
            self._obj.add_file('cds.fa', Acquire2().curl(files['cds']).gunzip())
            self._obj.register('cds', ['cds.fa'], lambda f: formats.Fasta(f['cds.fa'],
                                                                          seqType=formats.Sequence.DNATYPE))
        #fi
        
        if files['aa'] is not None:
            def id_map_func(inFile, outFile):
                fasta = formats.Fasta(inFile[0])
                with open(outFile, 'w') as ofd:
                    ofd.write('\t'.join(['gene', 'transcript', 'protein', 'symbol']) + '\n')
                    for seq in fasta:
                        fullName = fasta[seq].fullName
                        data = dict([(p.split(':')[0], ':'.join(p.split(':')[1:]))
                                     for p in fullName.split(' ') ])
                        ofd.write('\t'.join([data.get('gene', ''),
                                             data.get('transcript', ''),
                                             seq,
                                             data.get('gene_symbol', '') ]))
                        ofd.write('\n')
                    #efor
                #ewith
                return 0
              #edef
                
            aa_file  = Acquire2().curl(files['aa']).gunzip()
            ids_file = aa_file.func(id_map_func)
            
            self._obj.add_file('aa.fa', aa_file)
            self._obj.add_file('ids.tsv', ids_file)
            self._obj.register('aa', ['aa.fa'], lambda f: formats.Fasta(f['aa.fa'], seqType=formats.Sequence.PROTTYPE))
            self._obj.register('ids', ['ids.tsv'], lambda f: formats.MappingIndex(f['ids.tsv']))
        #fi
    #edef
        
    def _get_from_ensembl_index(self, grch37, release, organism):
        k = (grch37, release, organism)
        if True:#k not in self.ensembl_index:
            basedir = '/pub/grch37' if grch37 else '/pub'
            conn = ftplib.FTP("ftp.ensembl.org")
            conn.login()

            def get_gff3(conn, basedir, release, organism):
                # This is ugly, but the most stable way I was able to find the correct GFF3 file.
                # in GRCH37, it didnt match the obvious pattern...
                uri = [ u[0] for u in
                         sorted([ (line, len(line.split('.'))) 
                                 for line in
                                     conn.nlst("%s/release-%s/gff3/%s" % (basedir, release, organism)) if (len(line.split('.')) > 3) ], 
                                key=lambda x:x[1]) ]
                if len(uri) > 0:
                    return "ftp://ftp.ensembl.org/%s" % uri[0]
                #fi
                return None
            #edef

            def get_genome(conn, basedir, release, organism):
                uri = [ line for line in
                       conn.nlst("%s/release-%s/fasta/%s/dna" % (basedir, release, organism)) if 'dna.chromosome' in line ]
                return [ "ftp://ftp.ensembl.org/%s" % f for f in uri ]
            #edef

            def get_cds(conn, basedir, release, organism):
                uri = [ line for line in conn.nlst("%s/release-%s/fasta/%s/cds" % (basedir, release, organism)) if 'fa.gz' in line ]
                if len(uri) > 0:
                    return "ftp://ftp.ensembl.org/%s" % uri[0]
                #fi
                return None
            #edef

            def get_aa(conn, basedir, release, organism):
                uri = [ line for line in conn.nlst("%s/release-%s/fasta/%s/pep" % (basedir, release, organism)) if 'all.fa.gz' in line ]
                if len(uri) > 0:
                    return "ftp://ftp.ensembl.org/%s" % uri[0]
                #fi
                return None
            #edef
        
        
            self.ensembl_index[k] = {
                'gff' : get_gff3(conn, basedir, release, organism),
                'dna' : get_genome(conn, basedir, release, organism),
                'cds' : get_cds(conn, basedir, release, organism),
                'aa'  : get_aa(conn, basedir, release, organism)
            }
        #fi
        return self.ensembl_index[k]
    #edef
    
    def seq(self, ID):
        return self.gff.seq(ID, self.genome)
    #edef
#eclass