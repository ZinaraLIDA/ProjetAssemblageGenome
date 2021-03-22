import gzip

def readerFastq_writerTxt(fastaqFile, genome, parameter, outfile1, outfile2) :
    
    with gzip.open(fastaqFile, 'rt') as f, open(outfile1, 'w') as w1, \
    open(outfile2, 'w') as w2 :
        for line in f :
            if line[0] == '@' :
                line = f.readline()
                if line[0] == '@' :
                    line = f.readline()
                # Si le read appartient au génome, on le met dans le "genomeReads.txt",
                # sinon on le met dans "NoGenomeReads.txt"
                if Filter(line, genome, parameter, outfile1, outfile2) :
                    w1.write(line)
                if not Filter(line, genome, parameter, outfile1, outfile2) :
                    w2.write(line)

def Filter(read, genome, parameter, outfile1, outfile2) : 
    ...
    #return True or False

if __name__ == "__main__" :
    operation = readerFastq_writerTxt('reads-1M.fastq.gz', "sequenceSARS.fasta", parameter, \
    "genomeReads.txt", "NoGenomeReads.txt")
    