import gzip

import math
import mmh3
from bitarray import bitarray
    
def readerFastq_writerTxt(fastaqFile, genomeFile, p, length_kmer, outfile1, outfile2) :
    
    with gzip.open(fastaqFile, 'rt') as f2, open(genomeFile, 'r') as f1, \
    open(outfile1, 'w') as w1, open(outfile2, 'w') as w2 :
        genome = ''
        line1 = f1.readline()
        line1 = f1.readline()
        while line1 != "" :
            genome += line1.rstrip()
            line1 = f1.readline()
        list_kmers_genome = k_mers(genome, length_kmer)
        for line2 in f2 :
            if line2[0] == '@' :
                line2 = f2.readline()
                if line2[0] == '@' :
                    line2 = f2.readline()
                # Si le read appartient au génome, on le met dans le "genomeReads.txt",
                # sinon on le met dans "NoGenomeReads.txt"
                kmers_read = k_mers(line2, length_kmer)
                if Filter(kmers_read, p, list_kmers_genome) :
                    w1.write(line2 + '\n')
                if not Filter(kmers_read, p, list_kmers_genome) :
                    w2.write(line2 + '\n')

def k_mers(seq, k) :
    list_k_mers = []
    for i in range(len(seq) - k + 1) :
        list_k_mers.append(seq[i:i+k])
    return list_k_mers
         
class Filter :
    
    def __init__(self, kmers_read, p, list_kmers_genome) :
        self.__kmers_read = kmers_read
        self.__size = self.size_Filter(len(kmers_read), p)
        self.__hash_count = self.hash_count_Filter(self.__size, len(kmers_read))
        self.__bit_array = bitarray(self.__size)
        self.__bit_array.setall(0)
        self.__list_kmers_genome = list_kmers_genome
        for kmer_genome in self.__list_kmers_genome :
            self.add(kmer_genome)
    
    def add(self, string) :
        for seed in range(self.__hash_count) :
            result = mmh3.hash(string, seed)%self.__size
            self.__bit_array[result] = 1
        
    def lookup(self, string) :
        for seed in range(self.__hash_count) :
            result = mmh3.hash(string, seed)%self.__size
            if self.__bit_array[result] == 0 :
                return False
        return True
    
    def size_Filter(self, n, p) :
        m = -(n * math.log(p)) / (math.log(2)**2)
        return int(m)
    
    def hash_count_Filter(self, m, n) :
        k = (m/n) * math.log(2)
        return int(k)

    def __call__(self) :
        for kmer in self.__kmers_read :
            if  self.lookup(kmer)
            
if __name__ == "__main__" :
    operation = readerFastq_writerTxt('reads-1M.fastq.gz', \
    "sequenceSARS.fasta", 0.05, 25, "genomeReads.txt", "NoGenomeReads.txt")