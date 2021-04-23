import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse
import time
from bitarray import bitarray
import mmh3
import math

#======================================================================================================================

class Filter:

    # parameter (tailleKmer,probabilite faux positif)
    def __init__(self, readFilePATH , genomeFilePATH, parameter, outfile1, outfile2):
        self.readFilePATH = readFilePATH
        self.genomeFilePATH = genomeFilePATH
        self.parameter = parameter
        self.outfile1 = outfile1
        self.outfile2 = outfile2
        

    # Remplir les filtres de bloom avec le genome
    def FilterGenome(self) :
        """Builds a bloom filter for the genome and its reverse complement.
        The function calculates parameters such as the number of hash functions and the filter size.
        """
        genome = ""
        genomeReverseComplement = ""
        with open(self.genomeFilePATH, 'rt') as f:
            record = SeqIO.read(f, "fasta")
            genome = record.seq
            genomeReverseComplement = record.reverse_complement().seq
            nbKmer = len(genome) - self.parameter[0] + 1
            filterSize = int(-(nbKmer * math.log(self.parameter[1]))/((math.log(2))**2))
            self.filterTab = bitarray(filterSize)
            self.filterTabReverseComplement = bitarray(filterSize)
            self.filterSize = filterSize
            self.nbHash = round((filterSize/nbKmer) * math.log(2))
        
        taille_Kmer = self.parameter[0]
        positionStartKmer = 0

        while (positionStartKmer + taille_Kmer) < len(genome):
            Kmer = str(genome[positionStartKmer:positionStartKmer + taille_Kmer])
            self.addhashKmer(Kmer,self.nbHash,True)
            positionStartKmer+=1

        positionStartKmer = 0
        while (positionStartKmer + taille_Kmer) < len(genomeReverseComplement):
            Kmer = str(genomeReverseComplement[positionStartKmer:positionStartKmer + taille_Kmer])
            self.addhashKmer(Kmer,self.nbHash,False)
            positionStartKmer+=1

    def filterReads(self,read):
        """Check the presence of each kmer in the two bloom filters

        Args:
            read (Seq()): Read whose presence in the filter is determined

        Returns:
            (Boolean,Boolean): A tuple whose first parameter determines whether the read belongs to the genome, and the second parameter indicates whether it is on the positive or negative strand
        """
        readSequence = read
        taille_Read = len(readSequence)
        taille_Kmer = self.parameter[0]
        positionStartKmer=0
        nbTotalKmer=0
        nbKmerPositif=0
        nbKmerPositifReverseComp=0
        while (positionStartKmer + taille_Kmer) <= taille_Read:
            nbTotalKmer+=1
            Kmer = readSequence[positionStartKmer:positionStartKmer + taille_Kmer]
            if self.foundKmer(Kmer,True):
                nbKmerPositif+=1
            if self.foundKmer(Kmer,False):
                nbKmerPositifReverseComp+=1
            positionStartKmer+=20
        
        if nbKmerPositif/nbTotalKmer > 0.25:
            return (True,True)
        elif nbKmerPositifReverseComp/nbTotalKmer > 0.25:
            return (True,False)
        else:
            return (False,False)

    def trieReads(self):
        """Function that sorts the reads with the filter and writes to the corresponding file 
        """
        with gzip.open(self.readFilePATH, 'rt') as f, open(self.outfile1, 'w') as w1, open(self.outfile2, 'w') as w2 :
            for title, seq, qual in FastqGeneralIterator(f):
                result = self.filterReads(seq)
                if result[0]:
                    if result[1]:
                        w1.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                    else:
                        rev = str(Seq(seq).reverse_complement())
                        w1.write("@%s\n%s\n+\n%s\n" % (title, rev, qual))
                else:
                    w2.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))


    def foundKmer(self,Kmer,brin):
        """Function that looks for a kmer in the bloom filters according to the different iterations of the hash function

        Args:
            Kmer (str): kmer search in the filter
            brin (boolean): True if searching in the positive strand, false otherwise

        Returns:
            boolean: True if present in the filter, false otherwise
        """
        for i in range(self.nbHash):
            if brin:
                if self.filterTab[(mmh3.hash(Kmer,i) % self.filterSize)] == 0:
                    return False
            else:
                if self.filterTabReverseComplement[(mmh3.hash(Kmer,i) % self.filterSize)] == 0:
                    return False
        return True

    def addhashKmer(self,Kmer,n,brin):
        """Add a kmer in one of the filters

        Args:
            Kmer (str): kmer to be added in the filter
            n (int): Number of hash functions
            brin (boolean): True for the positive strand, false otherwise
        """
        for i in range(n):
            if brin:
                hashValue = (mmh3.hash(Kmer,i) % self.filterSize)
                self.filterTab[hashValue] = 1
            else:
                hashValue = (mmh3.hash(Kmer,i) % self.filterSize)
                self.filterTabReverseComplement[hashValue] = 1

#======================================================================================================================
# python mainFaster2.py ..\fichier-genome\reads-1M.fastq.gz ..\fichier-genome\genome_SARS_cov2.fasta 25 0.01 sortieReadGenome.fastq sortieReadNotGenom.fastq


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Filtre de bloom',description='Realise un filtre de bloom',usage="rech7.py /fastqReads /fastaGenom sizeKmer probabilite outfile1 outfile2")
    parser.add_argument('fastqReads',help='Read fastq',type=str)
    parser.add_argument('fastaGenom',help='Genom Fasta',type=str)
    parser.add_argument('sizeKmer',help='Size of the Kmer',type=str)
    parser.add_argument('probabilite',help='Probabilite',type=str)
    parser.add_argument('outfile1',help='File read which are in the genome',type=str)
    parser.add_argument('outfile2',help='File read which are not in the genome',type=str)

    args = parser.parse_args()

    fastqReadsPATH = args.fastqReads
    fastaGenomPATH = args.fastaGenom
    sizeKmer = args.sizeKmer
    probabilite = args.probabilite
    outfile1 = args.outfile1
    outfile2 = args.outfile2

    start_time = time.time()
    parameter = (int(sizeKmer),float(probabilite))
    filter = Filter(fastqReadsPATH,fastaGenomPATH,parameter,outfile1,outfile2)
    start_time_genom = time.time()
    filter.FilterGenome()
    finish_time_genom = time.time()
    filter.trieReads()

    totalTime = time.time() - start_time
    totalTimeForGenom = finish_time_genom - start_time_genom
    print('Total time in seconds for genom:', totalTimeForGenom)   
    print('Total time in seconds:', totalTime) 


    