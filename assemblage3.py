import time
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Assemblage:

    def __init__(self, readFilePATH , outfile, sizeKmer):
        self.readFilePATH = readFilePATH
        self.outfile = outfile
        self.sizeKmer = sizeKmer
        self.hashTable = dict()

    def indexerReads(self):
        """Function to index the kmer of all reads and deletes the weakly represented kmer
        """
        with open(self.readFilePATH, 'r') as f:
            for _, seq, _ in FastqGeneralIterator(f):
                self.readsToKmer(seq)
        self.cleanIndex()

    def cleanIndex(self):
        """ Function that removes all non-representative kmer from the hash table
        """
        for key in list(self.hashTable.keys()):
            _,occ = self.hashTable[key]
            if occ < 3:
                del self.hashTable[key]

    # dic : key = Kmer ; value = (kmerUsed,nbOccurence,Valide)
    def readsToKmer(self,read):
        """Function that stores all the kmer of a read

        Args:
            read ([str]): Read which should be indexed
        """
        readSequence = read
        taille_Read = len(readSequence)
        taille_Kmer = self.sizeKmer
        positionStartKmer=0
        while (positionStartKmer + taille_Kmer) <= taille_Read:
            Kmer = readSequence[positionStartKmer:positionStartKmer + taille_Kmer]
            if Kmer in self.hashTable:
                _,occ = self.hashTable[Kmer]
                self.hashTable[Kmer] = (False,occ+1)
            else:
                self.hashTable[Kmer] = (False,1)
            positionStartKmer+=1

    def doContigs(self):
        """ Create one or more contigs from the kmer indexer
        """

        contigsList = []
        contig = ""
        kmerStart = self.foundKmerNotUse()
        while kmerStart is not None:
            extendRight = self.extendRight(kmerStart)
            extendLeft = self.extendLeft(kmerStart)
            contig = extendLeft[:-len(kmerStart)] + extendRight
            if len(contig) > 300:
                contigsList.append(contig)
            kmerStart = self.foundKmerNotUse()
        #print("La liste de contig",contigsList)
        self.writeContigs(contigsList)

    def writeContigs(self,contigsList):
        with open(self.outfile, 'w') as f:
            idx = 0
            for contig in contigsList:
                record = SeqRecord(
                    Seq(contig),
                    id=str(idx),
                    name="None",
                    description="",
                )
                SeqIO.write(record, f, "fasta")
                idx += 1

    def extendRight(self,kmerStart):
        """Extends a kmer to the right from the possible kmer in the hashtable

        Args:
            kmerStart ([str]): kmer to extend

        Returns:
            [str]: The best possible contig from the kmerStart extended to the right
        """

        contig = kmerStart
        kmerSuiv = self.foundNextKmerRight(kmerStart)
        while kmerSuiv is not None:
            contig += kmerSuiv[-1]
            kmerSuiv = self.foundNextKmerRight(kmerSuiv)
        return contig

    def extendLeft(self,kmerStart):
        """Extends a kmer to the left from the possible kmer in the hashtable

        Args:
            kmerStart ([str]): kmer to extend

        Returns:
            [str]: The best possible contig from the kmerStart extended to the left
        """

        contig = kmerStart
        kmerSuiv = self.foundNextKmerLeft(kmerStart)
        while kmerSuiv is not None:
            contig = kmerSuiv[0] + contig
            kmerSuiv = self.foundNextKmerLeft(kmerSuiv)
        return contig
            

    def foundKmerNotUse(self):
        """Search in the table for a kmer not used

        Returns:
            [str]: The first unused kmer found
        """

        for key,value in self.hashTable.items():
            use,occ = value
            if use == False:
                self.hashTable[key] = (True,occ)
                return key
        return None


    def foundNextKmerRight(self,kmer):
        """Search in the hash table, if a kmer exists by testing with the four nucleotides at the last nucleotide on the right

        Args:
            kmer ([str]): kmer for which we are looking for a next right

        Returns:
            [str]: The next right kmer
        """

        kmerSup = None
        KmerSuivList = []
        for nuc in ['A','C','G','T']:
            kmerTest = kmer[1:] + nuc
            if kmerTest in self.hashTable:
                used,occ = self.hashTable.get(kmerTest)
                if not used:
                    KmerSuivList.append((kmerTest,occ))
        if KmerSuivList:
            KmerSuivList.sort(key = lambda x: x[1],reverse = True)
            kmerSup = KmerSuivList[0][0]
            used,occ = self.hashTable.get(kmerSup)
            self.hashTable[kmerSup] = (True,occ)
        return kmerSup
    
    def foundNextKmerLeft(self,kmer):
        """Search in the hash table, if a kmer exists by testing with the four nucleotides at the first nucleotide on the left

        Args:
            kmer ([str]): kmer for which we are looking for a next left

        Returns:
            [str]: The next left kmer
        """

        kmerSup = None
        KmerSuivList = []
        for nuc in ['A','C','G','T']:
            kmerTest = nuc + kmer[:-1]
            if kmerTest in self.hashTable:
                used,occ = self.hashTable.get(kmerTest)
                if not used:
                    KmerSuivList.append((kmerTest,occ))
        if KmerSuivList:
            KmerSuivList.sort(key = lambda x: x[1],reverse = True)
            kmerSup = KmerSuivList[0][0]
            used,occ = self.hashTable.get(kmerSup)
            self.hashTable[kmerSup] = (True,occ)
        return kmerSup


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Filtre de bloom',description='Realise un assemblage des reads',usage="assemblage.py /fastqReads outfile.fasta sizeKmer")
    parser.add_argument('fastqReads',help='Read fastq',type=str)
    parser.add_argument('outfileGenom',help='name of contigs file',type=str)
    parser.add_argument('sizeKmer',help='size Of Kmer',type=str)

    args = parser.parse_args()

    fastqReadsPATH = args.fastqReads
    outfile = args.outfileGenom
    sizeKmer = int(args.sizeKmer)

    start_time = time.time()

    assemb = Assemblage(fastqReadsPATH,outfile,sizeKmer)
    assemb.indexerReads()
    assemb.doContigs()
    
    totalTime = time.time() - start_time 
    print('Total time in seconds:', totalTime) 
