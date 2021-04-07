import gzip
import math
import mmh3
from bitarray import bitarray
    
def readerFastq_writerTxt(fastaqFile, genomeFile, parameters, outfile1, outfile2) :
    """[Ouvre un fichier fastaq compressé contenant des reads, écrit dans un fichier txt les reads appartenant au génome et dans un autre ceux qui ne le sont pas]
    Args:
        fastaqFile ([fichier fastaq]): [contient des reads]
        genomeFile ([fichier fasta]): [génome de référence]
        parameters ([tuple contenant 3 paramètres]): [1 : longueur de kmer, 2 : taux de faux positif, 3 : pourcentage seuil de décision]
        outfile1 ([fichier txt]): [reads appartenant au génome]
        outfile2 ([fichier txt]): [reads n'appartenant pas au génome]
    """       
    with gzip.open(fastaqFile, 'rt') as f2, open(genomeFile, 'r') as f1, \
    open(outfile1, 'w') as w1, open(outfile2, 'w') as w2 :
        genome = ''
        line1 = f1.readline()
        line1 = f1.readline()
        while line1 != "" :
            genome += line1.rstrip()
            line1 = f1.readline()
        fb = Filter(parameters, genome)
        for line2 in f2 :
            if line2[0] == '@' :
                line2 = f2.readline()
                if line2[0] == '@' :
                    line2 = f2.readline()
                if fb(line2) :
                    w1.write(line2)
                if not fb(line2) :
                    w2.write(line2)
         
class Filter :
    
    def __init__(self, parameters, genome) : # parameters = (kmer_length, fp, threshold)
        """[Réalise les opérations, retourne "True" si le read appartient au genome, "False" sinon]
        Args:
            parameters ([tuple]): [kmer_length : longueur de kmer, fp : pourcentage de faux positif, threshold : pourcentage seuil de décision si le read est trouvé ou pas]
            genome ([fasta]): [génome de référence]
        """        
        self.__parameters = parameters # les paramètres
        self.__kmers_genome_pos = self.k_mers(genome, parameters[0]) # les kmers du brin positif du génome
        self.__genome_neg = self.rev_comp(genome) # le reverse-complement du génome
        self.__kmers_genome_neg = self.k_mers(self.__genome_neg, parameters[0]) # les kmers du reverse-complement du génome
        self.__items_nb = len(self.__kmers_genome_pos) # le nombre d'éléments qu'on veut stocker dans le filtre
        self.__size = self.size_Filter(self.__items_nb, parameters[1]) # la taille du filtre
        self.__hash_nb = self.hash_count_Filter(self.__size, self.__items_nb) # le nombre de fonctions de hachage
        self.__bit_array_pos = bitarray(self.__size) # le filtre de bloom pour le brin positif du génome
        self.__bit_array_neg = bitarray(self.__size) # le filtre de bloom pour le brin négatif du génome
        self.__bit_array_pos.setall(0) # initialisation à 0 des valeurs du filtre pour le brin positif du génome
        self.__bit_array_neg.setall(0) # initialisation à 0 des valeurs du filtre pour le brin négatif du génome
        for kmer_pos in self.__kmers_genome_pos :
            self.add(kmer_pos,"pos") # remplissage du filtre pour le brin positif du génome par les kmers du read
        for kmer_neg in self.__kmers_genome_neg :
            self.add(kmer_neg,"neg") # remplissage du filtre pour le brin négatif du génome par les kmers du read
        
    def add(self, string, sens) :
        """[Ajoute un kmer dans le filtre]
        Args:
            string ([str]): [kmer du read]
            sens ([str]): ["pos" si brin positif du génome, "neg" si brin négatif du génome]
        """        
        if sens == "pos" : # Ajoute le kmer dans le filtre pour le brin positif du génome
            for seed in range(self.__hash_nb) :
                result = mmh3.hash(string, seed)%self.__size
                self.__bit_array_pos[result] = 1
        elif sens == "neg" : # Ajoute le kmer dans le filtre pour le brin négatif du génome
            for seed in range(self.__hash_nb) :
                result = mmh3.hash(string, seed)%self.__size
                self.__bit_array_neg[result] = 1
        
    def lookup(self, string, sens) :
        """[Vérifie si un kmer est présent ou pas dans le filtre]
        Args:
            string ([str]): [kmer du read]
            sens ([str]): ["pos" si brin positif du génome, "neg" si brin négatif du génome]
        Returns:
            [boolean]: ["True" si kmer probablement trouvé dans le filtre, "False" sinon]
        """        
        if sens == "pos" : # Cherche le kmer dans le filtre pour le brin positif du génome
            for seed in range(self.__hash_nb) :
                result = mmh3.hash(string, seed)%self.__size
                if self.__bit_array_pos[result] == 0 :
                    return False
            return True
        elif sens == "neg" : # Cherche le kmer dans le filtre pour le brin négatif du génome
            for seed in range(self.__hash_nb) :
                result = mmh3.hash(string, seed)%self.__size
                if self.__bit_array_neg[result] == 0 :
                    return False
            return True
    
    def size_Filter(self, n, p) :
        """[Calcule la taille du filtre]
        Args:
            n ([int]): [nombre d'éléments qu'on veut stocker dans le filtre]
            p ([float]): [pourcentage de faux positif]
        Returns:
            [int]: [taille du filtre]
        """        
        assert p > 0, "Le pourcentage de faux positif entrée en paramètre doit être supérieur à 0."
        m = -(n * math.log(p)) / (math.log(2)**2)
        return int(m)
    
    def hash_count_Filter(self, m, n) :
        """[Calcule le nombre de fonctions de hachage à utiliser]
        Args:
            m ([int]): [taille du filtre]
            n ([int]): [nombre d'éléments qu'on veut stocker dans le filtre]
        Returns:
            [int]: [nombre de fonctions de hachage]
        """        
        k = (m/n) * math.log(2)
        return int(k)

    def k_mers(self, seq, k) :
        """[Découpe une séquence en kmer]
        Args:
            seq ([str]): [séquence qu'on veut découper en kmer]
            k ([int]): [longueur de kmer]
        Returns:
            [list]: [liste des kmers]
        """        
        list_k_mers = []
        for i in range(len(seq) - k + 1) :
            list_k_mers.append(seq[i:i+k])
        return list_k_mers
    
    def percentage_kmer_read_found(self, read) :
        """[Calcule le pourcentage de kmers du read probablement trouvés dans le filtre]
        Args:
            read ([str]): [read dont on veut vérifier s'il est présent ou pas dans le génome]
        Returns:
            [tuple]: [1er élément : pourcentage de kmers du read trouvé dans le brin positif du génome, 2ème élément : pourcentage de kmers du read trouvé dans le brin négatif du génome]
        """        
        assert self.__parameters[0] < len(read), "Longueur du kmer entrée en paramètre trop long"
        kmers_read = self.k_mers(read, self.__parameters[0])
        nb_kmers_found_pos = 0
        nb_kmers_found_neg = 0
        for kmer in kmers_read :
            if self.lookup(kmer,"pos") :
                nb_kmers_found_pos += 1
            elif self.lookup(kmer,"neg") :
                nb_kmers_found_neg += 1
        p_pos = (nb_kmers_found_pos/len(kmers_read)) - self.__parameters[1]
        p_neg = (nb_kmers_found_neg/len(kmers_read)) - self.__parameters[1]
        return (p_pos,p_neg)

    def rev_comp(self, seq) :
        """[Trouve le reverse-complement d'une séquence d'ADN]
        Args:
            seq ([str]): [séquence dont on veut avoir le reverse-complement]
        Returns:
            [str]: [reverse-complement de la séquence]
        """        
        revcomp = ''
        for i in range(-1,-len(seq)-1,-1) :
            if seq[i] == "A" :
                revcomp += "T"
            elif seq[i] == "T" :
                revcomp += "A"
            elif seq[i] == "C" :
                revcomp += "G"
            elif seq[i] == "G" :
                revcomp += "C"
        return revcomp

    def __call__(self, read) :
        """[Appelle une instance de la classe, réalise les opérations sur read]
        Args:
            read ([str]): [read dont on veut vérifier s'il est présent ou pas dans le génome]
        Returns:
            [boolean]: ["True" si read appartient au génome, "False" sinon]
        """        
        if (self.percentage_kmer_read_found(read))[0] <= self.__parameters[2] and \
            (self.percentage_kmer_read_found(read))[1] <= self.__parameters[2] :
            return False
        else :
            return True

            
if __name__ == "__main__" :
    parameters = (50, 0.01, 0.65)
    operation = readerFastq_writerTxt('reads-1M.fastq.gz', \
    "sequenceSARS.fasta", parameters, "genomeReads.txt", "NoGenomeReads.txt")