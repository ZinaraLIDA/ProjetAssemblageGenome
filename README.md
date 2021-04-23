## RINGEVAL Allan - LIDAMAHASOLO Zinara

---

## Executer avec snakemake

> Il faut avoir créer le dossier Reads avec le fichier "reads.fastq" ainsi que le dossier
Genome avec le fichier "genome_SARS_cov2.fasta"

> Commande : snakemake -j 1 SARS.fasta

## Executer sans snakemake

### Pour le filtrage

> Exemple commande : python mainFaster2.py Reads/reads.fastq.gz Genome/genome_SARS_cov2.fasta 25 0.01 sortieReadGenome.fastq sortieReadNotGenom.fastq

---

> Syntax : mainFaster2.py [ReadsFile.fastq.gz] [GenomeFile.fasta] [SizeKmer] [FalseNegative] [outputReadGenome.fastq] [outputReadNotGenome.fastq]

### Pour l'assemblage

> Exemple commande : python assemblage3.py ReadFilterGenome.fastq outfile.fasta 35

---

> Syntax : assemblage3.py [sortieReadGenome.fastq] [outfile.fasta] [Sizekmer]


