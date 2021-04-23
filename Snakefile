
rule filtrage:
    input:
        reads = "Reads/reads-1M.fastq.gz",
        genome = "Genome/genome_SARS_cov2.fasta"
    params:
        sizeKmer = "25",
        falsePositiv = "0.01",
        sortieReadGenome = "SARS.fastq",
        sortieReadNotGenom = "OTHER.fastq"
    output:
        "SARS.fastq"
    shell:"python mainFaster2.py {input} {params.sizeKmer} {params.falsePositiv} {params.sortieReadGenome} {params.sortieReadNotGenom}"

rule assemblage:
    input:
        reads = "SARS.fastq"
    params:
        sizeKmer = "35",
        outfile = "SARS.fasta"
    output:
        "SARS.fasta"
    shell:"python assemblage3.py {input} {params.outfile} {params.sizeKmer}"
