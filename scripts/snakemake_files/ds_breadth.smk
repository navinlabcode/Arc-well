samples, = glob_wildcards("data/{sample}.bam")
samtools_path="/volumes/seq/code/3rd_party/samtools/samtools-1.10/samtools"

rule all:
    input:
        expand('covfile/{sample}.covhist.txt', sample=samples)

rule sambamba_markdup:
    input:
        "data/{sample}.bam",
    output:
        temp("marked/{sample}.bam")
    threads: 4
    shell:
        "/volumes/seq/code/3rd_party/git/sambamba-0.8.1/sambamba-0.8.1-linux-amd64-static markdup -r -t {threads} {input} {output}"

rule save_header:
    input: 
        "marked/{sample}.bam"
    output:
        temp("temp_header/{sample}.sam")
    shell:
        "{samtools_path} view {input} -H > {output}"

rule downsampling:
    input:
        "marked/{sample}.bam"
    output:
        temp("temp_downsampled/{sample}.sam")
    params:
        reads=500000
    shell:
        "{samtools_path} view {input} | shuf -n {params.reads} --random-source={input} >> {output}"

rule add_header:
    input: 
        header = 'temp_header/{sample}.sam',
        reads = 'temp_downsampled/{sample}.sam'
    output:
        temp('downsampled/{sample}.sam')
    shell:
        "cat {input.header} {input.reads} > {output}"

rule sort:
    input: 
        "downsampled/{sample}.sam"
    output:
        temp("sort/{sample}.bam")
    threads: 4
    shell:
        "{samtools_path} sort {input} -@ {threads} -o {output}"

rule coverage:
    input:
        "sort/{sample}.bam"
    output: 
        "covfile/{sample}.covhist.txt"
    params:
        read_size=50
    shell:
        "/volumes/lab/users/dminussi/software/bedtools2/bin/genomeCoverageBed -ibam {input} -fs {params.read_size} -g /volumes/seq/code/3rd_party/bedtools2-master/genomes/human.hg19.genome > {output} "
