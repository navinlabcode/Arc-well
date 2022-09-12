samples, = glob_wildcards("data/{sample}.bam")
samtools_path="/volumes/seq/code/3rd_party/samtools/samtools-1.10/samtools"

rule all:
    input:
        expand('downsampled_sort_1M/{sample}.bam', sample=samples),
        expand('downsampled_sort_750k/{sample}.bam', sample=samples),
        expand('downsampled_sort_500k/{sample}.bam', sample=samples),
        expand('downsampled_sort_250k/{sample}.bam', sample=samples),
        expand('downsampled_sort_125k/{sample}.bam', sample=samples),
        expand('downsampled_sort_75k/{sample}.bam', sample=samples),
        expand('downsampled_sort_50k/{sample}.bam', sample=samples)

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

# rule downsampling_1M:
#     input:
#         "marked/{sample}.bam"
#     output:
#         temp("temp_downsampled_1M/{sample}.sam")
#     params:
#         reads=1000000
#     shell:
#         "{samtools_path} view {input} | shuf -n {params.reads} --random-source={input} >> {output}"

# rule add_header_1M:
#     input: 
#         header = 'temp_header/{sample}.sam',
#         reads = 'temp_downsampled_1M/{sample}.sam'
#     output:
#         temp('downsampled_1M/{sample}.sam')
#     shell:
#         "cat {input.header} {input.reads} > {output}"

# rule sort_1M:
#     input: 
#         "downsampled_1M/{sample}.sam"
#     output:
#         "downsampled_sort_1M/{sample}.bam"
#     threads: 4
#     shell:
#         "{samtools_path} sort {input} -@ {threads} -o {output}"

# rule downsampling_750k:
#     input:
#         "marked/{sample}.bam"
#     output:
#         temp("temp_downsampled_750k/{sample}.sam")
#     params:
#         reads=750000
#     shell:
#         "{samtools_path} view {input} | shuf -n {params.reads} --random-source={input} >> {output}"

# rule add_header_750k:
#     input: 
#         header = 'temp_header/{sample}.sam',
#         reads = 'temp_downsampled_750k/{sample}.sam'
#     output:
#         temp('downsampled_750k/{sample}.sam')
#     shell:
#         "cat {input.header} {input.reads} > {output}"

# rule sort_750k:
#     input: 
#         "downsampled_750k/{sample}.sam"
#     output:
#         "downsampled_sort_750k/{sample}.bam"
#     threads: 4
#     shell:
#         "{samtools_path} sort {input} -@ {threads} -o {output}"

# rule downsampling_500k:
#     input:
#         "marked/{sample}.bam"
#     output:
#         temp("temp_downsampled_500k/{sample}.sam")
#     params:
#         reads=500000
#     shell:
#         "{samtools_path} view {input} | shuf -n {params.reads} --random-source={input} >> {output}"

# rule add_header_500k:
#     input: 
#         header = 'temp_header/{sample}.sam',
#         reads = 'temp_downsampled_500k/{sample}.sam'
#     output:
#         temp('downsampled_500k/{sample}.sam')
#     shell:
#         "cat {input.header} {input.reads} > {output}"

# rule sort_500k:
#     input: 
#         "downsampled_500k/{sample}.sam"
#     output:
#         "downsampled_sort_500k/{sample}.bam"
#     threads: 4
#     shell:
#         "{samtools_path} sort {input} -@ {threads} -o {output}"

# rule downsampling_250k:
#     input:
#         "marked/{sample}.bam"
#     output:
#         temp("temp_downsampled_250k/{sample}.sam")
#     params:
#         reads=250000
#     shell:
#         "{samtools_path} view {input} | shuf -n {params.reads} --random-source={input} >> {output}"

# rule add_header_250k:
#     input: 
#         header = 'temp_header/{sample}.sam',
#         reads = 'temp_downsampled_250k/{sample}.sam'
#     output:
#         temp('downsampled_250k/{sample}.sam')
#     shell:
#         "cat {input.header} {input.reads} > {output}"

# rule sort_250k:
#     input: 
#         "downsampled_250k/{sample}.sam"
#     output:
#         "downsampled_sort_250k/{sample}.bam"
#     threads: 4
#     shell:
#         "{samtools_path} sort {input} -@ {threads} -o {output}"

# rule downsampling_125k:
#     input:
#         "marked/{sample}.bam"
#     output:
#         temp("temp_downsampled_125k/{sample}.sam")
#     params:
#         reads=125000
#     shell:
#         "{samtools_path} view {input} | shuf -n {params.reads} --random-source={input} >> {output}"

# rule add_header_125k:
#     input: 
#         header = 'temp_header/{sample}.sam',
#         reads = 'temp_downsampled_125k/{sample}.sam'
#     output:
#         temp('downsampled_125k/{sample}.sam')
#     shell:
#         "cat {input.header} {input.reads} > {output}"

# rule sort_125k:
#     input: 
#         "downsampled_125k/{sample}.sam"
#     output:
#         "downsampled_sort_125k/{sample}.bam"
#     threads: 4
#     shell:
#         "{samtools_path} sort {input} -@ {threads} -o {output}"

rule downsampling_75k:
    input:
        "marked/{sample}.bam"
    output:
        temp("temp_downsampled_75k/{sample}.sam")
    params:
        reads=75000
    shell:
        "{samtools_path} view {input} | shuf -n {params.reads} --random-source={input} >> {output}"

rule add_header_75k:
    input: 
        header = 'temp_header/{sample}.sam',
        reads = 'temp_downsampled_75k/{sample}.sam'
    output:
        temp('downsampled_75k/{sample}.sam')
    shell:
        "cat {input.header} {input.reads} > {output}"

rule sort_75k:
    input: 
        "downsampled_75k/{sample}.sam"
    output:
        "downsampled_sort_75k/{sample}.bam"
    threads: 4
    shell:
        "{samtools_path} sort {input} -@ {threads} -o {output}"

rule downsampling_50k:
    input:
        "marked/{sample}.bam"
    output:
        temp("temp_downsampled_50k/{sample}.sam")
    params:
        reads=50000
    shell:
        "{samtools_path} view {input} | shuf -n {params.reads} --random-source={input} >> {output}"

rule add_header_50k:
    input: 
        header = 'temp_header/{sample}.sam',
        reads = 'temp_downsampled_50k/{sample}.sam'
    output:
        temp('downsampled_50k/{sample}.sam')
    shell:
        "cat {input.header} {input.reads} > {output}"

rule sort_50k:
    input: 
        "downsampled_50k/{sample}.sam"
    output:
        "downsampled_sort_50k/{sample}.bam"
    threads: 4
    shell:
        "{samtools_path} sort {input} -@ {threads} -o {output}"