configfile: 'config.json'

import pandas as pd

sample_info = pd.read_csv('fastq/filereport_read_run_PRJNA292073_tsv.txt', sep='\t')
chip = sample_info.loc[sample_info['sample_title'].str.startswith('ChIP-Seq'),]
chip = chip.loc[chip['fastq_ftp'].str.endswith('_2.fastq.gz'),]
srr = chip['run_accession'].values

rule all:
    input:
        expand('fastp_fq/{sample}_r1_trimmed.fq.gz', sample=srr),
        expand('fastp_fq/{sample}_r2_trimmed.fq.gz', sample=srr),
        expand('fastp_report/{sample}_fastp.html', sample=srr),
        expand('chromap_output/{sample}_aln.bed.gz', sample=srr),
        expand('chromap_output/{sample}_summary.txt', sample=srr),
        expand('macs2_output/{sample}/macs2.done', sample=srr),
        expand('macs2_output/{sample}/{sample}_treat_pileup.bw', sample=srr)

rule trim:
    input:
        r1='fastq/{sample}_1.fastq.gz',
        r2='fastq/{sample}_2.fastq.gz'
    output:
        r1='fastp_fq/{sample}_r1_trimmed.fq.gz',
        r2='fastp_fq/{sample}_r2_trimmed.fq.gz',
        qc='fastp_report/{sample}_fastp.html'
    log:
        'log/fastp/{sample}.log'
    threads: 16
    shell:
        ''' fastp -w {threads} -l 25 -x -g \
            --detect_adapter_for_pe \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            -h {output.qc} -j /dev/null > {log}
        '''

rule align:
    input:
        r1='fastp_fq/{sample}_r1_trimmed.fq.gz',
        r2='fastp_fq/{sample}_r2_trimmed.fq.gz'
    output:
        'chromap_output/{sample}_aln.bed.gz',
        qc='chromap_output/{sample}_summary.txt'
    log:
        'log/chromap/{sample}.log'
    threads: 16
    params:
        idx=config['idx'],
        fa=config['fa'],
        chromapother=config['chromapother']
    shell:
        ''' chromap -t {threads} -x {params.idx} -r {params.fa} \
            {params.chromapother} --summary {output.qc} \
            -1 {input.r1} -2 {input.r2} \
            -o chromap_output/{wildcards.sample}_aln.bed \
            > {log} && gzip chromap_output/{wildcards.sample}_aln.bed
        '''

rule peakCalling:
    input:
        'chromap_output/{sample}_aln.bed.gz'
    output:
        touch('macs2_output/{sample}/macs2.done')
    params:
        gsize=config['gsize'],
        broad=config['broad']
    log:
        'log/macs2/{sample}.log'
    shell:
        ''' macs2 callpeak -t {input} -f BEDPE -g {params.gsize} \
            {params.broad} -B --SPMR --nomodel -q 0.01 --keep-dup all \
            --outdir macs2_output/{wildcards.sample} \
            -n {wildcards.sample}
        '''

rule bdgToBw:
    input:
        'macs2_output/{sample}/macs2.done'
    output:
        'macs2_output/{sample}/{sample}_treat_pileup.bw'
    params:
        chromsize=config['chromsize']
    shell:
        ''' bdg2bw macs2_output/{wildcards.sample}/{wildcards.sample}_treat_pileup.bdg \
            {params.chromsize}
        '''
