# ChIP-seq Data Preprocessing

This repository contains detailed instructions that we often use in the lab to perform basic ChIP-seq preprocessing routines, which is an update to our original protocol published in STAR Protocols in 2020 (doi: [10.1016/j.xpro.2020.100187](https://doi.org/10.1016/j.xpro.2020.100187)). A minimum example data set is provided as fastq files under the directory `fastq`. The data contains 1 million sequence reads pairs from BATF ChIP and Input from mouse T helper cells. The data is taken from [ERP105570](https://www.ebi.ac.uk/ena/browser/view/ERP105570).

In order to run the pipeline, we need [chromap](https://github.com/haowenz/chromap) and [MACS2](https://pypi.org/project/MACS2/) installed in your system. On top of that, you also need the `bedGraphToBigWig`, `bedClip` and `fetchChromSizes` executables from [UCSC utilities](https://hgdownload.soe.ucsc.edu/admin/exe/). I have put them under the `misc` directory. Make sure they are in your `PATH` so that they can be called directly from the command line.


## Building Index For The Reference Genome

We need the genome of interest in `fasta` format, and we should have it with us. This workflow requires the file in placed under the `genome/{build}/fasta` directory.  If not, the genome sequences in `fasta` format can be downloaded from the [UCSC genome browser](https://genome.ucsc.edu) or [Ensembl](https://www.ensembl.org/index.html).

Use the mouse genome build `mm10` as an example, we can download the `mm10.fa.gz` from [this UCSC page](https://hgdownload.gi.ucsc.edu/goldenPath/mm10/bigZips/). Unzip the file and place the `mm10.fa` under `genome/mm10/fasta`. In addition, we also need the chromosome size file, which is a tab-delimited file containing the size in bp of each chromosome. This file is needed in the future, and we can get it by simply doing `fetchChromSizes mm10 > misc/mm10.chrom.sizes`.

Now we can build the genome index using `chromap` with 4 cores:

```console
chromap \
    -i -t 4 \
    -r genome/mm10/fasta/mm10.fa \
    -o genome/mm10/chromap_index/mm10_chromap.index
```

## Sequence Read Alignment

Once we have the genome index, we can align the reads to the reference genome:

```console
mkdir -p chromap_output

# map the chip sample
chromap \
    -t 4 \
    -x genome/mm10/chromap_index/mm10_chromap.index \
    -r genome/mm10/fasta/mm10.fa \
    --preset chip --trim-adapters \
    --summary chromap_output/BATF_ChIP_summary.csv \
    -1 fastq/BATF_ERR2213792_1M_r1.fq.gz \
    -2 fastq/BATF_ERR2213792_1M_r2.fq.gz \
    -o /dev/stdout | gzip > chromap_output/BATF_ChIP_fragments.bed.gz

# map the input control sample
chromap \
    -t 4 \
    -x genome/mm10/chromap_index/mm10_chromap.index \
    -r genome/mm10/fasta/mm10.fa \
    --preset chip --trim-adapters \
    --summary chromap_output/Input_summary.csv \
    -1 fastq/Input_ERR2213811_1M_r1.fq.gz \
    -2 fastq/Input_ERR2213811_1M_r2.fq.gz \
    -o /dev/stdout | gzip > chromap_output/Input_fragments.bed.gz
```

## Peak Calling

Once we have the mapped fragment files, we can perform the peak calling to identify the binding events. In this example, BATF is a transcription factor with punctate binding patterns. Therefore, we just stick to the normal narrowPeak mode:

```console
mkdir -p macs2_output

macs2 callpeak \
    -t chromap_output/BATF_ChIP_fragments.bed.gz \
    -c chromap_output/Input_fragments.bed.gz \
    -f BEDPE -g mm \
    -B --SPMR --nomodel -q 0.01 --keep-dup all \
    --outdir macs2_output \
    -n BATF
```

Once the peak calling is done, we should get a `BATF_treat_pileup.bdg` file under the `macs2_output` directory. We can generate a `bigWig` file from here to visualise the binding signal via a genome browser:

```
bdg2bw macs2_output/BATF_treat_pileup.bdg misc/mm10.chrom.sizes
```

If everything finishes without any problems, our directory should look like this:

```
.
├── chromap_output
│   ├── BATF_ChIP_fragments.bed.gz
│   ├── BATF_ChIP_summary.csv
│   ├── Input_fragments.bed.gz
│   └── Input_summary.csv
├── fastq
│   ├── BATF_ERR2213792_1M_r1.fq.gz
│   ├── BATF_ERR2213792_1M_r2.fq.gz
│   ├── Input_ERR2213811_1M_r1.fq.gz
│   └── Input_ERR2213811_1M_r2.fq.gz
├── genome
│   └── mm10
│       ├── chromap_index
│       │   └── mm10_chromap.index
│       └── fasta
│           └── mm10.fa
├── macs2_output
│   ├── BATF_control_lambda.bdg
│   ├── BATF_peaks.narrowPeak
│   ├── BATF_peaks.xls
│   ├── BATF_summits.bed
│   ├── BATF_treat_pileup.bdg
│   └── BATF_treat_pileup.bw
└── misc
    └── mm10.chrom.sizes

8 directories, 17 files
```

## Citation

Xu, W., Ye, Y., Sharrocks, A.D., Zhang, W., and Chen, X. (2020). Genome-wide Interrogation of Protein-DNA Interactions in Mammalian Cells Using ChIPmentation. STAR Protocols 1, 100187. https://doi.org/10.1016/j.xpro.2020.100187