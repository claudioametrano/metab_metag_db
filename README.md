
# Metabarcoding and metagenomic databases: taxonomic assignment and ecological metadata
This repository contains the materials for the "Databases in ecology and comparative genomics course": day 2.
It deals with the secondary databases developed to barcode diversity using amplicon sequencing (metabarcoding) and metagenomic data.


### Software required (on HPC)
- QIIME2 (container)
- ...
- FastQC
- MultiQC
### Software required (locally)
- a fasta file reader (MEGA, Aliview, Jalview, Bioedit ...)
- Terminal (Win, Linux and Mac) to ssh into the HPC and a client (e.g. Filezilla) for easy file transfer


### BEFORE WE START
Login to your account on the HPC and start an interactive session (as we won't run long, computationally heavy analyses that require to submit a job to the cluster queue):
```bash
ssh username@l2.gsc1.uni-graz.at

srun --mem=32G --ntasks=8 --cpus-per-task=1 --time=10:00:00 --pty bash
```

Download this repository:
```bash
git clone https://github.com/claudioametrano/metab_metag_db.git
```

Rename the folder containing the results, so you won't overwrite it running the analyses of this tutorial, and create a new results folder
```bash
mv results results_backup  
mkdir results
```

We then need to:
- Bbtain the qiime2 container and create the .sif file Singularity uses
- Start an interactive session in the container (you can also launch a single command  with `singularity run`) and check if QIIME2 is in there, 
```bash
singularity pull docker://quay.io/qiime2/amplicon:2024.10

singularity shell --bind "$(pwd)":/in --home "$(pwd)":/home/qiime2 amplicon_2024.10.sif

qiime --help
```
`qiime` tries to create a small cache under **$HOME** (`/home/qiime2/`).  
Inside a Singularity image the root filesystem is read-only, so we need to mount also the qiime home folder.

It is also possible to install software QiiME2 via Conda, in a dedicated conda environment, via a .yml recipe from QIIME website
```bash
conda env create -n qiime2 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
```


### 1- Nucleotide reference databases for metabarcoding and diversity assessment (an example of secondary database)

Since the introduction of the second-generation sequencing technologies in the mid-2000s (Illumina, formerly Solexa; Roche 454; ABI SOLiD), metabarcoding has unlocked the unprecedented possibility of barcoding diversity, by enabling the simultaneous recovery of thousands of taxonomic signatures in a single run. Yet, the power of this approach depends on comprehensive, accurately curated **reference databases** for the most widely used "unioversal" barcodes (e.g., COI, 16S, ITS), without which the sequence output cannot be reliably linked to taxa identities. Most of them are a reduced, curated (and often clustered) version of NCBI GenBank (or similar databases).

#### Why reference databases are crucial:
- **Taxonomic anchor:** Metabarcoding reads are just strings of bases until they can be matched to a reference sequence; curated databases turn anonymous fragments, grouped into OUTs (Operational Taxonomic Units), into named taxa, giving ecological meaning to the data.
- **Accuracy:** Curated reference libraries increase the accuracy of taxonomic assignment of metabarcoding data diminishing, mis-labelled, chimeric and artifact sequences. This is critical for decision making based on biodiversity trend, or control of biological matrices (e.g. yes, multi-flower honey, but from what flowers?! [plosone](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0134735]); [appl. pl. sci.](https://bsapubs.onlinelibrary.wiley.com/doi/full/10.3732/apps.1400066)).
- **Breadth of coverage:** Comprehensive, standardized region-specific reference set minimizes the microbial “dark-matter” problem ([Nature](https://www.nature.com/articles/nature12352)), making community profiles more complete, and analyses  more reproducible and comparable across studies.
- **Continuous update:** Because taxonomy and sequence diversity keep changing and growing, regularly updated databases provide a mechanism to incorporate new barcodes and revised names, ensuring metabarcoding results stay interpretable over time.

### Common Reference Databases for Metabarcoding

| Marker                            | Taxa Focus                       | Database                           | Link / Notes                                                                  |
| --------------------------------- | -------------------------------- | ---------------------------------- | ----------------------------------------------------------------------------- |
| **COI (rbcL, matk, ITS, 18S)**    | Animals (plants, fungi)          | **BOLD**: Barcode Of Life Database | [bold](https://www.boldsystems.org/)                                          |
| **COI** (additional mito markers) | Amimals (and other eukaryotes)   | **MIDORI**                         | [midori2](https://www.reference-midori.info)                                  |
| **ITS**                           | Fungi (plants, other eukaryotes) | **UNITE**                          | [unite](https://unite.ut.ee/)                                                 |
| **ITS**                           | plants                           | **PLANiTS**                        | https://academic.oup.com/database/article/doi/10.1093/database/baz155/5722079 |
| **18S rRNA, full rDNA**           | Eukaryotes                       | **PR2**                            | [pr2](https://pr2-database.org/)                                              |
| **16S/18S, 23S/28S rRNA**         | Bacteria, Archaea and Eukarya    | **SILVA**                          | [silva](https://www.arb-silva.de/)                                            |
| **16S rRNA**                      | Bacteria, Archaea                | **Greengenes**                     | https://www.nature.com/articles/s41587-023-01845-1                            |
These databases have usually a very simple structure, they are made by one or two files, containing:
- Reference sequences (usually in .fasta format)
- A taxonomy file with the taxonomy associated to each of the representative sequences

### Tools and pipelines commonly used in metabarcoding
After about two decades of metabarcoding there are plenty of tools and pipelines which were developed to analyze metabarcoding data, many of them composed by the same fundamental steps, also often sharing methods and piece of software (e.g. QIIME using DADA2 denoising algorithm)

| Pipeline / Platform | Core language & interface                    | Main approach (ASV vs OTU, multi-marker, etc.)                             | Reference                                                                                                                            |
| ------------------- | -------------------------------------------- | -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| **QIIME 2**         | Python + plugin framework (CLI, GUI, Galaxy) | Flexible; ASV (DADA2/Deblur) or OTU; integrates phylogenetics & statistics | [qiime2](https://qiime2.org/)                                                                                                        |
| **DADA2**           | R package                                    | Denoising → exact **A**mplicon **S**equence **V**ariants (ASVs)            | [dada2](https://benjjneb.github.io/dada2/)                                                                                           |
| **mothur**          | C++ binary with command shell                | OTU clustering (97 %), plus chimera checking                               | [https://mothur.org](https://mothur.org/?utm_source=chatgpt.com "mothur website")                                                    |
| **OBITools**        | Python scripts                               | Read filtering, taxonomic assignment, ecology metrics                      | [obitools](https://pythonhosted.org/OBITools/welcome.html?utm_source=chatgpt.com "Welcome to the OBITools - PythonHosted.org")       |
| **Anacapa Toolkit** | Conda/R + Snakemake                          | Multi-locus (COI, 12S, 18S) with custom reference building                 | [GitHub](https://github.com/limey-bean/Anacapa)                                                                                      |
| **MetaWorks**       | Snakemake workflow (Python/R)                | Multi-locus (ITS, COI) with VSEARCH + tax-assign                           | [GitHub](https://github.com/terrimporter/MetaWorks?utm_source=chatgpt.com "MetaWorks: A Multi-Marker Metabarcode Pipeline - GitHub") |
| **mBRAVE**          | Web cloud platform                           | OTU/ASV assignment against curated BOLD                                    | [mbrave.net](https://www.mbrave.net/?utm_source=chatgpt.com "mBRAVE - Metabarcoding at Scale")                                       |

### A typical metabarcoding experiment workflow
![metabarcoding](/images/metabarcoding_workflow.jpg)
modif. from [Pawlowsky et al. 2018](https://www.sciencedirect.com/science/article/pii/S0048969718316322)

#### 1- Experimental design:
- Question to answer using amplicon sequencing data
- Actual design: How many samples/replicates? How many markers/which organims target? How many libraries? What expected sequencing depth? How large is the budget?
- Metadata: by direct measurments? from public databases (especially for environmental dataset)? 

#### 2- Molecular biology laboratory procedures
![wetlab](/images/illumina.pdf)
from [Illumina](https://www.illumina.com/)
#### 3- Data analysis
This is an overview from [QIIME2](https://amplicon-docs.qiime2.org/en/latest/explanations/conceptual-overview.html) website, but most of these steps are similar no matter what pipeline you select
![qiime](/images/qiime_flow.png)
We are not going to produce our own data this time, we will instead start from metabarcoding data produced for a small marine project.

### ... Let's begin 
![miramare](/images/miramare.png)
The dataset is a toy version of an actual experiment.
We are going to use one of the library (16S rDNA) produced made to assess prokaryotic diversity in a coastal environment.

![16S](/images/16S_rDNA.png)
from [Fukuda et al. 2016](https://www.researchgate.net/publication/308040658_Molecular_Approaches_to_Studying_Microbial_Communities_Targeting_the_16S_Ribosomal_RNA_Gene)
### Needed files:
- Sequences (fastq) ->  ./data/raw_fatsq/16S_biochar_run2_10perc_sampled
- Metadata -> ./data/metadada.csv 
  let's take a look at them to understad the **experimental design**
- Reference database -> [SILVA](https://www.arb-silva.de)

 Let's download the last version of SILVA 99% similarity clustered version 
```bash 
cd data

wget https://www.arb-silva.de/fileadmin/silva_databases/current/ARB_files/SILVA_138.2_SSURef_NR99_03_07_24_opt.arb.gz

cd ..
```

#### **TASK 1**
> - Check on of the fastq file without decompressing them (It would be not convenient, as the software we use can deal with compressed archives)
> - Count the number of sequences per fastq file
> - Which kind of reads are these? (type, reads length)
>    (hint: use zless or zcat, zgrep and awk with length)
>

### Raw reads quality benchmark
Command line is great (he said), quick and versatile, but user friendly, interactive .html quality report are generated by software such as [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://github.com/MultiQC/MultiQC)

We will run them via existing containers: 
FastQC to benchmark each sample
```bash
singularity shell --bind "$(pwd)":/in  https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0

mkdir /in/results/fastqc_raw_out

fastqc /in/data/16S_biochar_run2_10perc_sampled/*.gz -o /in/results/fastqc_raw_out --threads 8 --nogroup

exit
```

MultiQC to aggregate results and highlight possible outlier
```bash
singularity shell --bind "$(pwd)":/in https://depot.galaxyproject.org/singularity/multiqc:1.26--pyhdfd78af_0

multiqc --verbose /in/results/fastqc_raw_out/ -o /in/results/multiqc_raw_out

exit
```
#### **TASK2**
> Now check the .html output of R1 and R2 fastq file from the some samples (download it from server by Filezilla and open them locally) and the aggregated report of MultuQC and try to answer the following questions:
> 1- Which of the warnings thrown by fasQC are of actual concern, and which are not? Why?
> 2- Do you notice any difference in term of quality between corresponding R1 and R2 files?
> 3- Does any sample show peculiar characteristics in term of quality, nucleotide composition, etc?
> 4- If any, what are the over-represented sequences? Are they of concern for subsequent analyses?
> 5- Does your sequence contains residual Illumina adapters/sequincing primers and marker's primers?
> 6- Could you explain why polyG sequences are common? Do they have biological meaning? (hint: look at bottom-right corner of /images/illumina.pdf)

![[P5_to_P7.png]]
from [here](https://teichlab.github.io/scg_lib_structs/methods_html/Illumina.html)
### **TASK3** 
>FastQC do not search for your custom primer (well, not by default) so now it is your turn to do so:
>Primer for 16S rRNA (V3-V4 region) 
>Forward: Pro341F (5’-CCTACGGGNBGCASCAG-3’)
>Reverse: Pro805R (5’-GACTACNVGGGTATCTAATCC-3’)]
>Do you notice anything unusual?
>[IUPAC nucleotide code](https://pmc.ncbi.nlm.nih.gov/articles/PMC2865858/)
>
> **Questions:**
> - Do all sequence have primers? Where in the sequence?
> - Is there any sequence that do not match your search? if so, why?
> (hint: use zgrep with regular expression (-E) and the regex e.g. "\[A | T]" when more than one character is possible, see grep --help. Note that "|" has here a different meaning in regex than it has as a pipe, here it means OR)


### Quality trimming, adapter and universal primer removal
Many software are available (fastp, Trimmomatic, cutadapt etc.) with various option and approaches to trimming, removing adapters and/or primers, even though these reads are of really have quality and denoising algoritm work usually well with raw reads, there is a little room for improvement (let's check after the trimming if it is worth it)

Fastp is one of the most versatile tool to pre-process short reads
# fastp dont support degenerate base, iether use another tool or ctadapt + fastp or simply cut the first x bases with DADA2 in denoising step
```bash
singularity shell --bind "$(pwd)":/in https://depot.galaxyproject.org/singularity/fastp:0.24.0--heae3180_1

mkdir /in/results/trimmed_fastq



IN_DIR="/in/data/16S_biochar_run2_10perc_sampled"
OUT_DIR="/in/results/trimmed_fastq"
for r1 in "$IN_DIR"/*_R1_001.fastq_10perc.fastq.gz; do
    # derive sample prefix by stripping the suffix
    sample=${r1%_R1_001.fastq_10perc.fastq.gz}
    r2=${sample}_R1_001.fastq_10perc.fastq.gz
    
    echo "▶ Trimming sample: $sample"

    fastp \
        -i  "$r1"  -I  "$r2" \
        -o  "$OUT_DIR/${sample}_R1.trimmed.fastq.gz" \
        -O  "$OUT_DIR/${sample}_R2.trimmed.fastq.gz" \
        --adapter_fasta=CCTACGGGNBGCASCAG \
        --adapter_sequence_r2=GACTACNVGGGTATCTAATCC \
        --detect_adapter_for_pe           \
        --trim_poly_g                      \
        --trim_poly_x                      \
        --cut_front --cut_tail             \
        --cut_window_size     4            \
        --cut_mean_quality    20           \
        --qualified_quality_phred 15       \
        --length_required  100     \
        --thread              8   \
        --html  "$OUT_DIR/${sample}.html" \

    echo "✓ Done with $sample"
done
  
```




### Trimmed reads quality 
Run again fastQC and MultiQC (in a different output folder!!) and check what happened.



Possible contaminant adapters removal and universal PCR primer removal
https://www.melbournebioinformatics.org.au/tutorials/tutorials/qiime2/qiime2/
https://github.com/otagoedna/edna_workshop_june2021/tree/master/docs
https://www.slideshare.net/evelienjongepier1/metabarcoding-qiime2-workshop-denoise-249464789#8
https://telatin.github.io/microbiome-bioinformatics/Metabarcoding-1/
https://amplicon-docs.qiime2.org/en/latest/tutorials/gut-to-soil.html
https://use.qiime2.org/en/latest/how-to-guides/merge-metadata.html
https://amplicon-docs.qiime2.org/en/latest/explanations/conceptual-overview.html
https://amplicon-docs.qiime2.org/en/latest/tutorials/gut-to-soil.html








### FINAL TASKS
BUILD A REPORT CONTAINING EVERY STEP YOU TOOK, REPORTING THE COMMANDS USED AND THEIR OUTPUT (meaningful examples are enough if the output is big!): 
#### **FINAL TASK C**

> 1) Pick a metabarcoding study from literature, with the following characteristics:
> - A reasonable amount of samples (metabarcoding surveys can be huge, even though per single sample data are usually quite manageable -> short reads amplicon sequencing).
> - Based on **one** barcode (if multiple barcodes libraries are used in the manuscript, you can select one, and only work on a subset of samples)
> - Clearly explained and reproducible methods 
> - NCBI SRA stored raw sequencing runs
> - Available metadata table 
> - It can be from whatever matrix, you have maximum freedom to select something which intrigues you.
> Some examples: Human gut/oral/skin/... microbiota,  eDNA from water, air (yes, aerobiology does exists), soil, aerosol, plant, fugal, animal microbiota, honey, herbal tea blend, ... 
> 1) Reproduce the main basic steps of a metabarcoding analysis following the main step we highlighted during the tutorial (alpha diversity, ordination. etc).

#### **FINAL TASK D**
> 1) Critically Compare the method and the results you obtained with the methods and results from the paper you selected.
> 
  - Did you use the same tool/pipeline the authors used to perform the analyses?
  - Did you get to their same conclusions?/changing the method noticeably/slightly affected the results? 
  - If that happened what were the main differences?