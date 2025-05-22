
# Metabarcoding and metagenomic databases: taxonomic assignment and ecological metadata
This repository contains the materials for the "Databases in ecology and comparative genomics course": day 2.
It deals with the secondary databases developed to barcode diversity using amplicon sequencing (metabarcoding) and metagenomic data.


### Software required (on HPC)
- QIIME2
- R
- ...
- 
### Software required (locally)
- a fasta file reader (MEGA, Aliview, Jalview, Bioedit ...)
- Terminal (Win, Linux and Mac) to ssh into the HPC and (optionally) GUI file client (e.g. Filezilla)
- ...

### BEFORE WE START
Login to your account on the HPC and start an interactive session (as we won't run long analyses that require to submit a job to the cluster queue):
```bash
ssh username@l2.gsc1.uni-graz.at

srun --mem=8G --ntasks=4 --cpus-per-task=1 --time=10:00:00 --pty bash
```

download this repository:
```bash
git clone https://github.com/claudioametrano/name-of-the-repository.git
```

Rename the folder containing the results, so you won't overwrite it running the analyses of this tutorial, and it will be available if needed:
```bash
mv results results_backup  
```

We will work on the HPC cluster (High Performace Computing) using Docker containers which have the specific software we need, then every command will need this before the actual command that launch that analysis
```bash
docker run --rm -u $(id -u):$(id -g) -v $(pwd):/in -w /in path/to/the/container/name-of-the-container command-to-launch
```
**docker run** ->	Starts a new Docker container.
**--rm** ->	Automatically removes the container once it finishes running. 
**-u \$(id -u):$(id -g)** ->	Runs the container as the current user (user ID and group ID), so the output files are not owned by root.
**-v $(pwd):/in** ->	Mounts the current directory (pwd) from your host into the container at /in. This makes your local files accessible to the container.
**-w /in** ->	Sets the working directory inside the container to /in (which maps to your current directory).

It is also possible to install software via Conda, in this case command of this tutorial can be launched as they are written, but you will need to install Anaconda/Miniconda and the following software:


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
| **16S rRNA**                      |                                  | **Greengenes**                     | https://www.nature.com/articles/s41587-023-01845-1                            |
These databases have usually a very simple structure, they are made by one or two files, containing:
- Reference sequences (usually in .fasta format)
- A taxonomy file with the taxonomy associated to each of the representative sequences

### Task 1


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

### 2- A typical metabarcoding data analysis workflow

![metabarcoding](/images/metabarcoding_workflow.jpg)
modif. from [Pawlowsky et al. 2018](https://www.sciencedirect.com/science/article/pii/S0048969718316322)


![wetlab](/images/wetlab.png)

We are not going to produce our own data this time, we will instead relay on 




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