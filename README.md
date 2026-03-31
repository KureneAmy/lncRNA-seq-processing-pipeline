# lncRNA-seq Processing Pipeline

The analysis pipeline for lncRNA sequencing data processes raw FASTQ files through sequential steps such as initial quality control, adapter trimming, read alignment, transcript/gene quantification, expression profiling, and result summarization, and it supports multiple samples. Also, we provide a fully containerized Singularity environment that bundles all required tools and dependencies, and with a single command, the entire workflow can be executed reproducibly on any compatible system.

# Part I Workflow

Here stands an throughout workflow of data analysis.

```mermaid
flowchart TD
    A[Raw FASTQ files] --> B[raw FastQC]
    B --> C[Adapter trimming<br/>Trim Galore]
    C --> D[FastQC<br/>Trimmed data QC]
    D --> E[Genome alignment<br/>STAR]
    E --> F[BAM processing<br/>samtools sort / index / filter]
    F --> G[Gene / transcript quantification<br/>featureCounts / RSEM]
    G --> H[RNA downstream analysis<br/>count matrix / GO enrichment]
    H --> I[lncRNA downstream analysis<br/>filtering / annotation / DEG / plotting]
    G -. optional .-> J[ERCC quantification]
    J -.-> H

    classDef input fill:#E8F1FF,stroke:#4A78C2,stroke-width:3px,color:#1F3B73;
    classDef qc fill:#FFF4D6,stroke:#D6A100,stroke-width:3px,color:#7A5A00;
    classDef align fill:#E6F7EA,stroke:#2E8B57,stroke-width:3px,color:#155B33;
    classDef quant fill:#F3E8FF,stroke:#8A4FD6,stroke-width:3px,color:#4A227A;
    classDef downstream fill:#f5d7ec,stroke:#d149a7,stroke-width:3px,color:#d65db1;
    classDef optional fill:#F2F2F2,stroke:#888888,stroke-width:3px,stroke-dasharray: 5 3,color:#444444;

    class A,C input;
    class B,D qc;
    class E,F align;
    class G quant;
    class H,I downstream;
    class J optional;
```

# Part II Requirements

1. **Recommended System Configuration**:

      * 8-core CPU
      * 24 GB RAM

2. **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

      * **Step 1: Install System Dependencies**

        ```bash
        # Update package lists and install dependencies
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
            libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        # Download and install Go
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Note: The script navigates to /mnt/share/software.
        # You can change this to your preferred directory for source code.
        cd /mnt/share/software

        # Download the Singularity CE source code
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure the build
        ./mconfig

        # Build Singularity (this can be time-consuming)
        cd builddir
        make

        # Install Singularity to the system
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version

        # Display help information
        singularity -h
        ```

3. **Snakemake**: Snakemake must be installed on your system and requires a Python 3 distribution.

   ```bash
   pip install snakemake
   ```

4. **Reference Data**: A directory containing the reference genome and gene annotation file is required. Depending on your analysis, additional reference files such as transcriptome FASTA or genome index files may also be needed.

   Example:

   ```bash
   mkdir -p reference
   cd reference

   # Genome FASTA
   wget [<genome_fasta_url>](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz)

   # Gene annotation file
   wget [<annotation_gtf_url>](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.primary_assembly.annotation.gtf.gz)

   # Make rRNA data
   awk '$3=="gene" && (/gene_type "rRNA_pseudogene"/ || /gene_type "rRNA"/) {print $1":"$4"-"$5}' gencode.vM38.primary_assembly.annotation.gtf > rRNA.annotation.txt
   singularity exec ../lncRNAseq.sif samtools faidx GRCm39.primary_assembly.genome.fa -r rRNA.annotation.txt > GRCm39.primary_assembly.rRNA.fa

   ### Build indices
   # Genome index
   singularity exec ../lncRNAseq.sif STAR --runThreadN 8 \
      --runMode genomeGenerate \
      --genomeDir ./star_genome_index \
      --genomeFastaFiles GRCm39.primary_assembly.genome.fa \
      --sjdbGTFfile gencode.vM38.primary_assembly.annotation.gtf \
      --sjdbOverhang 149
   # rRNA index
   mkdir -p star_rrna_index
   singularity exec ../lncRNAseq.sif STAR --runThreadN 8 \
      --runMode genomeGenerate \
      --genomeDir ./star_rrna_index \
      --genomeFastaFiles GRCm39.primary_assembly.rRNA.fa
   # RSEM reference
   mkdir -p rsem_index
   singularity exec ../lncRNAseq.sif rsem-prepare-reference --gtf gencode.vM38.primary_assembly.annotation.gtf \
      --num-threads 8 \
      GRCm39.primary_assembly.genome.fa \
      rsem_index/rsem_ref
   # ERCC index (optional)
   singularity exec ../lncRNAseq.sif STAR --runThreadN 8  \
      --runMode genomeGenerate  \
      --genomeDir ./star_ercc_index  \
      --genomeFastaFiles GRCm39.primary_assembly.genome.fa ./ercc/ERCC92.fa \
      --sjdbGTFfile gencode.vM38.primary_assembly.annotation.gtf ./ercc/ERCC92.gtf
   ```

5. **Data Preparation**: The data run by this pipeline should be organized as paired-end FASTQ files. If your data are stored in the SRA database, they can be downloaded and converted to FASTQ format first.

   Example:

   ```bash
   mkdir -p data/samples
   cd data/samples

   prefetch <SRR_ID_1>
   prefetch <SRR_ID_2>

   fastq-dump --gzip --split-files <SRR_ID_1>/<SRR_ID_1>.sra
   fastq-dump --gzip --split-files <SRR_ID_2>/<SRR_ID_2>.sra
   ```

6. **Required File Structure**

   ```bash
   root/
       ├── config.yaml
       ├── data
       │   ├── reference
       │   │   ├── genome.fa
       │   │   └── annotation.gtf
       │   └── samples
       │       ├── sample1_1.fastq.gz
       │       ├── sample1_2.fastq.gz
       │       ├── sample2_1.fastq.gz
       │       └── sample2_2.fastq.gz
       ├── lncRNA-seq.sif
       ├── scripts
       └── lncRNA-seq.smk
   ```

   - **lncRNA-seq.smk** — The main Snakemake workflow script.
   - **config.yaml** — Configuration file containing paths, parameters, and sample information.
     ⚠️ Must be located in the same directory as `lncRNA-seq.smk`.
   - **lncRNA-seq.sif** — Singularity container image with all required software and dependencies pre-installed.
   - **scripts/** — Auxiliary scripts used in the pipeline.

# Part III Running

* **Example code**

   * **Step 1: Edit `config.yaml`**

     ```bash
     # config.yaml
     # Please use ABSOLUTE paths and avoid adding "/" at the end of file folder's path
     
     
     # --- 1. Experimental Design ---
     # ERCC spike-in controls (optional)
     # If you have ERCC spike-ins, set to true and provide the path to other necessary files
     # If not using ERCC, set to false and other fields will be ignored
     ercc:
       use_ercc: true
       ercc_fasta: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/ref/ercc/ERCC92.fa"
       ercc_gtf: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/ref/ercc/ERCC92.gtf"
       ercc_conc: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/ref/ercc/ERCC92_conc.txt"
       ercc_star_index: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/ref/star_ercc_index"
       ercc_analysis_script: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/scripts/ercc_analysis.R"
     
     ercc_params:
       # Specify the path to the ERCC spike-in database
       dilution: 100
       ercc_ul: 2
       rna_pg_per_cell: 30
       total_rna_ug: 1
     
     # Scripts essential for whole pipeline (paths should be ABSOLUTE)
     scripts:
       # R script for differential expression analysis and GO enrichment
       rna_analysis: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/scripts/rna_seq_analysis.R"
       # R script for lncRNA downstream analysis
       lncRNA_analysis: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/scripts/lncRNA_analysis.R"
       # If using ERCC, provide the R script for ERCC analysis in Expermental Design section
     
     
     # --- 2. Run Parameters ---
     # Path to the Singularity container (lncRNAseq.sif)
     # Use ABSOLUTE path
     container: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/scripts/lncRNAseq.sif"
     
     # Output directory for all results
     # Use ABSOLUTE path
     output_dir: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/Output"
     
     # --- 3. Samples ---
     # Sample information with paths to FASTQ files
     # For single-end data: omit "R2" field
     # For paired-end data: include both "R1" and "R2" fields
     # Use ABSOLUTE paths for all FASTQ files
     # 
     # The "condition" field should match your experimental groups
     # and will be used in differential expression analysis
     samples:
       SRR7685881_control_Rep2:
         R1: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/data/SRR7685881_control_Rep2.fastq.gz"
         condition: "Control"
       
       SRR7685880_control_Rep1:
         R1: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/data/SRR7685880_control_Rep1.fastq.gz"
         condition: "Control"
         
       SRR7685879_treated_Rep2:
         R1: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/data/SRR7685879_treated_Rep2.fastq.gz"
         condition: "Treated"
     
       SRR7685878_treated_Rep1:
         R1: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/data/SRR7685878_treated_Rep1.fastq.gz"
         condition: "Treated"
         
     # Remember to keep the file name same with the sample name.
     # Example for single-end data:
     #  Sample_SE:
     #    R1: "/path/to/data/Sample_SE.fastq.gz"
     #    condition: "Control"/...
     #
     # Example for paired-end data:
     #  Sample_PE:
     #    R1: "/path/to/data/Sample_PE_R1.fastq.gz"
     #    R2: "/path/to/data/Sample_PE_R2.fastq.gz"
     #    condition: "Control"/...
     
     
     # --- 4. Reference Files ---
     # All paths should be ABSOLUTE paths
     # See README.md for instructions on preparing reference files
     ref:
       # GTF annotation file
       gtf: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/ref/gencode.vM38.primary_assembly.annotation.gtf"
       
       # Genome FASTA file
       fasta: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/ref/GRCm39.primary_assembly.genome.fa"
       
       # STAR rRNA index
       # Used for removing ribosomal RNA reads
       rrna_star_index: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/ref/star_rrna_index"
       
       # STAR genome index
       star_index: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/ref/star_genome_index"
       
       # RSEM reference prefix
       # This should be the directory containing the RSEM reference files
       # The actual reference name should be 'rsem_ref'
       rsem_index: "/mnt1/4.NAS2025/zhangam/easyomics_test/lncRNA-seq/ref/rsem_index"
     
     
     # --- 5. Analysis Parameters ---
     # Number of CPU threads to use for parallel processing
     # Recommended: 8 cores
     threads: 8
     
     params:
       # Adapter sequence for Trim Galore.
       # Set to "" (empty string) to use auto-detection.
       trim_galore_adapter_sequence: ""
       # --- Align Mode Parameters ---
       # featureCounts strand specificity:
       #   0: unstranded
       #   1: stranded (first-strand)
       #   2: reversely stranded (second-strand, most common for Illumina dUTP methods)
       # Default: 2 (reversely stranded)
       featurecounts_strandedness: 2
       
     
     # --- 6. Differential Expression Analysis Settings ---
     analysis:
       # Contrast for differential expression analysis
       # Format: [treatment_group, control_group]
       # The first group will be compared against the second (control)
       # These names must match the "condition" values in the samples section
       contrast: ["Treated", "Control"]
       
       # FDR (False Discovery Rate) threshold for significance
       # Genes with adjusted p-value < fdr_threshold are considered significant
       fdr_threshold: 0.05
       
       # Log2 fold change threshold for filtering
       # Genes with |log2FoldChange| > lfc_threshold are considered differentially expressed
       lfc_threshold: 1
       
       # Annotation database for GO enrichment analysis
       # For human: "org.Hs.eg.db"
       # For mouse: "org.Mm.eg.db"
       # For other species, see Bioconductor annotation packages
       r_annotation_db: "org.Mm.eg.db"
       
       # Organism name for GO enrichment
       # For human: "hsa"
       # For mouse: "mmu"
       organism: "mmu"
       
       # Whether to perform DESeq2 analysis and the downstream lncRNA-seq analysis
       # If you don't have at least 2 replicated samples for each group, deseq2 should be set to false.
       # If replicated samples are available, you can set to true to perform DESeq2 and GSEA analyses, etc.
       deseq2: true
     ```

   * **Step 2: Dry-run and dag-make**

     Here `/absolute/path/to/root/` represents the root directory.

     ```bash
     # Dry-run
     snakemake -np \
       -s lncRNA-seq.smk \
       --use-singularity \
       --singularity-args "--bind /path/to/root/"

     # Dag-make
     snakemake -s lncRNA-seq.smk \
               --use-singularity \
               --singularity-args "--bind /path/to/root/" \
               --dag | \
     dot -Tsvg > dag.svg
     ```

     Please try dry-run and dag-make first to check pipeline usability and generate a workflow diagram.

   * **Step 3: Run snakemake**

     Here `/path/to/root/` represents the root directory.

     ```bash
     snakemake -s lncRNA-seq.smk \
       --cores 8 \
       --use-singularity \
       --singularity-args "--bind /path/to/root/"
     ```

* **Command Parameters**

   **edit `config.yaml`**
   - `samples`:(required) A list describing all input samples, including sample names and raw FASTQ file paths. Each sample entry should contain:
     - `R1`: path to the Read 1 FASTQ file
     - `R2`: path to the Read 2 FASTQ file
   - `output_dir`:(required) Path to the output directory where all results will be saved.
   - `container`:(required) Path to the Singularity container image containing all required software and dependencies.
   - `reference.genome`:(required) Path to the reference genome FASTA file.
   - `reference.annotation`:(required) Path to the gene annotation file in GTF/GFF format.
   - `threads`:(optional) Number of threads to use. Default: 8.
   - `overwrite`:(optional) Overwrite existing files in output and temporary folders without warning. Default: false.

   **run snakemake**
   - `--use-singularity`: Enables execution of rules within a Singularity container to ensure a fully reproducible environment.
   - `--singularity-args`: Allows passing additional arguments to the Singularity runtime, such as `--bind`.
   - `--cores`: Specifies the maximum number of CPU cores that Snakemake can use in parallel when executing workflow rules.
   - `--bind`: Specifies the directories to be mounted within the Singularity container. Include all required paths such as raw data, scripts, container images, and references.

# Part IV Output

* **Output Structure**

   ```bash
   output_dir/
       ├── qc
       │   ├── multiqc_data
       │   └── multiqc_report.html
       ├── trimmed
       ├── alignments
       ├── counts
       ├── expression
       ├── logs
       └── reports
   ```

* **Output Interpretation**

   - **`multiqc_report.html`**: Open `multiqc_report.html` in a web browser to explore all quality-control metrics interactively.
   - **`trimmed/`**: Adapter-trimmed reads after preprocessing.
   - **`alignments/`**: Genome alignment files in BAM format.
   - **`counts/`**: Gene-level or transcript-level count matrices for downstream analysis.
   - **`expression/`**: Normalized expression tables, including TPM/FPKM or differential expression results if supported by the pipeline.
   - **`reports/`**: Final summary plots and tables generated by the workflow.

# Part V Reference

[1] Schertzer, M. D., Murvin, M. M. and Calabrese, J. M. (2020). Using RNA Sequencing and Spike-in RNAs to Measure Intracellular Abundance of lncRNAs and mRNAs. Bio-protocol 10(19): e3772.

[2] Bora, F. E. (2026). lncRNA Analysis Pipeline (Version 1.0.0) [Computer software]. https://github.com/Cingoz-Lab/lncRNA-Analysis
