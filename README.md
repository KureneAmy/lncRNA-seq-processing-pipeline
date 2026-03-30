# lncRNA-seq Processing Pipeline

The analysis pipeline for lncRNA sequencing data processes raw FASTQ files through sequential steps such as initial quality control, adapter trimming, read alignment, transcript/gene quantification, expression profiling, and result summarization, and it supports multiple samples. Also, we provide a fully containerized Singularity environment that bundles all required tools and dependencies, and with a single command, the entire workflow can be executed reproducibly on any compatible system.

# Part I Workflow

Here stands an throughout workflow of data analysis.

<img width="1731" height="655" alt="workflow" src="https://github.com/user-attachments/assets/placeholder-workflow-image" />

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
   wget <genome_fasta_url>

   # Gene annotation file
   wget <annotation_gtf_url>
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
   - **scripts/** — Auxiliary scripts used in the pipeline, if any.

# Part III Running

* **Example code**

   * **Step 1: Edit `config.yaml`**

     ```bash
     # Please use ABSOLUTE PATH for all file paths
     samples:
       sample1:
         R1: "/absolute/path/to/data/samples/sample1_1.fastq.gz"
         R2: "/absolute/path/to/data/samples/sample1_2.fastq.gz"
       sample2:
         R1: "/absolute/path/to/data/samples/sample2_1.fastq.gz"
         R2: "/absolute/path/to/data/samples/sample2_2.fastq.gz"

     # Output directory, no "/" at the end
     output_dir: "/absolute/path/to/output"

     # Singularity container
     container: "/absolute/path/to/lncRNA-seq.sif"

     # Reference files
     reference:
       genome: "/absolute/path/to/data/reference/genome.fa"
       annotation: "/absolute/path/to/data/reference/annotation.gtf"

     # Analysis settings
     threads: 8
     overwrite: true
     ```

   * **Step 2: Dry-run and dag-make**

     Here `/absolute/path/to/root/` represents the root directory.

     ```bash
     # Dry-run
     snakemake -np \
       -s lncRNA-seq.smk \
       --use-singularity \
       --singularity-args "--bind /absolute/path/to/root/"

     # Dag-make
     snakemake -s lncRNA-seq.smk \
               --use-singularity \
               --singularity-args "--bind /absolute/path/to/root/" \
               --dag | \
     dot -Tsvg > dag.svg
     ```

     Please try dry-run and dag-make first to check pipeline usability and generate a workflow diagram.

   * **Step 3: Run snakemake**

     Here `/absolute/path/to/root/` represents the root directory.

     ```bash
     snakemake -s lncRNA-seq.smk \
       --cores 8 \
       --use-singularity \
       --singularity-args "--bind /absolute/path/to/root/"
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
