# lncRNA-seq_pipeline.smk

configfile: "config.yaml"

import os
import re
from pathlib import Path

# ===== 配置 =====
ERCC = bool(config["ercc"].get("use_ercc", False))
DESEQ2 = bool(config["analysis"].get("deseq2", True))
REPORT_ENABLE = bool(config.get("report", {}).get("enable", True))
SAMPLES = list(config["samples"].keys())
CONTRAST = config["analysis"]["contrast"]
OUTPUT_DIR = config["output_dir"]
USE_RSEM = config["analysis"].get("use_rsem", False)
if USE_RSEM:
    RSEM_GRP = config["ref"].get("rsem_grp")
    RSEM_REF_NAME = str(RSEM_GRP)[:-4]
    if not RSEM_GRP:
        raise ValueError("Missing required config: ref.rsem_grp (absolute path to .grp file)")
    if not os.path.isabs(str(RSEM_GRP)):
        raise ValueError(f"ref.rsem_grp must be an absolute path, got: {RSEM_GRP}")
    if not str(RSEM_GRP).endswith(".grp"):
        raise ValueError(f"ref.rsem_grp must be an absolute path ending with '.grp', got: {RSEM_GRP}")
else:
    RSEM_GRP = None





# ===== 样本分类 =====
PAIRED_SAMPLES = [s for s in SAMPLES if config["samples"][s].get("R2")]
SINGLE_SAMPLES = [s for s in SAMPLES if not config["samples"][s].get("R2")]

def build_sample_constraint(samples):
    if not samples:
        return r"(?!)"
    return "|".join(re.escape(s) for s in samples)

# ===== 辅助函数 =====
def get_fastq_r1(wildcards):
    return config["samples"][wildcards.sample]["R1"]

def get_fastq_r2(wildcards):
    return config["samples"][wildcards.sample].get("R2", "")

def is_paired_end(wildcards):
    return bool(config["samples"][wildcards.sample].get("R2"))

def any_paired_end():
    return any(s.get("R2") for s in config["samples"].values())

def sample_prefix(wildcards):
    # QC 输出文件统一由 config 里的 prefix 决定，不依赖输入 fastq 的文件名。
    return str(config["samples"][wildcards.sample].get("prefix", wildcards.sample))

def trimmed_paired_r1(wildcards):
    return f"{OUTPUT_DIR}/qc/trimmed/{sample_prefix(wildcards)}_val_1.fq.gz"

def trimmed_paired_r2(wildcards):
    return f"{OUTPUT_DIR}/qc/trimmed/{sample_prefix(wildcards)}_val_2.fq.gz"

def trimmed_single_r1(wildcards):
    return f"{OUTPUT_DIR}/qc/trimmed/{sample_prefix(wildcards)}_trimmed.fq.gz"

def get_final_outputs(wildcards):
    outputs = []
    
    # QC 输出
    if PAIRED_SAMPLES:
        outputs.extend([
            f"{OUTPUT_DIR}/logs/{sample}_fastqc_raw_paired.log" for sample in PAIRED_SAMPLES
        ])
        outputs.extend([
            f"{OUTPUT_DIR}/logs/{sample}_fastqc_trimmed_paired.log" for sample in PAIRED_SAMPLES
        ])
    
    if SINGLE_SAMPLES:
        outputs.extend([
            f"{OUTPUT_DIR}/logs/{sample}_fastqc_raw_single.log" for sample in SINGLE_SAMPLES
        ])
        outputs.extend([
            f"{OUTPUT_DIR}/logs/{sample}_fastqc_trimmed_single.log" for sample in SINGLE_SAMPLES
        ])
    
    # BAM 处理输出
    outputs.extend([
        f"{OUTPUT_DIR}/results/{sample}.unique.dedup.bam" for sample in SAMPLES
    ])
    
    # 定量输出
    if USE_RSEM:
        outputs.extend([
            f"{OUTPUT_DIR}/intermediate/align/rsem_temp/{sample}.genes.results" for sample in SAMPLES
        ])
    
    # DEA 和 GO 分析输出
    outputs.extend([
        f"{OUTPUT_DIR}/results/{CONTRAST[0]}_vs_{CONTRAST[1]}.deseq2_results.csv",
        f"{OUTPUT_DIR}/results/{CONTRAST[0]}_vs_{CONTRAST[1]}.volcano_plot.png",
        f"{OUTPUT_DIR}/results/{CONTRAST[0]}_vs_{CONTRAST[1]}.GO_enrichment.csv",
        f"{OUTPUT_DIR}/results/{CONTRAST[0]}_vs_{CONTRAST[1]}.GO_dotplot.png",
        f"{OUTPUT_DIR}/results/featureCounts.gene_counts.symbol.tsv",
        f"{OUTPUT_DIR}/results/multiqc_report.html",
    ])
    
    # ERCC 输出
    if ERCC:
        outputs.extend([
            f"{OUTPUT_DIR}/results/{sample}_ERCC_standard_curve.pdf" for sample in SAMPLES
        ])
        outputs.extend([
            f"{OUTPUT_DIR}/results/{sample}_mpc.txt" for sample in SAMPLES
        ])

    if DESEQ2:
        outputs.append(f"{OUTPUT_DIR}/lncRNA_analysis_output/tables/QC_summary.csv")

    # 报告输出
    if REPORT_ENABLE:
        outputs.append(f"{OUTPUT_DIR}/report/report.html")
        outputs.append(f"{OUTPUT_DIR}/report/report.md")
    
    return outputs

# ===== 主规则 =====
rule all:
    input:
        get_final_outputs

# ===== 配对端 QC 规则 =====
rule fastqc_raw_paired:
    input:
        r1=lambda wildcards: get_fastq_r1(wildcards),
        r2=lambda wildcards: get_fastq_r2(wildcards),
    output:
        log = f"{OUTPUT_DIR}/logs/{{sample}}_fastqc_raw_paired.log",
    threads: 2
    container: config["container"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/qc/raw
        fastqc "{input.r1}" "{input.r2}" -t {threads} -o {OUTPUT_DIR}/qc/raw > {output.log} 2>&1
        """

rule trim_galore_paired:
    input:
        log = f"{OUTPUT_DIR}/logs/{{sample}}_fastqc_raw_paired.log",
    output:
        log = f"{OUTPUT_DIR}/logs/{{sample}}_trim_galore_paired.log",
    params:
        r1=lambda w: get_fastq_r1(w),
        r2=lambda w: get_fastq_r2(w),
        basename=lambda w: sample_prefix(w),
        adapter_flag=lambda w: f"--adapter {config['params']['trim_galore_adapter_sequence']}" 
                              if config['params'].get('trim_galore_adapter_sequence') else "",
    threads: 4
    container: config["container"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/qc/trimmed
        trim_galore --cores {threads} {params.adapter_flag} --gzip --paired \
            -o {OUTPUT_DIR}/qc/trimmed "{params.r1}" "{params.r2}" --basename "{params.basename}" > {output.log} 2>&1
        """

rule fastqc_trimmed_paired:
    input:
        log = f"{OUTPUT_DIR}/logs/{{sample}}_trim_galore_paired.log",
    output:
        log = f"{OUTPUT_DIR}/logs/{{sample}}_fastqc_trimmed_paired.log",
    params:
        r1=lambda w: trimmed_paired_r1(w),
        r2=lambda w: trimmed_paired_r2(w),
    threads: 2
    container: config["container"]
    shell:
        """
        fastqc "{params.r1}" "{params.r2}" -t {threads} -o {OUTPUT_DIR}/qc/trimmed > {output.log} 2>&1
        """

rule bowtie2_remove_rrna_paired:
    input:
        log = f"{OUTPUT_DIR}/logs/{{sample}}_fastqc_trimmed_paired.log",
    output:
        r1=f"{OUTPUT_DIR}/intermediate/rrna_removed/{{sample}}.unmapped.R1.fastq.gz",
        r2=f"{OUTPUT_DIR}/intermediate/rrna_removed/{{sample}}.unmapped.R2.fastq.gz",
    params:
        prefix=f"{OUTPUT_DIR}/intermediate/rrna_removed/{{sample}}.",
        rrna_bt2_index=config["ref"]["rrna_bt2_index"],
        r1=lambda w: trimmed_paired_r1(w),
        r2=lambda w: trimmed_paired_r2(w),
    threads: config["threads"]
    container: config["container"]
    log:
        f"{OUTPUT_DIR}/logs/{{sample}}_bowtie2_rrna_removal_paired.log"
    shell:
        """
        bowtie2 -x {params.rrna_bt2_index} \
            -1 "{params.r1}" -2 "{params.r2}" \
            -p {threads} \
            --very-sensitive-local \
            -k 1 \
            --no-unal \
            --un-conc-gz {params.prefix}unmapped.fastq.gz \
            -S /dev/null \
            > {log} 2>&1
        """

# ===== 单端 QC 规则 =====
rule fastqc_raw_single:
    input:
        r1=lambda wildcards: get_fastq_r1(wildcards),
    output:
        log = f"{OUTPUT_DIR}/logs/{{sample}}_fastqc_raw_single.log",
    threads: 2
    container: config["container"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/qc/raw
        fastqc "{input.r1}" -t {threads} -o {OUTPUT_DIR}/qc/raw > {output.log} 2>&1
        """

rule trim_galore_single:
    input:
        log = f"{OUTPUT_DIR}/logs/{{sample}}_fastqc_raw_single.log",
    output:
        log = f"{OUTPUT_DIR}/logs/{{sample}}_trim_galore_single.log",
    params:
        r1=lambda w: get_fastq_r1(w),
        basename=lambda w: sample_prefix(w),
        adapter_flag=lambda w: f"--adapter {config['params']['trim_galore_adapter_sequence']}" 
                              if config['params'].get('trim_galore_adapter_sequence') else "",
    threads: 4
    container: config["container"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/qc/trimmed
        trim_galore --cores {threads} {params.adapter_flag} --gzip \
            -o {OUTPUT_DIR}/qc/trimmed "{params.r1}" --basename "{params.basename}" > {output.log} 2>&1
        """

rule fastqc_trimmed_single:
    input:
        log = f"{OUTPUT_DIR}/logs/{{sample}}_trim_galore_single.log",
    output:
        log = f"{OUTPUT_DIR}/logs/{{sample}}_fastqc_trimmed_single.log",
    params:
        r1=lambda w: trimmed_single_r1(w),
    threads: 2
    container: config["container"]
    shell:
        """
        fastqc "{params.r1}" -t {threads} -o {OUTPUT_DIR}/qc/trimmed > {output.log} 2>&1
        """

rule bowtie2_remove_rrna_single:
    input:
        log = f"{OUTPUT_DIR}/logs/{{sample}}_fastqc_trimmed_single.log",
    output:
        r1=f"{OUTPUT_DIR}/intermediate/rrna_removed/{{sample}}.unmapped.fastq.gz",
    params:
        prefix=f"{OUTPUT_DIR}/intermediate/rrna_removed/{{sample}}.",
        rrna_bt2_index=config["ref"]["rrna_bt2_index"],
        r1=lambda w: trimmed_single_r1(w),
    threads: config["threads"]
    container: config["container"]
    log:
        f"{OUTPUT_DIR}/logs/{{sample}}_bowtie2_rrna_removal_single.log"
    shell:
        """
        bowtie2 -x {params.rrna_bt2_index} \
            -U "{params.r1}" \
            -p {threads} \
            --very-sensitive-local \
            -k 1 \
            --no-unal \
            --un-gz {params.prefix}unmapped.fastq.gz \
            -S /dev/null \
            > {log} 2>&1
        """

# ===== 对齐规则（配对端）=====
rule star_align_paired:
    input:
        r1=f"{OUTPUT_DIR}/intermediate/rrna_removed/{{sample}}.unmapped.R1.fastq.gz",
        r2=f"{OUTPUT_DIR}/intermediate/rrna_removed/{{sample}}.unmapped.R2.fastq.gz",
    output:
        genome_bam=f"{OUTPUT_DIR}/intermediate/align/genome_bam/{{sample}}.Aligned.out.bam",
        transcript_bam=f"{OUTPUT_DIR}/intermediate/align/transcript_bam/{{sample}}.Aligned.toTranscriptome.out.bam",
        log=f"{OUTPUT_DIR}/logs/star/{{sample}}.Log.final.out",
    params:
        prefix=f"{OUTPUT_DIR}/intermediate/align/genome_bam/{{sample}}.",
        genomeDir=lambda w: config["ercc"]["ercc_star_index"] if ERCC else config["ref"]["star_index"]
    threads: config["threads"]
    container: config["container"]
    wildcard_constraints:
        sample=build_sample_constraint(PAIRED_SAMPLES)
    shell:
        """
        mkdir -p {OUTPUT_DIR}/intermediate/align/genome_bam {OUTPUT_DIR}/intermediate/align/transcript_bam {OUTPUT_DIR}/logs/star
        STAR --genomeDir {params.genomeDir} --runThreadN {threads} \
            --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM Unsorted \
            --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --quantMode TranscriptomeSAM GeneCounts
        mv {params.prefix}Aligned.toTranscriptome.out.bam {output.transcript_bam}
        mv {params.prefix}Log.final.out {output.log}
        """

# ===== 对齐规则（单端）=====
rule star_align_single:
    input:
        r1=f"{OUTPUT_DIR}/intermediate/rrna_removed/{{sample}}.unmapped.fastq.gz",
    output:
        genome_bam=f"{OUTPUT_DIR}/intermediate/align/genome_bam/{{sample}}.Aligned.out.bam",
        transcript_bam=f"{OUTPUT_DIR}/intermediate/align/transcript_bam/{{sample}}.Aligned.toTranscriptome.out.bam",
        log=f"{OUTPUT_DIR}/logs/star/{{sample}}.Log.final.out",
    params:
        prefix=f"{OUTPUT_DIR}/intermediate/align/genome_bam/{{sample}}.",
        genomeDir=lambda w: config["ercc"]["ercc_star_index"] if ERCC else config["ref"]["star_index"]
    threads: config["threads"]
    container: config["container"]
    wildcard_constraints:
        sample=build_sample_constraint(SINGLE_SAMPLES)
    shell:
        """
        mkdir -p {OUTPUT_DIR}/intermediate/align/genome_bam {OUTPUT_DIR}/intermediate/align/transcript_bam {OUTPUT_DIR}/logs/star
        STAR --genomeDir {params.genomeDir} --runThreadN {threads} \
            --readFilesIn {input.r1} --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM Unsorted \
            --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --quantMode TranscriptomeSAM GeneCounts
        mv {params.prefix}Aligned.toTranscriptome.out.bam {output.transcript_bam}
        mv {params.prefix}Log.final.out {output.log}
        """

# ===== BAM 处理（共享）=====
rule process_genome_bam:
    input:
        f"{OUTPUT_DIR}/intermediate/align/genome_bam/{{sample}}.Aligned.out.bam"
    output:
        bam = f"{OUTPUT_DIR}/results/{{sample}}.unique.dedup.bam",
        bai = f"{OUTPUT_DIR}/results/{{sample}}.unique.dedup.bam.bai",
        bigwig = f"{OUTPUT_DIR}/results/{{sample}}.unique.dedup.bw",
        log = f"{OUTPUT_DIR}/logs/bam_processing/{{sample}}.log"
    threads: 4
    container: config["container"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/results {OUTPUT_DIR}/logs/bam_processing
        SAMPLE_NAME={wildcards.sample}
        SAMPLE_OUT_DIR={OUTPUT_DIR}/results
        INPUT_BAM={input}
        THREADS={threads}

        echo "==> Processing sample: $SAMPLE_NAME" > {output.log}

        samtools view -q 255 -bS "$INPUT_BAM" | \
        samtools sort -n -@ $THREADS -o "$SAMPLE_OUT_DIR/${{SAMPLE_NAME}}.namesort.bam" - 2>> {output.log}

        samtools fixmate -m -r \
            "$SAMPLE_OUT_DIR/${{SAMPLE_NAME}}.namesort.bam" \
            "$SAMPLE_OUT_DIR/${{SAMPLE_NAME}}.fixmate.bam" 2>> {output.log}

        samtools sort -@ $THREADS \
            -o "$SAMPLE_OUT_DIR/${{SAMPLE_NAME}}.fixmate.sorted.bam" \
            "$SAMPLE_OUT_DIR/${{SAMPLE_NAME}}.fixmate.bam" 2>> {output.log}

        samtools markdup -r -s -@ $THREADS \
            "$SAMPLE_OUT_DIR/${{SAMPLE_NAME}}.fixmate.sorted.bam" \
            "{output.bam}" 2>> {output.log}

        samtools index "{output.bam}" 2>> {output.log}
        bamCoverage -b {output.bam} -o {output.bigwig} -p $THREADS --normalizeUsing CPM 2>> {output.log}
        rm -f "$SAMPLE_OUT_DIR/${{SAMPLE_NAME}}.namesort.bam" \
            "$SAMPLE_OUT_DIR/${{SAMPLE_NAME}}.fixmate.bam" \
            "$SAMPLE_OUT_DIR/${{SAMPLE_NAME}}.fixmate.sorted.bam"

        echo "==> BAM processing completed for $SAMPLE_NAME" >> {output.log}
        """

# ===== 定量规则 =====
rule featurecounts:
    input: 
        expand(f"{OUTPUT_DIR}/results/{{sample}}.unique.dedup.bam", sample=SAMPLES)
    output:
        counts=f"{OUTPUT_DIR}/intermediate/featureCounts.raw_counts.txt",
        summary=f"{OUTPUT_DIR}/logs/featurecounts_summary.txt",
    params:
        strandedness=config["params"]["featurecounts_strandedness"],
        is_paired="-p" if any_paired_end() else "",
        gtf=config["ref"]["gtf"],
    threads: config["threads"]
    container: config["container"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/intermediate {OUTPUT_DIR}/logs
        featureCounts -a {params.gtf} -o {output.counts} -s {params.strandedness} \
            {params.is_paired} -T {threads} --verbose {input} > {output.summary} 2>&1
        """

rule format_featurecounts:
    input:
        f"{OUTPUT_DIR}/intermediate/featureCounts.raw_counts.txt"
    output:
        f"{OUTPUT_DIR}/intermediate/featureCounts.raw_counts.tsv"
    container: config["container"]
    shell:
        """
        cp {input} {output}
        """

# ===== RSEM 定量 =====
if USE_RSEM:
    rule rsem_calculate_expression:
        input:
            bam=f"{OUTPUT_DIR}/intermediate/align/transcript_bam/{{sample}}.Aligned.toTranscriptome.out.bam",
            ref_done=RSEM_GRP,
        output:
            f"{OUTPUT_DIR}/intermediate/align/rsem_temp/{{sample}}.genes.results"
        params:
            ref_name=RSEM_REF_NAME,
            prefix=f"{OUTPUT_DIR}/intermediate/align/rsem_temp/{{sample}}",
            is_paired=lambda w: "--paired-end" if is_paired_end(w) else ""
        threads: config["threads"]
        container: config["container"]
        log:
            f"{OUTPUT_DIR}/logs/rsem/{{sample}}.log"
        benchmark:
            f"{OUTPUT_DIR}/logs/rsem/{{sample}}.benchmark.txt"
        shell:
            """
            mkdir -p {OUTPUT_DIR}/intermediate/align/rsem_temp
            rsem-calculate-expression --bam {params.is_paired} -p {threads} {input.bam} {params.ref_name} {params.prefix} > {log} 2>&1
            """

# ===== ERCC 分析（可选）=====
if ERCC:
    rule ercc_featurecounts:
        input:
            bam = f"{OUTPUT_DIR}/results/{{sample}}.unique.dedup.bam",
        output:
            counts=f"{OUTPUT_DIR}/intermediate/{{sample}}.featureCounts.ercc_counts.txt",
            summary=f"{OUTPUT_DIR}/logs/{{sample}}.featurecounts_ercc_summary.txt",
        params:
            strandedness=config["params"]["featurecounts_strandedness"],
            is_paired="-p" if any_paired_end() else "",
            ercc_gtf=config["ercc"]["ercc_gtf"],
        threads: config["threads"]
        container: config["container"]
        shell:
            """
            mkdir -p {OUTPUT_DIR}/intermediate {OUTPUT_DIR}/logs
            featureCounts -a {params.ercc_gtf} -o {output.counts} -s {params.strandedness} \
                {params.is_paired} -T {threads} --verbose {input} > {output.summary} 2>&1
            """

    rule ercc_analysis:
        input:
            ercc_counts = f"{OUTPUT_DIR}/intermediate/{{sample}}.featureCounts.ercc_counts.txt",
            gene_counts = f"{OUTPUT_DIR}/intermediate/featureCounts.raw_counts.txt",
            ercc_conc = config["ercc"]["ercc_conc"],
            bam = f"{OUTPUT_DIR}/results/{{sample}}.unique.dedup.bam"
        output:
            plot = f"{OUTPUT_DIR}/results/{{sample}}_ERCC_standard_curve.pdf",
            table = f"{OUTPUT_DIR}/results/{{sample}}_mpc.txt"
        params:
            script = config["ercc"]["ercc_analysis_script"],
            dilution = config["ercc_params"]["dilution"],
            ercc_ul = config["ercc_params"]["ercc_ul"],
            rna_pg_per_cell = config["ercc_params"]["rna_pg_per_cell"],
            total_rna_ug = config["ercc_params"]["total_rna_ug"],
            output_prefix = f"{OUTPUT_DIR}/results/{{sample}}"
        log:
            f"{OUTPUT_DIR}/logs/{{sample}}.ERCC_analysis.log"
        container: config["container"]
        shell:
            """
            mkdir -p {OUTPUT_DIR}/results {OUTPUT_DIR}/logs
            TOTAL_READS=$(samtools view -c -F 260 {input.bam})
            Rscript {params.script} \
                -d {params.dilution} \
                -e {params.ercc_ul} \
                -p {params.rna_pg_per_cell} \
                -r {params.total_rna_ug} \
                -o {params.output_prefix} \
                {input.ercc_counts} \
                {input.gene_counts} \
                {input.ercc_conc} \
                $TOTAL_READS \
                2> {log}
            """

# ===== 下游分析 =====
rule RNA_downstream_analysis:
    input:
        config="config.yaml",
        counts_file=f"{OUTPUT_DIR}/intermediate/featureCounts.raw_counts.tsv",
        rsem_files=expand(f"{OUTPUT_DIR}/intermediate/align/rsem_temp/{{sample}}.genes.results", sample=SAMPLES) if USE_RSEM else []
    output:
        results=f"{OUTPUT_DIR}/results/{CONTRAST[0]}_vs_{CONTRAST[1]}.deseq2_results.csv",
        plot=f"{OUTPUT_DIR}/results/{CONTRAST[0]}_vs_{CONTRAST[1]}.volcano_plot.png",
        go_csv=f"{OUTPUT_DIR}/results/{CONTRAST[0]}_vs_{CONTRAST[1]}.GO_enrichment.csv",
        go_plot=f"{OUTPUT_DIR}/results/{CONTRAST[0]}_vs_{CONTRAST[1]}.GO_dotplot.png",
        counts_matrix=f"{OUTPUT_DIR}/results/featureCounts.gene_counts.symbol.tsv",
    params:
        script=config["scripts"]["rna_analysis"],
        mode="align",
        contrast=f"{CONTRAST[0]},{CONTRAST[1]}",
        use_rsem=USE_RSEM,
        rsem_files_str=lambda w, input: ",".join(input.rsem_files) if input.rsem_files else "",
        tpm_matrix=f"{OUTPUT_DIR}/results/RSEM.gene_tpm.symbol.tsv" if USE_RSEM else None,
    log:
        f"{OUTPUT_DIR}/logs/{CONTRAST[0]}_vs_{CONTRAST[1]}_RNA.log"
    container: config["container"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/results {OUTPUT_DIR}/logs
        if [ "{params.use_rsem}" = "True" ]; then
            Rscript {params.script} \
                --mode {params.mode} \
                --config {input.config} \
                --input_counts {input.counts_file} \
                --input_rsem_files "{params.rsem_files_str}" \
                --contrast {params.contrast} \
                --output_prefix "{OUTPUT_DIR}/results/{CONTRAST[0]}_vs_{CONTRAST[1]}" \
                --output_counts_matrix "{output.counts_matrix}" \
                --output_tpm_matrix "{params.tpm_matrix}" \
                2> {log}
        else
            Rscript {params.script} \
                --mode {params.mode} \
                --config {input.config} \
                --input_counts {input.counts_file} \
                --contrast {params.contrast} \
                --output_prefix "{OUTPUT_DIR}/results/{CONTRAST[0]}_vs_{CONTRAST[1]}" \
                --output_counts_matrix "{output.counts_matrix}" \
                2> {log}
        fi
        """
if DESEQ2:
    rule lncRNA_downstream_analysis:
        input:
            deseq_results = f"{OUTPUT_DIR}/results/{CONTRAST[0]}_vs_{CONTRAST[1]}.deseq2_results.csv",
            gtf = config["ref"]["gtf"]
        output:
            summary_table = f"{OUTPUT_DIR}/lncRNA_analysis_output/tables/QC_summary.csv"
        params:
            script = config["scripts"]["lncRNA_analysis"],
            outdir = f"{OUTPUT_DIR}/lncRNA_analysis_output",
            padj_thr = config["analysis"]["fdr_threshold"],
            lfc_thr = config["analysis"]["lfc_threshold"]
        log:
            f"{OUTPUT_DIR}/logs/lncRNA_analysis.log"
        container: config["container"]
        shell:
            """
            mkdir -p {OUTPUT_DIR}/lncRNA_analysis_output/tables {OUTPUT_DIR}/logs
            Rscript {params.script} \
                --deg {input.deseq_results} \
                --gtf {input.gtf} \
                --outdir {params.outdir} \
                --padj {params.padj_thr} \
                --lfc {params.lfc_thr} \
                2> {log}
            """

rule multiqc:
    input:
        star_logs=expand(f"{OUTPUT_DIR}/logs/star/{{sample}}.Log.final.out", sample=SAMPLES),
        bam_logs=expand(f"{OUTPUT_DIR}/logs/bam_processing/{{sample}}.log", sample=SAMPLES),
        rna_analysis_log=f"{OUTPUT_DIR}/logs/{CONTRAST[0]}_vs_{CONTRAST[1]}_RNA.log",
    output:
        f"{OUTPUT_DIR}/results/multiqc_report.html",
    container: config["container"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/results
        multiqc {OUTPUT_DIR} -f -o {OUTPUT_DIR}/results
        """

# ===== 报告生成 =====
if REPORT_ENABLE:
    def _report_inputs(wildcards):
        """Collect all upstream outputs that should exist before report generation."""
        inputs = {
            "multiqc": f"{OUTPUT_DIR}/results/multiqc_report.html",
            "config": "config.yaml",
        }
        if DESEQ2:
            inputs["deseq2_results"] = (
                f"{OUTPUT_DIR}/results/{CONTRAST[0]}_vs_{CONTRAST[1]}.deseq2_results.csv"
            )
            inputs["qc_summary"] = (
                f"{OUTPUT_DIR}/lncRNA_analysis_output/tables/QC_summary.csv"
            )
        if ERCC:
            inputs["ercc_mpc"] = expand(
                f"{OUTPUT_DIR}/results/{{sample}}_mpc.txt", sample=SAMPLES
            )
        if USE_RSEM:
            inputs["rsem_tpm"] = f"{OUTPUT_DIR}/results/RSEM.gene_tpm.symbol.tsv"
        return inputs

    # Absolute path to this script so the rule works regardless of working directory.
    _REPORT_SCRIPT = str(Path(workflow.snakefile).parent / "report" / "generate_report.py")

    rule generate_report:
        input:
            unpack(_report_inputs)
        output:
            html = f"{OUTPUT_DIR}/report/report.html",
            md   = f"{OUTPUT_DIR}/report/report.md",
        params:
            script     = _REPORT_SCRIPT,
            output_dir = OUTPUT_DIR,
            report_dir = f"{OUTPUT_DIR}/report",
        log:
            f"{OUTPUT_DIR}/logs/generate_report.log"
        shell:
            """
            mkdir -p {params.report_dir}
            python {params.script} \
                --config {input.config} \
                --outdir {params.output_dir} \
                --report {params.report_dir} \
                > {log} 2>&1
            """
