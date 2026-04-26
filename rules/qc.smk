# =============================================================================
# Módulo 1: Control de calidad y trimming
# Descripción: Evalúa y mejora la calidad de las lecturas crudas antes
#              de la alineación. Calcula también la cobertura sobre el panel.
# Herramientas: fastp, FastQC, mosdepth, MultiQC
# Flujo: fastp → FastQC → mosdepth → MultiQC
# =============================================================================

rule fastp:
    """
    Elimina adaptadores de secuenciación, descarta bases de baja calidad
    y filtra lecturas demasiado cortas.
    Parámetros aplicados:
    - Longitud mínima de lectura: definida en config.yaml (por defecto 50 pb)
    - Calidad mínima Phred: definida en config.yaml (por defecto Q20)
    - Detección automática de adaptadores para datos paired-end
    """
    input:
        r1 = get_fastq_r1,
        r2 = get_fastq_r2,
    output:
        r1     = "results/{sample}/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2     = "results/{sample}/trimmed/{sample}_R2.trimmed.fastq.gz",
        html   = "results/{sample}/qc/fastp/{sample}_fastp.html",   # reporte visual
        json   = "results/{sample}/qc/fastp/{sample}_fastp.json",   # métricas en formato JSON
    log:
        "logs/{sample}/fastp.log"
    conda:
        "../envs/qc.yaml"
    threads:
        config["resources"]["default"]["threads"]
    resources:
        mem_mb = config["resources"]["default"]["mem_mb"]
    params:
        min_len      = config["qc"]["fastp"]["min_length"],
        quality      = config["qc"]["fastp"]["quality"],
        adapter_flag = "--detect_adapter_for_pe" if config["qc"]["fastp"]["detect_adapter"] else "",
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --html {output.html} --json {output.json} \
            --length_required {params.min_len} \
            --qualified_quality_phred {params.quality} \
            {params.adapter_flag} \
            --thread {threads} \
            2> {log}
        """

rule fastqc_trimmed:
    """
    Genera reportes de calidad detallados sobre las lecturas ya filtradas
    por fastp. Permite verificar visualmente que el trimming fue correcto.
    """
    input:
        r1 = "results/{sample}/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2 = "results/{sample}/trimmed/{sample}_R2.trimmed.fastq.gz",
    output:
        html_r1 = "results/{sample}/qc/fastqc/{sample}_R1.trimmed_fastqc.html",
        zip_r1  = "results/{sample}/qc/fastqc/{sample}_R1.trimmed_fastqc.zip",
        html_r2 = "results/{sample}/qc/fastqc/{sample}_R2.trimmed_fastqc.html",
        zip_r2  = "results/{sample}/qc/fastqc/{sample}_R2.trimmed_fastqc.zip",
    log:
        "logs/{sample}/fastqc.log"
    conda:
        "../envs/qc.yaml"
    threads: 2
    params:
        outdir = "results/{sample}/qc/fastqc/"
    shell:
        """
        fastqc {input.r1} {input.r2} \
            --outdir {params.outdir} \
            --threads {threads} \
            2> {log}
        """

rule mosdepth_coverage:
    """
    Calcula la cobertura de secuenciación sobre las regiones del panel.
    Genera un resumen de cobertura media y el porcentaje de bases cubiertas
    a diferentes umbrales (20x, 100x, 500x).
    """
    input:
        bam  = "results/{sample}/bam/{sample}.recal.bam",
        bai  = "results/{sample}/bam/{sample}.recal.bai",
        bed  = config["panel"]["bed"],
    output:
        summary    = "results/{sample}/qc/coverage/{sample}.mosdepth.summary.txt",
        regions    = "results/{sample}/qc/coverage/{sample}.regions.bed.gz",
        thresholds = "results/{sample}/qc/coverage/{sample}.thresholds.bed.gz",
    log:
        "logs/{sample}/mosdepth.log"
    conda:
        "../envs/qc.yaml"
    threads:
        config["resources"]["default"]["threads"]
    params:
        prefix = "results/{sample}/qc/coverage/{sample}",
        thres  = "20,100,500",
    shell:
        """
        mkdir -p results/{wildcards.sample}/qc/coverage/
        mosdepth \
            --by {input.bed} \
            --thresholds {params.thres} \
            --threads {threads} \
            {params.prefix} {input.bam} \
            2> {log}
        """

rule multiqc:
    """
    Integra todos los reportes de QC de las tres muestras en un único
    informe HTML
    """
    input:
        fastp    = expand("results/{sample}/qc/fastp/{sample}_fastp.json", sample=SAMPLES),
        fastqc   = expand("results/{sample}/qc/fastqc/{sample}_R1.trimmed_fastqc.zip", sample=SAMPLES),
        coverage = expand("results/{sample}/qc/coverage/{sample}.mosdepth.summary.txt", sample=SAMPLES),
    output:
        "results/qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    conda:
        "../envs/qc.yaml"
    params:
        outdir = "results/qc/",
        dirs   = "results/",
    shell:
        """
        multiqc {params.dirs} \
            --outdir {params.outdir} \
            --filename multiqc_report.html \
            --force \
            2> {log}
        """
