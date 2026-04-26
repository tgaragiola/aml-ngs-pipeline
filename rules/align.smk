# =============================================================================
# Módulo 2: Alineamiento
# Descripción: Alinea las lecturas al genoma de referencia hg38, filtra las
#              regiones del panel de interés y marca los duplicados de PCR.
# Herramientas: BWA-MEM, samtools, GATK MarkDuplicates
# Flujo: BWA-MEM → subset por panel → MarkDuplicates
# =============================================================================

rule bwa_mem:
    """
    Alinea las lecturas paired-end al genoma de referencia humano hg38
    usando el algoritmo BWA-MEM. El resultado es un archivo BAM ordenado
    por coordenadas genómicas.
    El read group (@RG) es necesario para que GATK identifique correctamente
    la muestra en pasos posteriores.
    """
    input:
        r1  = "results/{sample}/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2  = "results/{sample}/trimmed/{sample}_R2.trimmed.fastq.gz",
        ref = config["reference"]["genome"],
    output:
        bam = temp("results/{sample}/bam/{sample}.raw.bam"),  # temporal: se elimina al finalizar
    log:
        "logs/{sample}/bwa_mem.log"
    conda:
        "../envs/align.yaml"
    threads:
        config["resources"]["bwa"]["threads"]
    resources:
        mem_mb = config["resources"]["bwa"]["mem_mb"]
    params:
        rg = lambda wildcards: (
            f"@RG\\tID:{wildcards.sample}\\t"
            f"SM:{wildcards.sample}\\t"
            f"PL:ILLUMINA\\t"
            f"LB:{wildcards.sample}_lib1\\t"
            f"PU:{wildcards.sample}"
        )
    shell:
        """
        # Elimina archivos temporales de corridas anteriores que puedan causar conflictos
        rm -f results/{wildcards.sample}/bam/*.tmp.*.bam 2>/dev/null || true

        # Alinea con BWA-MEM y ordena el BAM resultante con samtools
        bwa mem \
            -t {threads} \
            -R '{params.rg}' \
            {input.ref} \
            {input.r1} {input.r2} \
        | samtools sort \
            -@ {threads} \
            -o {output.bam} \
            -m 2G \
        2> {log}
        """

rule subset_panel:
    """
    Filtra el BAM para conservar únicamente las lecturas que mapean sobre
    las regiones definidas en el archivo BED del panel de interés.
    Esto reduce el tamaño del archivo en ~97%, acelerando los pasos siguientes
    sin perder información relevante para el análisis.
    Paso incorporado por recomendación de la tutora del trabajo.
    """
    input:
        bam = "results/{sample}/bam/{sample}.raw.bam",
        bed = config["panel"]["bed"],
    output:
        bam = temp("results/{sample}/bam/{sample}.panel.bam"),
        bai = temp("results/{sample}/bam/{sample}.panel.bam.bai"),
    log:
        "logs/{sample}/subset_panel.log"
    conda:
        "../envs/align.yaml"
    threads: 2
    shell:
        """
        # Filtra el BAM por las regiones del BED y genera el índice
        samtools view -b -L {input.bed} {input.bam} \
            -o {output.bam} \
            2> {log}
        samtools index {output.bam}
        """

rule mark_duplicates:
    """
    Identifica y marca las lecturas duplicadas generadas durante la
    amplificación por PCR de la librería. Los duplicados no se eliminan
    sino que se marcan para que GATK los ignore en el llamado de variantes.
    """
    input:
        bam = "results/{sample}/bam/{sample}.panel.bam",
    output:
        bam     = temp("results/{sample}/bam/{sample}.markdup.bam"),
        bai     = temp("results/{sample}/bam/{sample}.markdup.bai"),
        metrics = "results/{sample}/qc/{sample}.markdup_metrics.txt",  # estadísticas de duplicados
    log:
        "logs/{sample}/markdup.log"
    conda:
        "../envs/gatk.yaml"
    threads:
        config["resources"]["gatk"]["threads"]
    resources:
        mem_mb = config["resources"]["gatk"]["mem_mb"]
    shell:
        """
        gatk MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            2> {log}
        """
