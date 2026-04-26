# =============================================================================
# Módulo 3: Calibración de scores de calidad de base (BQSR)
# Descripción: Corrige los sesgos sistemáticos en los scores de calidad de base
#              asignados por el secuenciador, mejorando la precisión del llamado
#              de variantes posterior.
# Herramientas: GATK BaseRecalibrator, GATK ApplyBQSR
# Flujo: BaseRecalibrator (genera tabla) → ApplyBQSR (aplica corrección)
# Recursos externos necesarios: dbSNP, Mills & 1000G (variantes conocidas)
# =============================================================================

rule base_recalibrator:
    """
    Analiza los patrones de error del secuenciador comparando las bases
    observadas con variantes conocidas (dbSNP, Mills, 1000G).
    Genera una tabla de recalibración que será aplicada en el paso siguiente.
    Las variantes conocidas se usan como referencia para distinguir errores
    reales del secuenciador de variantes verdaderas.
    """
    input:
        bam    = "results/{sample}/bam/{sample}.markdup.bam",
        bai    = "results/{sample}/bam/{sample}.markdup.bai",
        ref    = config["reference"]["genome"],
        bed    = config["panel"]["bed"],
        dbsnp  = config["known_sites"]["dbsnp"],    # base de datos de polimorfismos conocidos
        mills  = config["known_sites"]["mills"],    # indels conocidos de Mills & 1000G
        g1000  = config["known_sites"]["g1000"],    # SNPs de alta confianza del 1000 Genomes
    output:
        table = "results/{sample}/bqsr/{sample}.recal.table",  # tabla de recalibración
    log:
        "logs/{sample}/base_recalibrator.log"
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = config["resources"]["gatk"]["mem_mb"]
    params:
        padding = config["panel"]["padding"],  # margen extra alrededor de las regiones del panel
    shell:
        """
        gatk BaseRecalibrator \
            --input {input.bam} \
            --reference {input.ref} \
            --intervals {input.bed} \
            --interval-padding {params.padding} \
            --known-sites {input.dbsnp} \
            --known-sites {input.mills} \
            --known-sites {input.g1000} \
            --output {output.table} \
            2> {log}
        """

rule apply_bqsr:
    """
    Aplica la tabla de recalibración generada por BaseRecalibrator al BAM,
    produciendo un nuevo BAM con los scores de calidad de base corregidos.
    Este BAM calibrado es el que se usa para el llamado de variantes con Mutect2.
    """
    input:
        bam   = "results/{sample}/bam/{sample}.markdup.bam",
        bai   = "results/{sample}/bam/{sample}.markdup.bai",
        ref   = config["reference"]["genome"],
        table = "results/{sample}/bqsr/{sample}.recal.table",
        bed   = config["panel"]["bed"],
    output:
        bam = "results/{sample}/bam/{sample}.recal.bam",   # BAM calibrado final
        bai = "results/{sample}/bam/{sample}.recal.bai",   # índice del BAM calibrado
    log:
        "logs/{sample}/apply_bqsr.log"
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = config["resources"]["gatk"]["mem_mb"]
    params:
        padding = config["panel"]["padding"],
    shell:
        """
        gatk ApplyBQSR \
            --input {input.bam} \
            --reference {input.ref} \
            --intervals {input.bed} \
            --interval-padding {params.padding} \
            --bqsr-recal-file {input.table} \
            --output {output.bam} \
            --create-output-bam-index true \
            2> {log}
        """
