# =============================================================================
# Módulo 4: Llamado de variantes somáticas - Modo tumor-only
# Descripción: Detecta variantes somáticas (SNVs e indels) en las muestras
#              tumorales sin muestra normal pareada, usando GATK Mutect2.
#              Incluye filtrado de artefactos y estimación de contaminación.
# Herramientas: GATK Mutect2, FilterMutectCalls, SelectVariants, bcftools
# Flujo: Mutect2 → LearnReadOrientation → GetPileupSummaries →
#        CalculateContamination → FilterMutectCalls → SelectVariants
# =============================================================================

rule mutect2_tumor_only:
    """
    Detecta variantes somáticas usando GATK Mutect2 en modo tumor-only.
    Utiliza un Panel of Normals (PoN) para filtrar artefactos técnicos
    y un recurso germinal (gnomAD) para distinguir variantes somáticas
    de polimorfismos germinales frecuentes en la población.
    """
    input:
        bam  = "results/{sample}/bam/{sample}.recal.bam",
        bai  = "results/{sample}/bam/{sample}.recal.bai",
        ref  = config["reference"]["genome"],
        bed  = config["panel"]["bed"],
        pon  = config["mutect2"]["pon"],                # Panel of Normals
        germ = config["mutect2"]["germline_resource"],  # gnomAD para filtrar germinales
    output:
        vcf   = temp("results/{sample}/vcf/{sample}.mutect2.vcf.gz"),
        tbi   = temp("results/{sample}/vcf/{sample}.mutect2.vcf.gz.tbi"),
        stats = "results/{sample}/vcf/{sample}.mutect2.vcf.gz.stats",
        f1r2  = "results/{sample}/vcf/{sample}.f1r2.tar.gz",  # datos para modelo de orientación
    log:
        "logs/{sample}/mutect2.log"
    conda:
        "../envs/gatk.yaml"
    threads:
        config["resources"]["gatk"]["threads"]
    resources:
        mem_mb = config["resources"]["gatk"]["mem_mb"]
    params:
        padding            = config["panel"]["padding"],
        af_not_in_resource = config["mutect2"]["af_of_alleles_not_in_resource"],
    shell:
        """
        gatk Mutect2 \
            --input {input.bam} \
            --reference {input.ref} \
            --intervals {input.bed} \
            --interval-padding {params.padding} \
            --panel-of-normals {input.pon} \
            --germline-resource {input.germ} \
            --af-of-alleles-not-in-resource {params.af_not_in_resource} \
            --f1r2-tar-gz {output.f1r2} \
            --output {output.vcf} \
            --native-pair-hmm-threads {threads} \
            --dont-use-soft-clipped-bases true \
            2> {log}
        """

rule learn_read_orientation:
    """
    Aprende el modelo de artefactos de orientación de lectura (OxoG/FFPE).
    Estos artefactos son errores técnicos que pueden confundirse con variantes
    somáticas reales. El modelo generado se usa en el filtrado posterior.
    """
    input:
        f1r2 = "results/{sample}/vcf/{sample}.f1r2.tar.gz",
    output:
        model = "results/{sample}/vcf/{sample}.read_orientation_model.tar.gz",
    log:
        "logs/{sample}/learn_read_orientation.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk LearnReadOrientationModel \
            --input {input.f1r2} \
            --output {output.model} \
            2> {log}
        """

rule get_pileup_summaries:
    """
    Calcula un resumen de las bases observadas en cada posición del genoma
    para estimar la fracción de contaminación cruzada entre muestras.
    """
    input:
        bam  = "results/{sample}/bam/{sample}.recal.bam",
        bai  = "results/{sample}/bam/{sample}.recal.bai",
        germ = config["mutect2"]["germline_resource"],
        bed  = config["panel"]["bed"],
    output:
        pileup = "results/{sample}/vcf/{sample}.pileup_summary.table",
    log:
        "logs/{sample}/pileup_summaries.log"
    conda:
        "../envs/gatk.yaml"
    params:
        padding = config["panel"]["padding"],
    shell:
        """
        gatk GetPileupSummaries \
            --input {input.bam} \
            --variant {input.germ} \
            --intervals {input.bed} \
            --interval-padding {params.padding} \
            --output {output.pileup} \
            2> {log}
        """

rule calculate_contamination:
    """
    Estima la fracción de contaminación cruzada usando el resumen de pileup.
    Esta información se usa en FilterMutectCalls para descartar variantes
    que podrían ser producto de contaminación y no mutaciones reales.
    """
    input:
        pileup = "results/{sample}/vcf/{sample}.pileup_summary.table",
    output:
        contam   = "results/{sample}/vcf/{sample}.contamination.table",
        segments = "results/{sample}/vcf/{sample}.tumor_segments.table",
    log:
        "logs/{sample}/calculate_contamination.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk CalculateContamination \
            --input {input.pileup} \
            --tumor-segmentation {output.segments} \
            --output {output.contam} \
            2> {log}
        """

rule filter_mutect_calls:
    """
    Aplica filtros estadísticos sobre las variantes llamadas por Mutect2
    para eliminar artefactos técnicos, variantes germinales y contaminación.
    Las variantes que superan todos los filtros reciben la etiqueta PASS.
    """
    input:
        vcf      = "results/{sample}/vcf/{sample}.mutect2.vcf.gz",
        tbi      = "results/{sample}/vcf/{sample}.mutect2.vcf.gz.tbi",
        stats    = "results/{sample}/vcf/{sample}.mutect2.vcf.gz.stats",
        ref      = config["reference"]["genome"],
        model    = "results/{sample}/vcf/{sample}.read_orientation_model.tar.gz",
        contam   = "results/{sample}/vcf/{sample}.contamination.table",
        segments = "results/{sample}/vcf/{sample}.tumor_segments.table",
    output:
        vcf   = temp("results/{sample}/vcf/{sample}.filtered.vcf.gz"),
        tbi   = temp("results/{sample}/vcf/{sample}.filtered.vcf.gz.tbi"),
        stats = "results/{sample}/vcf/{sample}.filtering.stats",
    log:
        "logs/{sample}/filter_mutect_calls.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk FilterMutectCalls \
            --variant {input.vcf} \
            --reference {input.ref} \
            --orientation-bias-artifact-priors {input.model} \
            --contamination-table {input.contam} \
            --tumor-segmentation {input.segments} \
            --filtering-stats {output.stats} \
            --output {output.vcf} \
            2> {log}
        """

rule select_pass_variants:
    """
    Selecciona únicamente las variantes con etiqueta PASS y aplica filtros
    adicionales de frecuencia alélica mínima (VAF) y número mínimo de lecturas
    de soporte para el alelo alternativo, definidos en config.yaml.
    """
    input:
        vcf = "results/{sample}/vcf/{sample}.filtered.vcf.gz",
        tbi = "results/{sample}/vcf/{sample}.filtered.vcf.gz.tbi",
        ref = config["reference"]["genome"],
    output:
        vcf = "results/{sample}/vcf/{sample}.pass.vcf.gz",
        tbi = "results/{sample}/vcf/{sample}.pass.vcf.gz.tbi",
    log:
        "logs/{sample}/select_pass.log"
    conda:
        "../envs/gatk.yaml"
    params:
        min_af = config["mutect2"]["min_af"],   # VAF mínimo (por defecto 1%)
        min_dp = config["mutect2"]["max_alt_dp"], # lecturas alternativas mínimas (por defecto 5)
    shell:
        """
        # Selecciona variantes PASS
        gatk SelectVariants \
            --variant {input.vcf} \
            --reference {input.ref} \
            --exclude-filtered true \
            --output {output.vcf} \
            2> {log}

        # Aplica filtro adicional por VAF mínimo y soporte mínimo de lecturas
        bcftools view \
            -i "FORMAT/AF[0:0] >= {params.min_af} && FORMAT/AD[0:1] >= {params.min_dp}" \
            {output.vcf} \
        | bgzip -c > {output.vcf}.tmp && mv {output.vcf}.tmp {output.vcf}

        tabix -f -p vcf {output.vcf}
        """
