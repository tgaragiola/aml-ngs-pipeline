# =============================================================================
# Módulo 5: Anotación funcional
# Descripción: Anota las variantes con información funcional, frecuencias
#              poblacionales y clasificaciones clínicas previas.
# Herramientas: Ensembl VEP, bcftools annotate
# Flujo: VEP (anotación funcional) → bcftools (anotación ClinVar)
# =============================================================================

rule vep_annotate:
    """
    Anota cada variante con información funcional usando Ensembl VEP v111.
    Para cada variante determina:
    - El efecto sobre la proteína (missense, frameshift, sinónima, etc.)
    - El transcripto canónico afectado
    - La frecuencia en la población general (gnomAD)
    - Predicciones in silico de patogenicidad (CADD, SpliceAI)
    - Nomenclatura HGVS (HGVSc y HGVSp)
    Usa una caché local para no depender de conexión a internet durante el análisis.
    """
    input:
        vcf       = "results/{sample}/vcf/{sample}.pass.vcf.gz",
        tbi       = "results/{sample}/vcf/{sample}.pass.vcf.gz.tbi",
        ref       = config["reference"]["genome"],
        cache_dir = config["vep"]["cache_dir"],
    output:
        vcf   = "results/{sample}/annotation/{sample}.vep.vcf.gz",
        tbi   = "results/{sample}/annotation/{sample}.vep.vcf.gz.tbi",
        stats = "results/{sample}/annotation/{sample}.vep_stats.html",  # reporte de estadísticas
    log:
        "logs/{sample}/vep.log"
    conda:
        "../envs/vep.yaml"
    threads:
        config["resources"]["vep"]["threads"]
    resources:
        mem_mb = config["resources"]["vep"]["mem_mb"]
    params:
        cache_version = config["vep"]["cache_version"],
        species       = config["vep"]["species"],
        assembly      = config["vep"]["assembly"],
        extra         = config["vep"]["extra_flags"],
        plugins       = " ".join([f"--plugin {p}" for p in config["vep"]["plugins"]]),
    shell:
        """
        vep \
            --input_file {input.vcf} \
            --output_file {output.vcf} \
            --format vcf \
            --vcf \
            --compress_output bgzip \
            --stats_file {output.stats} \
            --cache \
            --cache_version {params.cache_version} \
            --dir_cache {input.cache_dir} \
            --species {params.species} \
            --assembly {params.assembly} \
            --fasta {input.ref} \
            --fork {threads} \
            {params.plugins} \
            {params.extra} \
            2> {log}
        tabix -p vcf {output.vcf}
        """

rule annotate_clinvar:
    """
    Agrega al VCF la información de ClinVar para cada variante.
    ClinVar es una base de datos pública que contiene clasificaciones
    clínicas de variantes reportadas por laboratorios de todo el mundo.
    Los campos agregados incluyen:
    - CLNSIG: clasificación clínica (Pathogenic, Benign, VUS, etc.)
    - CLNDN: nombre de la enfermedad asociada
    - CLNREVSTAT: nivel de revisión (número de estrellas)
    - CLNHGVS: nomenclatura HGVS de ClinVar
    """
    input:
        vcf     = "results/{sample}/annotation/{sample}.vep.vcf.gz",
        tbi     = "results/{sample}/annotation/{sample}.vep.vcf.gz.tbi",
        clinvar = "resources/clinvar/clinvar.vcf.gz",
    output:
        vcf = "results/{sample}/annotation/{sample}.annotated.vcf.gz",
        tbi = "results/{sample}/annotation/{sample}.annotated.vcf.gz.tbi",
    log:
        "logs/{sample}/clinvar_annotate.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools annotate \
            --annotations {input.clinvar} \
            --columns INFO/CLNDN,INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNHGVS,INFO/CLNVC \
            --output {output.vcf} \
            --output-type z \
            {input.vcf} \
            2> {log}
        tabix -p vcf {output.vcf}
        """

