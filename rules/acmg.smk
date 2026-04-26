# =============================================================================
# Módulo 6: Clasificación ACMG/AMP
# Descripción: Clasifica cada variante anotada según los criterios estándar
#              ACMG/AMP 2015, adaptados al contexto clínico de la LMA.
# Herramientas: Script Python desarrollado en este trabajo (acmg_classifier.py)
# Flujo: VCF anotado → clasificación ACMG → tabla resumen
# =============================================================================

rule acmg_classify:
    """
    Clasifica cada variante según los criterios ACMG/AMP 2015 adaptados a LMA.
    El script evalúa los siguientes criterios para cada variante:
    - PVS1: variante de pérdida de función en gen supresor tumoral
    - PM2: frecuencia poblacional muy baja en gnomAD
    - PM4: inserción/deleción in-frame en región no repetitiva
    - PP3: predicciones in silico deletéreas (CADD, SpliceAI)
    - BA1: frecuencia alta en gnomAD (>5%) → benigna directa
    - BP4: predicciones in silico benignas
    - AML_HOTSPOT: posición recurrentemente mutada en LMA
    La clasificación final puede ser: Pathogenic, Likely pathogenic,
    Uncertain significance (VUS), Likely benign o Benign.
    Los umbrales de frecuencia y genes de interés se definen en config.yaml.
    """
    input:
        vcf = "results/{sample}/annotation/{sample}.annotated.vcf.gz",
        tbi = "results/{sample}/annotation/{sample}.annotated.vcf.gz.tbi",
    output:
        vcf = "results/{sample}/acmg/{sample}.acmg_classified.vcf.gz",
        tbi = "results/{sample}/acmg/{sample}.acmg_classified.vcf.gz.tbi",
        tsv = "results/{sample}/acmg/{sample}.acmg_classified.tsv",  # tabla de variantes clasificadas
    log:
        "logs/{sample}/acmg_classify.log"
    conda:
        "../envs/python.yaml"
    params:
        script     = "scripts/acmg_classifier.py",
        aml_genes  = ",".join(
            config["acmg"]["aml_genes"]["tier1"] +
            config["acmg"]["aml_genes"]["tier2"]
        ),
        freq_benign       = config["acmg"]["freq_thresholds"]["benign"],        # por defecto 0.05
        freq_likelybenign = config["acmg"]["freq_thresholds"]["likely_benign"], # por defecto 0.01
        freq_pathogenic   = config["acmg"]["freq_thresholds"]["pathogenic"],    # por defecto 0.001
        clinvar_stars_min = config["acmg"]["clinvar_stars_min"],                # estrellas mínimas ClinVar
    shell:
        """
        python {params.script} \
            --input {input.vcf} \
            --output-vcf {output.vcf} \
            --output-tsv {output.tsv} \
            --aml-genes "{params.aml_genes}" \
            --freq-benign {params.freq_benign} \
            --freq-likely-benign {params.freq_likelybenign} \
            --freq-pathogenic {params.freq_pathogenic} \
            --clinvar-stars-min {params.clinvar_stars_min} \
            2> {log}
        tabix -p vcf {output.vcf}
        """

rule merge_acmg_summary:
    """
    Combina las tablas de variantes clasificadas de todas las muestras
    en un único archivo TSV resumen para facilitar la comparación entre muestras
    y el análisis global de los resultados.
    """
    input:
        tsvs = expand("results/{sample}/acmg/{sample}.acmg_classified.tsv", sample=SAMPLES),
    output:
        tsv = "results/summary/all_variants_acmg.tsv",
    log:
        "logs/merge_acmg_summary.log"
    conda:
        "../envs/python.yaml"
    params:
        script = "scripts/merge_summary.py",
    shell:
        """
        python {params.script} \
            --inputs {input.tsvs} \
            --output {output.tsv} \
            2> {log}
        """
