# =============================================================================
# Módulo 7: Generación de reportes
# Descripción: Genera un reporte HTML por muestra con las variantes detectadas
#              y su clasificación ACMG, en formato clínicamente interpretable.
# Herramientas: Script Python desarrollado en este trabajo (generate_report.py)
# =============================================================================

rule generate_html_report:
    """
    Genera un reporte HTML simple y legible por muestra, que incluye:
    - Nombre de la muestra y fecha de análisis
    - Total de variantes detectadas
    - Tabla con todas las variantes ordenadas por clasificación ACMG
      (patogénicas primero, benignas al final)
    - Para cada variante: gen, consecuencia, VAF, profundidad,
      clasificación ACMG, criterios activados y Tier ELN 2022
    El reporte está diseñado para ser interpretado por un profesional
    clínico sin necesidad de conocimientos bioinformáticos.
    """
    input:
        acmg_tsv = "results/{sample}/acmg/{sample}.acmg_classified.tsv",
    output:
        html = "results/{sample}/report/{sample}_report.html",
    log:
        "logs/{sample}/report.log"
    conda:
        "../envs/python.yaml"
    params:
        sample_id = lambda wildcards: wildcards.sample,
    shell:
        """
        python scripts/generate_report.py \
            --sample {params.sample_id} \
            --acmg-tsv {input.acmg_tsv} \
            --output {output.html} \
            2> {log}
        """
