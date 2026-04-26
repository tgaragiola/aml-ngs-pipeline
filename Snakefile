# =============================================================================
# AML NGS Pipeline - Snakemake
# Descripción: Pipeline bioinformático para la identificación, anotación y 
#              clasificación de variantes somáticas en leucemia mieloide aguda (LMA)
# Modo: Tumor-only (sin muestra normal pareada)
# Panel: Dirigido a regiones de interés que puede definir el usuario mediante un archivo BED
#        (por defecto: 23 genes de relevancia clínica en LMA según ELN 2022)
# Referencia: hg38 (GRCh38)
# Autor: María Trinidad Garagiola
# =============================================================================
import pandas as pd
from pathlib import Path
# ── Configuración ─────────────────────────────────────────────────────────────
# En este archivo están definidas las rutas y parámetros del flujo de trabajo
configfile: "config/config.yaml"
# Cargar tabla de muestras
samples_df = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
SAMPLES = samples_df.index.tolist()
# ── Funciones helper ──────────────────────────────────────────────────────────
# Estas funciones obtienen las rutas a los archivos FASTQ de cada muestra
# a partir de la tabla de muestras, usando el nombre de la muestra como índice
def get_fastq_r1(wildcards):
    return samples_df.loc[wildcards.sample, "fastq_r1"]
def get_fastq_r2(wildcards):
    return samples_df.loc[wildcards.sample, "fastq_r2"]
def get_sample_id(wildcards):
    return wildcards.sample
# ── Output final (regla all) ──────────────────────────────────────────────────
# Esta regla define los archivos finales que debe generar el pipeline.
# Snakemake trabaja hacia atrás desde estos archivos para determinar
# qué pasos necesita ejecutar.
rule all:
    input:
        "results/qc/multiqc_report.html",
        expand("results/{sample}/acmg/{sample}.acmg_classified.vcf.gz", sample=SAMPLES),
        expand("results/{sample}/report/{sample}_report.html", sample=SAMPLES),
        "results/summary/all_variants_acmg.tsv"
# ── Incluir módulos ───────────────────────────────────────────────────────────
include: "rules/qc.smk"              # Control de calidad (fastp, FastQC, MultiQC, mosdepth)
include: "rules/align.smk"           # Alineación al genoma (BWA, subset por panel)
include: "rules/bqsr.smk"            # Calibración de calidad de base (GATK BQSR)
include: "rules/variant_calling.smk" # Llamado de variantes (GATK Mutect2)
include: "rules/annotation.smk"      # Anotación funcional (VEP, ClinVar)
include: "rules/acmg.smk"            # Clasificación ACMG
include: "rules/report.smk"          # Generación de reportes HTML y tabla resumen
