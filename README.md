# AML NGS Pipeline — Guía de instalación y uso
# Clasificación de variantes ACMG | Panel dirigido | hg38
Este flujo de trabajo fue diseñado para la identificación, anotación y clasificación de variantes somáticas en leucemia mieloide aguda (LMA), implementado con Snakemake. 
Desarrollado como parte de un Trabajo de Fin de Máster en Bioinformática.

# Requisitos previos

   - Python ≥ 3.10
   - Conda instalado
   - Snakemake ≥ 8.0
```bash
conda install -c conda-forge snakemake
```

# 1. Clonar el repositorio
```bash
git clone https://github.com/tgaragiola/aml-ngs-pipeline
cd aml-ngs-pipeline
```
# 2. Obtener datos de secuenciación

El pipeline acepta cualquier dato de secuenciación paired-end en formato FASTQ. Los datos deben corresponder a muestras de pacientes con LMA secuenciadas con panel dirigido o exoma completo (WES), con lecturas paired-end de 100-150 pb.
# Datos utilizados en este trabajo

Las muestras utilizadas para el desarrollo de este pipeline provienen del estudio CBF-AML del MD Anderson Cancer Center (BioProject PRJNA1445381), disponibles en SRA:

```bash
# Instalar sra-tools
conda install -c bioconda sra-tools

# Descargar las muestras utilizadas en este trabajo
fastq-dump --split-files --gzip SRR37858073 -O data/
fastq-dump --split-files --gzip SRR37858074 -O data/
fastq-dump --split-files --gzip SRR37858079 -O data/

# Renombrar para que coincidan con el formato del pipeline
mv data/SRR37858073_1.fastq.gz data/CBF_AML_010_R1.fastq.gz
mv data/SRR37858073_2.fastq.gz data/CBF_AML_010_R2.fastq.gz
mv data/SRR37858074_1.fastq.gz data/CBF_AML_009_R1.fastq.gz
mv data/SRR37858074_2.fastq.gz data/CBF_AML_009_R2.fastq.gz
mv data/SRR37858079_1.fastq.gz data/CBF_AML_008_R1.fastq.gz
mv data/SRR37858079_2.fastq.gz data/CBF_AML_008_R2.fastq.gz

```

# Para usar otros datos o datos propios:

Se puede usar cualquier dato de secuenciación paired-end. Solo se necesita actualizar config/samples.tsv con los nombres y rutas de las muestras:

sample        fastq_r1                    fastq_r2
MI_MUESTRA    data/MI_MUESTRA_R1.fastq.gz data/MI_MUESTRA_R2.fastq.gz

# 3. Descargar recursos genómicos

La forma más simple es ejecutar el script incluido, que descarga los datos de SRA y todos los recursos, y puede volver a ejecutarse si se interrumpe (saltea lo ya descargado):

```bash
mamba install -c conda-forge -c bioconda sra-tools bwa samtools bcftools htslib pigz
bash descargar_recursos.sh 2>&1 | tee descargas.log
```

El script realiza los siguientes pasos (espacio necesario: ~100 GB):

- Descarga las tres muestras de SRA y las renombra al formato del pipeline.
- Descarga el genoma de referencia hg38 y los recursos de GATK desde el repositorio público del Broad Institute (`gcp-public-data--broad-references`; el repositorio anterior, `genomics-public-data`, ya no permite el acceso libre) y genera el índice de BWA.
- Descarga el Panel of Normals y gnomAD (`gatk-best-practices/somatic-hg38`).
- Descarga ClinVar y **renombra los cromosomas** (`1` → `chr1`, `MT` → `chrM`) para que coincidan con el genoma de referencia. Sin este paso, la anotación con ClinVar no encuentra ninguna variante y no emite errores.

La caché de VEP (~25 GB) se descarga por separado; con `wget -c` la descarga puede retomarse si se corta:

```bash
mkdir -p resources/vep_cache && cd resources/vep_cache
wget -c https://ftp.ensembl.org/pub/release-111/variation/indexed_vep_cache/homo_sapiens_vep_111_GRCh38.tar.gz
tar -xzf homo_sapiens_vep_111_GRCh38.tar.gz && rm homo_sapiens_vep_111_GRCh38.tar.gz
cd ../..
```

# 4. Configurar el pipeline

```bash
# Editar config/config.yaml con las rutas a los recursos
nano config/config.yaml

# Editar config/samples.tsv con las muestras
nano config/samples.tsv

# Si se tiene un panel BED propio, se puede reemplazar:
cp /ruta/a/mi_panel.bed config/aml_panel.bed

```
## El pipeline incluye por defecto un panel de 23 genes de relevancia clínica en LMA según las guías ELN 2022 (config/aml_panel_full.bed). Pero se puede utilizar cualquier otro archivo BED con las regiones de interés deseadas.

## Notas de versión

- **Noviembre 2026:** se corrigió el script `scripts/acmg_classifier.py`. (1) Los campos de la anotación de VEP se leían según una lista fija que no coincidía con el orden real; ahora el orden se lee desde el encabezado del VCF. Con la versión anterior, la frecuencia de gnomAD y la puntuación CADD se leían de columnas equivocadas y las variantes de cambio de sentido se clasificaban erróneamente como benignas. (2) Los campos de ClinVar con comas (p. ej. `CLNREVSTAT`) se leían como tuplas y el número de estrellas resultaba siempre 0. (3) El criterio BP4 ya no se aplica cuando no hay predicción CADD disponible.

# 5. Ejecución
```bash
# Verificar que todo está en orden (dry-run)
snakemake --dry-run --cores 4 --use-conda

# Ejecución local
snakemake --cores 4 --use-conda

# Ejecución en segundo plano (recomendado para sesiones largas)
nohup snakemake --cores 4 --use-conda > logs/pipeline.log 2>&1 &
```
# 6. Estructura de outputs
```
results/
├── qc/
│   └── multiqc_report.html           # Reporte QC integrado de todas las muestras
├── {sample}/
│   ├── trimmed/                      # FASTQs filtrados por fastp
│   ├── bam/
│   │   └── {sample}.recal.bam        # BAM final calibrado (BQSR)
│   ├── qc/
│   │   ├── fastp/                    # Métricas de control de calidad
│   │   └── coverage/                 # Cobertura sobre el panel
│   ├── vcf/
│   │   └── {sample}.pass.vcf.gz      # Variantes filtradas (PASS)
│   ├── annotation/
│   │   └── {sample}.annotated.vcf.gz # VCF anotado (VEP + ClinVar)
│   ├── acmg/
│   │   ├── {sample}.acmg_classified.vcf.gz # VCF con clasificación ACMG
│   │   └── {sample}.acmg_classified.tsv     # Tabla de variantes clasificadas
│   └── report/
│       └── {sample}_report.html      # Reporte HTML por muestra
└── summary/
    └── all_variants_acmg.tsv         # Tabla resumen de todas las muestras
```

# 7. Interpretar la clasificación ACMG

| Clasificación | Significado |
|---|---|
| Pathogenic | Variante patogénica — requiere reporte clínico |
| Likely pathogenic | Probablemente patogénica — considerar acción |
| Uncertain significance (VUS) | Significado incierto — requiere seguimiento |
| Likely benign | Probablemente benigna |
| Benign | Benigna — descartada |
# Criterios específicos para LMA implementados en este pipeline: 
- AML_HOTSPOT: variante en posición recurrentemente mutada en LMA 
- AML_TIER1_GENE: gen de Tier 1 ELN 2022


# 8. Datos ùtiles:

- Modificando el config/aml_panel.bed con un archivo BED propio y actualizando la ruta en config/config.yaml se puede modificar el subset de genes.

- Editando las listas de genes AML_TIER1_GENES y AML_TIER2_GENES en scripts/acmg_classifier.py y modificando el archivo BED se pueden agregar mas genes de interés a la clasificación.

- El pipeline està diseñado para muestra modo tumor-only.
