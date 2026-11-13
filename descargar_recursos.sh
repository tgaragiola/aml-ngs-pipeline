#!/usr/bin/env bash
# =============================================================================
# Fase 2 - Descarga y preparación de datos y recursos para aml-ngs-pipeline
# Uso (dentro de ~/aml-ngs-pipeline, con el ambiente "snakemake" activo):
#     bash descargar_recursos.sh 2>&1 | tee descargas.log
# Se puede volver a ejecutar si se corta: saltea lo que ya está descargado.
# Uso: bash descargar_recursos.sh 2>&1 | tee descargas.log
# Espacio necesario: ~100 GB. Tiempo estimado: 4-8 h según la conexión.
# =============================================================================
set -euo pipefail
cd ~/aml-ngs-pipeline
THREADS=8

paso() { echo; echo "=============================================="; echo ">> $1  ($(date '+%H:%M'))"; echo "=============================================="; }

# Descarga un archivo; si es un índice .tbi que no existe en el servidor, lo genera con tabix
bajar() {  # bajar URL CARPETA
  local url=$1 dir=$2 f; f=$(basename "$url")
  [[ -s "$dir/$f" ]] && return 0
  if ! wget -c -q --show-progress -P "$dir" "$url"; then
    rm -f "$dir/$f"
    if [[ $f == *.tbi ]]; then
      echo "   $f no está en el servidor: se genera localmente con tabix"
      tabix -f -p vcf "$dir/${f%.tbi}"
    else
      echo "   ERROR al descargar $url"; return 1
    fi
  fi
}

mkdir -p data resources/{hg38,dbsnp,known_indels,pon,vep_cache,clinvar} tmp_sra

# -----------------------------------------------------------------------------
paso "1/6 Datos de secuenciación (SRA)"
declare -A MUESTRAS=( [SRR37858079]=CBF_AML_008 [SRR37858074]=CBF_AML_009 [SRR37858073]=CBF_AML_010 )
for SRR in "${!MUESTRAS[@]}"; do
  NOMBRE=${MUESTRAS[$SRR]}
  if [[ -s data/${NOMBRE}_R1.fastq.gz && -s data/${NOMBRE}_R2.fastq.gz ]]; then
    echo "   $NOMBRE ya descargada, se saltea."; continue
  fi
  echo "   Descargando $SRR -> $NOMBRE"
  prefetch --max-size 100G -O tmp_sra "$SRR"
  fasterq-dump --split-files -e $THREADS -t tmp_sra -O tmp_sra "tmp_sra/$SRR"
  pigz -p $THREADS "tmp_sra/${SRR}_1.fastq" "tmp_sra/${SRR}_2.fastq"
  mv "tmp_sra/${SRR}_1.fastq.gz" "data/${NOMBRE}_R1.fastq.gz"
  mv "tmp_sra/${SRR}_2.fastq.gz" "data/${NOMBRE}_R2.fastq.gz"
  rm -rf "tmp_sra/$SRR"
done
rm -rf tmp_sra

# -----------------------------------------------------------------------------
paso "2/6 Genoma de referencia hg38 + índice BWA (el índice tarda ~1 h)"
# Repositorio actualizado del Broad (el anterior, genomics-public-data, devuelve 403)
BASE="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
cd resources/hg38
for f in Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict; do
  bajar "$BASE/$f" .
done
[[ -s Homo_sapiens_assembly38.fasta.bwt ]] || bwa index Homo_sapiens_assembly38.fasta
ln -sf Homo_sapiens_assembly38.fasta hg38.fa
ln -sf Homo_sapiens_assembly38.fasta.fai hg38.fa.fai
ln -sf Homo_sapiens_assembly38.dict hg38.dict
for ext in amb ann bwt pac sa; do ln -sf Homo_sapiens_assembly38.fasta.$ext hg38.fa.$ext; done
cd ~/aml-ngs-pipeline

# -----------------------------------------------------------------------------
paso "3/6 Variantes conocidas para BQSR (GATK bundle)"
for f in Homo_sapiens_assembly38.dbsnp138.vcf.gz Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi; do
  bajar "$BASE/$f" resources/dbsnp
done
for f in Mills_and_1000G_gold_standard.indels.hg38.vcf.gz Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
         1000G_phase1.snps.high_confidence.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi; do
  bajar "$BASE/$f" resources/known_indels
done

# -----------------------------------------------------------------------------
paso "4/6 Panel of Normals y gnomAD (Mutect2)"
GATK_BASE="https://storage.googleapis.com/gatk-best-practices/somatic-hg38"
for f in 1000g_pon.hg38.vcf.gz 1000g_pon.hg38.vcf.gz.tbi; do
  bajar "$GATK_BASE/$f" resources/pon
done
for f in af-only-gnomad.hg38.vcf.gz af-only-gnomad.hg38.vcf.gz.tbi; do
  bajar "$GATK_BASE/$f" resources
done

# -----------------------------------------------------------------------------
paso "5/6 ClinVar (con cromosomas renombrados a 'chr1', 'chr2'... para que coincidan con hg38)"
cd resources/clinvar
if [[ ! -s clinvar.vcf.gz.tbi ]]; then
  bajar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz .
  mv clinvar.vcf.gz clinvar_original.vcf.gz
  # El VCF de ClinVar no declara los contigs en el encabezado, por lo que se
  # renombran directamente en las líneas de datos.
  zcat clinvar_original.vcf.gz \
    | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$1 = ($1=="MT") ? "chrM" : "chr"$1; print}' \
    | bgzip > clinvar.vcf.gz
  tabix -f -p vcf clinvar.vcf.gz
  rm -f clinvar_original.vcf.gz
fi
echo "   Comprobación (debe empezar con 'chr'): $(zcat clinvar.vcf.gz | grep -v '^#' | head -1 | cut -f1)"
cd ~/aml-ngs-pipeline

# -----------------------------------------------------------------------------
paso "6/6 Caché de Ensembl VEP v111 (~25 GB, la descarga más larga)"
if [[ ! -d resources/vep_cache/homo_sapiens/111_GRCh38 ]]; then
  cd resources/vep_cache
  wget -c https://ftp.ensembl.org/pub/release-111/variation/indexed_vep_cache/homo_sapiens_vep_111_GRCh38.tar.gz
  tar -xzf homo_sapiens_vep_111_GRCh38.tar.gz
  rm -f homo_sapiens_vep_111_GRCh38.tar.gz
  cd ~/aml-ngs-pipeline
fi

paso "LISTO. Resumen del espacio usado:"
du -sh data resources/* 2>/dev/null
echo
echo "Si todo terminó sin errores, ya se puede ejecutar el pipeline (ver README, sección 5)."
