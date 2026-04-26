#!/usr/bin/env python3
"""
acmg_classifier.py
==================
Clasificación de variantes según criterios ACMG/AMP 2015
con reglas específicas para Leucemia Mieloide Aguda (AML).

Referencia: Richards et al. Genetics in Medicine (2015) 17:405-424
            Li et al. J Mol Diagn (2017) - adaptación para variantes somáticas
"""

import argparse
import sys
import logging
from dataclasses import dataclass, field
from pathlib import Path

import pysam

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr
)
log = logging.getLogger(__name__)


# =============================================================================
# Estructura de datos
# =============================================================================

@dataclass
class ACMGEvidence:
    """Criterios ACMG/AMP activados para una variante."""
    # Criterios patogénicos muy fuertes
    PVS1: bool = False  # pérdida de función en gen donde LoF causa enfermedad

    # Criterios patogénicos fuertes
    PS1: bool = False   # mismo cambio aminoacídico que variante patogénica conocida
    PS3: bool = False   # estudios funcionales — usado como proxy para hotspots AML

    # Criterios patogénicos moderados
    PM1: bool = False   # localizado en dominio funcional crítico
    PM2: bool = False   # frecuencia muy baja en gnomAD
    PM4: bool = False   # indel in-frame en región no repetitiva
    PM5: bool = False   # missense diferente en misma posición que variante patogénica

    # Criterios patogénicos de soporte
    PP2: bool = False   # missense en gen con baja tasa de variantes benignas
    PP3: bool = False   # predicciones in silico deletéreas (CADD, SpliceAI)
    PP5: bool = False   # reportado como patogénico por fuente confiable

    # Criterios benignos
    BA1: bool = False   # frecuencia >5% en gnomAD — benigna directa
    BS1: bool = False   # frecuencia mayor de lo esperado para la enfermedad
    BP1: bool = False   # missense en gen donde solo LoF causa enfermedad
    BP3: bool = False   # indel in-frame en repetición
    BP4: bool = False   # predicciones in silico benignas
    BP6: bool = False   # reportado como benigno por fuente confiable
    BP7: bool = False   # sinónima sin efecto en splicing predicho

    # Criterios específicos para LMA (extensión del estándar ACMG)
    AML_HOTSPOT:    bool = False  # posición recurrentemente mutada en LMA
    AML_TIER1_GENE: bool = False  # gen de Tier 1 según ELN 2022

    def active_pathogenic(self) -> list:
        pvs = [c for c in ["PVS1"] if getattr(self, c)]
        ps  = [c for c in ["PS1","PS3"] if getattr(self, c)]
        pm  = [c for c in ["PM1","PM2","PM4","PM5"] if getattr(self, c)]
        pp  = [c for c in ["PP2","PP3","PP5"] if getattr(self, c)]
        aml = [c for c in ["AML_HOTSPOT","AML_TIER1_GENE"] if getattr(self, c)]
        return pvs + ps + pm + pp + aml

    def active_benign(self) -> list:
        ba  = [c for c in ["BA1"] if getattr(self, c)]
        bs  = [c for c in ["BS1"] if getattr(self, c)]
        bp  = [c for c in ["BP1","BP3","BP4","BP6","BP7"] if getattr(self, c)]
        return ba + bs + bp


@dataclass
class ClassifiedVariant:
    """Resultado de clasificación ACMG de una variante."""
    sample: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    consequence: str
    hgvs_c: str
    hgvs_p: str
    vaf: float
    depth: int
    alt_reads: int
    gnomad_af: float
    clinvar_sig: str
    clinvar_stars: int
    cadd_phred: float
    splice_ai: float
    evidence: ACMGEvidence = field(default_factory=ACMGEvidence)
    classification: str = "Uncertain significance"
    acmg_score: int = 0
    criteria_met: str = ""
    aml_tier: str = ""


# =============================================================================
# Hotspots AML conocidos
# =============================================================================

AML_HOTSPOTS = {
    # gen: lista de (hgvs_p_pattern, descripcion)
    "FLT3": [
        ("p.Asp835", "FLT3 D835 TKD mutation"),
        ("p.Ile836", "FLT3 I836 TKD mutation"),
        ("ITD",      "FLT3 Internal Tandem Duplication"),
    ],
    "NPM1": [
        ("p.Trp288Cysfs", "NPM1 type A insertion"),
        ("p.Trp288Leufs", "NPM1 type B insertion"),
        ("p.Trp288Argfs", "NPM1 type D insertion"),
        ("fs",            "NPM1 frameshift exon12"),
    ],
    "DNMT3A": [
        ("p.Arg882", "DNMT3A R882 hotspot"),
    ],
    "IDH1": [
        ("p.Arg132", "IDH1 R132 hotspot"),
    ],
    "IDH2": [
        ("p.Arg140", "IDH2 R140 hotspot"),
        ("p.Arg172", "IDH2 R172 hotspot"),
    ],
    "KIT": [
        ("p.Asp816", "KIT D816 hotspot"),
        ("p.Val559", "KIT V559 hotspot"),
    ],
    "NRAS": [
        ("p.Gly12", "NRAS G12 hotspot"),
        ("p.Gly13", "NRAS G13 hotspot"),
        ("p.Gln61", "NRAS Q61 hotspot"),
    ],
    "KRAS": [
        ("p.Gly12", "KRAS G12 hotspot"),
        ("p.Gly13", "KRAS G13 hotspot"),
    ],
    "TP53": [
        ("p.Arg248", "TP53 R248 hotspot"),
        ("p.Arg175", "TP53 R175 hotspot"),
        ("p.Arg273", "TP53 R273 hotspot"),
    ],
}

AML_TIER1_GENES = {
    "FLT3", "NPM1", "DNMT3A", "IDH1", "IDH2",
    "TET2", "RUNX1", "TP53", "ASXL1", "CEBPA"
}

AML_TIER2_GENES = {
    "NRAS", "KRAS", "PTPN11", "KIT", "SF3B1",
    "SRSF2", "U2AF1", "STAG2", "RAD21", "SMC1A",
    "SMC3", "EZH2", "KDM6A"
}

LOF_CONSEQUENCES = {
    "stop_gained", "frameshift_variant", "splice_acceptor_variant",
    "splice_donor_variant", "start_lost", "stop_lost",
    "transcript_ablation", "transcript_amplification"
}

MISSENSE_CONSEQUENCES = {"missense_variant"}

SYNONYMOUS_CONSEQUENCES = {
    "synonymous_variant", "stop_retained_variant"
}


# =============================================================================
# Parser de info VEP en el VCF
# =============================================================================

def parse_vep_info(info_dict: dict, vep_field: str = "CSQ") -> dict:
    """Extrae el primer transcript canónico de la anotación VEP."""
    csq_raw = info_dict.get(vep_field, "")
    if isinstance(csq_raw, tuple):
        csq_raw = ",".join(csq_raw)
    if not csq_raw:
        return {}

    # Header CSQ: Allele|Consequence|IMPACT|SYMBOL|Gene|...
    vep_keys = [
        "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene",
        "Feature", "BIOTYPE", "HGVSc", "HGVSp",
        "Existing_variation", "AF", "gnomADe_AF", "gnomADg_AF",
        "CADD_PHRED", "SpliceAI_pred_DS_AG", "CANONICAL",
        "SIFT", "PolyPhen", "DOMAINS",
    ]

    transcripts = csq_raw.split(",")
    for t in transcripts:
        fields = t.split("|")
        d = dict(zip(vep_keys, fields + [""] * max(0, len(vep_keys) - len(fields))))
        if d.get("CANONICAL") == "YES":
            return d
    # Si no hay canónico, tomar el primero
    fields = transcripts[0].split("|")
    return dict(zip(vep_keys, fields + [""] * max(0, len(vep_keys) - len(fields))))


def safe_float(val: str, default: float = 0.0) -> float:
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def get_clinvar_stars(revstat: str) -> int:
    """Convierte CLNREVSTAT de ClinVar a número de estrellas."""
    star_map = {
        "practice_guideline": 4,
        "reviewed_by_expert_panel": 3,
        "criteria_provided,_multiple_submitters,_no_conflicts": 2,
        "criteria_provided,_conflicting_interpretations": 1,
        "criteria_provided,_single_submitter": 1,
        "no_assertion_criteria_provided": 0,
        "no_assertion_provided": 0,
    }
    return star_map.get(revstat.lower().replace(" ", "_"), 0)


def is_aml_hotspot(gene: str, hgvs_p: str, consequence: str) -> bool:
    """Verifica si la variante cae en un hotspot conocido de AML."""
    if gene not in AML_HOTSPOTS:
        return False
    for pattern, _ in AML_HOTSPOTS[gene]:
        if pattern in hgvs_p or pattern in consequence:
            return True
    return False


# =============================================================================
# Motor de clasificación ACMG
# =============================================================================

class ACMGClassifier:
    def __init__(
        self,
        aml_genes: set,
        freq_benign: float = 0.05,
        freq_likely_benign: float = 0.01,
        freq_pathogenic: float = 0.001,
        clinvar_stars_min: int = 1,
    ):
        self.aml_genes = aml_genes
        self.freq_benign = freq_benign
        self.freq_likely_benign = freq_likely_benign
        self.freq_pathogenic = freq_pathogenic
        self.clinvar_stars_min = clinvar_stars_min

    def classify(self, variant: ClassifiedVariant) -> ClassifiedVariant:
        """Aplica todos los criterios ACMG/AMP y retorna la clasificación."""
        ev = variant.evidence
        gene = variant.gene
        csq  = variant.consequence
        gnomad_af = variant.gnomad_af
        hgvs_p = variant.hgvs_p
        cadd  = variant.cadd_phred
        splice_ai = variant.splice_ai
        clinvar_sig = variant.clinvar_sig.lower()
        clinvar_stars = variant.clinvar_stars

        # ── Criterios patogénicos ─────────────────────────────────────────────

        # PVS1: LoF en gen donde LoF es mecanismo conocido de enfermedad
        if any(c in csq for c in LOF_CONSEQUENCES):
            if gene in AML_TIER1_GENES or gene in AML_TIER2_GENES:
                ev.PVS1 = True

        # PS1: Misma variante aminoacídica que otra patogénica conocida
        if clinvar_stars >= self.clinvar_stars_min:
            if any(k in clinvar_sig for k in ["pathogenic", "likely_pathogenic"]):
                ev.PS1 = True

        # PS3: Estudios funcionales bien establecidos (hotspot AML = proxy)
        if is_aml_hotspot(gene, hgvs_p, csq):
            ev.PS3 = True
            ev.AML_HOTSPOT = True

        # PM1: Localizado en dominio funcional crítico sin variantes benignas
        if gene in AML_TIER1_GENES:
            ev.PM1 = True
            ev.AML_TIER1_GENE = True

        # PM2: Ausente o muy baja frecuencia en controles poblacionales
        if gnomad_af < self.freq_pathogenic:
            ev.PM2 = True

        # PM4: Indel en marco de lectura o cambio stop en región no repetitiva
        if "inframe_insertion" in csq or "inframe_deletion" in csq:
            ev.PM4 = True

        # PM5: Variante missense diferente en misma posición que una patogénica
        if "missense_variant" in csq and ev.PS1:
            ev.PM5 = True

        # PP2: Missense en gen con baja tasa de variantes missense benignas
        if "missense_variant" in csq and gene in AML_TIER1_GENES:
            ev.PP2 = True

        # PP3: Predicciones computacionales dañinas
        if cadd >= 20 or splice_ai >= 0.2:
            ev.PP3 = True

        # PP5: Reportado como patogénico por fuente confiable
        if ev.PS1 and clinvar_stars >= 2:
            ev.PP5 = True

        # ── Criterios benignos ────────────────────────────────────────────────

        # BA1: Frecuencia alélica >5% en gnomAD
        if gnomad_af > self.freq_benign:
            ev.BA1 = True

        # BS1: Frecuencia alélica mayor de lo esperado para la enfermedad
        elif gnomad_af > self.freq_likely_benign:
            ev.BS1 = True

        # BP1: Missense en gen donde solo LoF causa la enfermedad
        # (NO aplica en AML donde missense también es mecanismo)

        # BP4: Predicciones computacionales benignas
        if cadd < 10 and splice_ai < 0.05:
            ev.BP4 = True

        # BP6: Reportado como benigno por fuente confiable
        if clinvar_stars >= self.clinvar_stars_min:
            if any(k in clinvar_sig for k in ["benign", "likely_benign"]):
                ev.BP6 = True

        # BP7: Sinónima sin efecto en splicing predicho
        if any(c in csq for c in SYNONYMOUS_CONSEQUENCES) and splice_ai < 0.05:
            ev.BP7 = True

        # ── Scoring y clasificación final ─────────────────────────────────────
        variant.evidence = ev
        variant.classification, variant.acmg_score = self._score(ev)
        variant.criteria_met = ",".join(ev.active_pathogenic() + ev.active_benign())

        # AML tier
        if gene in AML_TIER1_GENES:
            variant.aml_tier = "Tier1_ELN2022"
        elif gene in AML_TIER2_GENES:
            variant.aml_tier = "Tier2_ELN2022"
        else:
            variant.aml_tier = "Other"

        return variant

    def _score(self, ev: ACMGEvidence) -> tuple[str, int]:
        """
        Implementa las reglas de combinación de criterios ACMG 2015.
        Retorna (clasificación, score numérico de patogenicidad).
        """
        pvs = ev.PVS1
        ps  = sum([ev.PS1, ev.PS3])
        pm  = sum([ev.PM1, ev.PM2, ev.PM4, ev.PM5])
        pp  = sum([ev.PP2, ev.PP3, ev.PP5])
        ba  = ev.BA1
        bs  = sum([ev.BS1])
        bp  = sum([ev.BP1, ev.BP3, ev.BP4, ev.BP6, ev.BP7])
        aml_boost = sum([ev.AML_HOTSPOT, ev.AML_TIER1_GENE])

        # Score numérico simplificado (para ordenar)
        score = pvs * 8 + ps * 4 + pm * 2 + pp * 1 + aml_boost * 2 - ba * 16 - bs * 4 - bp * 1

        # Benign stand-alone
        if ba:
            return "Benign", score

        # Pathogenic
        if (pvs and ps >= 1) or (pvs and pm >= 2) or (pvs and pm >= 1 and pp >= 1) or \
           (pvs and pp >= 2) or (ps >= 2) or (ps >= 1 and pm >= 3) or \
           (ps >= 1 and pm >= 2 and pp >= 2) or (ps >= 1 and pm >= 1 and pp >= 4):
            return "Pathogenic", score

        # Likely pathogenic
        if (pvs and pm == 1) or (ps == 1 and pm >= 1) or (ps == 1 and pp >= 2) or \
           (pm >= 3) or (pm == 2 and pp >= 2) or (pm == 1 and pp >= 4) or \
           (ev.AML_HOTSPOT and pm >= 1):
            return "Likely pathogenic", score

        # Likely benign
        if (bs >= 1 and bp >= 1) or (bp >= 2):
            return "Likely benign", score

        # Benign
        if (bs >= 2):
            return "Benign", score

        return "Uncertain significance", score


# =============================================================================
# Procesamiento del VCF
# =============================================================================

def process_vcf(
    input_vcf: str,
    output_vcf: str,
    output_tsv: str,
    classifier: ACMGClassifier,
    sample_name: str,
):
    vcf_in  = pysam.VariantFile(input_vcf)
    sample_id = list(vcf_in.header.samples)[0] if vcf_in.header.samples else sample_name

    # Añadir campos ACMG al header
    vcf_in.header.info.add("ACMG_CLASS",    1, "String", "ACMG classification")
    vcf_in.header.info.add("ACMG_CRITERIA", 1, "String", "ACMG criteria met")
    vcf_in.header.info.add("ACMG_SCORE",    1, "Integer","ACMG pathogenicity score")
    vcf_in.header.info.add("AML_TIER",      1, "String", "AML ELN 2022 gene tier")
    vcf_in.header.info.add("AML_HOTSPOT",   0, "Flag",   "AML recurrent hotspot")

    vcf_out = pysam.VariantFile(output_vcf + ".tmp.vcf", "w", header=vcf_in.header)
    tsv_rows = []

    tsv_header = [
        "sample","chrom","pos","ref","alt","gene","consequence",
        "hgvs_c","hgvs_p","vaf","depth","alt_reads",
        "gnomad_af","clinvar_sig","clinvar_stars",
        "cadd_phred","splice_ai","acmg_classification","acmg_score",
        "acmg_criteria","aml_tier"
    ]

    for rec in vcf_in.fetch():
        info = dict(rec.info)
        vep  = parse_vep_info(info)

        # Extraer VAF y profundidad del formato de muestra
        try:
            fmt = rec.samples[sample_id]
            ad  = fmt.get("AD", (0, 0))
            af  = fmt.get("AF", [0])[0] if fmt.get("AF") else (ad[1] / max(sum(ad), 1))
            dp  = sum(ad)
            alt_reads = ad[1]
        except Exception:
            af, dp, alt_reads = 0.0, 0, 0

        gnomad_af = safe_float(vep.get("gnomADg_AF") or vep.get("gnomADe_AF"))

        v = ClassifiedVariant(
            sample      = sample_name,
            chrom       = rec.chrom,
            pos         = rec.pos,
            ref         = rec.ref,
            alt         = rec.alts[0] if rec.alts else ".",
            gene        = vep.get("SYMBOL", "."),
            consequence = vep.get("Consequence", "."),
            hgvs_c      = vep.get("HGVSc", "."),
            hgvs_p      = vep.get("HGVSp", "."),
            vaf         = round(af, 4),
            depth       = dp,
            alt_reads   = alt_reads,
            gnomad_af   = gnomad_af,
            clinvar_sig = str(info.get("CLNSIG", ".")),
            clinvar_stars = get_clinvar_stars(str(info.get("CLNREVSTAT", ""))),
            cadd_phred  = safe_float(vep.get("CADD_PHRED")),
            splice_ai   = safe_float(vep.get("SpliceAI_pred_DS_AG")),
        )

        v = classifier.classify(v)

        # Escribir al VCF
        rec.info["ACMG_CLASS"]    = v.classification
        rec.info["ACMG_CRITERIA"] = v.criteria_met if v.criteria_met else "."
        rec.info["ACMG_SCORE"]    = v.acmg_score
        rec.info["AML_TIER"]      = v.aml_tier
        if v.evidence.AML_HOTSPOT:
            rec.info["AML_HOTSPOT"] = True
        vcf_out.write(rec)

        # Agregar a TSV
        tsv_rows.append([
            v.sample, v.chrom, v.pos, v.ref, v.alt,
            v.gene, v.consequence, v.hgvs_c, v.hgvs_p,
            v.vaf, v.depth, v.alt_reads,
            v.gnomad_af, v.clinvar_sig, v.clinvar_stars,
            v.cadd_phred, v.splice_ai,
            v.classification, v.acmg_score, v.criteria_met, v.aml_tier
        ])

    vcf_out.close()
    vcf_in.close()

    # Comprimir y ordenar VCF
    import subprocess
    subprocess.run(
        f"bcftools sort {output_vcf}.tmp.vcf | bgzip > {output_vcf}",
        shell=True, check=True
    )
    Path(f"{output_vcf}.tmp.vcf").unlink(missing_ok=True)

    # Escribir TSV
    with open(output_tsv, "w") as f:
        f.write("\t".join(tsv_header) + "\n")
        for row in tsv_rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    log.info(f"Clasificadas {len(tsv_rows)} variantes para {sample_name}")
    n_path = sum(1 for r in tsv_rows if "Pathogenic" in str(r[18]))
    log.info(f"  Patogénicas/Probablemente patogénicas: {n_path}")
    log.info(f"  VUS: {sum(1 for r in tsv_rows if r[18]=='Uncertain significance')}")


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="Clasificador ACMG para variantes de AML")
    p.add_argument("--input",              required=True)
    p.add_argument("--output-vcf",         required=True)
    p.add_argument("--output-tsv",         required=True)
    p.add_argument("--aml-genes",          default="")
    p.add_argument("--freq-benign",        type=float, default=0.05)
    p.add_argument("--freq-likely-benign", type=float, default=0.01)
    p.add_argument("--freq-pathogenic",    type=float, default=0.001)
    p.add_argument("--clinvar-stars-min",  type=int,   default=1)
    return p.parse_args()


def main():
    args = parse_args()

    aml_genes = set(args.aml_genes.split(",")) if args.aml_genes else set()
    aml_genes |= AML_TIER1_GENES | AML_TIER2_GENES

    classifier = ACMGClassifier(
        aml_genes          = aml_genes,
        freq_benign        = args.freq_benign,
        freq_likely_benign = args.freq_likely_benign,
        freq_pathogenic    = args.freq_pathogenic,
        clinvar_stars_min  = args.clinvar_stars_min,
    )

    sample_name = Path(args.input).stem.replace(".annotated", "").replace(".vcf", "")
    log.info(f"Iniciando clasificación ACMG para muestra: {sample_name}")

    process_vcf(
        input_vcf   = args.input,
        output_vcf  = args.output_vcf,
        output_tsv  = args.output_tsv,
        classifier  = classifier,
        sample_name = sample_name,
    )

    log.info("Clasificación completada.")


if __name__ == "__main__":
    main()
