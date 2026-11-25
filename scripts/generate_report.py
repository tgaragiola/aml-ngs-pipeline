#!/usr/bin/env python3
"""
generate_report.py
Genera un informe preliminar en HTML por muestra, con las variantes detectadas,
su clasificación y la información de calidad necesaria para interpretarlas.

Versión 2 (nov. 2026). Cambios respecto de la versión original:
  - Se eliminó la columna "Tier ELN 2022": la ELN no clasifica variantes
    individuales sino el riesgo de cada paciente.
  - Se agregó la nomenclatura HGVS (transcrito, ADNc y proteína).
  - Se agregó un resumen de cobertura del panel y, por gen, la cobertura
    mediana de sus exones (a partir de los resultados de mosdepth). Los genes
    con cobertura insuficiente se informan como "no evaluables", no como
    "no detectados".
  - Se agregó una tabla de genes de relevancia clínica en LMA (OMS/ICC 2022,
    ELN 2022 y terapias dirigidas) con su resultado y su cobertura.
  - Se agregaron observaciones automáticas por variante (VAF baja, baja
    profundidad, VAF compatible con origen germinal, variantes contiguas que
    podrían ser un único evento, genes con terapia dirigida).
  - Se agregó una sección de alcance y limitaciones del análisis.
El informe sigue siendo PRELIMINAR: requiere la revisión de un profesional.
"""

import argparse
import gzip
import html
from collections import defaultdict
from datetime import date
from statistics import median
from urllib.parse import unquote

import pandas as pd

# ---------------------------------------------------------------------------
# Tablas de referencia
# ---------------------------------------------------------------------------

CLASIFICACION_ES = {
    "Pathogenic": "Patogénica",
    "Likely pathogenic": "Probablemente patogénica",
    "Uncertain significance": "Significado incierto",
    "Likely benign": "Probablemente benigna",
    "Benign": "Benigna",
}
ORDEN = list(CLASIFICACION_ES.keys())
COLOR = {
    "Pathogenic": "#b03a2e",
    "Likely pathogenic": "#d35400",
    "Uncertain significance": "#7f6000",
    "Likely benign": "#1e8449",
    "Benign": "#1e8449",
}

CONSECUENCIA_ES = {
    "missense_variant": "Cambio de sentido",
    "frameshift_variant": "Cambio del marco de lectura",
    "stop_gained": "Codón de terminación prematuro",
    "inframe_insertion": "Inserción en marco",
    "inframe_deletion": "Deleción en marco",
    "splice_acceptor_variant": "Sitio aceptor de empalme",
    "splice_donor_variant": "Sitio donante de empalme",
    "splice_region_variant": "Región de empalme",
    "synonymous_variant": "Sinónima",
    "5_prime_UTR_variant": "Región 5' no traducida",
    "3_prime_UTR_variant": "Región 3' no traducida",
    "intron_variant": "Intrónica",
    "start_lost": "Pérdida del codón de inicio",
    "stop_lost": "Pérdida del codón de terminación",
}

# Genes con terapia dirigida aprobada en LMA (capítulo 4.5 del trabajo)
GENES_TERAPIA = {"FLT3", "IDH1", "IDH2", "NPM1"}

# Relevancia clínica resumida de los genes (OMS/ICC 2022, ELN 2022, ELN 2024)
RELEVANCIA = {
    "NPM1": "Define entidad (OMS/ICC); riesgo ELN; terapia dirigida",
    "CEBPA": "Define entidad (región bZIP); riesgo favorable ELN",
    "FLT3": "Riesgo ELN; terapia dirigida (FLT3-ITD no evaluable con este flujo)",
    "IDH1": "Terapia dirigida",
    "IDH2": "Terapia dirigida",
    "TP53": "Entidad propia en ICC; riesgo adverso ELN (VAF ≥10 %)",
    "ASXL1": "Relacionado con mielodisplasia; riesgo adverso ELN",
    "EZH2": "Relacionado con mielodisplasia; riesgo adverso ELN",
    "RUNX1": "Relacionado con mielodisplasia; riesgo adverso ELN",
    "SF3B1": "Relacionado con mielodisplasia; riesgo adverso ELN",
    "SRSF2": "Relacionado con mielodisplasia; riesgo adverso ELN",
    "STAG2": "Relacionado con mielodisplasia; riesgo adverso ELN",
    "U2AF1": "Relacionado con mielodisplasia; riesgo adverso ELN",
    "KIT": "Señalización; no modifica el riesgo ELN en LMA-CBF",
    "KRAS": "Señalización; pronóstico en tratamientos de menor intensidad (ELN 2024)",
    "NRAS": "Señalización; pronóstico en tratamientos de menor intensidad (ELN 2024)",
    "PTPN11": "Señalización",
    "DNMT3A": "Metilación del ADN; frecuente en hematopoyesis clonal",
    "TET2": "Metilación del ADN; frecuente en hematopoyesis clonal",
    "KDM6A": "Modificador de la cromatina",
    "RAD21": "Complejo de cohesinas",
    "SMC1A": "Complejo de cohesinas",
    "SMC3": "Complejo de cohesinas",
}

# Genes de la ELN 2022 que no forman parte del panel
GENES_ELN_FUERA_DEL_PANEL = ["BCOR", "ZRSR2", "DDX41"]

UMBRAL_COBERTURA = 100   # lecturas
UMBRAL_VAF_BAJA = 0.05


# ---------------------------------------------------------------------------
# Funciones auxiliares
# ---------------------------------------------------------------------------

def texto(valor) -> str:
    """Convierte un valor a texto seguro para HTML; vacío si falta."""
    if valor is None or (isinstance(valor, float) and pd.isna(valor)):
        return ""
    v = str(valor)
    return "" if v in (".", "nan") else html.escape(v)


def formatear_hgvs(hgvs_c, hgvs_p) -> str:
    """Devuelve 'ENST...: c.XXX<br>p.(XXX)' a partir de los campos de VEP."""
    c = unquote(texto(hgvs_c))
    p = unquote(texto(hgvs_p))
    partes = []
    if c:
        transcrito, _, cambio = c.partition(":")
        partes.append(f"{cambio} <span class='tr'>({transcrito})</span>"
                      if cambio else c)
    if p:
        cambio_p = p.partition(":")[2] or p
        if cambio_p.startswith("p.") and not cambio_p.startswith("p.("):
            cambio_p = f"p.({cambio_p[2:]})"
        partes.append(cambio_p)
    return "<br>".join(partes) if partes else "—"


def leer_cobertura(regions_path, summary_path):
    """Lee mosdepth (cobertura media por región del panel).

    Para cada gen calcula la cobertura mediana de sus exones y cuántos exones
    tienen una cobertura media menor al umbral. Se usa la mediana por exón y no
    el porcentaje de bases, porque el archivo de regiones incluye los exones
    completos con sus regiones no traducidas, que el panel puede no capturar.
    """
    if not regions_path:
        return None, None, {}
    por_gen = defaultdict(list)
    with gzip.open(regions_path, "rt") as fh:
        for linea in fh:
            campos = linea.rstrip("\n").split("\t")
            por_gen[campos[3].replace("_exon", "")].append(float(campos[4]))
    media = None
    if summary_path:
        with open(summary_path) as fh:
            for linea in fh:
                campos = linea.split("\t")
                if campos[0] == "total_region":
                    media = float(campos[3])
    todas = [c for v in por_gen.values() for c in v]
    resumen = (sum(c >= UMBRAL_COBERTURA for c in todas), len(todas))
    gen = {g: (median(v), sum(c < UMBRAL_COBERTURA for c in v), len(v))
           for g, v in por_gen.items()}
    return media, resumen, gen


def observaciones(fila, df) -> str:
    """Observaciones automáticas para orientar la revisión profesional."""
    obs = []
    vaf, prof = float(fila["vaf"]), int(fila["depth"])
    criterios = str(fila.get("acmg_criteria", ""))
    if vaf < UMBRAL_VAF_BAJA:
        obs.append("VAF &lt;5 %: confirmar con un método independiente")
    if prof < UMBRAL_COBERTURA:
        obs.append(f"Profundidad &lt;{UMBRAL_COBERTURA}×")
    if 0.40 <= vaf <= 0.60 or vaf >= 0.90:
        if "AML_HOTSPOT" in criterios:
            obs.append("VAF compatible con origen germinal, pero es una "
                       "posición con mutaciones somáticas recurrentes en LMA")
        else:
            obs.append("VAF compatible con origen germinal: evaluar con "
                       "muestra normal pareada")
    # Variantes contiguas sostenidas por las mismas lecturas
    vecinas = df[(df["chrom"] == fila["chrom"])
                 & ((df["pos"] - fila["pos"]).abs().between(1, 3))
                 & (df["alt_reads"] == fila["alt_reads"])
                 & (df["depth"] == fila["depth"])]
    for _, v in vecinas.iterrows():
        cambio = unquote(texto(v["hgvs_c"])).partition(":")[2] or "variante contigua"
        obs.append(f"Posible evento complejo junto con {cambio}: revisar el alineamiento")
    if fila["gene"] in GENES_TERAPIA:
        obs.append("Gen con terapia dirigida aprobada")
    return "<br>".join(obs) if obs else "—"


# ---------------------------------------------------------------------------
# Programa principal
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True, help="Nombre de la muestra")
    parser.add_argument("--acmg-tsv", required=True, help="Tabla de variantes clasificadas")
    parser.add_argument("--regions", default=None, help="mosdepth regions.bed.gz")
    parser.add_argument("--coverage-summary", default=None, help="mosdepth summary.txt")
    parser.add_argument("--output", required=True, help="Archivo HTML de salida")
    args = parser.parse_args()

    df = pd.read_csv(args.acmg_tsv, sep="\t")
    df["acmg_classification"] = pd.Categorical(
        df["acmg_classification"], categories=ORDEN, ordered=True)
    df = df.sort_values(["acmg_classification", "gene", "pos"]).reset_index(drop=True)

    media, resumen, cob_gen = leer_cobertura(args.regions, args.coverage_summary)

    # --- Tabla de variantes ---
    filas = ""
    for _, f in df.iterrows():
        clas = str(f["acmg_classification"])
        filas += f"""
        <tr>
          <td><i>{texto(f['gene'])}</i></td>
          <td>{formatear_hgvs(f.get('hgvs_c'), f.get('hgvs_p'))}</td>
          <td>{CONSECUENCIA_ES.get(f['consequence'], texto(f['consequence']))}</td>
          <td class="num">{f"{float(f['vaf']) * 100:.1f}".replace('.', ',')} %</td>
          <td class="num">{int(f['alt_reads'])}/{int(f['depth'])}</td>
          <td style="color:{COLOR.get(clas, '#333')};font-weight:bold">{CLASIFICACION_ES.get(clas, clas)}</td>
          <td class="crit">{texto(f.get('acmg_criteria')).replace(',', ', ')}</td>
          <td class="obs">{observaciones(f, df)}</td>
        </tr>"""
    if not filas:
        filas = "<tr><td colspan='8'>No se detectaron variantes que superen los filtros de calidad.</td></tr>"

    # --- Tabla de genes de relevancia clínica ---
    detectados = defaultdict(list)
    for _, f in df.iterrows():
        cambio = unquote(texto(f.get("hgvs_p"))).partition(":")[2] \
            or unquote(texto(f.get("hgvs_c"))).partition(":")[2]
        detectados[f["gene"]].append(cambio)
    filas_genes = ""
    for gen, relevancia in RELEVANCIA.items():
        cob = cob_gen.get(gen)
        insuficiente = cob is not None and cob[0] < UMBRAL_COBERTURA
        if gen in detectados:
            resultado = ", ".join(detectados[gen])
        elif insuficiente:
            resultado = "<b>No evaluable: cobertura insuficiente</b>"
        else:
            resultado = "No detectada"
        if cob is None:
            cob_txt, exones_txt = "—", "—"
        else:
            cob_txt = f"{cob[0]:.0f}×"
            exones_txt = f"{cob[1]}/{cob[2]}"
        alerta = " baja" if insuficiente else ""
        filas_genes += f"""
        <tr><td><i>{gen}</i></td><td>{relevancia}</td><td>{resultado}</td>
            <td class="num{alerta}">{cob_txt}</td><td class="num">{exones_txt}</td></tr>"""

    cobertura_txt = ""
    if media is not None:
        estado = "adecuada" if media >= UMBRAL_COBERTURA else "INSUFICIENTE"
        cobertura_txt = (f"Cobertura media del panel: <b>{media:.0f}×</b> ({estado}; mínimo "
                         f"{UMBRAL_COBERTURA}×) &nbsp;|&nbsp; Exones con cobertura media "
                         f"≥{UMBRAL_COBERTURA}×: <b>{resumen[0]}/{resumen[1]}</b>")

    pagina = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<title>Informe preliminar NGS - {html.escape(args.sample)}</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 32px 40px; color: #1f2d3a; font-size: 12px; }}
  h1 {{ font-size: 21px; margin: 0 0 4px 0; }}
  h2 {{ font-size: 14px; margin: 22px 0 6px 0; border-bottom: 2px solid #1f2d3a; padding-bottom: 3px; }}
  .aviso {{ background: #fdecea; border-left: 5px solid #b03a2e; padding: 7px 10px; margin: 10px 0;
            font-weight: bold; color: #7b241c; }}
  .meta {{ color: #4d5d6c; margin: 3px 0; }}
  table {{ border-collapse: collapse; width: 100%; }}
  th {{ background: #1f2d3a; color: #fff; padding: 6px 7px; text-align: left; font-size: 11px; }}
  td {{ padding: 5px 7px; border-bottom: 1px solid #d5dbe0; vertical-align: top; }}
  td.num {{ text-align: right; white-space: nowrap; }}
  td.crit {{ font-size: 10.5px; }}
  td.obs {{ font-size: 10.5px; color: #4d5d6c; }}
  td.baja {{ color: #b03a2e; font-weight: bold; }}
  .tr {{ color: #7a8793; font-size: 10px; }}
  ul {{ margin: 4px 0; padding-left: 18px; }} li {{ margin: 2px 0; }}
  .pie {{ margin-top: 18px; font-size: 10.5px; color: #4d5d6c; }}
</style>
</head>
<body>
<h1>Informe preliminar de variantes somáticas por NGS — LMA</h1>
<p class="meta">Muestra: <b>{html.escape(args.sample)}</b> &nbsp;|&nbsp; Fecha de análisis: {date.today():%d/%m/%Y}
 &nbsp;|&nbsp; Referencia: GRCh38 &nbsp;|&nbsp; Anotación: Ensembl VEP 111, gnomAD, ClinVar
 &nbsp;|&nbsp; Panel: {len(RELEVANCIA)} genes (exones)</p>
<p class="meta">{cobertura_txt}</p>
<div class="aviso">Informe preliminar generado automáticamente por un flujo de trabajo NO validado para uso
clínico. Debe ser revisado e interpretado por un profesional antes de cualquier decisión.</div>

<h2>1. Variantes detectadas ({len(df)})</h2>
<table>
<thead><tr><th>Gen</th><th>Variante (HGVS)</th><th>Consecuencia</th><th>VAF</th><th>Lecturas<br>alt./total</th>
<th>Clasificación<br>(ACMG/AMP adaptada)</th><th>Criterios</th><th>Observaciones</th></tr></thead>
<tbody>{filas}
</tbody>
</table>

<h2>2. Genes de relevancia clínica en LMA incluidos en el panel</h2>
<table>
<thead><tr><th>Gen</th><th>Relevancia (OMS/ICC 2022, ELN 2022/2024)</th><th>Resultado</th>
<th>Cobertura<br>mediana</th><th>Exones<br>&lt;{UMBRAL_COBERTURA}×</th></tr></thead>
<tbody>{filas_genes}
</tbody>
</table>
<p class="meta">Genes de la ELN 2022 no incluidos en el panel: {", ".join(GENES_ELN_FUERA_DEL_PANEL)}.
La cobertura mediana se calcula sobre los exones de cada gen; cuando es menor
a {UMBRAL_COBERTURA}× se resalta en rojo y la ausencia de variantes no puede interpretarse como un resultado negativo.
Los exones con regiones no traducidas extensas pueden no estar incluidos en el diseño de captura.</p>

<h2>3. Alcance y limitaciones</h2>
<ul>
<li>No detecta reordenamientos cromosómicos ni genes de fusión (p. ej., <i>RUNX1::RUNX1T1</i>,
<i>CBFB::MYH11</i>, <i>PML::RARA</i>): requieren citogenética, FISH o RT-PCR.</li>
<li>No evalúa de forma confiable FLT3-ITD ni variaciones en el número de copias.</li>
<li>Análisis sin muestra normal pareada: el origen somático o germinal de las variantes no puede confirmarse.</li>
<li>Límite de detección no establecido (flujo no validado); umbral de detección configurado: VAF ≥1 %.</li>
<li>La clasificación es una adaptación automatizada de ACMG/AMP 2015, diseñada para variantes germinales. La
interpretación somática (AMP/ASCO/CAP 2017) y la integración con la clasificación OMS/ICC y el riesgo ELN de
cada paciente requieren la revisión de un profesional.</li>
</ul>
<p class="pie">Generado con aml-ngs-pipeline (github.com/tgaragiola/aml-ngs-pipeline).</p>
</body>
</html>"""

    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(pagina)
    print(f"Informe guardado en: {args.output}")


if __name__ == "__main__":
    main()
