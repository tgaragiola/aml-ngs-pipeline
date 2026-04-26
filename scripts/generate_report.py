#!/usr/bin/env python3
"""
generate_report.py
Genera un reporte HTML por muestra con las variantes detectadas
y su clasificación ACMG, para facilitar la interpretación clínica.
"""

import argparse
import pandas as pd
from datetime import date

def main():
    # --- Argumentos de entrada ---
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample",   required=True, help="Nombre de la muestra")
    parser.add_argument("--acmg-tsv", required=True, help="Tabla de variantes clasificadas")
    parser.add_argument("--output",   required=True, help="Archivo HTML de salida")
    args = parser.parse_args()

    # --- Leer la tabla de variantes ---
    df = pd.read_csv(args.acmg_tsv, sep="\t")

    # --- Ordenar variantes: patogénicas primero, benignas al final ---
    orden = ["Pathogenic", "Likely pathogenic", "Uncertain significance",
             "Likely benign", "Benign"]
    df["acmg_classification"] = pd.Categorical(
        df["acmg_classification"], categories=orden, ordered=True
    )
    df = df.sort_values("acmg_classification")

    # --- Generar las filas de la tabla HTML ---
    # Cada variante se muestra en una fila con color según su clasificación:
    # rojo = patogénica, naranja = VUS, verde = benigna
    filas = ""
    for _, row in df.iterrows():
        clasificacion = row.get("acmg_classification", ".")

        if "Pathogenic" in str(clasificacion):
            color = "#c0392b"  # rojo
        elif "Uncertain" in str(clasificacion):
            color = "#e67e22"  # naranja
        else:
            color = "#27ae60"  # verde

        filas += f"""
        <tr>
            <td>{row.get('gene', '.')}</td>
            <td>{row.get('consequence', '.')}</td>
            <td>{float(row.get('vaf', 0)):.2%}</td>
            <td>{int(row.get('depth', 0))}</td>
            <td>{int(row.get('alt_reads', 0))}</td>
            <td style="color:{color}; font-weight:bold">{clasificacion}</td>
            <td>{row.get('acmg_criteria', '.')}</td>
            <td>{row.get('aml_tier', '.')}</td>
        </tr>"""

    # --- Construir el HTML completo ---
    html = f"""<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <title>Reporte NGS - {args.sample}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; color: #2c3e50; }}
        h1 {{ color: #2c3e50; }}
        p  {{ font-size: 13px; color: #7f8c8d; }}
        table {{ border-collapse: collapse; width: 100%; margin-top: 20px; }}
        th {{ background-color: #2c3e50; color: white; padding: 10px; text-align: left; font-size: 12px; }}
        td {{ padding: 8px 10px; border-bottom: 1px solid #ddd; font-size: 12px; }}
        tr:hover {{ background-color: #f5f5f5; }}
    </style>
</head>
<body>
    <h1>Reporte de variantes NGS — LMA</h1>
    <p>Muestra: <strong>{args.sample}</strong> &nbsp;|&nbsp;
       Fecha: {date.today()} &nbsp;|&nbsp;
       Clasificacion: ACMG/AMP 2015 · ELN 2022 · hg38</p>
    <p>Total de variantes detectadas: <strong>{len(df)}</strong></p>

    <table>
        <thead>
            <tr>
                <th>Gen</th>
                <th>Consecuencia</th>
                <th>VAF</th>
                <th>Profundidad</th>
                <th>Reads alt.</th>
                <th>Clasificacion ACMG</th>
                <th>Criterios</th>
                <th>Tier ELN 2022</th>
            </tr>
        </thead>
        <tbody>
            {filas}
        </tbody>
    </table>

    <p style="margin-top:30px; font-size:11px; color:#aaa">
        Reporte generado automaticamente. Requiere revision por profesional calificado.
    </p>
</body>
</html>"""

    # --- Guardar el archivo HTML ---
    with open(args.output, "w") as f:
        f.write(html)
    print(f"Reporte guardado en: {args.output}")

if __name__ == "__main__":
    main()
