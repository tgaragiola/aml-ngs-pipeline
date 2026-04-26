#!/usr/bin/env python3
"""
merge_summary.py
Combina las tablas de variantes clasificadas ACMG de todas las muestras
en una única tabla resumen ordenada por score de patogenicidad.
"""

import argparse
import pandas as pd
import sys

def main():
    # Argumentos de entrada
    p = argparse.ArgumentParser()
    p.add_argument("--inputs", nargs="+", required=True, 
                   help="Lista de archivos TSV de variantes clasificadas (uno por muestra)")
    p.add_argument("--output", required=True, 
                   help="Archivo TSV de salida con todas las variantes combinadas")
    args = p.parse_args()

    # Leer cada TSV y agregarlo a la lista
    dfs = []
    for f in args.inputs:
        try:
            df = pd.read_csv(f, sep="\t")
            dfs.append(df)
        except Exception as e:
            print(f"Warning: no se pudo leer {f}: {e}", file=sys.stderr)

    # Verificar que al menos un archivo fue leído correctamente
    if not dfs:
        print("Error: ningún TSV pudo ser leído.", file=sys.stderr)
        sys.exit(1)

    # Combinar todas las tablas en una sola
    merged = pd.concat(dfs, ignore_index=True)

    # Ordenar por score ACMG (mayor primero), luego por muestra y gen
    merged.sort_values(
        ["acmg_score", "sample", "gene"],
        ascending=[False, True, True],
        inplace=True
    )

    # Guardar la tabla combinada
    merged.to_csv(args.output, sep="\t", index=False)
    print(
        f"Resumen guardado: {args.output} "
        f"({len(merged)} variantes, {merged['sample'].nunique()} muestras)",
        file=sys.stderr
    )

if __name__ == "__main__":
    main()
