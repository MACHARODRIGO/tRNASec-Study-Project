#!/usr/bin/env python3
"""
fetch_rnacentral.py

Descarga secuencias de tRNA-Sec desde RNAcentral y las anota con datos de
modificaciones de MODOMICS (si están disponibles).

Uso:
    python fetch_rnacentral.py --query "tRNA-Sec Homo sapiens" \
        --modifications-file data/raw/modomics_mods.csv \
        --output data/processed/annotated_trna_sec.csv
"""

import argparse
import logging
import csv
import json
from pathlib import Path
from typing import List, Dict, Optional

import requests
import pandas as pd

# Configurar logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def fetch_rnacentral_entries(query: str, limit: int = 100) -> List[Dict[str, str]]:
    """
    Consulta la API de RNAcentral para un término de búsqueda y devuelve los primeros
    resultados como una lista de diccionarios.

    Parameters
    ----------
    query: str
        Cadena de búsqueda, p. ej., "tRNA-Sec".
    limit: int
        Número máximo de entradas a recuperar.

    Returns
    -------
    List[Dict[str, str]]
        Lista de registros con campos URS y secuencia.
    """
    url = "https://rnacentral.org/api/v1/rna"
    params = {
        "query": query,
        "limit": limit,
        "format": "json"
    }
    logging.info("Consultando RNAcentral con query: %s", query)
    response = requests.get(url, params=params)
    response.raise_for_status()
    data = response.json()

    entries = []
    for hit in data.get("results", []):
        entry = {
            "urs_id": hit["rna_id"],
            "sequence": hit.get("sequence", "")
        }
        entries.append(entry)
    logging.info("Se obtuvieron %d secuencias de RNAcentral", len(entries))
    return entries


def load_modifications(modifications_file: str) -> pd.DataFrame:
    """
    Carga las modificaciones de MODOMICS desde un CSV. El CSV debe tener
    columnas 'urs_id' y 'modifications'.

    Parameters
    ----------
    modifications_file: str
        Ruta al archivo CSV con las modificaciones.

    Returns
    -------
    pd.DataFrame
        DataFrame con las modificaciones.
    """
    logging.info("Cargando modificaciones desde %s", modifications_file)
    df_mods = pd.read_csv(modifications_file)
    if "urs_id" not in df_mods.columns or "modifications" not in df_mods.columns:
        raise ValueError(
            "El archivo de modificaciones debe contener las columnas 'urs_id' y 'modifications'")
    return df_mods


def annotate_sequences(entries: List[Dict[str, str]],
                       mod_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """
    Anota una lista de entradas de RNAcentral con información de modificaciones.

    Parameters
    ----------
    entries: List[Dict[str, str]]
        Lista de entradas de RNAcentral (urs_id y sequence).
    mod_df: pd.DataFrame, opcional
        DataFrame con las modificaciones (urs_id y modifications).

    Returns
    -------
    pd.DataFrame
        DataFrame combinado.
    """
    df = pd.DataFrame(entries)
    if mod_df is not None:
        df = df.merge(mod_df, on="urs_id", how="left")
    return df


def parse_args() -> argparse.Namespace:
    """Analiza los argumentos de línea de comandos."""
    parser = argparse.ArgumentParser(description="Descarga y anota tRNA‑Sec de RNAcentral.")
    parser.add_argument("--query", required=True,
                        help="Cadena de búsqueda para la API de RNAcentral.")
    parser.add_argument("--output", required=True,
                        help="Ruta donde guardar el CSV anotado.")
    parser.add_argument("--modifications-file",
                        help="Archivo CSV con datos de modificaciones de MODOMICS (opcional).")
    parser.add_argument("--limit", type=int, default=200,
                        help="Número de secuencias a descargar (por defecto: 200).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    entries = fetch_rnacentral_entries(args.query, args.limit)

    mod_df = None
    if args.modifications_file:
        mod_df = load_modifications(args.modifications_file)

    annotated_df = annotate_sequences(entries, mod_df)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    annotated_df.to_csv(output_path, index=False)
    logging.info("Archivo guardado en: %s", output_path)


if __name__ == "__main__":
    main()
