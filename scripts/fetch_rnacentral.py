"""
📄 fetch_rnacentral.py
Functions to query RNAcentral for tRNA-Sec sequences, enrich with genes,
fetch related publications, and export results.
"""

import pandas as pd
import requests
from tqdm import tqdm
from typing import Optional, List
from collections import Counter
import os
from ratelimit import limits, sleep_and_retry

# Límite de la API de RNAcentral
CALLS = 20
PERIOD = 1

@sleep_and_retry
@limits(calls=CALLS, period=PERIOD)
def get_sequence(urs):
    seq_response = requests.get(
        f'https://rnacentral.org/api/v1/rna/{urs}.fasta',
        timeout=15
    )
    seq_response.raise_for_status()
    sequence_lines = seq_response.text.splitlines()[1:]
    sequence = "".join(sequence_lines)
    return sequence


def fetch_rna_sequences(
    custom_query: str,
    max_results: Optional[int] = None,
    fields: List[str] = ["description", "species", "rna_type", "so_rna_type_name", "expert_db", "length"],
    filename: Optional[str] = None,
    batch_size: int = 100
) -> pd.DataFrame:
    base_url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral'
    query_params = {
        'query': custom_query,
        'fields': ','.join(fields),
        'format': 'json',
        'size': batch_size,
        'start': 0
    }

    print(f"📋 Query: {custom_query}")

    # Primera consulta para obtener el total
    response = requests.get(base_url, params=query_params, timeout=30)
    response.raise_for_status()
    data = response.json()
    total_found = int(data['hitCount'])
    print(f"⚡ Encontrados {total_found} resultados")

    # Calcular número total de resultados
    if max_results and max_results < total_found:
        total_to_fetch = max_results
        print(f"⚡ Obteniendo {max_results} de {total_found} resultados")
    else:
        total_to_fetch = total_found
        print(f"⚡ Obteniendo todos los {total_found} resultados")

    # Extraer resultados con paginación
    results = []
    pbar = tqdm(total=total_to_fetch, desc="Extrayendo secuencias")
    start_index = 0

    while start_index < total_to_fetch and len(results) < total_to_fetch:
        query_params['start'] = start_index
        current_batch_size = min(batch_size, total_to_fetch - len(results))
        query_params['size'] = current_batch_size

        response = requests.get(base_url, params=query_params, timeout=30)
        response.raise_for_status()
        data = response.json()

        for entry in data['entries']:
            if len(results) >= total_to_fetch:
                break

            full_id = entry['id']
            urs = full_id.split('_')[0]

            metadata = {}
            for field in fields:
                field_data = entry['fields'].get(field, [])
                if field == "so_rna_type_name" and isinstance(field_data, list) and field_data:
                    metadata[field] = field_data[-1]
                elif isinstance(field_data, list) and len(field_data) > 1:
                    metadata[field] = "; ".join(str(x) for x in field_data)
                elif field_data:
                    metadata[field] = field_data[0] if isinstance(field_data, list) else field_data
                else:
                    metadata[field] = ""

            try:
                sequence = get_sequence(urs)
                results.append({
                    "URS_ID": full_id,
                    "Sequence": sequence,
                    "Length": len(sequence),
                    **metadata
                })
                pbar.update(1)
            except requests.exceptions.RequestException as e:
                print(f"⚠️ Error obteniendo secuencia {urs}: {e}")
                continue

        start_index += len(data['entries'])

    pbar.close()
    df = pd.DataFrame(results)

    if not df.empty and filename:
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        df.to_csv(filename, index=False, encoding='utf-8')
        print(f"💾 Datos guardados en: {filename}")

    return df


def enrich_with_genes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Enriquecer DataFrame con genes y publicaciones desde la API RNAcentral.
    """
    enriched_data = []

    for urs in tqdm(df["URS_ID"], desc="Enriqueciendo con genes"):
        try:
            url = f"https://rnacentral.org/api/v1/rna/{urs}"
            response = requests.get(url, timeout=15)
            response.raise_for_status()
            data = response.json()

            enriched_data.append({
                "URS_ID": urs,
                "genes": "; ".join(data.get("genes", [])) if data.get("genes") else "",
                "publications": data.get("publications", 0)
            })
        except requests.exceptions.RequestException as e:
            print(f"⚠️ Error en {urs}: {e}")
            enriched_data.append({
                "URS_ID": urs,
                "genes": "",
                "publications": 0
            })

    df_enriched = pd.DataFrame(enriched_data)
    return df.merge(df_enriched, on="URS_ID", how="left")
    
@sleep_and_retry
@limits(calls=20, period=1)
def get_publications(genes):
    """
    Devuelve publicaciones asociadas a una lista de genes desde RNAcentral LitScan.
    """
    if not genes:
        return pd.DataFrame(columns=["titles", "pmcids", "links"])

    format_gene_query = ['job_id:"{}"'.format(g) for g in genes]
    url = f"https://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral-litscan?query=entry_type:Publication AND ({' OR '.join(format_gene_query)})&fields=title,pmcid&format=json"
    response = requests.get(url, timeout=15)
    response.raise_for_status()
    r = response.json()

    if not r.get("entries"):
        return pd.DataFrame(columns=["titles", "pmcids", "links"])

    titles = [hit['fields']['title'][0] for hit in r['entries'] if 'fields' in hit and 'title' in hit['fields']]
    pmcids = [hit['fields']['pmcid'][0] for hit in r['entries'] if 'fields' in hit and 'pmcid' in hit['fields']]
    links = [f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/" for pmcid in pmcids]

    return pd.DataFrame({"titles": titles, "pmcids": pmcids, "links": links})


def export_publications(df: pd.DataFrame, out_path: str):
    """
    Genera un CSV con todas las publicaciones asociadas a cada URS_ID.
    """
    pub_records = []

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Extrayendo publicaciones"):
        genes = row["genes"].split("; ") if row["genes"] else []
        urs = row["URS_ID"]

        if genes:
            pub_df = get_publications(genes)
            if not pub_df.empty:
                for _, pub in pub_df.iterrows():
                    pub_records.append({
                        "URS_ID": urs,
                        "gene": "; ".join(genes),
                        "title": pub["titles"],
                        "pmcid": pub["pmcids"],
                        "link": pub["links"]
                    })

    if pub_records:
        df_pubs = pd.DataFrame(pub_records)
        df_pubs.to_csv(out_path, index=False, encoding="utf-8")
        print(f"✅ Publicaciones guardadas en: {out_path}")
    else:
        print("ℹ️ No se encontraron publicaciones asociadas")

def fetch_all_trna_sec(max_results: Optional[int] = None, filename: Optional[str] = None):
    query = 'tRNA AND so_rna_type_name:"Selenocysteinyl_tRNA"'
    return fetch_rna_sequences(
        custom_query=query,
        max_results=max_results,
        fields=["description", "species", "rna_type", "so_rna_type_name", "expert_db", "length"],
        filename=filename,
        batch_size=50
    )


def run_full_pipeline(maxsize: Optional[int] = None, outdir: str = "data/raw") -> pd.DataFrame:
    """
    Run the complete pipeline:
    1. Fetch tRNA-Sec sequences from RNAcentral.
    2. Enrich with gene information.
    3. Export results and associated publications.

    Args:
        maxsize (Optional[int]): Maximum number of sequences to retrieve (None = all).
        outdir (str): Output directory to save CSV files.

    Returns:
        pd.DataFrame: Enriched DataFrame with sequences and gene annotations.
    """
    os.makedirs(outdir, exist_ok=True)
    base_name = "all" if maxsize is None else str(maxsize)

    out_csv = os.path.join(outdir, f"{base_name}_trna_sec.csv")
    pubs_csv = os.path.join(outdir, f"{base_name}_trna_sec_publications.csv")

    print("=== Fetching tRNA-Sec with full metadata ===")

    # Step 1: Basic download
    df_sample = fetch_all_trna_sec(max_results=maxsize)

    if df_sample.empty:
        print("❌ No valid sequences found")
        return pd.DataFrame()

    # Step 2: Enrichment with genes
    df_full = enrich_with_genes(df_sample)
    df_full.to_csv(out_csv, index=False, encoding="utf-8")
    print(f"✅ Final CSV with genes saved to: {out_csv}")

    # Step 3: Export publications
    export_publications(df_full, pubs_csv)

    return df_full