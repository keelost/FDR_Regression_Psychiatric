# file: build_gene_trait_table.py
# Purpose: From a significant-results CSV, deduplicate genes and build a table of
#          Ensembl ID, UniProt ID, full name, and cognition/psychiatry traits.
# Requirements: pip install requests pandas

import os
import time
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

SUMMARY_DIR = os.path.join(os.environ.get("FDRREG_RESULTS_DIR", "data/pipeline"), "01.extra.analysis", "11.summary_table")
INPUT_CSV = os.path.join(SUMMARY_DIR, "significant_fdrreg_results.csv")
OUTPUT_CSV = os.path.join(SUMMARY_DIR, "gene_trait_table.csv")
# Intermediate cache so a crash/interruption does not force a full restart.
CACHE_CSV = os.path.join(SUMMARY_DIR, "gene_trait_cache.csv")

MYGENE_URL = "https://mygene.info/v3/query"
OT_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

MYGENE_BATCH_SIZE = 1000     # mygene POST batch size.
OT_MAX_WORKERS = 8           # Concurrent Open Targets requests.
OT_MAX_RETRIES = 3           # Retries per gene on transient failure.

# Diseases are kept if they belong to any of these Open Targets therapeutic areas.
RELEVANT_THERAPEUTIC_AREAS = {
    "nervous system disease",
    "psychiatric disorder",
}

# Measurement / phenotype terms live outside the disease branches above,
# so they are matched by keyword. Extend as needed.
MEASUREMENT_KEYWORDS = [
    "cognit", "intelligence", "educational attainment", "memory", "iq",
    "reasoning", "brain", "cortical", "cerebellar", "cerebral",
    "tau", "risk-taking", "risk taking", "insomnia", "sleep",
    "alcohol", "smoking", "nicotine", "wellbeing", "well-being",
    "neurotic", "mood", "depress", "anxiety",
]


def load_unique_genes(path):
    """Load the CSV and deduplicate by Ensembl gene ID."""
    df = pd.read_csv(path)
    # Keep the first occurrence of each Ensembl ID.
    df_unique = df.drop_duplicates(subset="ENSEMBL_GENE_ID", keep="first")
    # Keep only the columns we need downstream.
    genes = df_unique[["ENSEMBL_GENE_ID", "gene_name"]].reset_index(drop=True)
    print(f"Total rows: {len(df)}; unique genes: {len(genes)}")
    return genes


def annotate_by_ensembl(ensembl_ids):
    """Map Ensembl gene IDs to UniProt ID and full name via mygene.info, in batches."""
    annotations = {}
    for start in range(0, len(ensembl_ids), MYGENE_BATCH_SIZE):
        batch = ensembl_ids[start:start + MYGENE_BATCH_SIZE]
        params = {
            "q": ",".join(batch),
            "scopes": "ensembl.gene",
            "fields": "uniprot.Swiss-Prot,name",
            "species": "human",
        }
        resp = requests.post(MYGENE_URL, data=params, timeout=120)
        resp.raise_for_status()
        records = resp.json()

        for rec in records:
            query_id = rec.get("query")
            if rec.get("notfound"):
                annotations.setdefault(query_id, {"uniprot_id": "", "full_name": ""})
                continue

            uniprot = rec.get("uniprot", {})
            swiss = uniprot.get("Swiss-Prot", "") if isinstance(uniprot, dict) else ""
            uniprot_id = ";".join(swiss) if isinstance(swiss, list) else swiss

            # If an Ensembl ID maps to several mygene hits, keep the first non-empty.
            if query_id not in annotations or not annotations[query_id]["uniprot_id"]:
                annotations[query_id] = {
                    "uniprot_id": uniprot_id,
                    "full_name": rec.get("name", ""),
                }
        print(f"Annotated {min(start + MYGENE_BATCH_SIZE, len(ensembl_ids))}"
              f"/{len(ensembl_ids)} genes")
    return annotations


def fetch_ot_traits(ensembl_id):
    """Fetch associated diseases with therapeutic areas from Open Targets, with retries."""
    query = """
    query targetDiseases($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        associatedDiseases {
          rows {
            score
            disease {
              id
              name
              therapeuticAreas { id name }
            }
          }
        }
      }
    }
    """
    variables = {"ensemblId": ensembl_id}
    last_exc = None
    for attempt in range(OT_MAX_RETRIES):
        try:
            resp = requests.post(
                OT_GRAPHQL_URL,
                json={"query": query, "variables": variables},
                timeout=60,
            )
            resp.raise_for_status()
            data = resp.json().get("data", {})
            target = data.get("target")
            if not target or not target.get("associatedDiseases"):
                return []
            rows = target["associatedDiseases"]["rows"]
            return sorted(rows, key=lambda r: r.get("score", 0), reverse=True)
        except requests.RequestException as exc:
            last_exc = exc
            time.sleep(2 ** attempt)  # Exponential backoff: 1s, 2s, 4s.
    raise last_exc


def filter_traits(rows, min_score=0.0):
    """
    Keep a trait if either:
      1) its disease belongs to a relevant therapeutic area, or
      2) its name matches a neuro/cognition/behavior measurement keyword.
    """
    kept = []
    for r in rows:
        if r.get("score", 0) < min_score:
            continue

        disease = r["disease"]
        name = disease["name"]
        name_lc = name.lower()

        areas = {ta["name"].lower() for ta in disease.get("therapeuticAreas", [])}
        in_area = bool(areas & {a.lower() for a in RELEVANT_THERAPEUTIC_AREAS})
        by_keyword = any(kw in name_lc for kw in MEASUREMENT_KEYWORDS)

        if in_area or by_keyword:
            kept.append(name)

    seen, unique = set(), []
    for n in kept:
        if n not in seen:
            seen.add(n)
            unique.append(n)
    return unique


def _process_gene(ensembl_id, symbol, annotations):
    """Worker: fetch and filter traits for one gene. Returns an output row dict."""
    info = annotations.get(ensembl_id, {"uniprot_id": "", "full_name": ""})
    traits = []
    try:
        ot_rows = fetch_ot_traits(ensembl_id)
        traits = filter_traits(ot_rows)
    except requests.RequestException as exc:
        print(f"Open Targets query failed for {symbol} ({ensembl_id}): {exc}")
    return {
        "Gene symbol": symbol,
        "Ensembl ID": ensembl_id,
        "UniProt ID": info["uniprot_id"],
        "Full name": info["full_name"],
        "Associated cognitive Function and psychiatry-related traits in OpenTarget database":
            ", ".join(traits),
    }


def build_table(genes):
    ensembl_ids = genes["ENSEMBL_GENE_ID"].tolist()
    annotations = annotate_by_ensembl(ensembl_ids)

    # Resume support: skip genes already present in the cache.
    done = set()
    cached_rows = []
    if os.path.exists(CACHE_CSV):
        cached = pd.read_csv(CACHE_CSV)
        done = set(cached["Ensembl ID"].tolist())
        cached_rows = cached.to_dict("records")
        print(f"Resuming: {len(done)} genes already cached")

    todo = genes[~genes["ENSEMBL_GENE_ID"].isin(done)]
    output_rows = list(cached_rows)

    completed = 0
    with ThreadPoolExecutor(max_workers=OT_MAX_WORKERS) as executor:
        futures = {
            executor.submit(
                _process_gene, row["ENSEMBL_GENE_ID"], row["gene_name"], annotations
            ): row["ENSEMBL_GENE_ID"]
            for _, row in todo.iterrows()
        }
        for future in as_completed(futures):
            output_rows.append(future.result())
            completed += 1
            # Flush to cache every 50 genes.
            if completed % 50 == 0:
                pd.DataFrame(output_rows).to_csv(CACHE_CSV, index=False)
                print(f"Processed {completed}/{len(todo)} genes")

    # Reorder output to match the original gene order.
    df = pd.DataFrame(output_rows)
    order = {eid: i for i, eid in enumerate(ensembl_ids)}
    df["_order"] = df["Ensembl ID"].map(order)
    df = df.sort_values("_order").drop(columns="_order").reset_index(drop=True)
    return df


if __name__ == "__main__":
    genes = load_unique_genes(INPUT_CSV)
    df = build_table(genes)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Saved {len(df)} genes to {OUTPUT_CSV}")
