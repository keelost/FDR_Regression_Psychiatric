#!/usr/bin/env python3

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


TARGETS = [
    "adhd2016",
    "adhd2019",
    "asd2015.pgc",
    "asd2019",
    "bd2012",
    "bd2018",
    "mdd2013",
    "mdd2019",
    "mddco",
    "sa.bpd.pgc",
    "sa.ipsych",
    "sa.scz.pgc",
    "scz2012",
    "scz2014",
    "scz.eas2019",
    "wal.scz2018",
]

DEFAULT_ROOT = Path(os.environ.get("FDRREG_RESULTS_DIR", "data/pipeline"))

DEFAULT_OUTPUT_DIR = (
    DEFAULT_ROOT
    / "01.extra.analysis"
    / "09.summary_tables"
)

RISK_LOCI_DIR = (
    DEFAULT_ROOT
    / "01.extra.analysis"
    / "05.risk_loci"
)


def normalize_column_name(column_name):
    """Normalize a column name for robust matching."""
    column_name = str(column_name).replace("\ufeff", "").strip().lower()
    return re.sub(r"[^a-z0-9]", "", column_name)


def normalize_metric_name(metric_name):
    """Normalize a metric name for robust matching."""
    metric_name = str(metric_name).replace("\ufeff", "").strip().lower()
    return re.sub(r"[^a-z0-9]", "", metric_name)


def read_result_file(path, file_type):
    """Read a comma-separated or whitespace-separated result file."""
    if file_type == "csv":
        df = pd.read_csv(path, sep=",", dtype=str)
    elif file_type == "whitespace":
        df = pd.read_csv(
            path,
            sep=r"\s+",
            engine="python",
            dtype=str,
        )
    else:
        raise ValueError(f"Unsupported file type: {file_type}")

    df.columns = [
        str(column).replace("\ufeff", "").strip()
        for column in df.columns
    ]

    return df


def find_column(df, aliases):
    """Find a column using a list of possible aliases."""
    normalized_columns = {
        normalize_column_name(column): column
        for column in df.columns
    }

    for alias in aliases:
        normalized_alias = normalize_column_name(alias)

        if normalized_alias in normalized_columns:
            return normalized_columns[normalized_alias]

    return None


def count_below_threshold(df, column_name, threshold):
    """Count valid numeric values strictly below a threshold."""
    values = pd.to_numeric(
        df[column_name],
        errors="coerce",
    )

    return int((values < threshold).sum())


def make_base_result(
    target,
    analysis_type,
    population,
    brain_region,
    file_path,
):
    """Create a common output record."""
    return {
        "target": target,
        "analysis_type": analysis_type,
        "population": population,
        "brain_region": brain_region,
        "file_path": file_path,
        "n_records": None,
        "qval_metric_or_column": None,
        "qval_lt_0.01_count": None,
        "theoretical_fdr_metric_or_column": None,
        "theoretical_fdr_lt_0.01_count": None,
        "qval_plink_clumped_file": None,
        "theoretical_fdr_plink_clumped_file": None,
        "status": None,
        "message": None,
    }


def summarize_standard_file(
    path,
    target,
    analysis_type,
    brain_region,
    root_dir,
    threshold,
    file_type,
    qval_aliases,
    theoretical_fdr_aliases,
):
    """
    Summarize a MAGMA, MetaXcan, or S-MultiXcan result file.
    """
    if path.exists():
        try:
            relative_path = str(path.relative_to(root_dir))
        except ValueError:
            relative_path = str(path)
    else:
        relative_path = str(path)

    result = make_base_result(
        target=target,
        analysis_type=analysis_type,
        population="",
        brain_region=brain_region,
        file_path=relative_path,
    )

    if not path.exists():
        result["status"] = "missing_file"
        result["message"] = f"File does not exist: {path}"
        return result

    try:
        df = read_result_file(path, file_type)
    except Exception as exc:
        result["status"] = "read_error"
        result["message"] = str(exc)
        return result

    qval_column = find_column(df, qval_aliases)
    theoretical_fdr_column = find_column(
        df,
        theoretical_fdr_aliases,
    )

    result["n_records"] = len(df)
    result["qval_metric_or_column"] = qval_column
    result[
        "theoretical_fdr_metric_or_column"
    ] = theoretical_fdr_column

    missing_columns = []

    if qval_column is None:
        missing_columns.append(
            "Missing q-value column; expected one of: "
            + ", ".join(qval_aliases)
        )
    else:
        result["qval_lt_0.01_count"] = count_below_threshold(
            df=df,
            column_name=qval_column,
            threshold=threshold,
        )

    if theoretical_fdr_column is None:
        missing_columns.append(
            "Missing theoretical FDR column; expected one of: "
            + ", ".join(theoretical_fdr_aliases)
        )
    else:
        result[
            "theoretical_fdr_lt_0.01_count"
        ] = count_below_threshold(
            df=df,
            column_name=theoretical_fdr_column,
            threshold=threshold,
        )

    if missing_columns:
        result["status"] = "missing_column"
        result["message"] = "; ".join(missing_columns)
    else:
        result["status"] = "ok"
        result["message"] = ""

    return result


def find_risk_loci_value(
    df,
    metric_aliases,
    threshold,
):
    """
    Extract n_risk_loci for a metric at the specified threshold.

    Returns
    -------
    value : int or None
        The n_risk_loci value.
    matched_metric : str or None
        The metric name found in the input table.
    status : str
        Extraction status.
    message : str
        Additional information.
    """
    required_columns = {
        "metric",
        "threshold",
        "n_risk_loci",
    }

    missing_columns = required_columns.difference(df.columns)

    if missing_columns:
        return (
            None,
            None,
            "missing_column",
            "Missing columns: " + ", ".join(sorted(missing_columns)),
        )

    normalized_aliases = {
        normalize_metric_name(alias)
        for alias in metric_aliases
    }

    normalized_metrics = df["metric"].map(normalize_metric_name)

    threshold_values = pd.to_numeric(
        df["threshold"],
        errors="coerce",
    )

    threshold_mask = np.isclose(
        threshold_values,
        threshold,
        rtol=0.0,
        atol=1e-12,
        equal_nan=False,
    )

    metric_mask = normalized_metrics.isin(normalized_aliases)

    selected = df.loc[
        threshold_mask & metric_mask
    ].copy()

    if selected.empty:
        return (
            None,
            None,
            "missing_metric_threshold",
            (
                f"No row found for metrics {metric_aliases} "
                f"at threshold {threshold}"
            ),
        )

    selected["n_risk_loci_numeric"] = pd.to_numeric(
        selected["n_risk_loci"],
        errors="coerce",
    )

    selected = selected.dropna(
        subset=["n_risk_loci_numeric"]
    )

    if selected.empty:
        return (
            None,
            None,
            "invalid_n_risk_loci",
            "n_risk_loci could not be converted to numeric values",
        )

    matched_metrics = sorted(
        selected["metric"].astype(str).unique()
    )

    # One row is expected for each population and metric.
    # If duplicate rows exist, do not silently sum them.
    if len(selected) > 1:
        return (
            None,
            ",".join(matched_metrics),
            "duplicate_rows",
            (
                f"Found {len(selected)} rows for metric "
                f"and threshold {threshold}"
            ),
        )

    value = int(selected.iloc[0]["n_risk_loci_numeric"])
    matched_metric = str(selected.iloc[0]["metric"])

    return value, matched_metric, "ok", ""


def find_plink_clumped_file(
    target_risk_dir,
    target,
    metric_aliases,
    threshold,
):
    """Find a corresponding PLINK .clumped file."""
    threshold_label = f"{threshold:g}".replace(".", "p")

    candidates = []

    for metric in metric_aliases:
        patterns = [
            f"{target}.{metric}.thr{threshold_label}.clumped",
            f"{target}.{metric}*thr{threshold_label}*.clumped",
            f"*{metric}*thr{threshold_label}*.clumped",
        ]

        for pattern in patterns:
            candidates.extend(
                target_risk_dir.glob(pattern)
            )

    unique_candidates = sorted(
        set(candidates),
        key=lambda path: str(path),
    )

    if not unique_candidates:
        return ""

    return ";".join(str(path) for path in unique_candidates)


def summarize_risk_loci(
    summary_path,
    target,
    root_dir,
    threshold,
):
    """
    Extract risk-locus counts from a precomputed risk_loci_summary.csv.

    One output row is generated for each population.
    """
    if summary_path.exists():
        try:
            relative_path = str(summary_path.relative_to(root_dir))
        except ValueError:
            relative_path = str(summary_path)
    else:
        relative_path = str(summary_path)

    if not summary_path.exists():
        result = make_base_result(
            target=target,
            analysis_type="Risk_loci",
            population="",
            brain_region="ALL",
            file_path=relative_path,
        )
        result["status"] = "missing_file"
        result["message"] = (
            f"File does not exist: {summary_path}"
        )
        return [result]

    try:
        df = pd.read_csv(summary_path)
    except Exception as exc:
        result = make_base_result(
            target=target,
            analysis_type="Risk_loci",
            population="",
            brain_region="ALL",
            file_path=relative_path,
        )
        result["status"] = "read_error"
        result["message"] = str(exc)
        return [result]

    df.columns = [
        str(column).replace("\ufeff", "").strip()
        for column in df.columns
    ]

    if "population" not in df.columns:
        df["population"] = ""

    populations = sorted(
        df["population"]
        .fillna("")
        .astype(str)
        .unique()
    )

    results = []

    for population in populations:
        population_df = df.loc[
            df["population"].fillna("").astype(str) == population
        ].copy()

        result = make_base_result(
            target=target,
            analysis_type="Risk_loci",
            population=population,
            brain_region="ALL",
            file_path=relative_path,
        )

        result["n_records"] = len(population_df)

        qval_value, qval_metric, qval_status, qval_message = (
            find_risk_loci_value(
                df=population_df,
                metric_aliases=[
                    "qval_bh",
                    "qval",
                ],
                threshold=threshold,
            )
        )

        fdr_value, fdr_metric, fdr_status, fdr_message = (
            find_risk_loci_value(
                df=population_df,
                metric_aliases=[
                    "fdr_theoretical",
                    "FDR_theoretical",
                    "FDR.theoretical",
                    "FDR.the",
                ],
                threshold=threshold,
            )
        )

        result["qval_metric_or_column"] = qval_metric
        result["qval_lt_0.01_count"] = qval_value

        result[
            "theoretical_fdr_metric_or_column"
        ] = fdr_metric
        result[
            "theoretical_fdr_lt_0.01_count"
        ] = fdr_value

        target_risk_dir = summary_path.parent

        result["qval_plink_clumped_file"] = (
            find_plink_clumped_file(
                target_risk_dir=target_risk_dir,
                target=target,
                metric_aliases=[
                    "qval_bh",
                    "qval",
                ],
                threshold=threshold,
            )
        )

        result["theoretical_fdr_plink_clumped_file"] = (
            find_plink_clumped_file(
                target_risk_dir=target_risk_dir,
                target=target,
                metric_aliases=[
                    "fdr_theoretical",
                    "FDR_theoretical",
                    "FDR.the",
                ],
                threshold=threshold,
            )
        )

        statuses = [qval_status, fdr_status]
        messages = [
            message
            for message in [qval_message, fdr_message]
            if message
        ]

        if all(status == "ok" for status in statuses):
            result["status"] = "ok"
            result["message"] = ""
        else:
            result["status"] = ";".join(
                status
                for status in statuses
                if status != "ok"
            )
            result["message"] = "; ".join(messages)

        results.append(result)

    return results


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize risk loci, MAGMA, MetaXcan, and "
            "S-MultiXcan results at a specified threshold."
        )
    )

    parser.add_argument(
        "--root-dir",
        type=Path,
        default=DEFAULT_ROOT,
        help="Root directory containing all target directories.",
    )

    parser.add_argument(
        "--risk-loci-dir",
        type=Path,
        default=RISK_LOCI_DIR,
        help="Directory containing risk-loci summary files.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for the output CSV.",
    )

    parser.add_argument(
        "--output-name",
        default="fdr_counts_summary.csv",
        help="Output CSV file name.",
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=0.01,
        help="Threshold used for all analyses.",
    )

    args = parser.parse_args()

    root_dir = args.root_dir.resolve()
    risk_loci_dir = args.risk_loci_dir.resolve()
    output_dir = args.output_dir.resolve()

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    results = []

    for target in TARGETS:
        target_dir = root_dir / target

        # Risk-locus summary.
        risk_loci_summary_file = (
            risk_loci_dir
            / target
            / f"{target}.risk_loci_summary.csv"
        )

        results.extend(
            summarize_risk_loci(
                summary_path=risk_loci_summary_file,
                target=target,
                root_dir=root_dir,
                threshold=args.threshold,
            )
        )

        # MAGMA gene-level results.
        magma_file = (
            target_dir
            / "05.magma_fdrreg"
            / f"{target}.gene.bio.fdrreg.txt"
        )

        results.append(
            summarize_standard_file(
                path=magma_file,
                target=target,
                analysis_type="MAGMA",
                brain_region="ALL",
                root_dir=root_dir,
                threshold=args.threshold,
                file_type="whitespace",
                qval_aliases=[
                    "qval",
                ],
                theoretical_fdr_aliases=[
                    "FDR.the",
                    "FDR_theoretical",
                    "FDR.theoretical",
                    "fdr_theoretical",
                ],
            )
        )

        # MetaXcan, one row for each brain region.
        metaxcan_dir = (
            target_dir
            / "07.metaxcan_fdrreg"
            / "01.fdrreg_results"
        )

        metaxcan_files = sorted(
            metaxcan_dir.glob(
                "*/*.gene.bio.fdrreg.txt"
            )
        )

        if metaxcan_files:
            for metaxcan_file in metaxcan_files:
                results.append(
                    summarize_standard_file(
                        path=metaxcan_file,
                        target=target,
                        analysis_type="MetaXcan",
                        brain_region=metaxcan_file.parent.name,
                        root_dir=root_dir,
                        threshold=args.threshold,
                        file_type="whitespace",
                        qval_aliases=[
                            "qval",
                        ],
                        theoretical_fdr_aliases=[
                            "FDR_theoretical",
                            "FDR.theoretical",
                            "fdr_theoretical",
                        ],
                    )
                )
        else:
            missing_metaxcan_file = (
                metaxcan_dir / "NO_METAXCAN_FILE"
            )

            results.append(
                summarize_standard_file(
                    path=missing_metaxcan_file,
                    target=target,
                    analysis_type="MetaXcan",
                    brain_region="ALL",
                    root_dir=root_dir,
                    threshold=args.threshold,
                    file_type="whitespace",
                    qval_aliases=[
                        "qval",
                    ],
                    theoretical_fdr_aliases=[
                        "FDR_theoretical",
                        "FDR.theoretical",
                        "fdr_theoretical",
                    ],
                )
            )

        # S-MultiXcan gene-level results.
        smultixcan_file = (
            target_dir
            / "09.smultixcan_fdrreg"
            / f"{target}.gene.bio.fdrreg.txt"
        )

        results.append(
            summarize_standard_file(
                path=smultixcan_file,
                target=target,
                analysis_type="S-MultiXcan",
                brain_region="ALL",
                root_dir=root_dir,
                threshold=args.threshold,
                file_type="whitespace",
                qval_aliases=[
                    "qval",
                ],
                theoretical_fdr_aliases=[
                    "FDR.the",
                    "FDR_theoretical",
                    "FDR.theoretical",
                    "fdr_theoretical",
                ],
            )
        )

    summary_df = pd.DataFrame(results)

    sort_columns = [
        "target",
        "analysis_type",
        "population",
        "brain_region",
    ]

    summary_df = summary_df.sort_values(
        sort_columns,
        kind="stable",
    ).reset_index(drop=True)

    output_file = output_dir / args.output_name

    summary_df.to_csv(
        output_file,
        index=False,
    )

    print(f"Summary written to: {output_file}")
    print(f"Number of summary rows: {len(summary_df)}")

    problem_df = summary_df.loc[
        summary_df["status"] != "ok"
    ]

    if not problem_df.empty:
        print(
            "\nWarning: some inputs could not be fully processed:",
            file=sys.stderr,
        )

        print(
            problem_df[
                [
                    "target",
                    "analysis_type",
                    "population",
                    "brain_region",
                    "status",
                    "message",
                ]
            ].to_string(index=False),
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
