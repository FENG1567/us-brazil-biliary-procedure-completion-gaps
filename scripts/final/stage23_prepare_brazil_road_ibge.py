#!/usr/bin/env python3
"""Build Brazil road-time and IBGE municipality covariates."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd
import requests


IBGE_BASE = "https://servicodados.ibge.gov.br/api/v3/agregados"
OSRM_BASE = "https://router.project-osrm.org/table/v1/driving"
USER_AGENT = "biliary-gap-research/1.0 (academic research; contact 18940835066@163.com)"


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, index=False)
    os.replace(tmp, path)


def fetch_json(session: requests.Session, url: str, path: Path) -> object:
    if path.exists() and path.stat().st_size > 100:
        return json.loads(path.read_text(encoding="utf-8"))
    response = session.get(url, timeout=120)
    response.raise_for_status()
    payload = response.json()
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    os.replace(tmp, path)
    return payload


def flatten_ibge(payload: object, value_name: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for variable in payload:
        unit = variable.get("unidade")
        for result in variable.get("resultados", []):
            for series in result.get("series", []):
                locality = series.get("localidade", {})
                for period, value in series.get("serie", {}).items():
                    rows.append(
                        {
                            "municipality_code": str(locality.get("id", "")),
                            "municipality_name": locality.get("nome"),
                            "period": int(period),
                            value_name: pd.to_numeric(value, errors="coerce"),
                            f"{value_name}_unit": unit,
                        }
                    )
    return pd.DataFrame(rows)


def macroregion(code: str) -> str:
    try:
        state = int(str(code)[:2])
    except ValueError:
        return "Unknown"
    if 11 <= state <= 17:
        return "North"
    if 21 <= state <= 29:
        return "Northeast"
    if 31 <= state <= 35:
        return "Southeast"
    if 41 <= state <= 43:
        return "South"
    if 50 <= state <= 53:
        return "Central-West"
    return "Unknown"


def safe_qcut(series: pd.Series) -> pd.Series:
    valid = series.dropna()
    result = pd.Series(pd.NA, index=series.index, dtype="Int64")
    if valid.nunique() < 4:
        return result
    ranked = valid.rank(method="first")
    result.loc[valid.index] = pd.qcut(ranked, 4, labels=[1, 2, 3, 4]).astype("Int64")
    return result


def load_pairs(con: duckdb.DuckDBPyConnection, source: Path) -> pd.DataFrame:
    query = f"""
        SELECT DISTINCT
          CAST(origin_muni AS VARCHAR) AS origin_muni,
          CAST(destination_muni AS VARCHAR) AS destination_muni,
          CAST(origin_longitude AS DOUBLE) AS origin_longitude,
          CAST(origin_latitude AS DOUBLE) AS origin_latitude,
          CAST(destination_longitude AS DOUBLE) AS destination_longitude,
          CAST(destination_latitude AS DOUBLE) AS destination_latitude,
          CAST(centroid_distance_km AS DOUBLE) AS centroid_distance_km
        FROM read_parquet('{source.as_posix()}')
        WHERE cross_system_primary_eligible
          AND prior_window_complete
          AND origin_muni IS NOT NULL
          AND destination_muni IS NOT NULL
          AND isfinite(CAST(origin_longitude AS DOUBLE))
          AND isfinite(CAST(origin_latitude AS DOUBLE))
          AND isfinite(CAST(destination_longitude AS DOUBLE))
          AND isfinite(CAST(destination_latitude AS DOUBLE))
        ORDER BY origin_muni, destination_muni
    """
    return con.execute(query).fetchdf()


def route_batches(
    session: requests.Session,
    pairs: pd.DataFrame,
    cache_path: Path,
    max_sources: int,
    sleep_seconds: float,
    max_pairs: int | None,
) -> pd.DataFrame:
    key_cols = ["origin_muni", "destination_muni"]
    columns = key_cols + [
        "road_duration_seconds",
        "road_distance_m",
        "origin_snap_m",
        "destination_snap_m",
        "route_status",
        "retrieved_utc",
        "osrm_service",
    ]
    if cache_path.exists():
        cache = pd.read_csv(cache_path, dtype={"origin_muni": str, "destination_muni": str})
    else:
        cache = pd.DataFrame(columns=columns)

    local = pairs[pairs.origin_muni == pairs.destination_muni].copy()
    local_rows = pd.DataFrame(
        {
            "origin_muni": local.origin_muni,
            "destination_muni": local.destination_muni,
            "road_duration_seconds": 0.0,
            "road_distance_m": 0.0,
            "origin_snap_m": 0.0,
            "destination_snap_m": 0.0,
            "route_status": "same_municipality",
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "osrm_service": "not_queried_same_municipality",
        }
    )
    cache = pd.concat([cache, local_rows], ignore_index=True)
    cache = cache.drop_duplicates(key_cols, keep="first")

    external = pairs[pairs.origin_muni != pairs.destination_muni].copy()
    already = set(zip(cache.origin_muni.astype(str), cache.destination_muni.astype(str)))
    pending = external[
        ~external.apply(lambda row: (str(row.origin_muni), str(row.destination_muni)) in already, axis=1)
    ].copy()
    if max_pairs is not None:
        pending = pending.head(max_pairs)

    coord_by_muni: dict[str, tuple[float, float]] = {}
    for _, row in pairs.iterrows():
        coord_by_muni[str(row.origin_muni)] = (float(row.origin_longitude), float(row.origin_latitude))
        coord_by_muni[str(row.destination_muni)] = (
            float(row.destination_longitude),
            float(row.destination_latitude),
        )

    # Pack several origins into one rectangular OSRM table request.  The
    # request returns all origin x destination combinations, but only observed
    # edges are retained.  Limits keep each request below 100 coordinates and
    # 1,200 matrix cells, which is substantially faster than a square matrix
    # containing one source and one destination for every observed edge.
    origin_groups: list[tuple[str, list[str]]] = []
    for origin, group in pending.groupby("origin_muni", sort=False):
        destinations_for_origin = list(dict.fromkeys(group.destination_muni.astype(str)))
        for start in range(0, len(destinations_for_origin), 80):
            origin_groups.append((str(origin), destinations_for_origin[start : start + 80]))

    batches: list[list[tuple[str, list[str]]]] = []
    current: list[tuple[str, list[str]]] = []
    current_origins: set[str] = set()
    current_destinations: set[str] = set()
    for origin, destinations_for_origin in origin_groups:
        proposed_origins = current_origins | {origin}
        proposed_destinations = current_destinations | set(destinations_for_origin)
        too_large = (
            len(proposed_origins) > max_sources
            or len(proposed_origins) + len(proposed_destinations) > 100
            or len(proposed_origins) * len(proposed_destinations) > 1200
        )
        if current and too_large:
            batches.append(current)
            current = []
            current_origins = set()
            current_destinations = set()
        current.append((origin, destinations_for_origin))
        current_origins.add(origin)
        current_destinations.update(destinations_for_origin)
    if current:
        batches.append(current)

    new_rows: list[dict[str, object]] = []
    total_batches = len(batches)
    for batch_number, batch_groups in enumerate(batches, start=1):
        origin_order = list(dict.fromkeys(origin for origin, _ in batch_groups))
        destination_order = list(
            dict.fromkeys(destination for _, destination_list in batch_groups for destination in destination_list)
        )
        coordinates = [
            f"{coord_by_muni[municipality][0]:.6f},{coord_by_muni[municipality][1]:.6f}"
            for municipality in origin_order + destination_order
        ]
        sources = [str(index) for index in range(len(origin_order))]
        destinations = [str(len(origin_order) + index) for index in range(len(destination_order))]
        url = OSRM_BASE + "/" + ";".join(coordinates)
        params = {
            "annotations": "duration,distance",
            "sources": ";".join(sources),
            "destinations": ";".join(destinations),
        }
        payload = None
        error = "unknown_error"
        for attempt in range(1, 5):
            try:
                response = session.get(url, params=params, timeout=180)
                if response.status_code == 200:
                    candidate = response.json()
                    if candidate.get("code") == "Ok":
                        payload = candidate
                        break
                    error = f"osrm_{candidate.get('code', 'unknown')}"
                else:
                    error = f"http_{response.status_code}"
            except Exception as exc:  # network errors are recorded, not silently dropped
                error = type(exc).__name__
            time.sleep(min(30.0, attempt * 3.0))

        retrieved = datetime.now(timezone.utc).isoformat()
        if payload is None:
            for origin, destination_list in batch_groups:
                for destination in destination_list:
                    new_rows.append(
                        {
                            "origin_muni": origin,
                            "destination_muni": destination,
                            "road_duration_seconds": math.nan,
                            "road_distance_m": math.nan,
                            "origin_snap_m": math.nan,
                            "destination_snap_m": math.nan,
                            "route_status": error,
                            "retrieved_utc": retrieved,
                            "osrm_service": OSRM_BASE,
                        }
                    )
        else:
            duration = payload.get("durations", [])
            distance = payload.get("distances", [])
            source_meta = payload.get("sources", [])
            destination_meta = payload.get("destinations", [])
            origin_index = {municipality: index for index, municipality in enumerate(origin_order)}
            destination_index = {municipality: index for index, municipality in enumerate(destination_order)}
            for origin, destination_list in batch_groups:
                source_index = origin_index[origin]
                for destination in destination_list:
                    target_index = destination_index[destination]
                    dsec = (
                        duration[source_index][target_index]
                        if source_index < len(duration) and target_index < len(duration[source_index])
                        else None
                    )
                    dmet = (
                        distance[source_index][target_index]
                        if source_index < len(distance) and target_index < len(distance[source_index])
                        else None
                    )
                    osnap = (
                        source_meta[source_index].get("distance") if source_index < len(source_meta) else None
                    )
                    dsnap = (
                        destination_meta[target_index].get("distance")
                        if target_index < len(destination_meta)
                        else None
                    )
                    status = "ok" if dsec is not None and dmet is not None else "no_route"
                    new_rows.append(
                        {
                            "origin_muni": origin,
                            "destination_muni": destination,
                            "road_duration_seconds": dsec,
                            "road_distance_m": dmet,
                            "origin_snap_m": osnap,
                            "destination_snap_m": dsnap,
                            "route_status": status,
                            "retrieved_utc": retrieved,
                            "osrm_service": OSRM_BASE,
                        }
                    )

        if batch_number % 10 == 0 or batch_number == total_batches:
            cache = pd.concat([cache, pd.DataFrame(new_rows)], ignore_index=True)
            new_rows = []
            cache = cache.drop_duplicates(key_cols, keep="last")
            atomic_csv(cache[columns], cache_path)
            print(
                f"ROUTING_PROGRESS batch={batch_number}/{total_batches} cache_rows={len(cache)}",
                flush=True,
            )
        time.sleep(sleep_seconds)

    if new_rows:
        cache = pd.concat([cache, pd.DataFrame(new_rows)], ignore_index=True)
        cache = cache.drop_duplicates(key_cols, keep="last")
        atomic_csv(cache[columns], cache_path)
    return cache[columns]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("project")
    parser.add_argument("outdir")
    parser.add_argument("--batch-size", type=int, default=20, help="Maximum origins per OSRM table request")
    parser.add_argument("--sleep-seconds", type=float, default=1.05)
    parser.add_argument("--max-pairs", type=int)
    args = parser.parse_args()

    project = Path(args.project).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    source = project / "derived/analysis_ready_v1.0/stage9/datasus_care_flow_network_geography_v1.0.parquet"
    if not source.exists():
        raise SystemExit(f"Missing source: {source}")

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    population_url = (
        f"{IBGE_BASE}/9514/periodos/2022/variaveis/93"
        "?localidades=N6[all]&classificacao=2[6794]|287[100362]|286[113635]"
    )
    gdp_url = f"{IBGE_BASE}/5938/periodos/2022/variaveis/37?localidades=N6[all]"
    population_raw = outdir / "ibge_9514_population_2022_raw.json"
    gdp_raw = outdir / "ibge_5938_gdp_2022_raw.json"
    population = flatten_ibge(fetch_json(session, population_url, population_raw), "population_2022")
    gdp = flatten_ibge(fetch_json(session, gdp_url, gdp_raw), "gdp_2022_thousand_brl")
    municipality = population.merge(
        gdp[["municipality_code", "gdp_2022_thousand_brl"]],
        on="municipality_code",
        how="outer",
        validate="one_to_one",
    )
    municipality["gdp_per_capita_2022_brl"] = (
        municipality.gdp_2022_thousand_brl * 1000.0 / municipality.population_2022
    )
    # SIH-SUS stores the six-digit municipality code without the final IBGE
    # check digit; preserve both identifiers and join on the documented prefix.
    municipality["municipality_code_sih6"] = municipality.municipality_code.astype(str).str[:6]
    if municipality.municipality_code_sih6.duplicated().any():
        raise RuntimeError("IBGE-to-SIH six-digit municipality mapping is not one-to-one")
    municipality["gdp_pc_quartile"] = safe_qcut(municipality.gdp_per_capita_2022_brl)
    municipality["population_quartile"] = safe_qcut(municipality.population_2022)
    municipality["macroregion"] = municipality.municipality_code.map(macroregion)
    municipality["ibge_population_source"] = population_url
    municipality["ibge_gdp_source"] = gdp_url
    municipality_path = outdir / "ibge_municipality_equity_2022.csv"
    atomic_csv(municipality, municipality_path)

    con = duckdb.connect()
    pairs = load_pairs(con, source)
    pairs_path = outdir / "municipality_pairs_input.csv"
    atomic_csv(pairs, pairs_path)
    cache_path = outdir / "osrm_municipality_pair_cache.csv"
    routes = route_batches(
        session,
        pairs,
        cache_path,
        args.batch_size,
        args.sleep_seconds,
        args.max_pairs,
    )
    merged_routes = pairs.merge(routes, on=["origin_muni", "destination_muni"], how="left", validate="one_to_one")
    merged_routes["road_duration_minutes"] = merged_routes.road_duration_seconds / 60.0
    merged_routes["road_distance_km"] = merged_routes.road_distance_m / 1000.0
    merged_routes["road_to_centroid_distance_ratio"] = (
        merged_routes.road_distance_km / merged_routes.centroid_distance_km.where(merged_routes.centroid_distance_km > 0)
    )
    external = merged_routes.origin_muni != merged_routes.destination_muni
    merged_routes["route_plausibility_flag"] = False
    merged_routes.loc[
        external
        & (
            (merged_routes.route_status != "ok")
            | (merged_routes.road_distance_km < 0.95 * merged_routes.centroid_distance_km)
            | (merged_routes.road_to_centroid_distance_ratio > 10)
            | (merged_routes.road_duration_minutes <= 0)
            | (merged_routes.road_duration_minutes > 72 * 60)
            | (merged_routes.origin_snap_m > 20000)
            | (merged_routes.destination_snap_m > 20000)
        ),
        "route_plausibility_flag",
    ] = True
    route_metrics_path = outdir / "municipality_pair_road_metrics.csv"
    atomic_csv(merged_routes, route_metrics_path)

    output_parquet = outdir / "datasus_care_flow_network_road_ibge_v1.0.parquet"
    source_sql = source.as_posix().replace("'", "''")
    route_sql = route_metrics_path.as_posix().replace("'", "''")
    muni_sql = municipality_path.as_posix().replace("'", "''")
    output_sql = output_parquet.as_posix().replace("'", "''")
    con.execute(
        f"""
        COPY (
          SELECT d.*,
                 r.road_duration_minutes,
                 r.road_distance_km,
                 r.origin_snap_m,
                 r.destination_snap_m,
                 r.route_status,
                 r.road_to_centroid_distance_ratio,
                 r.route_plausibility_flag,
                 m.population_2022 AS origin_population_2022,
                 m.gdp_2022_thousand_brl AS origin_gdp_2022_thousand_brl,
                 m.gdp_per_capita_2022_brl AS origin_gdp_per_capita_2022_brl,
                 m.gdp_pc_quartile AS origin_gdp_pc_quartile,
                 m.population_quartile AS origin_population_quartile,
                 m.macroregion AS origin_macroregion
          FROM read_parquet('{source_sql}') d
          LEFT JOIN read_csv_auto('{route_sql}', header=TRUE, all_varchar=FALSE) r
            ON CAST(d.origin_muni AS VARCHAR)=CAST(r.origin_muni AS VARCHAR)
           AND CAST(d.destination_muni AS VARCHAR)=CAST(r.destination_muni AS VARCHAR)
          LEFT JOIN read_csv_auto('{muni_sql}', header=TRUE, all_varchar=FALSE) m
            ON CAST(d.origin_muni AS VARCHAR)=CAST(m.municipality_code_sih6 AS VARCHAR)
        ) TO '{output_sql}' (FORMAT PARQUET, COMPRESSION ZSTD)
        """
    )

    total_pairs = len(merged_routes)
    external_rows = merged_routes[external]
    qc = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_parquet": str(source),
        "pair_count": total_pairs,
        "same_municipality_pairs": int((~external).sum()),
        "external_pairs": int(external.sum()),
        "external_routes_ok": int((external_rows.route_status == "ok").sum()),
        "external_routes_missing": int((external_rows.route_status != "ok").sum()),
        "plausibility_flags": int(external_rows.route_plausibility_flag.fillna(True).sum()),
        "median_external_road_minutes": float(external_rows.road_duration_minutes.median()),
        "median_external_road_km": float(external_rows.road_distance_km.median()),
        "median_road_to_centroid_ratio": float(external_rows.road_to_centroid_distance_ratio.median()),
        "ibge_population_rows": int(population.population_2022.notna().sum()),
        "ibge_gdp_rows": int(gdp.gdp_2022_thousand_brl.notna().sum()),
        "max_pairs_option": args.max_pairs,
        "road_definition": "OSRM driving route between official IBGE 2022 municipality centroids; cached access-date estimate, not patient-address travel time",
        "interpretation": "geographic access metric; not a causal referral exposure",
    }
    qc_path = outdir / "stage23_road_ibge_qc.json"
    qc_path.write_text(json.dumps(qc, ensure_ascii=False, indent=2), encoding="utf-8")

    print("Road/IBGE preparation complete.", flush=True)


if __name__ == "__main__":
    main()
