# Data Availability and Access Routes

The NRD is available through the Healthcare Cost and Utilization Project under a data-use agreement and cannot be redistributed. SIH-SUS, CNES, IBGE, and OSRM data are available from their source organizations subject to source terms. The public submission package contains no patient-level records.

| Data or material | Repository content | Access route | Limitation |
|---|---|---|---|
| US NRD 2018–2022 | Concept-level codebook, scripts, protocols, and permitted aggregate outputs | [HCUP controlled access](https://hcup-us.ahrq.gov/nrdoverview.jsp) under a data-use agreement | Patient-level NRD records cannot be redistributed |
| Fine-grained NRD hospital-year outputs | Not included; generation code is provided | Regenerate locally after obtaining approved NRD access | Exact small cells and row-level hospital-year outputs are not publicly redistributed |
| Brazil SIH-SUS/RD/SP 2018–2022 | Concept-level codebook, scripts, and aggregate pathway outputs | [DATASUS](https://datasus.saude.gov.br/) subject to source terms | No patient-level extract is included in the repository |
| CNES hospital capability | Code and aggregate-derived results | [CNES](https://cnes.datasus.gov.br/) subject to source terms | No raw linked hospital extract is redistributed |
| IBGE municipality population, GDP, and centroids | Code and aggregate geographic outputs | [IBGE](https://www.ibge.gov.br/) subject to source terms | Source version and retrieval date should be recorded by the analyst |
| OSRM road routing | Routing code and aggregate route summaries | [Project OSRM](https://project-osrm.org/) or a permitted OSRM service | Routes represent municipality-centroid travel, not patient-address travel time |
| Figure source data | Included in `data/source_data/` | Open within this repository | Aggregate or summary records only; no patient-level data |
| Codebooks and protocols | Included in `data/codebooks/` and `protocols/` | Open within this repository | Concept definitions and analysis specifications, not raw records |

The repository intentionally does not claim that restricted analyses are fully executable without approved access to the underlying source data. Public aggregate files are provided to support inspection of the reported figures and summary results.
