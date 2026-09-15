# GSMN_chlamy — iYH2021

Genome-scale metabolic model (GSMN) of the light-tolerant natural isolate
***Chlamydomonas reinhardtii* CC-4414**, reconstructed from a *de novo*
genome assembly and integrated with high-light (HL) and low-light (LL)
quantitative proteomics via the iMAT algorithm.

This repository accompanies the manuscript:

> **Integrative Proteomics and Genome-Scale Modeling Elucidate Metabolic
> Flux Shifts in *Chlamydomonas reinhardtii* CC-4414 under Light Stress**
> International Journal of Molecular Sciences (manuscript ID: ijms-4546231)

It contains the model files, the full reconstruction and analysis pipeline,
all intermediate and input data, and the scripts needed to reproduce every
computational result reported in the manuscript.

---

## 1. Repository structure

```
GSMN_chlamy/
├── model/
│   ├── iYH2021.mat                  # Final model, MATLAB/COBRA format
│   ├── iYH2021.xlsx                 # Final model, RAVEN-style spreadsheet export
│   ├── model_An.mat                 # Intermediate model (post-BBH orthology
│   │                                 #   replacement, pre-gap-filling) — see Section 4
│   └── draft_gapfilled_from_carveme.xml  # Post-gap-filling draft (SBML)
│
├── ComplementaryData/
│   ├── iCre1355_auto.xml             # Template model (photoautotrophic variant)
│   ├── iCre1355_hetero.xml           # Template model (heterotrophic variant)
│   ├── iCre1355_mixo.xml             # Template model (mixotrophic variant)
│   ├── iCre1355_Reactions_GPR.xlsx   # Template reaction/GPR table
│   ├── gsmn_vs_pro.tsv               # BLASTP: CC-4414 genome vs. JGI v5.5 (orthology)
│   ├── raw_pro.xlsx                  # Raw proteomics abundance (UniProt-indexed)
│   ├── BBH_mapping.xlsx              # BBH mapping, proteomics → CC-4414 genes
│   ├── GeneExpression.xlsx           # Gene-level HL/LL abundance (post-BBH mapping)
│   ├── gapfill_added_reactions.csv   # 337 reactions added during gap-filling
│   └── updated_4414.xlsx             # Updated GPR/metabolite/EC/subsystem tables
│
├── annotation/ , bbh/                # Supporting genome annotation and BLAST outputs
│
├── net_construct.m                   # Step 1: template-based reconstruction + gap-filling
├── proteomics_gsmn_EN.m              # Step 2: iMAT integration (English-translated,
│                                      #   annotated, with corrected essentialRxns list)
├── FVA_HL_LL_addon.m                 # Step 3: flux variability analysis (targeted +
│                                      #   full-network options)
├── biomass_objective_HL_LL.m         # Step 4: supplementary theoretical growth-rate
│                                      #   analysis (biomass objective, light-anchored)
├── supplementary_analysis.m          # Step 5: exchange-reaction sensitivity screen +
│                                      #   acetate/CO2 2-D sensitivity analysis
├── proteomics_gsmn_EN_with_FVA.m     # Steps 2+3 combined, single ready-to-run file
├── proteomics_gsmn_EN_with_FVA_and_biomass.m  # Steps 2+3+4 combined
│
├── HL_vs_LL_flux_diff.xlsx           # Output: full point-solution flux comparison
├── HL_vs_LL_flux_diff_filtered.xlsx  # Output: same, transport reactions excluded
├── FVA_HL_vs_LL_targeted.xlsx        # Output: FVA on the 21 reactions discussed in
│                                      #   the manuscript (Table 4)
├── FVA_HL_vs_LL_full.xlsx            # Output: FVA on all 210 shared HL/LL reactions
├── Acetate_CO2_growth_sensitivity.xlsx  # Output: 2-D acetate/CO2 sensitivity matrix
│
└── README.md
```

---

## 2. Software requirements

| Software | Version used | Notes |
|---|---|---|
| MATLAB | R2023a or later (The MathWorks, Inc.) | |
| COBRA Toolbox | v3.x | `initCobraToolbox` must be run first |
| Gurobi Optimizer | v13.0.1 (Gurobi Optimization, LLC) | Set as the active COBRA solver: `changeCobraSolver('gurobi','all')` |
| RAVEN Toolbox | v2.0 | Used during gap-filling (Step 1) |
| CarveMe | v1.6.6 | Used to generate the candidate gap-filling network (Step 1); run outside MATLAB |

A free academic Gurobi license is required. GLPK can be substituted for
exploratory use but was **not** used to generate the results in the
manuscript, and — as we found during peer review — different solvers can
select different optimal vertices for the underlying MILP/LP problems, so
exact reproduction of point-solution values requires Gurobi specifically.

---

## 3. Reconstruction and analysis pipeline

Run the scripts in the following order. Each step's outputs feed into the
next.

### Step 1 — Template-based reconstruction and gap-filling (`net_construct.m`)

- Loads the **iCre1355** template model ([Imam et al. 2015]; auto/mixo/hetero
  variants provided in `ComplementaryData/`).
- Identifies orthologous genes between CC-4414 and the JGI v5.5 reference
  genome by **best bidirectional hits (BBH)** via BLASTP
  (`gsmn_vs_pro.tsv`), using a **sequence identity ≥ 40% and E-value
  ≤ 1×10⁻¹⁰** cutoff (bit score is used only as a secondary tie-breaker,
  not the primary filter).
- Replaces template gene IDs with CC-4414 gene IDs wherever a BBH ortholog
  is found; genes without an ortholog (primarily organellar genes not
  called by the nuclear-genome-focused AUGUSTUS pipeline) retain their
  original iCre1355 ID under an `OLD_iCre_` prefix.
- Gap-fills the network using a candidate reaction set generated by
  **CarveMe v1.6.6** from the CC-4414 genome; candidate reactions are
  manually screened against reaction equation, EC number, gene
  associations, subsystem, and consistency with known *C. reinhardtii*
  metabolism before being added.
- **Output:** 337 reactions added, 3 reactions removed relative to
  iCre1355 (`ATPM_NGAM`, `INSH`, `TAT`; the intermediate model after
  orthology replacement but before gap-filling is saved as
  `model_An.mat`, provided for provenance/audit purposes).
- The biomass objective function (`Biomass_Chlamy_mixo`, mixotrophic
  variant, consistent with TAP-medium growth) is inherited **unchanged**
  from iCre1355, since strain-specific biomass composition data are not
  currently available for CC-4414.

### Step 2 — Proteomics integration via iMAT (`proteomics_gsmn_EN.m`)

- Proteomics data: ProteomeXchange PXD040080 / JPOST JPST002038
  (Suwannachuen et al. 2023), covering CC-4414 grown in TAP medium at
  25 °C, ambient CO₂, continuous illumination, liquid shake culture
  (180 rpm); HL = 1500 µmol photons m⁻² s⁻¹ (sampled 2 days after the
  shift to high light), LL = 50 µmol photons m⁻² s⁻¹; n = 3 biological
  replicates pooled prior to LC-MS/MS, n = 3 technical injections.
- Protein abundance is mapped to reactions via GPR (maximum abundance
  among subunits for multi-gene reactions — see the caveat on this
  convention in the manuscript, Section 3/4.6), log₂-transformed
  (`log2(x+1)`), and classified into high/low expression using the
  **25th/75th percentile** of each condition's distribution
  (following Zur et al., the original iMAT reference):
  - **HL:** 25th percentile = 17.12, 75th percentile = 19.88 (log₂ units)
  - **LL:** 25th percentile = 16.97, 75th percentile = 18.53 (log₂ units)
  - (computed from the *n* = 172 proteins with proteomic support;
    reactions lacking support are assigned the midpoint value)
- Central pathways (photosynthesis, Calvin–Benson cycle, glycolysis, TCA
  cycle) and the biomass/photon-handling reactions are protected with
  wide bounds (−1000 to 1000 mmol gDW⁻¹ h⁻¹) before iMAT is run, using the
  **corrected** `essentialRxns` reaction-ID list (the originally coded
  list contained 7 IDs — `PSII`, `ATPS`, `RBPC`, `PRK`, `SUCDH`, `FBA`,
  `FBP` — that do not exist under those exact names in iYH2021 and were
  silently unprotected; this has been corrected here).
- Light uptake (`EX_photonVis_e`) is constrained to **−1000** for the HL
  model and **−50** for the LL model. **Important:** these are model flux
  units, not a validated physical unit conversion — see the manuscript
  (Section 4.6) for the full explanation of how these values were chosen
  and why they should not be read as direct µmol-photon equivalents.
- FBA is performed on each context-specific model using
  **`PRISM_design_growth`** (a photon spectral-decomposition reaction) as
  the objective, **not** a biomass reaction — because the context-specific
  networks retain only ~11–14% of the full reaction set and cannot
  support biomass synthesis. Results reflect photon-handling capacity and
  flux redistribution, not growth.
- **Output:** `HL_vs_LL_flux_diff.xlsx` (all reactions),
  `HL_vs_LL_flux_diff_filtered.xlsx` (transport reactions excluded).

### Step 3 — Flux variability analysis (`FVA_HL_LL_addon.m`)

- Runs `fluxVariability()` (COBRA Toolbox, Gurobi, 100% of optimum) on the
  HL and LL context-specific models from Step 2.
- **Option A (default):** the 21 reactions discussed in the manuscript
  (Table 4) → `FVA_HL_vs_LL_targeted.xlsx`.
- **Option B (commented out; uncomment to run):** all 210 reactions
  shared between the HL and LL models → `FVA_HL_vs_LL_full.xlsx`.
- Console output reports the number of fully unconstrained reactions and
  the number of robust (non-overlapping) reactions.

### Step 4 — Supplementary theoretical growth-rate analysis (`biomass_objective_HL_LL.m`)

- Because the context-specific models cannot support biomass synthesis
  (Step 2), this step uses the **unrestricted** genome-scale network with
  `Biomass_Chlamy_mixo` as the objective, and anchors light availability
  by setting `EX_photonVis_e` to an upper bound of 1000 (HL) or 50 (LL).
- **Output:** growth rates (HL = 0.151589 h⁻¹, LL = 0.121333 h⁻¹) and FVA
  under this objective → `FVA_HL_vs_LL_biomass_objective.xlsx`.
- This is a distinct analysis from Step 2 and should not be conflated
  with it (see manuscript Section 3 for discussion).

### Step 5 — Exchange-reaction and acetate/CO₂ sensitivity screen (`supplementary_analysis.m`)

- Identifies which exchange reactions most strongly limit growth on the
  unrestricted network by individually doubling each uptake bound.
- Finds that predicted growth **plateaus once photon uptake exceeds
  approximately −67** (model flux units), and that **acetate (`EX_ac_e`)
  and CO₂ (`EX_co2_e`)**, not light, become the binding constraints beyond
  that point.
- Runs a 2-D sensitivity analysis over acetate and CO₂ uptake bounds
  (−2 to −4.5 mmol gDW⁻¹ h⁻¹ each, photon uptake fixed at −80) →
  `Acetate_CO2_growth_sensitivity.xlsx`.
- **Note:** the default acetate/CO₂ bounds (−2/−2) are inherited from the
  iCre1355 template and are **not** experimentally measured, strain-specific
  uptake rates for CC-4414; absolute growth-rate values from Steps 4–5
  should be treated as illustrative rather than validated predictions.

---

## 4. Reproducing the published results

```matlab
% 0. Setup
initCobraToolbox
changeCobraSolver('gurobi', 'all');

% 1. Load the base model (already reconstructed; see net_construct.m
%    if rebuilding from scratch)
model = readCbModel('model/iYH2021.mat');

% 2. Proteomics integration + iMAT (produces model_HL, model_LL)
run('proteomics_gsmn_EN.m')

% 3. Flux variability analysis (Table 4 in the manuscript)
run('FVA_HL_LL_addon.m')

% 4. Supplementary growth-rate analysis (Section 2.4 in the manuscript)
run('biomass_objective_HL_LL.m')

% 5. Exchange-reaction / acetate-CO2 sensitivity screen
run('supplementary_analysis.m')
```

Alternatively, `proteomics_gsmn_EN_with_FVA.m` and
`proteomics_gsmn_EN_with_FVA_and_biomass.m` combine the steps above into
single ready-to-run scripts.

---

## 5. Data availability

| Data | Source |
|---|---|
| Whole-genome sequencing (CC-4414) | NCBI SRA, accession SRX4395102 (Run: SRR7526648) |
| Proteomics (HL/LL) | ProteomeXchange PXD040080; JPOST JPST002038 |
| Source proteomics publication | Suwannachuen et al. (2023), *Int. J. Mol. Sci.* 24(9):8374 |
| iCre1355 template model | Imam et al. (2015), *PLOS Comput. Biol.* |

**Note on genome assembly accession:** SRX4395102/SRR7526648 are SRA
accessions for the raw sequencing reads only. No GenBank assembly
(`GCA_...`) accession currently exists for the assembled CC-4414 genome;
if a citable assembly accession is required, the assembled contigs should
be submitted separately as a new WGS/Genome submission.

---

## 6. Known limitations and caveats (see manuscript for full discussion)

- The context-specific (iMAT) models retain only a minority of the full
  network and cannot be used to assess growth/biomass; growth is
  addressed separately via the unrestricted-network analysis (Step 4).
- Flux variability analysis shows that most individual flux differences
  between HL and LL are **not** uniquely determined by the model and
  proteomics data; only a small subset (notably `PGK`, plus the
  objective-linked photon reactions) are robust to solution
  non-uniqueness. Non-robust differences should be treated as
  model-generated hypotheses, not confirmed findings.
- The maximum-abundance convention used for multi-subunit complexes may
  overestimate effective complex activity in some cases.
- Model photon-exchange flux units are **not** a validated 1:1 conversion
  of physical irradiance (µmol photons m⁻² s⁻¹).

---

## 7. Citation

If you use this model, data, or code, please cite:

> [Authors]. Integrative Proteomics and Genome-Scale Modeling Elucidate
> Metabolic Flux Shifts in *Chlamydomonas reinhardtii* CC-4414 under Light
> Stress. *International Journal of Molecular Sciences*, 2026.
> (full citation to be updated upon publication)

and this repository via its Zenodo DOI: `[DOI to be added upon release]`

---

## 8. License

[Specify license, e.g., MIT / CC-BY 4.0 — to be confirmed by the authors]

## 9. Contact

Pramote Chumnanpuen (pramote.c@ku.th) and Wanwipa Vongsangnak
(wanwipa.v@ku.ac.th), Kasetsart University.
