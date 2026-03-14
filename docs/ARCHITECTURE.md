# PepAI — Architecture Documentation

## Overview

PepAI is built as a single-file R Shiny application (`app.R`) organised into six logical sections:

```
app.R
 ├── Section A  — FASTA Parsing & Protein Analysis
 ├── Section B  — AI Agent Modules (6 agents)
 ├── Section C  — BindingDB Ligand Database
 ├── Section D  — Agent Graph Definitions (visNetwork)
 ├── Section E  — UI (dashboardPage layout)
 └── Section F  — Server (reactive logic + render functions)
```

---

## Data Flow

```
User Input (FASTA)
        │
        ▼
 ┌─────────────┐
 │  Section A  │  parse_fasta() → validate_fasta() → analyse_protein()
 │  FASTA/AA   │  Outputs: aa_freq, mw_kda, pI, net_charge,
 └──────┬──────┘           frac_hydro, frac_aromatic, target_class
        │
        ▼
 ┌─────────────┐
 │   Agent 1   │  agent1_protein_lm(sequence, prot_feat)
 │   ESM-2     │  Outputs: embedding[6], secondary_structure,
 └──────┬──────┘           druggability, confidence
        │
        ├──────────────────────────┐
        ▼                          ▼
 ┌─────────────┐           ┌─────────────┐
 │   Agent 2   │           │   Agent 3   │
 │   Pocket    │           │   Ligand    │
 │  Detection  │           │    GNN      │
 └──────┬──────┘           └──────┬──────┘
        │  best_pocket             │  top-K ligand hits
        │                          │
        └────────────┬─────────────┘
                     │
                     ▼
             protein_to_pharmacophore()
             query_ph = {hbd, hba, aromatic, hydro}
                     │
                     ▼
             ┌─────────────┐
             │   Agent 4   │  agent4_generate(prot_feat, query_ph,
             │  ProteinMPNN│                  hits, pocket_out, n_pep)
             └──────┬──────┘  Outputs: peptide data.frame (40–150 rows)
                    │         columns: sequence, ph_match, mpnn_score,
                    │                  binding_score, solubility, toxicity...
                    ▼
             ┌─────────────┐
             │   Agent 5   │  agent5_docking(peptides, pocket_out)
             │  AutoDock   │  Adds: pose_rmsd, vdw_energy, electrostatic,
             └──────┬──────┘        pocket_volume_used, docking_confidence
                    │
                    ▼
             ┌─────────────┐
             │   Agent 6   │  agent6_optimize(peptides, n_iter)
             │  Bayesian   │  Outputs: iteration history data.frame
             │  Optimizer  │  columns: best_score, mean_ph, n_candidates,
             └─────────────┘           mutation_type
```

---

## Reactive Graph (Shiny Server)

```
input$run_btn
    │
    ├──► fasta_parsed   ──► fasta_valid
    │         │
    │         ▼
    ├──► prot_feat  ──────────────────────────────────┐
    │         │                                        │
    │         ▼                                        │
    ├──► lm_out (Agent 1)                              │
    │         │                                        │
    │         ▼                                        │
    ├──► pocket_out (Agent 2)                          │
    │         │                                        │
    ├──► query_ph ◄────────────────────────────────────┘
    │         │
    ├──► effective_class
    │         │
    ├──► ligand_hits (Agent 3) ◄──── query_ph + effective_class
    │         │
    ├──► peptides_rv ◄────── Agent 4 + Agent 5
    │         │
    └──► optim_rv ◄───────── Agent 6
              │
              ▼
         pep_filt (reactive filter on peptides_rv)
```

---

## Key Functions Reference

### Section A — FASTA & Protein Analysis

| Function | Inputs | Outputs |
|----------|--------|---------|
| `parse_fasta(txt)` | Raw FASTA string | `list(ok, msg, header, sequence)` |
| `validate_fasta(txt)` | Raw FASTA string | `list(ok, msg)` |
| `analyse_protein(sequence)` | AA sequence string | 13-field list with all protein properties |

### Section B — Agent Functions

| Function | Key Inputs | Key Outputs |
|----------|-----------|-------------|
| `agent1_protein_lm(sequence, prot_feat)` | AA sequence, protein features | embedding, druggability, secondary structure |
| `agent2_pocket(prot_feat, lm_out)` | Protein features, LM output | pocket list, best_pocket, n_pockets |
| `agent3_ligand(query_ph, protein_class, top_n)` | Pharmacophore query, class | Ranked BindingDB hits data.frame |
| `agent4_generate(prot_feat, query_ph, hits, pocket_out, n_pep, seed)` | All upstream outputs | Peptide candidates data.frame |
| `agent5_docking(peptides, pocket_out)` | Peptide data.frame, pocket | Same data.frame + docking columns |
| `agent6_optimize(peptides, n_iter, seed)` | Peptides, iterations | Optimization history data.frame |

### Helper Functions

| Function | Description |
|----------|-------------|
| `protein_to_pharmacophore(pf)` | Converts protein fractions to integer pharmacophore counts |
| `tanimoto_score(q, hbd, hba, arom, hydro)` | Tanimoto similarity between two pharmacophore vectors |

---

## Replacing a Simulated Agent

Every agent returns a well-defined list/data.frame. To swap in a real model, maintain the same output signature:

### Agent 1 → Real ESM-2
```r
agent1_protein_lm <- function(sequence, prot_feat) {
  # Call real ESM-2 API
  result <- reticulate::py_run_string(sprintf("
    import esm
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    # ... encode sequence ...
  "))
  
  # Must return same structure:
  list(
    model         = "ESM-2 (650M, real)",
    embedding     = result$embedding,        # named numeric vector, length 6+
    secondary_structure = result$ss_fractions,
    druggability  = result$druggability,     # numeric 0-1
    confidence    = result$confidence,       # numeric 0-1
    recommended_class = prot_feat$target_class
  )
}
```

### Agent 4 → Real ProteinMPNN
```r
agent4_generate <- function(prot_feat, query_ph, hits, pocket_out, n_pep=40, seed=42) {
  # Call ProteinMPNN API / local installation
  # Must return data.frame with columns:
  # peptide, sequence, length, mw_kda, net_charge, mean_hydro,
  # ph_match, mpnn_score, binding_score, docking_score,
  # solubility, toxicity, n_hbd, n_hba, n_aromatic, n_hydrophobic, stage
}
```

---

## UI Architecture

```
dashboardPage
 ├── dashboardHeader      — PepAI logo + wordmark
 ├── dashboardSidebar     — 8 menuItems
 └── dashboardBody
      ├── CSS styles (inline, ~80 lines)
      └── tabItems
           ├── target      — FASTA input + protein stats
           ├── agents      — visNetwork graph + 6 agent cards
           ├── ligands     — BindingDB table + similarity plots
           ├── peptides    — Filtered table + best candidate card
           ├── optim       — Trajectory + iteration log + mini-plots
           ├── analytics   — 5 analytical plots
           ├── about       — Scientific framework + agent table
           └── developers  — Team cards + InteractIQ panel
```

---

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `shiny` | ≥1.7 | Web application framework |
| `shinydashboard` | ≥0.7 | Dashboard layout components |
| `visNetwork` | ≥2.1 | Interactive agent/pipeline graphs |
| `ggplot2` | ≥3.4 | Static plot generation |
| `plotly` | ≥4.10 | Interactive plot conversion |
| `DT` | ≥0.28 | Interactive data tables |
| `dplyr` | ≥1.1 | Data manipulation |
| `stringr` | ≥1.5 | String operations (FASTA parsing) |
| `scales` | ≥1.2 | Axis formatting |
| `htmlwidgets` | ≥1.6 | JS `onRender` callbacks |
