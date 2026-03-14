
# PepAI
### Pharmacophore-Guided Peptide Discovery Framework

![R](https://img.shields.io/badge/R-%3E%3D4.0-276DC3?logo=r)
![Shiny](https://img.shields.io/badge/Shiny-Interactive-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Research%20Prototype-orange)

---

## Overview

**PepAI** is a modular computational framework for **peptide drug discovery from protein sequence data**.

The platform integrates multiple analytical stages into a single pipeline that transforms a **target protein FASTA sequence** into **ranked candidate therapeutic peptides**.

PepAI demonstrates how a **multi-agent architecture** can coordinate tasks such as:

- protein sequence analysis  
- binding pocket inference  
- ligand knowledge mining  
- pharmacophore-guided peptide generation  
- docking-like interaction scoring  
- iterative candidate optimization  

The system is implemented as an **interactive R Shiny application**, allowing users to explore the full discovery workflow through a visual dashboard.

---

## Important Note

PepAI currently implements **computational proxy modules** that simulate the functional roles of modern AI systems.

The architecture is designed so that each module can be replaced with production-level models such as:

- ESM-2 protein language models  
- AlphaFold pocket prediction  
- ProteinMPNN sequence design  
- AutoDock molecular docking  

The current version focuses on **demonstrating the architecture of an autonomous peptide discovery pipeline** rather than running these models directly.

---

## Pipeline Overview

FASTA Input → Sequence Feature Analysis → Pocket Inference → Ligand Similarity Search → Pharmacophore Extraction → Peptide Generation → Docking-like Interaction Scoring → Optimization Loop → Ranked Peptide Candidates

---

## The 6-Agent Architecture

| Agent | Role | Methodology | Status |
|------|------|------|------|
| Agent 1 | Sequence Encoder | Physicochemical descriptors + statistical feature embeddings | Implemented |
| Agent 2 | Pocket Inference | Sequence-derived pocket descriptors | Implemented |
| Agent 3 | Ligand Knowledge Mining | BindingDB subset + pharmacophore similarity (Tanimoto) | Implemented |
| Agent 4 | Peptide Generator | Pharmacophore-guided residue sampling | Implemented |
| Agent 5 | Docking Evaluator | Heuristic docking-like interaction scoring | Implemented |
| Agent 6 | Optimization Loop | Iterative mutation and score-based refinement | Implemented |

---

## Key Features

- End-to-end **FASTA → peptide candidate** pipeline  
- Interactive **Shiny dashboard interface**  
- Visualization of **agent architecture and workflow**  
- Ligand similarity search using **BindingDB subset**  
- Pharmacophore-guided peptide design  
- Docking-like interaction scoring  
- Iterative candidate optimization  
- Real-time analytics and visualization tools  

---

## Installation

### Requirements

- R ≥ 4.0  
- RStudio (recommended)

### Run the Application

```r
# Clone repository
git clone https://github.com/InteractIQ/PepAI.git
cd PepAI

# Run
shiny::runApp("app.R")
```

Required packages will install automatically on first launch.

---

## Dependencies

```
shiny
shinydashboard
visNetwork
ggplot2
plotly
DT
dplyr
stringr
scales
htmlwidgets
```

---

## Scientific Methodology

### Pharmacophore Representation

Peptide candidates are represented using a simplified pharmacophore vector:

φ = (HBD, HBA, Aromatic, Hydrophobic)

Similarity between vectors is calculated using the **Tanimoto coefficient**.

### Optimization Loop

The final stage iteratively refines candidate peptides through mutation strategies including:

- residue substitution  
- charge modification  
- aromatic enrichment  
- terminal extension  

Each iteration evaluates candidate performance using docking-like scoring.

---

## Project Structure

```
PepAI
│
├── app.R
├── www/
│   └── pepai_logo.png
│
├── R/
│   ├── agents.R
│   ├── fasta_utils.R
│   └── bindingdb.R
│
├── docs/
│   ├── ARCHITECTURE.md
│   ├── SCIENCE.md
│   └── API_INTEGRATION.md
│
├── tests/
│   └── test_agents.R
│
└── README.md
```

---

## Roadmap

### Current Version
✔ modular multi-agent pipeline  
✔ Shiny interactive interface  
✔ pharmacophore-guided peptide generation  

### Planned Features

- real **ESM-2 embedding integration**
- **AlphaFold pocket prediction**
- **ProteinMPNN peptide generation**
- **AutoDock docking integration**
- molecular dynamics pre-screening

---

## Team

PepAI is developed by the **InteractIQ** team.

- Rik Ganguly — Bioinformatics Lead  
- Gautam Ahuja — Computational Scientist  
- Siddhant Poudyal — Research Scientist  

---

## License

MIT License

---

## Citation

```
Ganguly R., Ahuja G., Poudyal S.
PepAI: Pharmacophore-Guided Peptide Discovery Framework.
GitHub repository.
```

---

PepAI is a **research prototype** demonstrating a modular AI-inspired pipeline for peptide drug discovery.
