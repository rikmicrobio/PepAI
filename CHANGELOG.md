# Changelog

All notable changes to PepAI are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [1.0.0] — 2024

### Added
- Full 6-agent autonomous peptide discovery pipeline
- Agent 1: ESM-2 proxy — sequence embeddings, druggability scoring, secondary structure prediction
- Agent 2: AlphaFold2/DeepSite proxy — binding pocket detection, volume, key residues
- Agent 3: GNN/ChemBERTa proxy — BindingDB search with Tanimoto + GNN similarity
- Agent 4: ProteinMPNN proxy — pharmacophore-biased, pocket-anchored peptide generation
- Agent 5: AutoDock-GPU proxy — docking scoring, pose RMSD, energy breakdown
- Agent 6: Bayesian RL optimizer — iterative sequence refinement with 6 mutation strategies
- Interactive R Shiny dashboard with 8 tabs
- visNetwork agent architecture graph with hover tooltips
- Interactive peptide table with 7 real-time filter controls
- Best candidate residue-coloured breakdown view
- Bayesian optimization trajectory plot (native plotly — no RGB conversion bugs)
- BindingDB 40-compound curated subset (6 target classes)
- Pharmacophore match vs binding score, MPNN vs docking, solubility/toxicity analytics
- Amino acid frequency plots (input protein + top-10 peptides)
- Developers page with InteractIQ team contacts and LinkedIn profiles
- Auto-install of all required R packages on first launch
- BCR-ABL tyrosine kinase built-in example (UniProt P00519)
- MIT License

### Renamed
- Rebranded from PepDez-Agent → **PepAI**

---

## [Unreleased]

### Planned
- Real ESM-2 API integration
- SMILES / SDF export of peptide candidates
- CSV download buttons for all result tables
- Docker image
- REST API endpoint
- AlphaFold2 real pocket prediction
- ProteinMPNN real generation
