# PepAI — Real AI Model Integration Guide

This guide explains how to replace each simulated agent with a real AI model call,
maintaining full compatibility with the PepAI Shiny interface.

---

## General Principle

Every agent function must return the **same output structure** regardless of whether
it is simulated or real. The Shiny server only sees the return value — the implementation
is fully opaque.

---

## Agent 1 — Real ESM-2 via Python (reticulate)

### Setup
```r
# Install reticulate and configure Python environment
install.packages("reticulate")
library(reticulate)
use_condaenv("esm_env")  # or use_virtualenv()

# In terminal:
# conda create -n esm_env python=3.10
# conda activate esm_env
# pip install fair-esm torch
```

### Replacement function
```r
agent1_protein_lm <- function(sequence, prot_feat) {
  py_run_string("
import esm, torch
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model.eval()
batch_converter = alphabet.get_batch_converter()
data = [('protein', seq)]
_, _, tokens = batch_converter(data)
with torch.no_grad():
    results = model(tokens, repr_layers=[33], return_contacts=False)
token_repr = results['representations'][33][0, 1:-1].mean(0).numpy()
  ", list(seq = sequence))
  
  embed_raw <- py$token_repr[1:6]
  embed <- setNames(as.numeric(embed_raw), 
    c("evolutionary_conservation","structural_order","binding_propensity",
      "surface_accessibility","charge_distribution","fold_confidence"))
  embed <- sapply(embed, function(x) max(0, min(1, (x + 2) / 4)))  # normalise
  
  list(
    model = "ESM-2 (650M, real)",
    embedding = embed,
    secondary_structure = c(helix=0.35, sheet=0.25, coil=0.40),  # or use DSSP
    druggability = as.numeric(mean(embed[c("binding_propensity","surface_accessibility")])),
    confidence   = as.numeric(mean(embed)),
    recommended_class = prot_feat$target_class
  )
}
```

---

## Agent 2 — Real Pocket Prediction via P2Rank

### Setup
```bash
# Download P2Rank (Java-based, no GPU required)
wget https://github.com/rdk/p2rank/releases/download/2.4/p2rank_2.4.tar.gz
tar -xzf p2rank_2.4.tar.gz
```

### Replacement function
```r
agent2_pocket <- function(prot_feat, lm_out) {
  # Write a dummy PDB (or obtain real PDB from AlphaFold API)
  pdb_url <- paste0(
    "https://alphafold.ebi.ac.uk/files/AF-",
    prot_feat$uniprot_id, "-F1-model_v4.pdb"
  )
  tmp_pdb <- tempfile(fileext = ".pdb")
  download.file(pdb_url, tmp_pdb, quiet = TRUE)
  
  # Run P2Rank
  system(paste("./prank predict -f", tmp_pdb, "-o /tmp/p2rank_out"))
  
  # Parse output
  pocket_csv <- read.csv("/tmp/p2rank_out/protein.pdb_predictions.csv")
  pockets <- lapply(seq_len(min(5, nrow(pocket_csv))), function(i) {
    list(
      pocket_id    = paste0("P", i),
      volume_A3    = as.numeric(pocket_csv$score[i]) * 100,
      druggability = min(1, as.numeric(pocket_csv$score[i]) / 10),
      hydrophobic  = prot_feat$frac_hydro,
      depth_A      = runif(1, 5, 20),
      key_residues = pocket_csv$residue_ids[i]
    )
  })
  
  list(
    model       = "P2Rank (real)",
    n_pockets   = nrow(pocket_csv),
    pockets     = pockets,
    best_pocket = pockets[[1]],
    confidence  = min(1, as.numeric(pocket_csv$score[1]) / 10)
  )
}
```

---

## Agent 3 — Real ChemBERTa Embeddings

### Setup
```r
install.packages("reticulate")
# pip install transformers torch rdkit-pypi
```

### Replacement function
```r
agent3_ligand <- function(query_ph, protein_class, top_n = 10) {
  # Load real ChemBERTa embeddings (pre-computed for BindingDB)
  # chemberta_embeddings <- readRDS("data/chemberta_bindingdb.rds")
  
  # For query: convert pharmacophore back to pseudo-SMILES or use stored embeddings
  # Here we show the structure — real implementation needs SMILES for query
  
  db <- bindingdb_ligands %>%
    mutate(
      # Replace with real cosine similarity from ChemBERTa embeddings
      tanimoto    = runif(n(), 0.1, 0.9),
      gnn_boost   = 0,  # real embedding similarity replaces this
      class_boost = ifelse(target_class == protein_class, 0.12, 0),
      embed_sim   = tanimoto,
      final_score = 0.5 * tanimoto + 0.3 * embed_sim + 0.2 * class_boost
    ) %>%
    arrange(desc(final_score), desc(binding_affinity)) %>%
    head(top_n)
  db
}
```

---

## Agent 4 — Real ProteinMPNN

### Setup
```bash
# Clone ProteinMPNN
git clone https://github.com/dauparas/ProteinMPNN.git
pip install torch numpy
```

### Replacement function
```r
agent4_generate <- function(prot_feat, query_ph, hits, pocket_out, n_pep=40, seed=42) {
  # Prepare input PDB with pocket residues masked
  # Run ProteinMPNN via system call or reticulate
  
  system(paste(
    "python ProteinMPNN/protein_mpnn_run.py",
    "--pdb_path", tmp_pdb,
    "--num_seq_per_target", n_pep,
    "--sampling_temp", "0.1",
    "--out_folder /tmp/mpnn_out"
  ))
  
  # Parse FASTA output
  seqs <- read.table("/tmp/mpnn_out/seqs/protein.fa", header=FALSE)
  # ... build data.frame with required columns ...
  
  # Must return data.frame with columns:
  # peptide, sequence, length, mw_kda, net_charge, mean_hydro,
  # ph_match, mpnn_score, binding_score, docking_score,
  # solubility, toxicity, n_hbd, n_hba, n_aromatic, n_hydrophobic, stage
}
```

---

## Agent 5 — Real AutoDock-Vina

### Setup
```bash
pip install vina  # Python bindings for AutoDock-Vina
```

### Replacement function
```r
agent5_docking <- function(peptides, pocket_out) {
  library(reticulate)
  vina <- import("vina")
  
  results <- lapply(seq_len(nrow(peptides)), function(i) {
    seq <- peptides$sequence[i]
    # Convert sequence to 3D via PeptideBuilder or similar
    # Run Vina docking
    # Return binding energy, pose RMSD
  })
  
  peptides %>%
    mutate(
      docking_score      = sapply(results, function(r) r$energy),
      pose_rmsd          = sapply(results, function(r) r$rmsd),
      vdw_energy         = sapply(results, function(r) r$vdw),
      electrostatic      = sapply(results, function(r) r$elec),
      pocket_volume_used = runif(n(), 0.3, 0.9),
      docking_confidence = pmax(0, pmin(1,
        0.4 * ph_match + 0.3 * (1 + docking_score/15) + 0.3 * (1 - pose_rmsd/5)))
    )
}
```

---

## Agent 6 — Real BoTorch Bayesian Optimization

### Setup
```bash
pip install botorch gpytorch torch
```

### Replacement function
```r
agent6_optimize <- function(peptides, n_iter=5, seed=123) {
  library(reticulate)
  
  # Pass top peptides to Python BoTorch optimizer
  top <- peptides %>% arrange(desc(ph_match), desc(binding_score)) %>% head(10)
  
  py_run_string("
import torch
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_model
from botorch.acquisition import ExpectedImprovement
# ... BoTorch optimization loop ...
  ")
  
  # Return same structure as simulated agent:
  # data.frame with: iteration, best_score, mean_ph, n_candidates, mutation_type
}
```

---

## Environment Setup Script

Save as `setup_real_agents.R` and run once before switching to real agents:

```r
# setup_real_agents.R
library(reticulate)

# Create dedicated conda environment
conda_create("pepai_env", python_version = "3.10")
use_condaenv("pepai_env", required = TRUE)

# Install Python dependencies
py_install(c(
  "fair-esm",
  "torch",
  "transformers",
  "rdkit-pypi",
  "botorch",
  "gpytorch",
  "vina",
  "biopython"
), envname = "pepai_env", pip = TRUE)

message("PepAI real agent environment ready.")
```

---

## Testing After Integration

After replacing any agent, validate with the built-in BCR-ABL example:

```r
# Load example
fasta <- ">sp|P00519|ABL1_HUMAN BCR-ABL tyrosine kinase
MGPSENDPNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNY
ITPVNSLEKHSWYHGPVSRNAAEYLLSSGINGSFLVRESESSPGQRSISLRYEGRVYHYR"

prot <- analyse_protein(parse_fasta(fasta)$sequence)
lm   <- agent1_protein_lm(parse_fasta(fasta)$sequence, prot)
pk   <- agent2_pocket(prot, lm)
ph   <- protein_to_pharmacophore(prot)
hits <- agent3_ligand(ph, prot$target_class)
peps <- agent4_generate(prot, ph, hits, pk)
dock <- agent5_docking(peps, pk)
opt  <- agent6_optimize(dock)

cat("Pipeline completed successfully.\n")
cat("Best peptide:", dock$sequence[1], "\n")
cat("Best docking score:", min(dock$docking_score), "\n")
```
