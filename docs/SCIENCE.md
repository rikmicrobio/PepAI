# PepAI — Scientific Methodology

## Background

Peptide therapeutics represent one of the fastest-growing segments of the pharmaceutical industry,
combining the specificity of biologics with the chemical tractability of small molecules. However,
the vast sequence space (20^N for a peptide of length N) makes exhaustive screening computationally
infeasible. PepAI addresses this through a **pharmacophore-guided, AI-augmented search strategy**.

---

## Core Concepts

### 1. Pharmacophore Representation

A pharmacophore is a spatial arrangement of chemical features essential for binding. PepAI uses a
4-feature integer vector:

```
φ = (HBD, HBA, Aromatic, Hydrophobic)
```

Where:
- **HBD**: Hydrogen bond donor count (S, T, Y, N, Q, K, R, H residues)
- **HBA**: Hydrogen bond acceptor count (D, E, N, Q, S, T residues)
- **Aromatic**: Aromatic ring count (F, Y, W, H residues)
- **Hydrophobic**: Hydrophobic residue count (A, V, L, I, M, F, W, P residues)

### 2. Tanimoto Similarity

Pharmacophore matching uses the Tanimoto (Jaccard) coefficient:

```
T(A, B) = Σ min(Aᵢ, Bᵢ) / Σ max(Aᵢ, Bᵢ)
```

Range: [0, 1], where 1 = perfect pharmacophore overlap.

### 3. GNN-Enhanced Similarity (Agent 3)

Agent 3 combines Tanimoto with a simulated graph neural network boost:

```
final_score = 0.5 × tanimoto + 0.3 × embed_sim + 0.2 × class_boost
embed_sim   = tanimoto × U(0.85, 1.15) + GNN_boost
GNN_boost   ~ U(0, 0.18)
class_boost  = 0.12  if target_class matches, else 0
```

In production, `GNN_boost` would be replaced with a real learned embedding from
ChemBERTa or Mol2Vec applied to the SMILES string.

### 4. Druggability Scoring (Agent 1)

Agent 1 computes a composite druggability score from ESM-2 embedding axes:

```
D = 0.30 × binding_propensity
  + 0.25 × surface_accessibility
  + 0.25 × evolutionary_conservation
  + 0.20 × (1 − structural_order)
  + ε,   ε ~ N(0, 0.04²)
```

Where each axis is derived from the protein's physicochemical features and the
(simulated) ESM-2 reduced embedding.

### 5. Peptide Generation Strategy (Agent 4)

Agent 4 builds a **biased amino acid pool** using three layers of context:

1. **Pharmacophore context**: Over-represents residue classes matching the query pharmacophore
   (e.g., if HBD > 3, adds 3× copies of HBD residues to the pool)
2. **Charge context**: Adds charged residues if the protein is highly charged (frac_charged > 0.2)
3. **Pocket context**: Augments with residues matching key pocket residues from Agent 2
   (PHE/TRP/TYR → aromatic bias; ARG/LYS → positive charge bias; ASP/GLU → negative charge bias)

Peptide length range is protein-mass-adaptive:
```
len_min = max(5, round(mw_kDa / 8))
len_max = min(20, max(14, round(mw_kDa / 4)))
```

### 6. Docking Scoring (Agent 5)

Agent 5 computes docking confidence as a weighted function:

```
conf = 0.4 × ph_match
     + 0.3 × (1 + docking_score/15)
     + 0.3 × (1 − pose_rmsd/5)
```

Clamped to [0, 1]. In production, `docking_score` would be AutoDock-GPU's
free energy of binding (kcal/mol).

### 7. Bayesian Optimization Loop (Agent 6)

The optimization simulates a RL-guided mutation search:

```
score(iter) = score_best × (1 + δ × iter/T)
δ ~ U(0.02, 0.12) × (1 − 0.4 × iter/T)
```

Mutation strategies modelled:
- **Point mutation**: Single residue substitution
- **Alanine scan**: Systematic Ala replacement to identify critical residues
- **N-term extension**: Prepend 1–2 residues to N-terminus
- **C-term trim**: Remove terminal residue
- **Charge swap**: Swap K↔R, D↔E, or add/remove charged residue
- **Aromatic substitution**: F↔Y↔W↔H swap

---

## Target Class Detection

Agent 1 infers a protein target class from sequence composition:

| Class | Rule |
|-------|------|
| Kinase | aromatic > 0.14 AND hydrophobic > 0.45 |
| Protease | charged_pos > 0.12 AND hba > 0.20 |
| ACE | hbd_donors > 0.35 AND polar > 0.25 |
| GPCR-Opioid | aromatic > 0.10 AND mean_hydro > 0.5 |
| COX | charged_neg > 0.12 AND net_charge < −5 |
| General | fallback |

---

## BindingDB Subset

PepAI includes a curated 40-compound subset from BindingDB spanning 6 target classes:

| Class | Compounds |
|-------|-----------|
| Kinase | Imatinib, Dasatinib, Erlotinib, Gefitinib, Sorafenib, Vemurafenib, Lapatinib, Sunitinib, Nilotinib, Bosutinib |
| COX | Aspirin, Ibuprofen, Naproxen, Celecoxib, Diclofenac |
| ACE | Lisinopril, Captopril, Enalapril, Ramipril, Perindopril |
| Protease | Ritonavir, Lopinavir, Saquinavir, Indinavir, Atazanavir |
| GPCR-Opioid | Morphine, Fentanyl, Naloxone, Buprenorphine, Methadone |
| General | Metformin, Sitagliptin, Saxagliptin, Alogliptin, Vildagliptin, Vancomycin, Daptomycin, Linezolid, Rifampicin, Ciprofloxacin |

---

## Limitations & Caveats

1. **All AI scores are simulated**. Agent outputs are deterministic proxies, not real model predictions.
2. **No 3D structure**. Pocket detection is based on sequence statistics, not crystallographic data.
3. **Simplified ADMET**. Solubility and toxicity are pharmacophore-match proxies, not QSAR models.
4. **Linear peptides only**. Cyclic, stapled, or N-methylated peptides are not yet supported.
5. **No stereochemistry**. All generated peptides are assumed L-amino acids.

---

## References

1. Rives, A. et al. (2021). Biological structure and function emerge from scaling unsupervised
   learning to 250 million protein sequences. *PNAS* 118(15), e2016239118.

2. Jumper, J. et al. (2021). Highly accurate protein structure prediction with AlphaFold.
   *Nature* 596, 583–589.

3. Jiménez, J. et al. (2017). DeepSite: protein-binding site predictor using 3D-convolutional
   neural networks. *Bioinformatics* 33(19), 3036–3042.

4. Dauparas, J. et al. (2022). Robust deep learning–based protein sequence design using
   ProteinMPNN. *Science* 378(6615), 49–56.

5. Trott, O. & Olson, A.J. (2010). AutoDock Vina: improving the speed and accuracy of docking
   with a new scoring function, efficient optimization, and multithreading.
   *Journal of Computational Chemistry* 31(2), 455–461.

6. Balandat, M. et al. (2020). BoTorch: A Framework for Efficient Monte-Carlo Bayesian
   Optimization. *Advances in Neural Information Processing Systems* 33.

7. Wang, Y. et al. (2020). Improving protein–ligand docking using enhanced sampling of
   ligand conformational states. *J. Chemical Information and Modeling* 60(4), 2246–2257.
