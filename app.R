# ============================================================
#  PepAI — AI-Augmented Peptide Discovery Pipeline
#  Architecture: 6-Agent autonomous drug design framework
#    Agent 1 → Protein Language Model (ESM-2 proxy)
#    Agent 2 → Pocket Prediction (AlphaFold/DeepSite proxy)
#    Agent 3 → Ligand Knowledge (GNN/ChemBERTa proxy)
#    Agent 4 → Peptide Generator (ProteinMPNN proxy)
#    Agent 5 → Docking Evaluator (AutoDock-GPU proxy)
#    Agent 6 → Bayesian Optimizer (RL feedback loop)
#  Developers: Rik Ganguly, Gautam Ahuja, Siddhant Poudel
#  Contact: rikgangulybioinfo@gmail.com
# ============================================================

required_pkgs <- c("shiny","shinydashboard","visNetwork",
                   "ggplot2","DT","dplyr","plotly","stringr","scales","htmlwidgets")
missing_pkgs  <- required_pkgs[!required_pkgs %in% installed.packages()[,"Package"]]
if (length(missing_pkgs) > 0) {
  message("Installing: ", paste(missing_pkgs, collapse=", "))
  install.packages(missing_pkgs, dependencies=TRUE)
}
library(shiny); library(shinydashboard); library(visNetwork)
library(ggplot2); library(DT); library(dplyr); library(plotly)
library(stringr); library(scales); library(htmlwidgets)

# Serve the www/ folder (contains pepai_logo.png)
if (dir.exists("www")) addResourcePath("www", "www")

# ============================================================
# A — FASTA PARSING & PROTEIN ANALYSIS
# ============================================================
parse_fasta <- function(txt) {
  lines  <- str_split(str_trim(txt), "\n")[[1]]
  lines  <- str_trim(lines)
  hi     <- which(str_starts(lines, ">"))
  if (length(hi) == 0)
    return(list(ok=FALSE, msg="No FASTA header (>) found.", header="", sequence=""))
  header    <- str_remove(lines[hi[1]], "^>")
  seq_lines <- lines[(hi[1]+1):ifelse(length(hi)>1, hi[2]-1, length(lines))]
  sequence  <- paste(str_remove_all(toupper(seq_lines), "[^ACDEFGHIKLMNPQRSTVWY]"), collapse="")
  if (nchar(sequence) < 10)
    return(list(ok=FALSE, msg="Sequence too short (< 10 AA).", header=header, sequence=sequence))
  list(ok=TRUE, msg=paste0("Parsed OK — ", nchar(sequence), " amino acids."),
       header=header, sequence=sequence)
}

validate_fasta <- function(txt) {
  txt <- str_trim(txt)
  if (nchar(txt) < 5)   return(list(ok=FALSE, msg="Input empty or too short."))
  if (!str_detect(txt,">")) return(list(ok=FALSE, msg="Missing FASTA header. Start with '>'."))
  p <- parse_fasta(txt)
  if (!p$ok) return(list(ok=FALSE, msg=p$msg))
  bad <- str_remove_all(p$sequence, "[ACDEFGHIKLMNPQRSTVWY]")
  if (nchar(bad)>0) return(list(ok=FALSE, msg=paste0("Unknown characters: ", bad)))
  list(ok=TRUE, msg=p$msg)
}

AA_MW    <- c(A=89,R=174,N=132,D=133,C=121,E=147,Q=146,G=75,H=155,I=131,
              L=131,K=146,M=149,F=165,P=115,S=105,T=119,W=204,Y=181,V=117)
AA_HYDRO <- c(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,E=-3.5,Q=-3.5,G=-0.4,
              H=-3.2,I=4.5,L=3.8,K=-3.9,M=1.9,F=2.8,P=-1.6,S=-0.8,
              T=-0.7,W=-0.9,Y=-1.3,V=4.2)
AA_CHARGE<- c(A=0,R=1,N=0,D=-1,C=0,E=-1,Q=0,G=0,H=0.1,I=0,
              L=0,K=1,M=0,F=0,P=0,S=0,T=0,W=0,Y=0,V=0)

analyse_protein <- function(sequence) {
  aas  <- strsplit(sequence,"")[[1]]
  n    <- length(aas)
  freq <- table(factor(aas, levels=names(AA_MW))) / n

  mw_total   <- sum(sapply(aas, function(a) AA_MW[a]),    na.rm=TRUE)
  mean_hydro <- mean(sapply(aas, function(a) AA_HYDRO[a]),na.rm=TRUE)
  net_charge <- sum(sapply(aas, function(a) AA_CHARGE[a]),na.rm=TRUE)

  hydrophobic  <- sum(freq[c("A","V","L","I","M","F","W","P")], na.rm=TRUE)
  polar        <- sum(freq[c("S","T","Y","N","Q")],             na.rm=TRUE)
  charged_pos  <- sum(freq[c("R","K","H")],                     na.rm=TRUE)
  charged_neg  <- sum(freq[c("D","E")],                         na.rm=TRUE)
  aromatic     <- sum(freq[c("F","Y","W","H")],                 na.rm=TRUE)
  hbd_donors   <- sum(freq[c("S","T","Y","N","Q","K","R","H")], na.rm=TRUE)
  hba_accept   <- sum(freq[c("D","E","N","Q","S","T")],         na.rm=TRUE)
  pI           <- 7 + net_charge * 0.5

  target_class <- dplyr::case_when(
    aromatic>0.14 & hydrophobic>0.45 ~ "Kinase",
    charged_pos>0.12 & hba_accept>0.20 ~ "Protease",
    hbd_donors>0.35 & polar>0.25  ~ "ACE",
    aromatic>0.10 & mean_hydro>0.5 ~ "GPCR-Opioid",
    charged_neg>0.12 & net_charge < -5 ~ "COX",
    TRUE ~ "General"
  )

  list(n_aa=n, mw_kda=round(mw_total/1000,2),
       mean_hydro=round(mean_hydro,3), net_charge=round(net_charge,1),
       pI=round(pI,1), frac_hydro=round(hydrophobic,3),
       frac_polar=round(polar,3), frac_charged=round(charged_pos+charged_neg,3),
       frac_aromatic=round(aromatic,3), hbd_donors=round(hbd_donors,3),
       hba_acceptors=round(hba_accept,3), target_class=target_class,
       aa_freq=as.data.frame(freq))
}

# ============================================================
# B — AI AGENT MODULES (simulated)
# ============================================================

# Agent 1: Protein Language Model (ESM-2 proxy)
# Generates sequence embeddings + confidence scores
agent1_protein_lm <- function(sequence, prot_feat) {
  set.seed(sum(utf8ToInt(substring(sequence,1,min(20,nchar(sequence))))))
  n <- nchar(sequence)

  # Simulated ESM-2 embedding (32-dim reduced to 6 interpretable axes)
  embed <- c(
    evolutionary_conservation = round(runif(1,0.55,0.98),3),
    structural_order          = round(1 - prot_feat$frac_hydro * 0.4 + runif(1,-0.1,0.1),3),
    binding_propensity        = round(prot_feat$frac_aromatic * 2 + runif(1,0,0.3),3),
    surface_accessibility     = round(prot_feat$frac_polar + runif(1,0,0.2),3),
    charge_distribution       = round(abs(prot_feat$net_charge)/10 + runif(1,0,0.1),3),
    fold_confidence           = round(runif(1,0.60,0.97),3)
  )
  embed <- sapply(embed, function(x) max(0, min(1, x)))

  # Predicted secondary structure fractions
  ss <- c(
    helix  = round(runif(1, 0.25, 0.55), 2),
    sheet  = round(runif(1, 0.10, 0.35), 2),
    coil   = round(runif(1, 0.15, 0.40), 2)
  )
  ss <- round(ss / sum(ss), 2)

  # Druggability score (0–1)
  druggability <- round(
    0.3 * embed["binding_propensity"] +
    0.25 * embed["surface_accessibility"] +
    0.25 * embed["evolutionary_conservation"] +
    0.2  * (1 - embed["structural_order"]) +
    rnorm(1, 0, 0.04), 3)
  druggability <- max(0, min(1, druggability))

  list(
    model         = "ESM-2 (650M, simulated)",
    embedding     = embed,
    secondary_structure = ss,
    druggability  = druggability,
    confidence    = round(mean(embed), 3),
    recommended_class = prot_feat$target_class
  )
}

# Agent 2: Pocket Prediction (AlphaFold2 / DeepSite proxy)
agent2_pocket <- function(prot_feat, lm_out) {
  set.seed(42)
  n_pockets <- sample(2:5, 1)

  pockets <- lapply(seq_len(n_pockets), function(i) {
    volume      <- round(runif(1, 200, 900), 0)
    druggability<- round(lm_out$druggability * runif(1, 0.7, 1.2), 3)
    druggability<- max(0, min(1, druggability))
    hydrophobic_ratio <- round(prot_feat$frac_hydro + runif(1,-0.1,0.15), 3)
    list(
      pocket_id     = paste0("P", i),
      volume_A3     = volume,
      druggability  = druggability,
      hydrophobic   = max(0, min(1, hydrophobic_ratio)),
      depth_A       = round(runif(1, 5, 20), 1),
      key_residues  = paste(sample(c("ARG","LYS","ASP","GLU","PHE","TRP","TYR","HIS",
                                     "LEU","ILE","VAL"), sample(3:6,1)), collapse=", ")
    )
  })

  # Rank by druggability
  pockets <- pockets[order(sapply(pockets, function(p) p$druggability), decreasing=TRUE)]
  best    <- pockets[[1]]

  list(
    model        = "DeepSite / AlphaFold2 (simulated)",
    n_pockets    = n_pockets,
    pockets      = pockets,
    best_pocket  = best,
    confidence   = round(best$druggability * 0.9 + runif(1,0,0.1), 3)
  )
}

# Agent 3: Ligand Knowledge Agent (GNN / ChemBERTa proxy)
agent3_ligand <- function(query_ph, protein_class, top_n=10) {
  # Tanimoto + learned embedding boost (simulated GNN score)
  db <- bindingdb_ligands %>%
    rowwise() %>%
    mutate(
      tanimoto   = {
        qv <- c(query_ph$hbd, query_ph$hba, query_ph$aromatic, query_ph$hydro)
        rv <- c(ph_hbd, ph_hba, ph_arom, ph_hydro)
        inter <- sum(pmin(qv,rv)); union <- sum(pmax(qv,rv))
        if (union==0) 0 else round(inter/union, 4)
      },
      gnn_boost  = round(runif(1, 0, 0.18), 4),   # simulated learned embedding score
      class_boost= ifelse(target_class==protein_class, 0.12, 0),
      embed_sim  = round(tanimoto * runif(1, 0.85, 1.15) + gnn_boost, 4)
    ) %>%
    ungroup() %>%
    mutate(
      embed_sim = pmin(embed_sim, 1),
      final_score = round(0.5*tanimoto + 0.3*embed_sim + 0.2*class_boost, 4)
    ) %>%
    arrange(desc(final_score), desc(binding_affinity)) %>%
    head(top_n)
  db
}

# Agent 4: Peptide Generator (ProteinMPNN / diffusion proxy)
AA_ROLES <- list(
  hydrophobic  = c("A","V","L","I","M","F","W","P"),
  hbd_donor    = c("S","T","Y","N","Q","K","R","H"),
  hba_acceptor = c("D","E","N","Q","S","T"),
  aromatic     = c("F","Y","W","H"),
  charged_pos  = c("K","R","H"),
  charged_neg  = c("D","E"),
  all          = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
)

protein_to_pharmacophore <- function(pf) {
  list(hbd=round(pf$hbd_donors*12), hba=round(pf$hba_acceptors*14),
       aromatic=round(pf$frac_aromatic*10), hydro=round(pf$frac_hydro*18))
}

tanimoto_score <- function(q, hbd, hba, arom, hydro) {
  qv <- c(q$hbd,q$hba,q$aromatic,q$hydro); rv <- c(hbd,hba,arom,hydro)
  inter <- sum(pmin(qv,rv)); union <- sum(pmax(qv,rv))
  if (union==0) return(0); round(inter/union,4)
}

agent4_generate <- function(prot_feat, query_ph, hits, pocket_out,
                             n_pep=40, seed=42) {
  set.seed(seed)
  pool <- character(0)
  if (query_ph$hbd    > 3) pool <- c(pool, rep(AA_ROLES$hbd_donor,    3))
  if (query_ph$hba    > 3) pool <- c(pool, rep(AA_ROLES$hba_acceptor, 3))
  if (query_ph$aromatic>2) pool <- c(pool, rep(AA_ROLES$aromatic,     4))
  if (query_ph$hydro  > 4) pool <- c(pool, rep(AA_ROLES$hydrophobic,  3))
  if (prot_feat$frac_charged>0.2)
    pool <- c(pool, rep(c(AA_ROLES$charged_pos, AA_ROLES$charged_neg), 2))
  if (length(pool)<8) pool <- AA_ROLES$all

  # Key residues from best pocket bias the pool further
  pocket_key <- pocket_out$best_pocket$key_residues
  if (grepl("PHE|TRP|TYR",pocket_key)) pool <- c(pool, rep(AA_ROLES$aromatic,2))
  if (grepl("ARG|LYS",pocket_key))     pool <- c(pool, rep(AA_ROLES$charged_pos,2))
  if (grepl("ASP|GLU",pocket_key))     pool <- c(pool, rep(AA_ROLES$charged_neg,2))

  len_min <- max(5,  round(prot_feat$mw_kda/8))
  len_max <- min(20, max(14, round(prot_feat$mw_kda/4)))
  best_aff<- max(hits$binding_affinity, na.rm=TRUE)

  peptides <- lapply(seq_len(n_pep), function(i) {
    len  <- sample(len_min:len_max, 1)
    seq  <- paste(sample(pool, len, replace=TRUE), collapse="")
    aa_v <- strsplit(seq,"")[[1]]
    n_hbd  <- sum(aa_v %in% AA_ROLES$hbd_donor)
    n_hba  <- sum(aa_v %in% AA_ROLES$hba_acceptor)
    n_arom <- sum(aa_v %in% AA_ROLES$aromatic)
    n_hydro<- sum(aa_v %in% AA_ROLES$hydrophobic)
    mw_pep <- round(sum(sapply(aa_v,function(a) AA_MW[a]),na.rm=TRUE)/1000,2)
    charge <- sum(sapply(aa_v,function(a) AA_CHARGE[a]),na.rm=TRUE)
    hydro  <- round(mean(sapply(aa_v,function(a) AA_HYDRO[a]),na.rm=TRUE),2)
    ph_m   <- tanimoto_score(query_ph,n_hbd,n_hba,n_arom,n_hydro)

    # Simulated ProteinMPNN score
    mpnn_score <- round(ph_m * runif(1,0.8,1.1) + rnorm(1,0,0.05),3)
    mpnn_score <- max(0,min(1,mpnn_score))

    bind_sc    <- round(best_aff * ph_m * runif(1,0.65,1.05), 2)
    dock_sc    <- round(-5.5 - ph_m*4.5 + rnorm(1,0,0.9), 2)
    sol        <- round(max(0,min(1,0.25 + ph_m*0.55 + rnorm(1,0,0.08))),3)
    tox        <- round(max(0,min(1,0.65 - ph_m*0.45 + rnorm(1,0,0.08))),3)
    stage      <- sample(
      c("Generation","Binding Prediction","Structure Modeling",
        "Docking","Optimization","Final Candidate"),
      1, prob=c(0.20,0.20,0.15,0.15,0.15,0.15)*(0.5+ph_m))

    data.frame(peptide=paste0("PEP-",sprintf("%03d",i)),
               sequence=seq, length=len, mw_kda=mw_pep,
               net_charge=round(charge,1), mean_hydro=hydro,
               ph_match=round(ph_m,3), mpnn_score=mpnn_score,
               binding_score=bind_sc, docking_score=dock_sc,
               solubility=sol, toxicity=tox,
               n_hbd=n_hbd, n_hba=n_hba, n_aromatic=n_arom,
               n_hydrophobic=n_hydro, stage=stage)
  })
  bind_rows(peptides) %>% arrange(desc(ph_match), desc(mpnn_score))
}

# Agent 5: Docking Evaluator (AutoDock-GPU proxy)
agent5_docking <- function(peptides, pocket_out) {
  set.seed(7)
  pv <- pocket_out$best_pocket$volume_A3
  peptides %>%
    mutate(
      pose_rmsd     = round(runif(n(), 0.5, 3.5), 2),
      vdw_energy    = round(-3 - ph_match*2 + rnorm(n(),0,0.5), 2),
      electrostatic = round(-1 - ph_match*1.5 + rnorm(n(),0,0.4), 2),
      pocket_volume_used = round(pmin(length * 110, pv) / pv, 2),
      docking_confidence = round(pmax(0, pmin(1,
                           0.4*ph_match + 0.3*(1+docking_score/15) +
                           0.3*(1-pose_rmsd/5))), 3)
    )
}

# Agent 6: Bayesian Optimizer (RL feedback loop, simulated)
agent6_optimize <- function(peptides, n_iter=5, seed=123) {
  set.seed(seed)
  top <- peptides %>% arrange(desc(ph_match), desc(binding_score)) %>% head(5)

  history <- lapply(seq_len(n_iter), function(iter) {
    improvement <- runif(1, 0.02, 0.12) * (1 - iter/n_iter * 0.4)
    data.frame(
      iteration   = iter,
      best_score  = round(max(top$binding_score) * (1 + improvement * iter/n_iter), 3),
      mean_ph     = round(mean(top$ph_match)     * (1 + improvement * 0.5), 3),
      n_candidates= max(1, round(nrow(peptides) * (1 - 0.12*iter))),
      mutation_type = sample(c("Point mutation","Alanine scan","N-term ext.",
                               "C-term trim","Charge swap","Aromatic sub."), 1)
    )
  })
  bind_rows(history)
}

# ============================================================
# C — BINDINGDB LIGAND DATABASE
# ============================================================
bindingdb_ligands <- data.frame(
  name=c(
    "Imatinib","Dasatinib","Erlotinib","Gefitinib","Sorafenib",
    "Vemurafenib","Lapatinib","Sunitinib","Nilotinib","Bosutinib",
    "Aspirin","Ibuprofen","Naproxen","Celecoxib","Diclofenac",
    "Lisinopril","Captopril","Enalapril","Ramipril","Perindopril",
    "Ritonavir","Lopinavir","Saquinavir","Indinavir","Atazanavir",
    "Morphine","Fentanyl","Naloxone","Buprenorphine","Methadone",
    "Metformin","Sitagliptin","Saxagliptin","Alogliptin","Vildagliptin",
    "Vancomycin","Daptomycin","Linezolid","Rifampicin","Ciprofloxacin"
  ),
  target_class=c(rep("Kinase",10),rep("COX",5),rep("ACE",5),
                 rep("Protease",5),rep("GPCR-Opioid",5),
                 rep("General",5),rep("General",5)),
  binding_affinity=c(8.9,8.5,7.8,7.6,8.1,9.0,8.3,7.9,8.7,8.0,
                     6.2,5.8,6.1,7.4,6.9,9.2,8.7,8.5,8.8,8.6,
                     9.5,9.1,8.8,9.3,8.9,8.4,9.0,7.5,9.2,8.1,
                     5.5,7.8,8.2,7.9,8.0,9.8,8.5,7.2,8.9,7.6),
  ph_hbd  =c(3,4,3,3,4,3,4,4,3,3,2,2,2,3,2,4,3,4,4,4,4,4,4,4,4,3,2,3,3,2,2,3,3,3,3,5,4,3,4,3),
  ph_hba  =c(5,6,5,5,5,4,6,5,5,5,3,2,3,4,3,5,4,5,5,5,6,6,5,6,6,4,3,4,4,3,2,5,4,4,4,7,5,4,6,4),
  ph_arom =c(4,5,4,4,4,4,5,4,4,4,2,1,2,3,2,2,1,2,2,2,4,4,4,4,4,3,3,3,3,2,1,3,2,2,2,3,2,2,3,2),
  ph_hydro=c(5,5,4,4,5,4,5,4,5,4,2,3,3,3,3,3,2,3,3,3,5,5,5,5,5,4,4,3,4,3,1,3,2,2,2,2,3,2,3,2),
  mw_est  =c(493,488,393,446,464,489,581,398,529,530,180,206,230,381,296,
             405,217,348,416,368,720,628,670,614,704,285,336,327,467,309,
             130,407,315,339,304,1449,1621,337,822,331),
  stringsAsFactors=FALSE
)

# ============================================================
# D — PIPELINE & AGENT GRAPH DEFINITIONS
# ============================================================
make_agent_nodes <- function() {
  data.frame(
    id=1:6,
    label=c(
      "Agent 1\nProtein LM",
      "Agent 2\nPocket Pred.",
      "Agent 3\nLigand GNN",
      "Agent 4\nPeptide Gen.",
      "Agent 5\nDocking Eval.",
      "Agent 6\nBayesian Opt."
    ),
    group=c("lm","pocket","ligand","generator","docking","optimizer"),
    title=c(
      "<b>Agent 1 — Protein Language Model</b><br>ESM-2 (650M) proxy: sequence embeddings, secondary structure prediction, druggability scoring.",
      "<b>Agent 2 — Pocket Prediction</b><br>AlphaFold2 / DeepSite proxy: binding pocket detection, volume estimation, key residue identification.",
      "<b>Agent 3 — Ligand Knowledge Agent</b><br>GNN / ChemBERTa proxy: learned molecular embeddings + Tanimoto similarity against BindingDB.",
      "<b>Agent 4 — Peptide Generator</b><br>ProteinMPNN proxy: pharmacophore-biased sequence generation anchored to pocket key residues.",
      "<b>Agent 5 — Docking Evaluator</b><br>AutoDock-GPU proxy: pose RMSD, VdW/electrostatic energies, pocket volume utilisation.",
      "<b>Agent 6 — Bayesian Optimizer</b><br>RL feedback loop: iterative point mutations, alanine scanning, sequence refinement."
    ),
    stringsAsFactors=FALSE
  )
}

make_agent_edges <- function() {
  data.frame(
    from=c(1,2,3,1,4,5),
    to  =c(2,4,4,3,5,6),
    arrows="to",
    label=c("embeddings","pocket profile","ligand embeds","seq context","scores","refine"),
    stringsAsFactors=FALSE
  )
}

agent_group_cols <- list(
  lm       =list(background="#CECBF6",border="#534AB7",highlight=list(background="#AFA9EC",border="#3C3489")),
  pocket   =list(background="#9FE1CB",border="#0F6E56",highlight=list(background="#5DCAA5",border="#085041")),
  ligand   =list(background="#B5D4F4",border="#185FA5",highlight=list(background="#85B7EB",border="#0C447C")),
  generator=list(background="#FAC775",border="#854F0B",highlight=list(background="#EF9F27",border="#633806")),
  docking  =list(background="#F5C4B3",border="#993C1D",highlight=list(background="#F0997B",border="#712B13")),
  optimizer=list(background="#C0DD97",border="#3B6D11",highlight=list(background="#97C459",border="#27500A"))
)

EXAMPLE_FASTA <- ">sp|P00519|ABL1_HUMAN BCR-ABL tyrosine kinase | Homo sapiens
MGPSENDPNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNY
ITPVNSLEKHSWYHGPVSRNAAEYLLSSGINGSFLVRESESSPGQRSISLRYEGRVYHYR"

# ============================================================
# E — UI
# ============================================================
ui <- dashboardPage(
  skin="blue",
  dashboardHeader(
    title=tags$span(
      style="display:flex;align-items:center;gap:8px;",
      tags$img(src="pepai_logo.png",
               style="height:38px;width:auto;border-radius:6px;object-fit:contain;",
               onerror="this.style.display='none'"),
      tags$span(style="font-weight:900;font-size:20px;letter-spacing:.5px;",
        tags$span(style="color:#fff;","Pep"),
        tags$span(style="color:#f5a623;","AI"))
    )
  ),

  dashboardSidebar(
    tags$head(tags$style(HTML("
      .skin-blue .sidebar-menu > li.active > a,
      .skin-blue .sidebar-menu > li:hover > a { border-left:4px solid #e74c3c; }
      .sidebar-menu .treeview-menu > li > a { padding-left:30px; }
    "))),
    sidebarMenu(
      menuItem("Target Protein",     tabName="target",    icon=icon("dna")),
      menuItem("AI Agent View",      tabName="agents",    icon=icon("robot")),
      menuItem("BindingDB Hits",     tabName="ligands",   icon=icon("database")),
      menuItem("Generated Peptides", tabName="peptides",  icon=icon("flask")),
      menuItem("Optimization Loop",  tabName="optim",     icon=icon("rotate")),
      menuItem("Analytics",          tabName="analytics", icon=icon("chart-bar")),
      menuItem("About / Science",    tabName="about",     icon=icon("book-open")),
      menuItem("Developers",         tabName="developers",icon=icon("users"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      /* ── Global ── */
      textarea#fasta_input { font-family:'Courier New',monospace; font-size:12px;
        background:#1e1e2e; color:#cdd6f4; border:1px solid #444;
        border-radius:6px; padding:10px; line-height:1.5; }
      .content-wrapper,.right-side { background:#f0f2f5 !important; }

      /* ── Cards ── */
      .ph-card  { background:#fff; border-radius:8px; padding:10px 14px;
        margin-bottom:8px; border-left:4px solid #3c8dbc;
        box-shadow:0 1px 4px rgba(0,0,0,0.07); }
      .ph-label { font-weight:700; font-size:10px; color:#888;
        text-transform:uppercase; letter-spacing:.5px; }
      .ph-value { font-size:22px; font-weight:700; color:#2c3e50; }
      .aa-pill  { display:inline-block; padding:2px 8px; border-radius:10px;
        font-size:12px; font-weight:700; color:#fff; margin:2px;
        font-family:'Courier New',monospace; }
      .section-info { font-size:13px; color:#555; margin-bottom:12px; line-height:1.6; }

      /* ── AI Agent tab ── */
      .agent-card { background:#fff; border-radius:12px; padding:18px;
        box-shadow:0 2px 10px rgba(0,0,0,0.09); margin-bottom:16px;
        border-top:4px solid #3c8dbc; }
      .agent-card h4 { margin:0 0 4px; font-size:15px; font-weight:700; color:#2c3e50; }
      .agent-subtitle { font-size:11px; color:#888; margin-bottom:10px;
        text-transform:uppercase; letter-spacing:.4px; }
      .agent-metric { display:inline-block; background:#f0f4f8; border-radius:8px;
        padding:6px 12px; margin:3px; font-size:12px; color:#2c3e50; font-weight:600; }
      .agent-metric span { display:block; font-size:10px; color:#888;
        font-weight:400; text-transform:uppercase; }
      .ai-badge { display:inline-block; padding:2px 9px; border-radius:12px;
        font-size:10px; font-weight:700; color:#fff; margin-right:5px;
        background:#8e44ad; letter-spacing:.3px; }
      .feedback-arrow { text-align:center; font-size:22px; color:#e74c3c;
        margin:4px 0; line-height:1; }
      .optim-row { background:#fff; border-radius:8px; padding:10px 14px;
        margin-bottom:6px; box-shadow:0 1px 4px rgba(0,0,0,0.06);
        display:flex; align-items:center; gap:14px; }
      .optim-iter { font-size:11px; font-weight:700; color:#888;
        text-transform:uppercase; min-width:55px; }
      .optim-bar { flex:1; height:8px; background:#eee; border-radius:4px; overflow:hidden; }
      .optim-fill { height:100%; border-radius:4px;
        background:linear-gradient(90deg,#27ae60,#2ecc71); }
    "))),

    tabItems(

      # ── TAB 1: TARGET INPUT ────────────────────────────────
      tabItem(tabName="target",
        fluidRow(
          box(
            title=tagList(icon("dna")," Paste Target Protein FASTA Sequence"),
            width=12, status="danger", solidHeader=TRUE,
            p(class="section-info",
              "Paste a standard UniProt / NCBI FASTA. The 6-agent AI pipeline will:",
              tags$ol(style="margin:6px 0;",
                tags$li("Run a Protein Language Model (ESM-2 proxy) for embeddings & druggability"),
                tags$li("Predict binding pockets (AlphaFold2 / DeepSite proxy)"),
                tags$li("Search BindingDB with a GNN-enhanced ligand similarity agent"),
                tags$li("Generate peptides anchored to the pocket (ProteinMPNN proxy)"),
                tags$li("Evaluate poses with a docking agent (AutoDock-GPU proxy)"),
                tags$li("Iteratively optimise with a Bayesian RL agent")
              )
            ),
            textAreaInput("fasta_input", label=NULL, value=EXAMPLE_FASTA,
                          width="100%", rows=6,
                          placeholder=">ProteinName\nMGPSENDPNLFVALYD..."),
            fluidRow(
              column(3, numericInput("top_k","Top-K Ligand Hits",value=10,min=3,max=40)),
              column(3, numericInput("n_pep","Peptides to Generate",value=40,min=10,max=150,step=10)),
              column(3, numericInput("n_iter","Optimizer Iterations",value=5,min=2,max=15)),
              column(3, selectInput("class_override","Override Target Class",
                choices=c("Auto-detect","Kinase","COX","ACE","Protease","GPCR-Opioid","General"),
                selected="Auto-detect"))
            ),
            fluidRow(column(12,
              actionButton("run_btn","  Run AI Agent Pipeline \u25B6",
                icon=icon("robot"), class="btn btn-danger btn-lg",
                style="margin-top:8px;font-weight:700;"),
              actionButton("load_example","Load BCR-ABL Example",
                icon=icon("vial"), class="btn btn-default btn-sm",
                style="margin-top:8px;margin-left:10px;")
            ))
          )
        ),
        fluidRow(
          uiOutput("fasta_status_box"),
          uiOutput("protein_stats_box")
        ),
        fluidRow(uiOutput("aa_freq_box"))
      ),

      # ── TAB 3: AI AGENT VIEW ───────────────────────────────
      tabItem(tabName="agents",
        fluidRow(column(12,
          div(style="background:linear-gradient(135deg,#1a1a2e,#16213e,#0f3460);
                     border-radius:12px;padding:20px 24px;margin-bottom:18px;
                     color:#fff;",
            h3(style="margin:0 0 4px;font-size:20px;",
               icon("robot"),"  PepAI Architecture"),
            p(style="margin:0;opacity:.75;font-size:13px;",
              "6 autonomous AI agents operating in a feedback loop for ",
              strong(style="color:#e74c3c;","de novo peptide drug design"))
          )
        )),
        fluidRow(
          column(12, visNetworkOutput("agent_graph", height="340px"))
        ),
        br(),
        fluidRow(
          column(6, uiOutput("agent1_card")),
          column(6, uiOutput("agent2_card"))
        ),
        fluidRow(
          column(6, uiOutput("agent3_card")),
          column(6, uiOutput("agent4_card"))
        ),
        fluidRow(
          column(6, uiOutput("agent5_card")),
          column(6, uiOutput("agent6_card"))
        )
      ),

      # ── TAB 4: BINDINGDB HITS ──────────────────────────────
      tabItem(tabName="ligands",
        fluidRow(
          valueBoxOutput("vb_n_hits"),
          valueBoxOutput("vb_best_aff"),
          valueBoxOutput("vb_top_class")
        ),
        fluidRow(box(
          title=tagList(tags$span(class="ai-badge","Agent 3"),
                        " BindingDB Ligands — GNN-Enhanced Similarity"),
          width=12, status="primary", solidHeader=TRUE,
          p(class="section-info",
            "Ligands ranked by combined Tanimoto fingerprint + simulated GNN embedding similarity."),
          DTOutput("ligand_table")
        )),
        fluidRow(
          box(title="Similarity Distribution",width=6,status="info",solidHeader=TRUE,
              plotlyOutput("sim_hist",height="260px")),
          box(title="Binding Affinity vs Similarity",width=6,status="success",solidHeader=TRUE,
              plotlyOutput("sim_aff",height="260px"))
        )
      ),

      # ── TAB 5: GENERATED PEPTIDES ──────────────────────────
      tabItem(tabName="peptides",
        fluidRow(
          valueBoxOutput("vb_total_pep"),
          valueBoxOutput("vb_final_pep"),
          valueBoxOutput("vb_top_bind")
        ),
        fluidRow(
          box(title="Filters",width=3,status="warning",solidHeader=TRUE,
            sliderInput("ph_min","Min Pharmacophore Match",min=0,max=1,value=0.25,step=0.01),
            sliderInput("mpnn_min","Min MPNN Score",min=0,max=1,value=0.2,step=0.01),
            sliderInput("bind_r","Binding Score",min=0,max=12,value=c(3,12),step=0.1),
            sliderInput("dock_r","Docking Score",min=-16,max=-2,value=c(-16,-3),step=0.1),
            sliderInput("sol_min","Min Solubility",min=0,max=1,value=0.2,step=0.01),
            sliderInput("tox_max","Max Toxicity",min=0,max=1,value=0.7,step=0.01),
            checkboxGroupInput("stage_sel","Stage",
              choices=c("Generation","Binding Prediction","Structure Modeling",
                        "Docking","Optimization","Final Candidate"),
              selected=c("Generation","Binding Prediction","Structure Modeling",
                         "Docking","Optimization","Final Candidate")),
            actionButton("reset_pep","Reset",class="btn-sm btn-default",icon=icon("undo"))
          ),
          box(
            title=tagList(tags$span(class="ai-badge","Agent 4+5"),
                          " Pocket-Anchored Peptide Candidates"),
            width=9, status="primary", solidHeader=TRUE,
            DTOutput("peptide_table")
          )
        ),
        fluidRow(box(
          title="Best Candidate — Residue Breakdown",
          width=12, status="success", solidHeader=TRUE,
          uiOutput("best_card")
        ))
      ),

      # ── TAB 6: OPTIMIZATION LOOP ───────────────────────────
      tabItem(tabName="optim",
        fluidRow(column(12,
          div(style="background:linear-gradient(135deg,#1a1a2e,#16213e);
                     border-radius:12px;padding:16px 22px;margin-bottom:16px;color:#fff;",
            h4(style="margin:0 0 3px;", icon("rotate"),"  Agent 6 — Bayesian Optimization Loop"),
            p(style="margin:0;opacity:.7;font-size:12px;",
              "Iterative RL-guided sequence refinement: generate → dock → score → mutate → repeat")
          )
        )),
        fluidRow(
          box(title="Optimization Score Trajectory",width=8,status="primary",solidHeader=TRUE,
              plotlyOutput("optim_plot",height="320px")),
          box(title="Iteration Log",width=4,status="info",solidHeader=TRUE,
              uiOutput("optim_table"))
        ),
        fluidRow(
          box(title="Pharmacophore Match Improvement",width=6,status="success",solidHeader=TRUE,
              plotlyOutput("optim_ph",height="260px")),
          box(title="Candidate Pool Reduction",width=6,status="warning",solidHeader=TRUE,
              plotlyOutput("optim_pool",height="260px"))
        )
      ),

      # ── TAB 7: ANALYTICS ───────────────────────────────────
      tabItem(tabName="analytics",
        fluidRow(
          box(title="Pharmacophore Match vs Binding Score",
              width=6,status="primary",solidHeader=TRUE,
              plotlyOutput("ph_bind",height="300px")),
          box(title="MPNN Score vs Docking Score",
              width=6,status="info",solidHeader=TRUE,
              plotlyOutput("mpnn_dock",height="300px"))
        ),
        fluidRow(
          box(title="Candidates by Stage",
              width=6,status="success",solidHeader=TRUE,
              plotlyOutput("bar_stage",height="270px")),
          box(title="Solubility vs Toxicity",
              width=6,status="danger",solidHeader=TRUE,
              plotlyOutput("sol_tox",height="270px"))
        ),
        fluidRow(
          box(title="Residue Frequency in Top Peptides",
              width=12,status="warning",solidHeader=TRUE,
              plotlyOutput("aa_freq_plot",height="230px"))
        )
      ),

      # ── TAB 8: ABOUT / SCIENCE ─────────────────────────────
      tabItem(tabName="about",
        fluidRow(
          box(
            title=tagList(icon("book-open")," Scientific Framework"),
            width=8, status="info", solidHeader=TRUE,
            h4("What PepAI Is"),
            p("PepAI is a ", strong("pharmacophore-guided computational peptide discovery framework"),
              " with an AI-agent architecture overlay. The current implementation uses",
              " simulated AI modules; each can be replaced with a real model:"),
            tags$table(class="table table-bordered table-sm",style="font-size:12px;",
              tags$thead(tags$tr(
                tags$th("Agent"), tags$th("Role"), tags$th("Real Model"), tags$th("Status")
              )),
              tags$tbody(
                tags$tr(tags$td("Agent 1"),tags$td("Protein Language Model"),
                        tags$td("ESM-2 650M (Meta)"),tags$td(tags$span(style="color:#e67e22;","Simulated"))),
                tags$tr(tags$td("Agent 2"),tags$td("Pocket Prediction"),
                        tags$td("AlphaFold2 + DeepSite"),tags$td(tags$span(style="color:#e67e22;","Simulated"))),
                tags$tr(tags$td("Agent 3"),tags$td("Ligand GNN"),
                        tags$td("ChemBERTa / Mol2Vec"),tags$td(tags$span(style="color:#e67e22;","Simulated"))),
                tags$tr(tags$td("Agent 4"),tags$td("Peptide Generator"),
                        tags$td("ProteinMPNN / RFdiffusion"),tags$td(tags$span(style="color:#e67e22;","Simulated"))),
                tags$tr(tags$td("Agent 5"),tags$td("Docking Evaluator"),
                        tags$td("AutoDock-GPU / Glide"),tags$td(tags$span(style="color:#e67e22;","Simulated"))),
                tags$tr(tags$td("Agent 6"),tags$td("Bayesian Optimizer"),
                        tags$td("BoTorch / REINFORCE"),tags$td(tags$span(style="color:#e67e22;","Simulated")))
              )
            ),
            hr(),
            h4("How to Describe This in a Paper"),
            tags$blockquote(style="border-left:4px solid #3c8dbc;padding:10px 16px;
                                   background:#f0f7ff;border-radius:0 8px 8px 0;font-style:italic;",
              "PepAI is a pharmacophore-guided computational peptide discovery pipeline,",
              " designed as a modular framework for integrating AI-driven agents",
              " (protein language models, pocket predictors, generative design models)",
              " into an autonomous drug discovery loop."
            ),
            hr(),
            h4("Feedback Loop"),
            p("The Bayesian Optimizer (Agent 6) creates a closed feedback loop:"),
            tags$code(style="display:block;background:#1e1e2e;color:#cdd6f4;padding:12px;border-radius:8px;",
              "generate peptides → dock → score → mutate top-K → repeat"),
            p("This is the key feature that distinguishes an AI-agent system from a static pipeline.")
          ),
          box(
            title=tagList(icon("circle-info")," Quick Reference"),
            width=4, status="warning", solidHeader=TRUE,
            h5("Example Protein"),
            p(tags$code("BCR-ABL tyrosine kinase"), br(), "UniProt: P00519"),
            hr(),
            h5("Key Outputs"),
            tags$ul(
              tags$li("ESM-2 druggability score"),
              tags$li("Pocket volume & key residues"),
              tags$li("GNN-enhanced ligand hits"),
              tags$li("MPNN-scored peptides"),
              tags$li("Docking confidence"),
              tags$li("Bayesian optimization trace")
            ),
            hr(),
            h5("Built With"),
            p("R · Shiny · visNetwork · plotly · DT"),
            hr(),
            p(em(style="font-size:11px;color:#888;",
                 "All AI scores are simulated for demonstration.",
                 " Replace agent functions with real API calls for production."))
          )
        )
      ),

      # ── TAB 9: DEVELOPERS ──────────────────────────────────
      tabItem(tabName="developers",

        # ── Hero banner ──────────────────────────────────────
        fluidRow(column(12,
          div(style="background:linear-gradient(135deg,#0d1b2a,#1b2a4a,#1a3a5c);
                     border-radius:16px;padding:32px 36px;margin-bottom:24px;
                     color:#fff;text-align:center;position:relative;overflow:hidden;",
            div(style="position:absolute;top:-40px;right:-40px;width:200px;height:200px;
                       border-radius:50%;background:rgba(245,166,35,0.07);pointer-events:none;"),
            div(style="position:absolute;bottom:-60px;left:-30px;width:250px;height:250px;
                       border-radius:50%;background:rgba(60,141,188,0.07);pointer-events:none;"),
            tags$img(src="pepai_logo.png",
                     style="height:80px;width:auto;margin-bottom:16px;border-radius:10px;
                            filter:drop-shadow(0 4px 12px rgba(0,0,0,0.4));",
                     onerror="this.style.display='none'"),
            h2(style="margin:0 0 6px;font-size:30px;font-weight:900;letter-spacing:1px;",
              tags$span(style="color:#fff;","Pep"),
              tags$span(style="color:#f5a623;","AI")),
            p(style="margin:0 0 14px;opacity:.7;font-size:14px;","AI-Augmented Peptide Discovery Pipeline"),
            div(style="display:inline-block;background:rgba(245,166,35,0.18);
                       border:1px solid rgba(245,166,35,0.4);border-radius:20px;
                       padding:6px 20px;",
              tags$span(style="font-size:12px;font-weight:700;color:#f5a623;letter-spacing:.5px;",
                icon("users"),"  Co-Founders — InteractIQ")),
            p(style="margin:10px 0 0;font-size:12px;opacity:.55;",
              "AI-first approach for drug discovery problems & service provider")
          )
        )),

        # ── Developer cards ───────────────────────────────────
        fluidRow(

          # ── Rik Ganguly ──
          column(4,
            div(style="background:#fff;border-radius:16px;overflow:hidden;
                       box-shadow:0 6px 28px rgba(0,0,0,0.10);margin-bottom:16px;",
              # Card header
              div(style="background:linear-gradient(135deg,#1a5276,#2980b9);
                         padding:24px 20px;text-align:center;",
                div(style="width:72px;height:72px;border-radius:50%;
                           background:rgba(255,255,255,0.2);border:3px solid rgba(255,255,255,0.5);
                           display:flex;align-items:center;justify-content:center;
                           margin:0 auto 12px;font-size:26px;color:#fff;font-weight:900;",
                    "RG"),
                h4(style="margin:0 0 3px;color:#fff;font-size:18px;font-weight:800;","Rik Ganguly"),
                p(style="margin:0 0 8px;color:rgba(255,255,255,0.75);font-size:11px;
                         text-transform:uppercase;letter-spacing:.7px;",
                  "Co-Founder · Lead Developer & Bioinformatician"),
                div(style="display:inline-block;background:rgba(245,166,35,0.25);
                           border-radius:12px;padding:3px 12px;",
                  tags$span(style="font-size:10px;font-weight:700;color:#f5a623;",
                    icon("building"),"  InteractIQ"))
              ),
              # Contact info body
              div(style="padding:18px 20px;",
                p(style="font-size:10px;font-weight:800;color:#3c8dbc;text-transform:uppercase;
                         letter-spacing:.8px;margin:0 0 12px;border-bottom:1px solid #eef2f7;
                         padding-bottom:8px;",
                  icon("address-card"),"  Contact Info"),
                # LinkedIn
                div(style="display:flex;align-items:flex-start;gap:10px;margin-bottom:10px;",
                  div(style="width:28px;height:28px;border-radius:6px;background:#0077b5;
                             display:flex;align-items:center;justify-content:center;flex-shrink:0;",
                      tags$span(style="color:#fff;font-size:12px;font-weight:900;","in")),
                  div(
                    p(style="margin:0 0 2px;font-size:10px;font-weight:700;color:#888;
                             text-transform:uppercase;","Rik's Profile"),
                    tags$a(href="https://linkedin.com/in/rik-ganguly-70a7a197",target="_blank",
                      style="font-size:11px;color:#0077b5;word-break:break-all;",
                      "linkedin.com/in/rik-ganguly-70a7a197")
                  )
                ),
                # Phone
                div(style="display:flex;align-items:flex-start;gap:10px;margin-bottom:10px;",
                  div(style="width:28px;height:28px;border-radius:6px;background:#27ae60;
                             display:flex;align-items:center;justify-content:center;flex-shrink:0;",
                      icon("phone",style="color:#fff;font-size:11px;")),
                  div(
                    p(style="margin:0 0 2px;font-size:10px;font-weight:700;color:#888;
                             text-transform:uppercase;","Phone"),
                    p(style="margin:0;font-size:12px;color:#2c3e50;font-weight:600;",
                      "+91-6350433960"),
                    p(style="margin:0;font-size:10px;color:#aaa;","Work")
                  )
                ),
                # Email
                div(style="display:flex;align-items:flex-start;gap:10px;",
                  div(style="width:28px;height:28px;border-radius:6px;background:#e74c3c;
                             display:flex;align-items:center;justify-content:center;flex-shrink:0;",
                      icon("envelope",style="color:#fff;font-size:11px;")),
                  div(
                    p(style="margin:0 0 2px;font-size:10px;font-weight:700;color:#888;
                             text-transform:uppercase;","Email"),
                    tags$a(href="mailto:rikgangulybioinfo@gmail.com",
                      style="font-size:11px;color:#e74c3c;word-break:break-all;",
                      "rikgangulybioinfo@gmail.com")
                  )
                )
              )
            )
          ),

          # ── Gautam Ahuja ──
          column(4,
            div(style="background:#fff;border-radius:16px;overflow:hidden;
                       box-shadow:0 6px 28px rgba(0,0,0,0.10);margin-bottom:16px;",
              div(style="background:linear-gradient(135deg,#1e8449,#27ae60);
                         padding:24px 20px;text-align:center;",
                div(style="width:72px;height:72px;border-radius:50%;
                           background:rgba(255,255,255,0.2);border:3px solid rgba(255,255,255,0.5);
                           display:flex;align-items:center;justify-content:center;
                           margin:0 auto 12px;font-size:26px;color:#fff;font-weight:900;",
                    "GA"),
                h4(style="margin:0 0 3px;color:#fff;font-size:18px;font-weight:800;","Gautam Ahuja"),
                p(style="margin:0 0 8px;color:rgba(255,255,255,0.75);font-size:11px;
                         text-transform:uppercase;letter-spacing:.7px;",
                  "Co-Founder · Developer & Computational Scientist"),
                div(style="display:inline-block;background:rgba(245,166,35,0.25);
                           border-radius:12px;padding:3px 12px;",
                  tags$span(style="font-size:10px;font-weight:700;color:#f5a623;",
                    icon("building"),"  InteractIQ"))
              ),
              div(style="padding:18px 20px;",
                p(style="font-size:10px;font-weight:800;color:#27ae60;text-transform:uppercase;
                         letter-spacing:.8px;margin:0 0 12px;border-bottom:1px solid #eef2f7;
                         padding-bottom:8px;",
                  icon("address-card"),"  Contact Info"),
                # LinkedIn
                div(style="display:flex;align-items:flex-start;gap:10px;margin-bottom:10px;",
                  div(style="width:28px;height:28px;border-radius:6px;background:#0077b5;
                             display:flex;align-items:center;justify-content:center;flex-shrink:0;",
                      tags$span(style="color:#fff;font-size:12px;font-weight:900;","in")),
                  div(
                    p(style="margin:0 0 2px;font-size:10px;font-weight:700;color:#888;
                             text-transform:uppercase;","Gautam's Profile"),
                    tags$a(href="https://linkedin.com/in/gautam8387",target="_blank",
                      style="font-size:11px;color:#0077b5;word-break:break-all;",
                      "linkedin.com/in/gautam8387")
                  )
                ),
                # Website
                div(style="display:flex;align-items:flex-start;gap:10px;margin-bottom:10px;",
                  div(style="width:28px;height:28px;border-radius:6px;background:#8e44ad;
                             display:flex;align-items:center;justify-content:center;flex-shrink:0;",
                      icon("globe",style="color:#fff;font-size:11px;")),
                  div(
                    p(style="margin:0 0 2px;font-size:10px;font-weight:700;color:#888;
                             text-transform:uppercase;","Website"),
                    tags$a(href="https://gautam8387.github.io/",target="_blank",
                      style="font-size:11px;color:#8e44ad;word-break:break-all;",
                      "gautam8387.github.io/"),
                    p(style="margin:0;font-size:10px;color:#aaa;","Personal")
                  )
                ),
                # Email
                div(style="display:flex;align-items:flex-start;gap:10px;",
                  div(style="width:28px;height:28px;border-radius:6px;background:#e74c3c;
                             display:flex;align-items:center;justify-content:center;flex-shrink:0;",
                      icon("envelope",style="color:#fff;font-size:11px;")),
                  div(
                    p(style="margin:0 0 2px;font-size:10px;font-weight:700;color:#888;
                             text-transform:uppercase;","Email"),
                    tags$a(href="mailto:goutamahuja8387@gmail.com",
                      style="font-size:11px;color:#e74c3c;word-break:break-all;",
                      "goutamahuja8387@gmail.com")
                  )
                )
              )
            )
          ),

          # ── Siddhant Poudel ──
          column(4,
            div(style="background:#fff;border-radius:16px;overflow:hidden;
                       box-shadow:0 6px 28px rgba(0,0,0,0.10);margin-bottom:16px;",
              div(style="background:linear-gradient(135deg,#922b21,#e74c3c);
                         padding:24px 20px;text-align:center;",
                div(style="width:72px;height:72px;border-radius:50%;
                           background:rgba(255,255,255,0.2);border:3px solid rgba(255,255,255,0.5);
                           display:flex;align-items:center;justify-content:center;
                           margin:0 auto 12px;font-size:26px;color:#fff;font-weight:900;",
                    "SP"),
                h4(style="margin:0 0 3px;color:#fff;font-size:18px;font-weight:800;","Siddhant Poudel"),
                p(style="margin:0 0 8px;color:rgba(255,255,255,0.75);font-size:11px;
                         text-transform:uppercase;letter-spacing:.7px;",
                  "Co-Founder · Developer & Research Scientist"),
                div(style="display:inline-block;background:rgba(245,166,35,0.25);
                           border-radius:12px;padding:3px 12px;",
                  tags$span(style="font-size:10px;font-weight:700;color:#f5a623;",
                    icon("building"),"  InteractIQ"))
              ),
              div(style="padding:18px 20px;",
                p(style="font-size:10px;font-weight:800;color:#e74c3c;text-transform:uppercase;
                         letter-spacing:.8px;margin:0 0 12px;border-bottom:1px solid #eef2f7;
                         padding-bottom:8px;",
                  icon("address-card"),"  Contact Info"),
                # LinkedIn
                div(style="display:flex;align-items:flex-start;gap:10px;margin-bottom:10px;",
                  div(style="width:28px;height:28px;border-radius:6px;background:#0077b5;
                             display:flex;align-items:center;justify-content:center;flex-shrink:0;",
                      tags$span(style="color:#fff;font-size:12px;font-weight:900;","in")),
                  div(
                    p(style="margin:0 0 2px;font-size:10px;font-weight:700;color:#888;
                             text-transform:uppercase;","Siddhant's Profile"),
                    tags$a(href="https://linkedin.com/in/siddhant-poudyal-3394381b2",target="_blank",
                      style="font-size:11px;color:#0077b5;word-break:break-all;",
                      "linkedin.com/in/siddhant-poudyal-3394381b2")
                  )
                ),
                # Email
                div(style="display:flex;align-items:flex-start;gap:10px;",
                  div(style="width:28px;height:28px;border-radius:6px;background:#e74c3c;
                             display:flex;align-items:center;justify-content:center;flex-shrink:0;",
                      icon("envelope",style="color:#fff;font-size:11px;")),
                  div(
                    p(style="margin:0 0 2px;font-size:10px;font-weight:700;color:#888;
                             text-transform:uppercase;","Email"),
                    tags$a(href="mailto:sidblitz1602@gmail.com",
                      style="font-size:11px;color:#e74c3c;word-break:break-all;",
                      "sidblitz1602@gmail.com")
                  )
                )
              )
            )
          )
        ),

        # ── InteractIQ + About panel ──────────────────────────
        fluidRow(column(12,
          div(style="background:#fff;border-radius:16px;padding:0;
                     box-shadow:0 4px 20px rgba(0,0,0,0.08);overflow:hidden;",

            # InteractIQ banner
            div(style="background:linear-gradient(90deg,#0d1b2a,#1b3a5c);
                       padding:18px 28px;display:flex;align-items:center;gap:16px;",
              div(style="width:44px;height:44px;border-radius:10px;
                         background:linear-gradient(135deg,#f5a623,#e67e22);
                         display:flex;align-items:center;justify-content:center;flex-shrink:0;",
                  tags$span(style="font-size:18px;font-weight:900;color:#fff;","IQ")),
              div(
                h4(style="margin:0 0 2px;color:#fff;font-size:16px;font-weight:800;","InteractIQ"),
                p(style="margin:0;color:rgba(255,255,255,0.6);font-size:12px;",
                  "AI-first approach for drug discovery problems & service provider")
              )
            ),

            div(style="padding:24px 28px;",
              fluidRow(
                column(6,
                  div(style="border-left:4px solid #3c8dbc;padding:14px 18px;
                             background:#f0f7ff;border-radius:0 10px 10px 0;margin-bottom:14px;",
                    h5(style="margin:0 0 7px;color:#1a5276;font-weight:800;",
                       icon("flask"),"  About PepAI"),
                    p(style="margin:0;font-size:13px;color:#444;line-height:1.7;",
                      "PepAI is an open-source pharmacophore-guided AI peptide discovery framework.",
                      " It uses a 6-agent autonomous loop to design, evaluate, and refine",
                      " peptide drug candidates against user-supplied protein targets.")
                  ),
                  div(style="border-left:4px solid #27ae60;padding:14px 18px;
                             background:#f0fff8;border-radius:0 10px 10px 0;",
                    h5(style="margin:0 0 7px;color:#1e8449;font-weight:800;",
                       icon("envelope"),"  General Contact"),
                    p(style="margin:0 0 6px;font-size:13px;color:#444;",
                      "For queries, collaborations, or feedback:"),
                    tags$a(href="mailto:rikgangulybioinfo@gmail.com",
                      style="color:#2980b9;font-weight:700;font-size:13px;",
                      "rikgangulybioinfo@gmail.com")
                  )
                ),
                column(6,
                  div(style="border-left:4px solid #f5a623;padding:14px 18px;
                             background:#fffbf0;border-radius:0 10px 10px 0;margin-bottom:14px;",
                    h5(style="margin:0 0 10px;color:#d68910;font-weight:800;",
                       icon("code"),"  Technology Stack"),
                    div(style="display:flex;flex-wrap:wrap;gap:6px;",
                      lapply(c("R","Shiny","shinydashboard","visNetwork","plotly","DT","ggplot2","dplyr"),
                        function(t) tags$span(style="background:#fff3cd;color:#856404;border-radius:10px;
                                                     padding:3px 10px;font-size:11px;font-weight:700;",t))
                    )
                  ),
                  div(style="border-left:4px solid #8e44ad;padding:14px 18px;
                             background:#f9f0ff;border-radius:0 10px 10px 0;",
                    h5(style="margin:0 0 10px;color:#6c3483;font-weight:800;",
                       icon("robot"),"  Simulated AI Models"),
                    div(style="display:flex;flex-wrap:wrap;gap:6px;",
                      lapply(c("ESM-2 650M","AlphaFold2","DeepSite","ChemBERTa","ProteinMPNN","AutoDock-GPU","BoTorch"),
                        function(t) tags$span(style="background:#f0e6ff;color:#6c3483;border-radius:10px;
                                                     padding:3px 10px;font-size:11px;font-weight:700;",t))
                    )
                  )
                )
              )
            )
          )
        ))
      )
    )
  )
)

# ============================================================
# F — SERVER
# ============================================================
server <- function(input, output, session) {

  observeEvent(input$load_example, {
    updateTextAreaInput(session,"fasta_input",value=EXAMPLE_FASTA)
    updateSelectInput(session,"class_override",selected="Kinase")
  })
  observeEvent(input$reset_pep, {
    updateSliderInput(session,"ph_min",  value=0.25)
    updateSliderInput(session,"mpnn_min",value=0.2)
    updateSliderInput(session,"bind_r",  value=c(3,12))
    updateSliderInput(session,"dock_r",  value=c(-16,-3))
    updateSliderInput(session,"sol_min", value=0.2)
    updateSliderInput(session,"tox_max", value=0.7)
  })

  # ── Core reactives ───────────────────────────────────────
  fasta_parsed <- eventReactive(input$run_btn, { req(input$fasta_input); parse_fasta(input$fasta_input) })
  fasta_valid  <- eventReactive(input$run_btn, { req(input$fasta_input); validate_fasta(input$fasta_input) })
  prot_feat    <- eventReactive(input$run_btn, {
    req(fasta_parsed()); fp <- fasta_parsed()
    if (!fp$ok) return(NULL)
    analyse_protein(fp$sequence)
  })
  lm_out <- eventReactive(input$run_btn, {
    req(prot_feat(), fasta_parsed())
    agent1_protein_lm(fasta_parsed()$sequence, prot_feat())
  })
  pocket_out <- eventReactive(input$run_btn, {
    req(prot_feat(), lm_out()); agent2_pocket(prot_feat(), lm_out())
  })
  query_ph <- eventReactive(input$run_btn, {
    req(prot_feat()); protein_to_pharmacophore(prot_feat())
  })
  effective_class <- eventReactive(input$run_btn, {
    req(prot_feat())
    if (input$class_override != "Auto-detect") input$class_override
    else prot_feat()$target_class
  })
  ligand_hits <- eventReactive(input$run_btn, {
    req(query_ph(), effective_class())
    agent3_ligand(query_ph(), effective_class(), top_n=input$top_k)
  })
  peptides_rv <- eventReactive(input$run_btn, {
    req(prot_feat(), query_ph(), ligand_hits(), pocket_out())
    raw <- agent4_generate(prot_feat(), query_ph(), ligand_hits(),
                           pocket_out(), n_pep=input$n_pep, seed=7)
    agent5_docking(raw, pocket_out())
  })
  optim_rv <- eventReactive(input$run_btn, {
    req(peptides_rv()); agent6_optimize(peptides_rv(), n_iter=input$n_iter, seed=123)
  })

  pep_filt <- reactive({
    req(peptides_rv())
    peptides_rv() %>%
      filter(ph_match>=input$ph_min, mpnn_score>=input$mpnn_min,
             binding_score>=input$bind_r[1], binding_score<=input$bind_r[2],
             docking_score>=input$dock_r[1], docking_score<=input$dock_r[2],
             solubility>=input$sol_min, toxicity<=input$tox_max,
             stage %in% input$stage_sel)
  })

  # ── FASTA STATUS ─────────────────────────────────────────
  output$fasta_status_box <- renderUI({
    req(fasta_valid())
    v   <- fasta_valid()
    col <- if (v$ok) "success" else "danger"
    ic  <- if (v$ok) "check-circle" else "times-circle"
    fp  <- fasta_parsed()
    box(title=tagList(icon(ic)," FASTA Validation"), width=4, status=col, solidHeader=TRUE,
      p(v$msg),
      if (v$ok && !is.null(fp)) tagList(
        p(strong("Header: "), fp$header),
        p(strong("Length: "), paste0(nchar(fp$sequence)," AA"))
      )
    )
  })

  output$protein_stats_box <- renderUI({
    req(prot_feat(), lm_out(), pocket_out())
    pf <- prot_feat(); lm <- lm_out(); pk <- pocket_out()
    box(title=tagList(tags$span(class="ai-badge","Agent 1+2"),
                      " Protein Profile + Pocket Prediction"),
        width=8, status="warning", solidHeader=TRUE,
      fluidRow(
        column(3,div(class="ph-card",div(class="ph-label","MW (kDa)"),div(class="ph-value",pf$mw_kda))),
        column(3,div(class="ph-card",div(class="ph-label","Length (AA)"),div(class="ph-value",pf$n_aa))),
        column(3,div(class="ph-card",div(class="ph-label","Druggability"),
                     div(class="ph-value",style="color:#8e44ad;",lm$druggability))),
        column(3,div(class="ph-card",div(class="ph-label","Pockets Found"),
                     div(class="ph-value",pk$n_pockets)))
      ),
      fluidRow(
        column(3,div(class="ph-card",div(class="ph-label","pI (est.)"),div(class="ph-value",pf$pI))),
        column(3,div(class="ph-card",div(class="ph-label","Net Charge"),div(class="ph-value",pf$net_charge))),
        column(3,div(class="ph-card",div(class="ph-label","Best Pocket Vol. (Å³)"),
                     div(class="ph-value",pk$best_pocket$volume_A3))),
        column(3,div(class="ph-card",div(class="ph-label","Inferred Class"),
                     div(class="ph-value",style="font-size:15px;",pf$target_class)))
      ),
      p(style="font-size:11px;color:#888;margin:6px 0 0;",
        strong("Key pocket residues: "), pk$best_pocket$key_residues)
    )
  })

  output$aa_freq_box <- renderUI({
    req(prot_feat())
    box(title="Residue Composition of Input Protein",width=12,status="primary",solidHeader=TRUE,
        plotlyOutput("aa_freq_prot",height="220px"))
  })

  output$aa_freq_prot <- renderPlotly({
    req(prot_feat())
    df <- prot_feat()$aa_freq; colnames(df) <- c("AA","Freq"); df$Freq <- as.numeric(df$Freq)
    df$role <- dplyr::case_when(
      df$AA %in% c("F","Y","W","H") ~ "Aromatic",
      df$AA %in% c("A","V","L","I","M","P") ~ "Hydrophobic",
      df$AA %in% c("D","E") ~ "Neg. Charged",
      df$AA %in% c("K","R") ~ "Pos. Charged",
      TRUE ~ "Polar/Other"
    )
    p <- ggplot(df,aes(x=AA,y=Freq,fill=role,text=paste0(AA,": ",round(Freq*100,1),"%")))+
      geom_col()+scale_y_continuous(labels=scales::percent)+
      labs(x=NULL,y="Fraction",fill="Role")+
      theme_minimal(base_size=12)+theme(legend.position="bottom",legend.text=element_text(size=9))
    ggplotly(p,tooltip="text")%>%layout(legend=list(orientation="h",y=-0.35))
  })

  # ── AGENT GRAPH ───────────────────────────────────────────
  output$agent_graph <- renderVisNetwork({
    nodes <- make_agent_nodes(); edges <- make_agent_edges()
    visNetwork(nodes,edges) %>%
      visNodes(shape="box",
               font=list(size=15,face="Arial Bold",color="#2c3e50",bold=list(size=15)),
               margin=list(top=14,bottom=14,left=18,right=18),
               borderWidth=2.5, shadow=list(enabled=TRUE,size=8),
               widthConstraint=list(minimum=140,maximum=170),
               heightConstraint=list(minimum=58)) %>%
      visEdges(color=list(color="#bdc3c7",highlight="#e74c3c"),
               smooth=list(type="curvedCW",roundness=0.2),
               width=2.5, font=list(size=10,color="#888"),
               arrows=list(to=list(enabled=TRUE,scaleFactor=1.1))) %>%
      visGroups(groupname="lm",       color=agent_group_cols$lm) %>%
      visGroups(groupname="pocket",   color=agent_group_cols$pocket) %>%
      visGroups(groupname="ligand",   color=agent_group_cols$ligand) %>%
      visGroups(groupname="generator",color=agent_group_cols$generator) %>%
      visGroups(groupname="docking",  color=agent_group_cols$docking) %>%
      visGroups(groupname="optimizer",color=agent_group_cols$optimizer) %>%
      visHierarchicalLayout(direction="LR",levelSeparation=210,nodeSpacing=120) %>%
      visPhysics(enabled=FALSE) %>%
      visInteraction(dragNodes=FALSE,zoomView=TRUE) %>%
      visEvents(stabilized="function(){this.fit({animation:{duration:500}});}")
  })

  # ── AGENT CARDS ───────────────────────────────────────────
  mk_agent_card <- function(n, color, icon_name, title, subtitle, body_ui) {
    div(class="agent-card", style=paste0("border-top-color:",color,";"),
      div(style="display:flex;align-items:center;gap:10px;margin-bottom:10px;",
        div(style=paste0("background:",color,";width:36px;height:36px;border-radius:50%;",
            "display:flex;align-items:center;justify-content:center;color:#fff;font-size:16px;"),
            icon(icon_name)),
        div(h4(paste0("Agent ",n," — ",title)),
            div(class="agent-subtitle", subtitle))
      ),
      body_ui
    )
  }

  output$agent1_card <- renderUI({
    if (is.null(lm_out())) {
      return(mk_agent_card(1,"#8e44ad","brain","Protein Language Model","ESM-2 650M (simulated)",
        p(class="section-info","Run the pipeline to activate.")))
    }
    lm <- lm_out()
    mk_agent_card(1,"#8e44ad","brain","Protein Language Model","ESM-2 650M (simulated)",
      tagList(
        fluidRow(
          column(4,div(class="agent-metric",span("Druggability"),
                       paste0(round(lm$druggability*100,1),"%"))),
          column(4,div(class="agent-metric",span("Confidence"),
                       paste0(round(lm$confidence*100,1),"%"))),
          column(4,div(class="agent-metric",span("Class"),lm$recommended_class))
        ),
        div(style="margin-top:10px;",
          p(style="font-size:11px;color:#888;margin-bottom:4px;",strong("Embedding axes:")),
          div(style="display:flex;flex-wrap:wrap;gap:4px;",
            lapply(names(lm$embedding), function(k)
              div(class="agent-metric",span(k), round(lm$embedding[[k]],3))
            )
          )
        ),
        div(style="margin-top:8px;",
          p(style="font-size:11px;color:#888;margin-bottom:4px;",strong("Predicted 2° structure:")),
          lapply(names(lm$secondary_structure), function(k)
            tags$span(class="agent-metric",span(k), scales::percent(lm$secondary_structure[[k]]))
          )
        )
      )
    )
  })

  output$agent2_card <- renderUI({
    if (is.null(pocket_out())) {
      return(mk_agent_card(2,"#16a085","search","Pocket Prediction","AlphaFold2 / DeepSite (simulated)",
        p(class="section-info","Run the pipeline to activate.")))
    }
    pk <- pocket_out()
    bp <- pk$best_pocket
    mk_agent_card(2,"#16a085","search","Pocket Prediction","AlphaFold2 / DeepSite (simulated)",
      tagList(
        fluidRow(
          column(3,div(class="agent-metric",span("Pockets"),pk$n_pockets)),
          column(3,div(class="agent-metric",span("Best Vol."),paste0(bp$volume_A3," Å³"))),
          column(3,div(class="agent-metric",span("Druggability"),bp$druggability)),
          column(3,div(class="agent-metric",span("Depth"),paste0(bp$depth_A," Å")))
        ),
        p(style="font-size:11px;color:#888;margin:8px 0 2px;",strong("Key pocket residues:")),
        p(style="font-size:12px;font-family:monospace;color:#2c3e50;", bp$key_residues),
        p(style="font-size:11px;color:#888;margin:6px 0 2px;",strong("All pockets:")),
        div(style="display:flex;flex-wrap:wrap;gap:4px;",
          lapply(pk$pockets, function(p)
            div(class="agent-metric",
              span(p$pocket_id), paste0("Vol:",p$volume_A3," D:",p$druggability))
          )
        )
      )
    )
  })

  output$agent3_card <- renderUI({
    if (is.null(ligand_hits())) {
      return(mk_agent_card(3,"#2980b9","database","Ligand Knowledge Agent","GNN / ChemBERTa (simulated)",
        p(class="section-info","Run the pipeline to activate.")))
    }
    hits <- ligand_hits()
    mk_agent_card(3,"#2980b9","database","Ligand Knowledge Agent","GNN / ChemBERTa (simulated)",
      tagList(
        fluidRow(
          column(4,div(class="agent-metric",span("Hits"),nrow(hits))),
          column(4,div(class="agent-metric",span("Best pKi"),max(hits$binding_affinity))),
          column(4,div(class="agent-metric",span("Top Class"),
                       names(sort(table(hits$target_class),decreasing=TRUE))[1]))
        ),
        p(style="font-size:11px;color:#888;margin:8px 0 2px;",strong("Top 3 ligands:")),
        div(style="display:flex;flex-wrap:wrap;gap:4px;",
          lapply(seq_len(min(3,nrow(hits))), function(i)
            div(class="agent-metric",
              span(hits$name[i]),
              paste0("Sim:",round(hits$embed_sim[i],3)," pKi:",hits$binding_affinity[i]))
          )
        )
      )
    )
  })

  output$agent4_card <- renderUI({
    if (is.null(peptides_rv())) {
      return(mk_agent_card(4,"#d35400","flask","Peptide Generator","ProteinMPNN / RFdiffusion (simulated)",
        p(class="section-info","Run the pipeline to activate.")))
    }
    peps <- peptides_rv()
    best <- peps[1,]
    mk_agent_card(4,"#d35400","flask","Peptide Generator","ProteinMPNN / RFdiffusion (simulated)",
      tagList(
        fluidRow(
          column(3,div(class="agent-metric",span("Generated"),nrow(peps))),
          column(3,div(class="agent-metric",span("Best MPNN"),max(peps$mpnn_score))),
          column(3,div(class="agent-metric",span("Best Ph.Match"),max(peps$ph_match))),
          column(3,div(class="agent-metric",span("Best Binding"),max(peps$binding_score)))
        ),
        p(style="font-size:11px;color:#888;margin:8px 0 2px;",strong("Top candidate:")),
        p(style="font-family:monospace;font-size:13px;letter-spacing:1px;color:#2c3e50;",
          best$sequence),
        p(style="font-size:11px;color:#888;",
          paste0("Length: ",best$length," AA  |  MW: ",best$mw_kda," kDa  |  Charge: ",best$net_charge))
      )
    )
  })

  output$agent5_card <- renderUI({
    if (is.null(peptides_rv())) {
      return(mk_agent_card(5,"#c0392b","atom","Docking Evaluator","AutoDock-GPU / Glide (simulated)",
        p(class="section-info","Run the pipeline to activate.")))
    }
    peps <- peptides_rv()
    mk_agent_card(5,"#c0392b","atom","Docking Evaluator","AutoDock-GPU / Glide (simulated)",
      tagList(
        fluidRow(
          column(3,div(class="agent-metric",span("Docked"),nrow(peps))),
          column(3,div(class="agent-metric",span("Best Dock"),min(peps$docking_score))),
          column(3,div(class="agent-metric",span("Mean RMSD"),round(mean(peps$pose_rmsd),2))),
          column(3,div(class="agent-metric",span("Mean Conf."),round(mean(peps$docking_confidence),3)))
        ),
        p(style="font-size:11px;color:#888;margin:8px 0 2px;",
          strong("Energy breakdown (mean across candidates):")),
        fluidRow(
          column(6,div(class="agent-metric",span("VdW Energy"),round(mean(peps$vdw_energy),2))),
          column(6,div(class="agent-metric",span("Electrostatic"),round(mean(peps$electrostatic),2)))
        )
      )
    )
  })

  output$agent6_card <- renderUI({
    if (is.null(optim_rv())) {
      return(mk_agent_card(6,"#27ae60","rotate","Bayesian Optimizer","BoTorch / REINFORCE (simulated)",
        p(class="section-info","Run the pipeline to activate.")))
    }
    opt <- optim_rv()
    mk_agent_card(6,"#27ae60","rotate","Bayesian Optimizer","BoTorch / REINFORCE (simulated)",
      tagList(
        fluidRow(
          column(3,div(class="agent-metric",span("Iterations"),nrow(opt))),
          column(3,div(class="agent-metric",span("Final Score"),max(opt$best_score))),
          column(3,div(class="agent-metric",span("Ph. Improvement"),
                       paste0("+",round((max(opt$mean_ph)-min(opt$mean_ph))*100,1),"%"))),
          column(3,div(class="agent-metric",span("Last Mutation"),tail(opt$mutation_type,1)))
        ),
        p(style="font-size:11px;color:#888;margin:8px 0 2px;",strong("Mutation strategy log:")),
        div(style="display:flex;flex-wrap:wrap;gap:4px;",
          lapply(opt$mutation_type, function(m)
            tags$span(class="agent-metric",style="font-size:11px;", m)
          )
        )
      )
    )
  })

  # ── BINDINGDB VALUE BOXES ─────────────────────────────────
  output$vb_n_hits    <- renderValueBox({
    req(ligand_hits())
    valueBox(nrow(ligand_hits()),"BindingDB Hits",icon=icon("database"),color="purple")
  })
  output$vb_best_aff  <- renderValueBox({
    req(ligand_hits())
    valueBox(max(ligand_hits()$binding_affinity),"Best pKi",icon=icon("star"),color="yellow")
  })
  output$vb_top_class <- renderValueBox({
    req(ligand_hits())
    valueBox(names(sort(table(ligand_hits()$target_class),decreasing=TRUE))[1],
             "Dominant Class",icon=icon("bullseye"),color="green")
  })

  output$ligand_table <- renderDT({
    req(ligand_hits())
    df <- ligand_hits() %>%
      select(name,target_class,binding_affinity,tanimoto,gnn_boost,embed_sim,
             ph_hbd,ph_hba,ph_arom,ph_hydro,mw_est) %>%
      rename("Ligand"=name,"Target Class"=target_class,"pKi"=binding_affinity,
             "Tanimoto"=tanimoto,"GNN Boost"=gnn_boost,"Embed Sim."=embed_sim,
             "HBD"=ph_hbd,"HBA"=ph_hba,"Arom."=ph_arom,"Hydro."=ph_hydro,"MW (Da)"=mw_est)
    datatable(df,options=list(pageLength=10,scrollX=TRUE),rownames=FALSE) %>%
      formatStyle("Embed Sim.",background=styleColorBar(c(0,1),"#B5D4F4"),
                  backgroundSize="90% 70%",backgroundRepeat="no-repeat",backgroundPosition="left") %>%
      formatRound(c("Tanimoto","GNN Boost","Embed Sim.","pKi"),3)
  })

  output$sim_hist <- renderPlotly({
    req(ligand_hits())
    p <- ggplot(ligand_hits(),aes(x=embed_sim,fill=target_class))+
      geom_histogram(binwidth=0.04,color="white")+
      labs(x="Embedding Similarity",y="Count",fill="Class")+
      theme_minimal(base_size=12)+theme(legend.position="bottom")
    ggplotly(p)
  })
  output$sim_aff <- renderPlotly({
    req(ligand_hits())
    p <- ggplot(ligand_hits(),aes(x=embed_sim,y=binding_affinity,
                                  color=target_class,text=name))+
      geom_point(size=3,alpha=0.85)+
      labs(x="Embedding Similarity",y="pKi")+
      theme_minimal(base_size=12)+theme(legend.position="none")
    ggplotly(p,tooltip=c("text","x","y"))
  })

  # ── PEPTIDE VALUE BOXES ───────────────────────────────────
  output$vb_total_pep <- renderValueBox({
    req(pep_filt())
    valueBox(nrow(pep_filt()),"Candidates",icon=icon("vial"),color="blue")
  })
  output$vb_final_pep <- renderValueBox({
    req(pep_filt())
    valueBox(sum(pep_filt()$stage=="Final Candidate"),
             "Final Candidates",icon=icon("star"),color="yellow")
  })
  output$vb_top_bind  <- renderValueBox({
    req(pep_filt())
    val <- if(nrow(pep_filt())>0) max(pep_filt()$binding_score) else "—"
    valueBox(val,"Best Binding Score",icon=icon("chart-line"),color="green")
  })

  output$peptide_table <- renderDT({
    req(pep_filt())
    df <- pep_filt() %>%
      select(peptide,sequence,length,mw_kda,net_charge,ph_match,mpnn_score,
             binding_score,docking_score,docking_confidence,solubility,toxicity,stage) %>%
      rename("ID"=peptide,"Sequence"=sequence,"Len"=length,"MW"=mw_kda,
             "Charge"=net_charge,"Ph.Match"=ph_match,"MPNN"=mpnn_score,
             "Binding"=binding_score,"Docking"=docking_score,
             "Dock.Conf."=docking_confidence,"Sol."=solubility,"Tox."=toxicity,"Stage"=stage)
    datatable(df,options=list(pageLength=12,scrollX=TRUE),rownames=FALSE,filter="top") %>%
      formatStyle("Ph.Match",background=styleColorBar(c(0,1),"#9FE1CB"),
                  backgroundSize="90% 70%",backgroundRepeat="no-repeat",backgroundPosition="left") %>%
      formatStyle("MPNN",background=styleColorBar(c(0,1),"#CECBF6"),
                  backgroundSize="90% 70%",backgroundRepeat="no-repeat",backgroundPosition="left") %>%
      formatStyle("Tox.",background=styleColorBar(c(0,1),"#F5C4B3"),
                  backgroundSize="90% 70%",backgroundRepeat="no-repeat",backgroundPosition="left") %>%
      formatRound(c("MW","Ph.Match","MPNN","Binding","Docking","Dock.Conf.","Sol.","Tox."),3)
  })

  output$best_card <- renderUI({
    req(peptides_rv())
    best <- peptides_rv() %>% arrange(desc(ph_match),desc(mpnn_score)) %>% head(1)
    if (nrow(best)==0) return(p("Run the pipeline first."))
    aa_col <- function(a) {
      if (a %in% c("F","Y","W","H")) "#534AB7"
      else if (a %in% c("A","V","L","I","M","P")) "#E8593C"
      else if (a %in% c("K","R")) "#185FA5"
      else if (a %in% c("D","E")) "#993C1D"
      else "#0F6E56"
    }
    pills <- lapply(strsplit(best$sequence,"")[[1]], function(a)
      tags$span(class="aa-pill",style=paste0("background:",aa_col(a),";"),a))
    tagList(
      fluidRow(
        column(2,div(class="ph-card",div(class="ph-label","Peptide"),div(class="ph-value",style="font-size:16px;",best$peptide))),
        column(2,div(class="ph-card",div(class="ph-label","Ph. Match"),div(class="ph-value",best$ph_match))),
        column(2,div(class="ph-card",div(class="ph-label","MPNN Score"),div(class="ph-value",best$mpnn_score))),
        column(2,div(class="ph-card",div(class="ph-label","Binding"),div(class="ph-value",best$binding_score))),
        column(2,div(class="ph-card",div(class="ph-label","Docking"),div(class="ph-value",best$docking_score))),
        column(2,div(class="ph-card",div(class="ph-label","Dock. Conf."),div(class="ph-value",best$docking_confidence)))
      ),
      div(style="margin-top:10px;",
        p(strong("Residue key: "),
          tags$span(class="aa-pill",style="background:#534AB7;","Aromatic"),
          tags$span(class="aa-pill",style="background:#E8593C;","Hydrophobic"),
          tags$span(class="aa-pill",style="background:#185FA5;","K/R (+)"),
          tags$span(class="aa-pill",style="background:#993C1D;","D/E (−)"),
          tags$span(class="aa-pill",style="background:#0F6E56;","Polar")),
        div(style="margin-top:8px;line-height:2.2;",pills)
      ),
      fluidRow(
        column(3,div(class="ph-card",div(class="ph-label","H-Bond Donors"),  div(class="ph-value",best$n_hbd))),
        column(3,div(class="ph-card",div(class="ph-label","H-Bond Acceptors"),div(class="ph-value",best$n_hba))),
        column(3,div(class="ph-card",div(class="ph-label","Aromatic Res."),  div(class="ph-value",best$n_aromatic))),
        column(3,div(class="ph-card",div(class="ph-label","Hydrophobic Res."),div(class="ph-value",best$n_hydrophobic)))
      )
    )
  })

  # ── OPTIMIZATION LOOP PLOTS ───────────────────────────────
  output$optim_plot <- renderPlotly({
    if (is.null(optim_rv())) {
      p <- ggplot(data.frame(x=1:5, y=c(5,6,7,8,9)),aes(x=x,y=y))+
        annotate("text",x=3,y=7,label="Run the pipeline to see optimization trajectory",
                 color="#888",size=5,hjust=0.5)+
        theme_void()
      return(ggplotly(p))
    }
    opt <- optim_rv()
    # Ensure best_score is numeric and valid
    opt$best_score <- as.numeric(opt$best_score)
    if (all(is.na(opt$best_score))) {
      p <- ggplot()+annotate("text",x=0,y=0,label="No score data available",color="#888",size=5)+theme_void()
      return(ggplotly(p))
    }
    # Use plotly directly to avoid ggplotly RGB conversion issues with shape=21
    hover_txt <- paste0("Iteration: ", opt$iteration,
                        "<br>Best Score: ", round(opt$best_score, 3),
                        "<br>Mutation: ", opt$mutation_type)
    plot_ly(opt, x=~iteration) %>%
      add_ribbons(ymin=~best_score*0.97, ymax=~best_score,
                  fillcolor="rgba(39,174,96,0.12)", line=list(color="transparent"),
                  name="Range", hoverinfo="none", showlegend=FALSE) %>%
      add_lines(y=~best_score,
                line=list(color="#27ae60", width=2.5),
                name="Best Score", hoverinfo="none", showlegend=FALSE) %>%
      add_markers(y=~best_score,
                  marker=list(color="#27ae60", size=11,
                               line=list(color="#fff", width=2.5)),
                  text=hover_txt, hoverinfo="text",
                  name="Score", showlegend=FALSE) %>%
      layout(
        xaxis=list(title="Optimization Iteration", tickvals=opt$iteration,
                   showgrid=TRUE, gridcolor="#f0f0f0"),
        yaxis=list(title="Best Binding Score", showgrid=TRUE, gridcolor="#f0f0f0"),
        hovermode="x unified",
        plot_bgcolor="#fff", paper_bgcolor="#fff",
        margin=list(l=50,r=20,t=20,b=50)
      )
  })

  output$optim_table <- renderUI({
    if (is.null(optim_rv()))
      return(div(style="text-align:center;padding:20px;color:#bbb;font-size:13px;",
                 "Run the pipeline to see the iteration log."))
    opt <- optim_rv()
    div(
      lapply(seq_len(nrow(opt)), function(i) {
        row <- opt[i,]
        pct <- min(1,(row$best_score - min(opt$best_score)) /
                     max(0.01,max(opt$best_score)-min(opt$best_score)))
        div(class="optim-row",
          div(class="optim-iter", paste0("Iter ",row$iteration)),
          div(class="optim-bar",
            div(class="optim-fill",style=paste0("width:",round(pct*100),"%;")),
          ),
          div(style="min-width:40px;font-size:12px;font-weight:700;color:#27ae60;",
              row$best_score),
          div(style="font-size:10px;color:#888;",row$mutation_type)
        )
      })
    )
  })

  output$optim_ph <- renderPlotly({
    if (is.null(optim_rv())) return(plot_ly() %>% layout(xaxis=list(visible=FALSE),yaxis=list(visible=FALSE)))
    p <- ggplot(optim_rv(),aes(x=iteration,y=mean_ph))+
      geom_col(fill="#9FE1CB",color="#0F6E56",linewidth=0.5)+
      labs(x="Iteration",y="Mean Pharmacophore Match")+
      theme_minimal(base_size=12)
    ggplotly(p)
  })

  output$optim_pool <- renderPlotly({
    if (is.null(optim_rv())) return(plot_ly() %>% layout(xaxis=list(visible=FALSE),yaxis=list(visible=FALSE)))
    p <- ggplot(optim_rv(),aes(x=iteration,y=n_candidates))+
      geom_area(fill="#FAC775",alpha=0.5,color="#854F0B",linewidth=0.8)+
      geom_point(color="#854F0B",size=3)+
      labs(x="Iteration",y="Remaining Candidates")+
      theme_minimal(base_size=12)
    ggplotly(p)
  })

  # ── ANALYTICS PLOTS ───────────────────────────────────────
  output$ph_bind <- renderPlotly({
    req(pep_filt())
    p <- ggplot(pep_filt(),aes(x=ph_match,y=binding_score,color=stage,
                               text=paste0(peptide,"\n",sequence)))+
      geom_point(alpha=0.75,size=2.5)+
      geom_smooth(method="lm",se=FALSE,color="#333",linewidth=0.4)+
      labs(x="Pharmacophore Match",y="Binding Score")+
      theme_minimal(base_size=12)+theme(legend.position="none")
    ggplotly(p,tooltip="text")
  })
  output$mpnn_dock <- renderPlotly({
    req(pep_filt())
    p <- ggplot(pep_filt(),aes(x=mpnn_score,y=docking_score,color=docking_confidence,
                               text=paste0(peptide,"\nMPNN: ",mpnn_score,"\nDock: ",docking_score)))+
      geom_point(alpha=0.8,size=2.5)+
      scale_color_gradient(low="#FAC775",high="#0F6E56")+
      labs(x="MPNN Score",y="Docking Score",color="Confidence")+
      theme_minimal(base_size=12)
    ggplotly(p,tooltip="text")
  })
  output$bar_stage <- renderPlotly({
    req(pep_filt())
    df <- pep_filt() %>% count(stage)
    p  <- ggplot(df,aes(x=reorder(stage,n),y=n,fill=stage,text=paste0(stage,": ",n)))+
      geom_col(show.legend=FALSE)+coord_flip()+
      labs(x=NULL,y="Count")+theme_minimal(base_size=12)
    ggplotly(p,tooltip="text")
  })
  output$sol_tox <- renderPlotly({
    req(pep_filt())
    p <- ggplot(pep_filt(),aes(x=solubility,y=toxicity,color=stage,
                               text=paste0(peptide,"\nSol:",solubility," Tox:",toxicity)))+
      geom_point(alpha=0.7,size=2.5)+
      geom_hline(yintercept=0.5,linetype="dashed",color="red",linewidth=0.4)+
      geom_vline(xintercept=0.3,linetype="dashed",color="steelblue",linewidth=0.4)+
      labs(x="Solubility",y="Toxicity")+
      theme_minimal(base_size=12)+theme(legend.position="none")
    ggplotly(p,tooltip="text")
  })
  output$aa_freq_plot <- renderPlotly({
    req(pep_filt()); if(nrow(pep_filt())==0) return(NULL)
    top10 <- pep_filt() %>% arrange(desc(ph_match)) %>% head(10)
    all_aa<- strsplit(paste(top10$sequence,collapse=""),"")[[1]]
    df    <- as.data.frame(table(AA=all_aa)) %>% arrange(desc(Freq))
    df$role <- dplyr::case_when(
      df$AA %in% c("F","Y","W","H") ~ "Aromatic",
      df$AA %in% c("A","V","L","I","M","P") ~ "Hydrophobic",
      df$AA %in% c("D","E") ~ "Neg. Charged",
      df$AA %in% c("K","R") ~ "Pos. Charged",
      TRUE ~ "Polar/Other")
    p <- ggplot(df,aes(x=reorder(AA,-Freq),y=Freq,fill=role,text=paste0(AA,": ",Freq)))+
      geom_col()+labs(x="Amino Acid",y="Count in Top-10",fill="Role")+
      theme_minimal(base_size=12)+theme(legend.position="bottom")
    ggplotly(p,tooltip="text")%>%layout(legend=list(orientation="h",y=-0.35))
  })
}

# ============================================================
# G — LAUNCH
# ============================================================
shinyApp(ui=ui, server=server)
