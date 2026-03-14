# tests/test_agents.R
# Unit tests for PepAI agent functions
# Run with: source("tests/test_agents.R")

source("app.R", local = TRUE)

cat("==============================\n")
cat("  PepAI Agent Function Tests  \n")
cat("==============================\n\n")

pass <- 0; fail <- 0

check <- function(name, expr) {
  tryCatch({
    result <- eval(expr)
    if (isTRUE(result)) {
      cat(sprintf("  [PASS] %s\n", name)); pass <<- pass + 1
    } else {
      cat(sprintf("  [FAIL] %s — returned FALSE\n", name)); fail <<- fail + 1
    }
  }, error = function(e) {
    cat(sprintf("  [ERROR] %s — %s\n", name, conditionMessage(e))); fail <<- fail + 1
  })
}

# ── FASTA Parsing ──────────────────────────────────────────────
cat("--- FASTA Parsing ---\n")
valid_fasta <- ">TestProtein\nACDEFGHIKLMNPQRSTVWY"
check("parse_fasta returns ok=TRUE for valid FASTA",
      quote(parse_fasta(valid_fasta)$ok == TRUE))
check("parse_fasta extracts correct sequence length",
      quote(nchar(parse_fasta(valid_fasta)$sequence) == 20))
check("parse_fasta fails on missing header",
      quote(parse_fasta("ACDEFGHIKLMNPQRSTVWY")$ok == FALSE))
check("parse_fasta fails on short sequence",
      quote(parse_fasta(">Test\nACDE")$ok == FALSE))
check("validate_fasta returns ok=TRUE for valid input",
      quote(validate_fasta(valid_fasta)$ok == TRUE))
check("validate_fasta returns ok=FALSE for empty input",
      quote(validate_fasta("")$ok == FALSE))

# ── Protein Analysis ───────────────────────────────────────────
cat("\n--- Protein Analysis ---\n")
test_seq <- "ACDEFGHIKLMNPQRSTVWY"
pf <- analyse_protein(test_seq)
check("analyse_protein returns a list",
      quote(is.list(pf)))
check("n_aa is correct",
      quote(pf$n_aa == 20))
check("mw_kda is positive",
      quote(pf$mw_kda > 0))
check("pI is numeric",
      quote(is.numeric(pf$pI)))
check("frac_hydro in [0,1]",
      quote(pf$frac_hydro >= 0 && pf$frac_hydro <= 1))
check("target_class is non-empty string",
      quote(nchar(pf$target_class) > 0))
check("aa_freq is a data.frame",
      quote(is.data.frame(pf$aa_freq)))

# ── Tanimoto Score ─────────────────────────────────────────────
cat("\n--- Tanimoto Score ---\n")
q <- list(hbd=3, hba=4, aromatic=2, hydro=5)
check("tanimoto_score returns value in [0,1]",
      quote({t <- tanimoto_score(q,3,4,2,5); t >= 0 && t <= 1}))
check("tanimoto_score = 1 for identical vectors",
      quote(tanimoto_score(q, 3, 4, 2, 5) == 1))
check("tanimoto_score = 0 for zero union",
      quote(tanimoto_score(list(hbd=0,hba=0,aromatic=0,hydro=0),0,0,0,0) == 0))

# ── Agent 1 ────────────────────────────────────────────────────
cat("\n--- Agent 1 (Protein LM) ---\n")
lm <- agent1_protein_lm(test_seq, pf)
check("agent1 returns a list",
      quote(is.list(lm)))
check("embedding has 6 named elements",
      quote(length(lm$embedding) == 6))
check("druggability in [0,1]",
      quote(lm$druggability >= 0 && lm$druggability <= 1))
check("secondary_structure sums to ~1",
      quote(abs(sum(lm$secondary_structure) - 1) < 0.05))

# ── Agent 2 ────────────────────────────────────────────────────
cat("\n--- Agent 2 (Pocket Prediction) ---\n")
pk <- agent2_pocket(pf, lm)
check("agent2 returns a list",
      quote(is.list(pk)))
check("n_pockets >= 1",
      quote(pk$n_pockets >= 1))
check("best_pocket has volume_A3",
      quote(!is.null(pk$best_pocket$volume_A3)))
check("best_pocket druggability in [0,1]",
      quote(pk$best_pocket$druggability >= 0 && pk$best_pocket$druggability <= 1))

# ── Agent 3 ────────────────────────────────────────────────────
cat("\n--- Agent 3 (Ligand GNN) ---\n")
ph <- protein_to_pharmacophore(pf)
hits <- agent3_ligand(ph, pf$target_class, top_n = 5)
check("agent3 returns a data.frame",
      quote(is.data.frame(hits)))
check("agent3 returns <= top_n rows",
      quote(nrow(hits) <= 5))
check("embed_sim column present",
      quote("embed_sim" %in% colnames(hits)))

# ── Agent 4 ────────────────────────────────────────────────────
cat("\n--- Agent 4 (Peptide Generator) ---\n")
peps <- agent4_generate(pf, ph, hits, pk, n_pep = 10, seed = 42)
check("agent4 returns a data.frame",
      quote(is.data.frame(peps)))
check("correct number of peptides",
      quote(nrow(peps) == 10))
check("sequence column present",
      quote("sequence" %in% colnames(peps)))
check("ph_match in [0,1]",
      quote(all(peps$ph_match >= 0 & peps$ph_match <= 1)))

# ── Agent 5 ────────────────────────────────────────────────────
cat("\n--- Agent 5 (Docking) ---\n")
dock <- agent5_docking(peps, pk)
check("agent5 returns a data.frame",
      quote(is.data.frame(dock)))
check("docking_confidence column added",
      quote("docking_confidence" %in% colnames(dock)))
check("pose_rmsd column added",
      quote("pose_rmsd" %in% colnames(dock)))

# ── Agent 6 ────────────────────────────────────────────────────
cat("\n--- Agent 6 (Bayesian Optimizer) ---\n")
opt <- agent6_optimize(dock, n_iter = 3, seed = 42)
check("agent6 returns a data.frame",
      quote(is.data.frame(opt)))
check("correct number of iterations",
      quote(nrow(opt) == 3))
check("best_score column present",
      quote("best_score" %in% colnames(opt)))
check("mutation_type column present",
      quote("mutation_type" %in% colnames(opt)))

# ── Summary ────────────────────────────────────────────────────
cat(sprintf("\n==============================\n"))
cat(sprintf("  Results: %d passed, %d failed\n", pass, fail))
cat(sprintf("==============================\n"))
if (fail == 0) cat("  All tests PASSED ✓\n") else cat("  Some tests FAILED ✗\n")
