# Contributing to PepAI

Thank you for your interest in contributing to PepAI! This document outlines the process for contributing code, documentation, and ideas.

## Ways to Contribute

- 🐛 **Bug reports** — Open an issue with a reproducible example
- 💡 **Feature requests** — Open an issue describing the use case
- 🔬 **Real agent integrations** — Replace a simulated agent with a real model
- 📖 **Documentation** — Improve scientific explanations or usage guides
- 🎨 **UI enhancements** — Improve the Shiny dashboard design
- 🧪 **Tests** — Add unit tests for agent functions

## Development Workflow

### 1. Fork & Clone
```bash
git clone https://github.com/YOUR_USERNAME/PepAI.git
cd PepAI
```

### 2. Create a Feature Branch
```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/bug-description
```

### 3. Make Your Changes
- Follow tidyverse R style conventions
- Add inline comments for complex logic
- Update `docs/` if changing architecture

### 4. Test Locally
```r
# Always test with the BCR-ABL example before submitting
shiny::runApp("app.R")
# Click "Load BCR-ABL Example" → "Run AI Agent Pipeline ▶"
# Verify all 8 tabs render without errors
```

### 5. Commit & Push
```bash
git add .
git commit -m "feat: add real ESM-2 agent integration"
git push origin feature/your-feature-name
```

### 6. Open a Pull Request
- Write a clear description of what changed and why
- Reference any related issues (`Closes #42`)
- Include a screenshot if changing the UI

## Code Style

```r
# Good — tidyverse style, clear variable names
analyse_protein <- function(sequence) {
  aas <- strsplit(sequence, "")[[1]]
  n   <- length(aas)
  
  mw_total <- sum(sapply(aas, function(a) AA_MW[a]), na.rm = TRUE)
  # ...
}

# Avoid — no spacing, cryptic names
analyse_protein<-function(s){aas<-strsplit(s,"")[[1]];n<-length(aas)...}
```

## Commit Message Format
```
type: short description (max 72 chars)

Optional longer explanation.

Closes #issue_number
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`

## Questions?

Open an issue or reach out directly:
- Rik Ganguly — rikgangulybioinfo@gmail.com
- Gautam Ahuja — goutamahuja8387@gmail.com
- Siddhant Poudel — sidblitz1602@gmail.com
