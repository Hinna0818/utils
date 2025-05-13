# 🔧 Bioinformatics Utility Functions (`utils`)
Some light-weight functions in R or Python.


## 🧰 Ready-to-go functions
### 🧬 Gene & ID Conversion
- `convert_mouse_to_human()`：MM_Symbol → HS_Symbol
- `convert_human_to_mouse()`：HS_Symbol → SS_Symbol
- Support NCBI Entrez ID and Bioconductor database（`org.Mm.eg.db`, `org.Hs.eg.db`）

### 🧮 Epidemiology Model Metrics
- `ORCI()`：calculate OR(odds ratio) and 95% CI of logistic or cox model
- `PAR()`：calculate PAR of the model
- `PFP`：calculate PFP(preventable fraction for the study population) of the model
