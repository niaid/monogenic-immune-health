AGE analysis reproducibility notes

This document records exact commands used in-session to reproduce the age-group differential-protein analysis (Python OLS + R limma) on Windows using Miniconda (user install) and a conda R env.

1) Python per-protein OLS (environment used earlier)

# create env and install Python deps
conda create -y -n ihm-py python=3.10
conda activate ihm-py
pip install numpy pandas scipy statsmodels seaborn matplotlib scikit-learn

# run Python script (produces perprotein_ols_results.csv and age_diff_heatmap_ols.png)
python scripts\age_grouping_examples.py

2) Miniconda user install workaround (if user homedir contains non-ASCII)

# Download Miniconda installer (manual or via PowerShell)
# Example PowerShell download (save to %TEMP%):
# Invoke-WebRequest -Uri "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe" -OutFile "$env:TEMP\Miniconda3-latest-Windows-x86_64.exe"

# Install to ASCII-only location (no admin):
# Start-Process -FilePath "$env:TEMP\Miniconda3-latest-Windows-x86_64.exe" -ArgumentList "/InstallationType=JustMe","/AddToPath=0","/RegisterPython=0","/S","/D=C:\\Users\\Public\\Miniconda3" -Wait

# Add conda to PATH for current session (PowerShell):
$env:Path = "C:\Users\Public\Miniconda3\Scripts;C:\Users\Public\Miniconda3\Library\bin;C:\Users\Public\Miniconda3\condabin;" + $env:Path

3) Accept conda Terms of Service (required before some installs)

# Example (PowerShell):
C:\Users\Public\Miniconda3\Scripts\conda.exe tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
C:\Users\Public\Miniconda3\Scripts\conda.exe tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
C:\Users\Public\Miniconda3\Scripts\conda.exe tos accept --override-channels --channel https://repo.anaconda.com/pkgs/msys2

4) Create conda R env and install r-base (conda-forge)

C:\Users\Public\Miniconda3\Scripts\conda.exe create -y -n copilot-r -c conda-forge r-base

Notes: conda will download large packages (gcc/gfortran, mkl, etc.). This may take many minutes and hundreds of MB.

5) R package installation (safe TMP workaround used here)

# Create a safe tmp directory (ASCII path)
mkdir C:\Temp

# Create a small R script: C:\Temp\install_limma.R with contents:
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cran.rstudio.com/")
# BiocManager::install(c("limma"), ask=FALSE)

# Run inside conda env, forcing safe TMP/TEMP (PowerShell -> conda run uses cmd /c):
C:\Users\Public\Miniconda3\Scripts\conda.exe run -n copilot-r cmd /c "set TMP=C:\\Temp&& set TEMP=C:\\Temp&& Rscript C:\\Temp\\install_limma.R"

# Install pheatmap (CRAN) similarly:
C:\Users\Public\Miniconda3\Scripts\conda.exe run -n copilot-r cmd /c "set TMP=C:\\Temp&& set TEMP=C:\\Temp&& Rscript -e \"install.packages('pheatmap', repos='https://cran.rstudio.com/')\""

6) Run the R analysis script (limma pipeline)

# Run with safe TMP/TEMP to avoid non-ASCII homedir issues:
C:\Users\Public\Miniconda3\Scripts\conda.exe run -n copilot-r cmd /c "set TMP=C:\\Temp&& set TEMP=C:\\Temp&& Rscript scripts\\age_grouping_examples.R"

Expected outputs (written to outputs/):
- outputs/limma_results_trend.csv
- outputs/age_diff_heatmap_limma.png
- outputs/perprotein_ols_results.csv (from Python) and outputs/age_diff_heatmap_ols.png (from Python)

7) Cleanup (optional)
- Remove Miniconda installer from %TEMP% and temporary R scripts in C:\Temp
- Optionally remove conda env: conda env remove -n copilot-r

8) Troubleshooting notes
- pandas.qcut: if it errors, check df['age'] for NaN or identical values; use pd.cut with manual bins as fallback.
- R package install failing to create temp directories often stems from non-ASCII homedir. Use safe TMP (C:\Temp) or install Miniconda to ASCII-only path.
- If network/push fails, run git push manually from a machine with network access.

Contact
- For questions about these exact commands or to adjust for your system, reply and I'll help adapt them.
