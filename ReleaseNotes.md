# Release notes

## Nigthly
- Function `logpdf` computes now the full *log-density* of the iGMRF as defined in Eq. (3.13) by Rue & Held (2002). In the previous version, only the log-density up to an additive constant were returned. 
- The full conditional logpdf values do not rely on Distributions.jl anymore to improve performance.