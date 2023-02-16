---
layout: post
title:  "Competing Risks"
date:   2022-06-15 
tags: [survival analysis, R, literature]
---
 
Main references:

- [Competing Risk Regression Models for Epidemiologic Data](https://academic.oup.com/aje/article/170/2/244/111339)
- [Competing risks in epidemiology: possibilities and pitfalls](https://academic.oup.com/ije/article/41/3/861/829598)

For clinical research and cardiology:

- [Introduction to the Analysis of Survival Data in the Presence of Competing Risks](https://www.ahajournals.org/doi/10.1161/circulationaha.115.017719)
- [Importance of Considering Competing Risks in Time-to-Event Analyses](https://www.ahajournals.org/doi/10.1161/CIRCOUTCOMES.118.004580)
- [Competing risks and the clinical community: irrelevance or ignorance?](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3575691/)

More advanced:

- [More details on Fine and Grey](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7216972/)
- [Estimating adjusted probability](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8209)

Multi-state models:

- [Tutorial in biostatistics: competing risks and multi-state models](https://onlinelibrary.wiley.com/doi/10.1002/sim.2712)

R packages and tutorials:

- [Multi-state models and competing risks](https://cran.r-project.org/web/packages/survival/vignettes/compete.pdf), describing how to evaluate competing risks within the `survival` package
- [SemiCompRisks](https://journal.r-project.org/archive/2019/RJ-2019-038/RJ-2019-038.pdf)
- `cmprsk` package, including basic function for Fine and Grey model