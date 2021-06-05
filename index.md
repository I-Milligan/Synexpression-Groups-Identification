## Introduction

Synexpression Groups (SGs) consist of genes that are involved in the same biological processes with shared complex regulation of expression that are in proximity spatially. However, the ability to predict gene co-expression from frequent remote chromatin interactions is difficult, even though there is a correlation between chromatin interactions and gene co-expression. There is increasing evidence that suggest these SGs are vital components of the complex multicellular matrix of all biological species.  Nonetheless, there has been little systemic genome-wide research into SGs; the exact means to which they are regulated, their cluster frequency and how widespread they occur within the genome remains unclear.


The two most common measures for gene expression data are Pearson’s correlation coefficient and the Euclidean distance. 

### Integrated Development Environment (IDE) and Versions Used
- For Python Code development **Jupyter Notebooks** running **Python 3.8.5 64-bit** was utilized.
- For R Code development **RStudio** version **1.4.1106** was utilized.

Pearson's Correlatiouubn Coefficient was initially used to measure of the strength of a linear association between two expressions
### Pearson's Correlation
[Link to Code Used](https://i-milligan.github.io/Synexpression-Groups-Identification/pearson)

Euclidean Distance is sensitive to scaling and differences in average expression level.  It was used to validate the findings from the Pearson's Correlation Coefficient findings.
### Euclidean Distance
[Link to Code Used](https://i-milligan.github.io/Synexpression-Groups-Identification/euclidean)

Once synexpression clusters were identified, they were put through the Biomart database to determine each gene’s chromosomal location. The probability of genes being located on the same chromosome was calculated utilizing the Fisher’s exact test (p-value of ≤ 0.05) to determine if there are non-random associations between two variables.
### Fisher's Exact Test
[Link to Code Used](fisher.md)

### Rank changes between Pearson's Correlation and Euclidean distance
[Link to Code Usded](rankchange.md)
