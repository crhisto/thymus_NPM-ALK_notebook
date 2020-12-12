
## R Notebook with scripts for running the bulk deconvolution and simulations.

You can find the following files in the repository:
1. **R Notebook with examples:** [thymus_NPM-ALK_notebook.Rmd](thymus_NPM-ALK_notebook.Rmd)
2. **Html with examples:** [thymus_NPM-ALK_notebook.html](https://htmlpreview.github.io/?https://github.com/crhisto/thymus_NPM-ALK_notebook/blob/master/thymus_NPM-ALK_notebook.html)
2. **Additional R functions:** [generic_functions.R](scripts/generic_functions.R)
3. **R file with the Bulk data object:** [eset.thymus.bulk.sparse.RData.zip](data/eset.thymus.bulk.sparse.RData.zip)

## Libraries 

In order to execute the complete pipeline, at the beginning of the R Notebook the following modified libraries must be installed:
1. *SCDC*:    
   - **Repository**: https://github.com/crhisto/SCDC
   - **Modifications**: 
     - Compatibility with sparse matrices using: `dgCMatrix` objects in R.
     - Routines with parallelization
     - Dynamic threshold for markers selection
     - Improvements in logs and so on.
   - **Original repository**: https://github.com/meichendong/SCDC
         
2. *Biobase*: 
   - **Repository**: https://github.com/crhisto/Biobase.
   - **Modifications**: Compatibility with sparse matrices using: `dgCMatrix` objects in R.
   - **Original repository**: https://github.com/Bioconductor/Biobase
            
3. *Xbioc*:   
   - **Repository**: https://github.com/crhisto/xbioc. Original repository: 
   - **Modifications**: Compatibility with sparse matrices using: `dgCMatrix` objects in R.
   - **Original repository**: https://github.com/renozao/xbioc
            
Publication: [Requirement of DNMT1 to orchestrate epigenomic reprogramming for NPM-ALK–driven lymphomagenesis](https://www.doi.org/10.26508/lsa.202000794)
