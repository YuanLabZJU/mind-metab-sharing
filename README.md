# mind-metab-sharing
Code repository for the project 'Circulating metabolomic profile of the MIND diet and its relation to cognition in middle-aged and older adults'. 


## Descriptions

This code repository documents the code scripts for the project 'Circulating metabolomic profile of the MIND diet and its relation to cognition in middle-aged and older adults'. 

## Dataset required for analysis

- UK Biobank: Data are available on request from https://www.ukbiobank.ac.uk.

- Whitehall II: Data are available on request for bone fide investigators from https://www.ucl.ac.uk/epidemiology-health-care/research/epidemiology-and-public-health/research/whitehall-ii.

## File structure

- /Analysis: Code scripts for data analysis

    - readin.R: Read in datasets
    
    - mind-mb.R: Associations between aMIND and metabolites
    
    - mind-mets.R: Construction of the MIND-MetS
    
    - mediation.R: Mediation on cognitive function (WHII) and dementia (UKB)

- /Supplementary Tables.xlsx: Meta data for figure repruduction of figures

- /Figures: Code scripts for data visulization

    - Figure 1
    
        - Fig1B.R: Associations between aMIND and metabolites
        
        - Fig1C.R: Associations between food groups and metabolites
        
    - Figure 2
    
        - Fig2A.R: Weight shrinkage in the elastic net
        
        - Fig2B.R: Weights of metabolites in the elastic net model
        
        - Fig2CDE.R: Correlations between aMIND and MIND-MetS in the UKB-Training, UKB-Validation, and WHII.

## Citation

Please cite this paper with:

    Chen, H.; Shen, J.; Tao, Y.; et al. Circulating metabolomic profile of the MIND diet and its relation to cognition in middle-aged and older adults. iMetaOmics. in press
