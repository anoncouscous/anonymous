# HLAlleleBias

Finding peptide targets that bind to HLAs and trigger an immune response is challenging, and peptide-HLA (pHLA) binding prediction is one of the crucial steps in the development of personalized peptide vaccines. Machine Learning (ML) pHLA binding prediction tools are trained on vast amounts of pHLA binding data. ML predictions are effective in guiding the search for therapeutic peptide targets. However, the use of datasets with imbalanced allele content raises concerns about biased performance toward certain geographic populations. We examine the bias of two ML-based pan-allele pHLA binding affinity predictors. We aim to draw attention to the potential therapeutic consequences of this bias, and we challenge the use of the term "pan-allele" to describe models trained with currently available public datasets.

In this repo you can find all the scripts used to perform the given analysis.

## Datasets and methods

### 1. Mapping HLA alleles to geographic populations and classifying them by income level 

<details>
  <summary>1.1 Allele Frequency Net Database (AFND)  scraping</summary>
  
To get the frequencies of HLA alleles across different populations we query the AFND database to get the most recent allele frequencies. More details [here](https://github.com/RomFas/HLAlleleBias/blob/main/AFND_population_frequencies/README.md).


Code for running this is in the `AFND_population_frequencies` folder. This code allows you to:
- Download the AFND HLA allele frequencies per locus for all available populations
- Combine the AFND frequencies across different loci
- Visualize and analyze the AFND data

The result of this step is a population allele frequency map `AFND_data_locus_all.csv` later used for calculating the population coverage as the source of ground truth allele frequency labels.
</details>

<details>
  <summary>1.2 Classifying geographic populations by their level of income</summary>


Code for running this is in the `WorldBank_Income_levels` folder.  We download the country income levels from the World Bank [here](https://datahelpdesk.worldbank.org/knowledgebase/articles/906519-world-bank-country-and-lending-groups) (current classification by income in XLSX format) and process them.


The result of this step is a `population_income_map.csv` which maps the AFND populations to countries and respective income levels
  
</details>

### 2. Examining data bias

Here we calculate the population coverage of different training datasets (NetMHCpan4.1 and MHCFlurry2.0 - binding affinity and mass spectrometry data). The related scripts are in the `Datasets_population_coverage` folder.

<details>
  <summary>2.1 Download and prepare the datasets</summary>

Find the detailed instructions and notebooks for downloading the datasets inside the `./datasets/` folder [here](https://github.com/RomFas/HLAlleleBias/blob/main/Datasets_population_coverage/datasets/README.md). 

- Prepare the NetMHCpan4.1 datasets 
- Prepare the MHCFlurry2.0 datastes
- Visualize the dataset content 
    
</details>


<details>
  <summary>2.2 Calculate the population coverage</summary>

Use the `./CalculateCoverage_IEDB.ipynb` notebook to calculate the population coverage of the datasets. This notebook will guide you through the following steps:

- Downloading and setting up the IEDB population coverage tool from [here](http://tools.iedb.org/population/download/)
- Running the population coverage for a single population and a single dataset
- Running the population coverage for all datasets and all populations
    
</details>

    

<details>
  <summary>2.3. Visualize the population coverage of the datasets</summary>
  
Run the `./Visualize_results.ipynb` notebook to see how different datasets cover the given populations and how this coverage relates to different income levels and ancestries. 


</details>

### 3. Examining algorithmic bias

To examine the algorithmic bias we use a recently published independent dataset collected by [Pyke et al.](https://doi.org/10.1016/j.mcpro.2023.100506). This dataset contains MS data for some rare alleles that were not sampled before (i.e., A*02:52, B*15:13). More details [here](https://github.com/RomFas/HLAlleleBias/blob/main/Algorithmic_bias_analysis/README.md).

The code for running this analysis is in the folder  `Algorithmic_bias_analysis`. This code allows you to:
- `Prepare_dataset.R`: Filter and prepare the dataset from Pyke et al. used in the analysis. Namely, peptides with any chemical modifications are removed, and negatives are sampled from the human proteome. 
- `Motif_generation.R`: Generate motifs for calculating FOOP scores (see paper for more details)
- `Data_analysis.R`: Calculate PPV/FOOP scores and generate plots 

###

## Additional Resources

### Databases
<details>

  
  - [AFND](http://www.allelefrequencies.net/) Allele frequency net database

  Currated frequencies across regions/ethnicities!! (what we need)
  A couple of resources to access this data
  - [immunotation](https://github.com/imkeller/immunotation) R package
  - [hladownload](https://github.com/ramonamerong/hladownload) python package    

- [IPD-IMGT](https://www.ebi.ac.uk/ipd/imgt/hla/)  Immuno Polymorphism Database,  international ImMunoGeneTics information system 

  Contains raw sample data - if we need more details.

- [IEDB](https://www.iedb.org) The immune epitope database 

  Contains more raw data (on the binding affinity / mass spec side).
  Contains the [population coverage](http://tools.iedb.org/population/) module
  
- [hla.alleles.org](http://hla.alleles.org)
</details>

### Literature
<details>
  
- Works that mention/address the possible bias in the therapeutics.
  
  - [Bui et al.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-153 ), 2005: Predicting population coverage of T-cell epitope-based diagnostics and vaccines
 
    _A disproportionate amount of MHC polymorphism occurs in positions constituting the peptide-binding region, and as a result, MHC molecules exhibit a widely varying binding specificity. In the design of peptide-based vaccines and diagnostics, the issue of population coverage in relation to MHC polymorphism is further complicated by the fact that different HLA types are expressed at dramatically different frequencies in different ethnicities. Thus, without careful consideration, a vaccine or diagnostic with ethnically biased population coverage could result._
    
  - [Oyarzun et al.](https://www.sciencedirect.com/science/article/pii/S0264410X15000845?casa_token=6a7cmf4BWAMAAAAA:V_V4as7_oKy-DKHZoLtMIq7z3s9LTctfMaFN73q3ClZ89b9kbRWSrq8_-3GZHadcsM6bPYD2D_E#bib0305), 2015: A bioinformatics tool for epitope-based vaccine design that accounts for human ethnic diversity: Application to emerging infectious diseases

    _Predivac-2.0 is a novel approach in epitope-based vaccine design, particularly suited to be applied to virus-related emerging infectious diseases, because the geographic distributions of the viruses are well defined and ethnic populations in need of vaccination can be determined (“ethnicity-oriented approach”). Predivac-2.0 is accessible through the website http://predivac.biosci.uq.edu.au/._
    
   - [Sarkizova et al.](https://www.nature.com/articles/s41587-019-0322-9.pdf?proof=t%C2%A0) - HLAAthena

   - [Pyke et al.](https://pubmed.ncbi.nlm.nih.gov/34126241/) 2021: Precision Neoantigen Discovery Using Large-scale Immunopeptidomes and Composite Modeling of MHC Peptide Presentation

</details>
  
### Figure generation
