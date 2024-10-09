# hek293aav
Improving HEK293-based AAV-production using GSMMs, and a multi-omics approach.

HEK293 cells are a versatile cell line extensively used in the production of recombinant proteins and viral vectors, notably AAV [5]. Despite their high transfection efficiency and adaptability to various culture conditions, challenges remain in achieving sufficient yields of active viral particles. This study presents a comprehensive multi-omics analysis of two HEK293 strains under good manufacturing practice conditions, focusing on the metabolic and cellular responses during AAV production. The investigation included lipidomic, exometabolomic, and transcriptomic profiling across different conditions and time points. Gsmms were reconstructed for these strains to elucidate metabolic shifts and identify potential bottlenecks in AAV production. Notably, the study revealed significant differences between a high-producing (HP) and a low-producing (LP) HEK293 strains, highlighting pseudohypoxia in the LP strain. Key findings include the identification of hypoxia-inducible factor 1-alpha (HIF1alpha) as a critical regulator in the LP strain, linking pseudohypoxia to poor AAV productivity. Inhibition of HIF1alpha resulted in immediate cessation of cell growth and a 2-fold increase in viral capsid production, albeit with a decreased number of viral genomes, impacting the full-to-empty particle ratio. This suggests that while HIF1alpha inhibition enhances capsid assembly, it simultaneously hampers nucleotide synthesis via the pentose phosphate pathway, necessary for genome packaging.

## Reconstruction approach for HEK293F strain:
Data from Dietmair, et al [1] were used to reconstruct HEK293F-specific GSMMs and validate them with experimental data based on the most recent reconstruction of the human metabolism, Human-GEM [2].
Reconstructions were performed using CORDA [3], according to the authors suggestions (reconstruction_according_to_corda.ipynb) and using an adapted approach (reconstruction_adapted_corda.ipynb). Both reconstructions were compared with each other (reconstruction.comparison.ipynb), simulated growth rates from Quek, et al [4], and experimental growth rates from Dietmair, et al. [1].

Datasets:
  1) [human1 gsmm](https://github.com/SysBioChalmers/Human-GEM/blob/main/model/Human-GEM.xml) [2]
  2) GSE36094_processed_data_filtered_by_detection-pvalue.txt [1]

## Reconstruction of HEK293 strains for AAV production:

A) Biomass composition:
  1) Biomass composition data are stored in hek_data.xlsx
  2) Lipidome data are stored in Resilts_Lipidomics_cells.xlsx
     For conversion of lipids, refmet.csv was used, obtained from https://www.lipidmaps.org/
  3) Plots for biomass composition and lipidome analysis were generated in plots_hek.R

B) Cellular Density:
  1) Data are stored in diameter_aggregates.xlsx
  2) Plot for cellular density comparison was generated in plots_hek.R

C) Calculation of exchange rates:
  1) Measurements from metabolomic analyses are stored in Results_Metabolomics_Media_1.xlsx, gln_depletion.xlsx, and lac_ammonia_glx_fitting.xlsx
  2) specific exchange rates were calculated in exchange_rates_calculation_final.R and stored in specific_exchange_rates.csv
  3) Limiting compounds were identified based on the results from exchange_rates_calculation_final.R
  4) Experimental growth rates were calculated based on biomass data in the script exchange_rates_calculation_final.R and stored as growth_results.csv

D) GSMM reconstruction based on transcriptomic data:
  1) Transcriptomic data were obtained from [6] and are stored in tpm_P.csv for the HP strain and tpm_I.csv for the LP strain. PCA was performed on transcriptomes in the script plots_hek.R
  2) Reconstruction of GSMMs based on the Human-GEM [2] was performed in P_strain_reconstruction.py, and I_strain_reconstruction.py, resulting in a binary reaction matrix (reaction_matrix.csv or reaction_matrix_LP.csv, and reaction_matrix_HP.csv for each individual strain)
  3) Reconstructed models were compared using Logistic PCA [7] in the script plots_hek.R
  4) Exchange rates were added as constraints in the corresponding models to obtain simulated growth rates via FBA. Simulated growth rates were calculated in growth_reconstruction_final.ipynb, and stored in growth_results_w_pred.csv. Comparison between experimental and simulated growth rates was performed in the script plots_hek.R
  5) AAV production reactions were added using the capsid sequences of AAV8 from NCBI, stored in the file sequences_1.fasta. Production envelopes were created in growth_reconstruction_final.ipynb
  6) pFBA, and reaction analysis was performed in growth_reconstruction_final.ipynb

E) Validation experiments
  1) Titers from validation experiments after adding media supplementation, and HIF1alpha inhibitor, are stored in aav_inh_results_final.xlsx
  2) Plots from validation experiments were generated in validation_experiments.ipynb


## References
<a id="1">[1]</a> 
Dietmair, Stefanie and Hodson, Mark P and Quek, Lake-Ee and Timmins, Nicholas E and Gray, Peter and Nielsen, Lars K
A multi-omics analysis of recombinant protein production in Hek293 cells
PLOS One. (2012)

<a id="2">[2]</a> 
Robinson, Jonathan L., et al. 
An atlas of human metabolism.
Science signaling 13.624 (2020)

<a id="3">[3]</a> 
Schultz, Andre and Qutub, Amina A
Reconstruction of tissue-specific metabolic networks using CORDA
PLoS computational biology. 12.3 (2016)

<a id="4">[4]</a> 
Quek, Lake-Ee and Dietmair, Stefanie and Hanscho, Michael and Mart{\'\i}nez, Ver{\'o}nica S and Borth, Nicole and Nielsen, Lars K
Reducing Recon 2 for steady-state flux analysis of HEK cell culture
Elsevier. 184 (2014)

<a id="5">[5]</a> 
Bulcha, Jote T and Wang, Yi and Ma, Hong and Tai, Phillip WL and Gao, Guangping
Viral vector platforms within the gene therapy landscape
Nature Publishing Group. 6 (2021)

<a id="6">[6]</a> 
Pistek, Martina and Kahlig, Carolin-Isabel and Hackl, Matthias and Unterthurner, Sabine and Kraus, Barbara and Grabherr, Reingard and Grillari, Johannes and Hernandez Bort, Juan A
Comprehensive mRNA-sequencing-based characterization of three HEK-293 cell lines during an rAAV production process for gene therapy applications
Wiley Online Library. 18 (2023)

<a id="7">[7]</a> 
Zehetner, Leopold and Szeliova, Diana and Kraus, Barbara and Hernandez Bort, Juan A and Zanghellini, Juergen
Logistic PCA explains differences between genome-scale metabolic models in terms of metabolic pathways
PLOS Computational Biology. 20 (2024)
