# hek293aav
Improving HEK293-based AAV-production using GSMMs, and a multi-omics approach

## Reconstruction approach for HEK293F strain:
Data from Dietmair, et al [1] were used to reconstruct HEK293F-specific GSMMs and validate them with experimental data based on the most recent reconstruction of the human metabolism, Human-GEM [2].
Reconstructions were performed using CORDA [3], according to the authors suggestions (reconstruction_according_to_corda.ipynb) and using an adapted approach (reconstruction_adapted_corda.ipynb). Both reconstructions were compared with each other (reconstruction.comparison.ipynb), simulated growth rates from Quek, et al [4], and experimental growth rates from Dietmair, et al. [1].

Datasets:
  1) [human1 gsmm](https://github.com/SysBioChalmers/Human-GEM/blob/main/model/Human-GEM.xml)
  2) GSE36094_processed_data_filtered_by_detection-pvalue.txt

## Reconstruction of HEK293 strains for AAV production:




## References
<a id="1">[1]</a> 
Dietmair, Stefanie and Hodson, Mark P and Quek, Lake-Ee and Timmins, Nicholas E and Gray, Peter and Nielsen, Lars K
A multi-omics analysis of recombinant protein production in Hek293 cells
PLOS One. 2012; https://doi.org/10.1371/journal.pone.0043394

<a id="2">[2]</a> 
Robinson, Jonathan L., et al. 
An atlas of human metabolism.
Science signaling 13.624 (2020): eaaz1482.

