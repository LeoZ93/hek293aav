This folder contains Jupyter notebooks and data to reconstruct GSMMs of HEK293F cells during exponential growth phase.

Transcriptomic data were obtained from Dietmair, et al. and are stored in GSE36094_processed_data_filtered_by_detection-pvalue.txt
GSMM reconstruction was based on the most recent global reconstruction of the human metabolism, Human-Gem, obtained from Robinson, et al.
HEK293F-specific GSMMs were generated using CORDA in two approaches:
  - first approach was according to the authors suggestions (reconstruction_according_to_corda.ipynb)
  - second approach was adapted from the original suggestions (reconstruction_adapted_corda.ipynb)
Results from both approaches were used to reconstruct growth rates and compare to experimental data (reconstruction_comparison.ipynb)
