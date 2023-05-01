import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import pandas as pd
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

args = sys.argv

input_dir_prefix = args[1]
input_dir_suffix = args[2]
samples          = args[3:]



for sample in samples:
# get number of cells cellrangers called per sample
  sample_summary = pd.read_csv(input_dir_prefix + sample +  '/outs/metrics_summary.csv')
  n_cells = sample_summary['Estimated Number of Cells'].str.replace(',','')
  n_cells = int(n_cells)
  if n_cells <= 500:
      multiplet_rate = 0.004
  elif n_cells > 500 and n_cells < 2000:
      multiplet_rate = 0.008 
  elif n_cells >= 2000 and n_cells < 3000:
      multiplet_rate = 0.016 
  elif n_cells >= 3000 and n_cells < 4000:
      multiplet_rate = 0.024
  elif n_cells >= 4000 and n_cells < 5000:
      multiplet_rate = 0.032
  elif n_cells >= 5000 and n_cells < 6000:
      multiplet_rate = 0.040
  elif n_cells >= 6000 and n_cells < 7000:
      multiplet_rate = 0.048
  elif n_cells >= 7000 and n_cells < 8000:
      multiplet_rate = 0.056
  elif n_cells >= 8000 and n_cells < 9000:
      multiplet_rate = 0.064
  else:
     multiplet_rate = 0.072

for sample in samples:
  input_dir = input_dir_prefix + sample + input_dir_suffix
  print("reading from: ", input_dir)
  counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
  print("Processing ", sample)
  scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=multiplet_rate)
  doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                          min_cells=3,
                                                          min_gene_variability_pctl=85,
                                                          n_prin_comps=30)
  np.savetxt(sample + "_scrublet.score", doublet_scores)
  np.savetxt(sample + "_scrublet.logic", predicted_doublets)

