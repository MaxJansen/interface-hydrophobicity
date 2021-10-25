# Steps to assess MaSIF PPI Search performance on interfaces of varying hydrophobicity
Run scripts to annotate hydrophobicity on MaSIF test set.
Will not work for all.
Run scripts to precompute MaSIF test set and compute descriptors.

Use /scripts/06-05-2021-rank-default-testing-hydro.ipynb to match and see the intersection of these sets.

From this list of proteins, select the 100 most hydrophobic ones, 100 the middle ones and the 100 least hydrophobic ones.

Run comparison step.
The slurm commands will yield and output in ```/work/upcorreia/users/jansen/masif_hydrophobicity_test/comparison/masif_ppi_search/masif_descriptors_nn```.

In there, get the number of correctly ranked proteins in the top 1,
top 10 and top 100. From the different slurm outputs with:
```
grep -i "number" slurm-877*.out > all_hydro_slurms.txt
```

Then use this file in a local notebook to make barplots. The plots are based on the average results for the three categories of hydrophobicity.

Make similar plots based on the log output from the compute descriptors step.

Go to:
```
/work/upcorreia/users/jansen/masif_hydrophobicity_test/data/masif_ppi_search/descriptors/sc05/all_feat
```

Run:
```
grep -i ": roc auc" log.txt > roc_auc_for_all.txt
```

Similarly, plot these results locally with a notebook.
