# Bioinformatic analysis of metabarcoding data for soil invertebrates

## 3. Diversity measures and statistical analyses

### 3.1. Alpha diversity
Alpha diversity refers to the diversity within a single sample or community, typically measured as the number and relative abundance of taxa present. Before computing diversity metrics we must deal with the fact that different sequencing depths in metabarcoding datasets can bias computation, due to unequal sampling effort. So metabarcoding datasets must be normalized somehow before calculating alpha diversity. A quick way for doing it (even if probably not the best way) is to randomly subsample the samples to the same depth (rarefaction). Usually the depth of the smallest sample is used. Let's have a look to sequencing depth variation in our samples.

```bash
qiime feature-table summarize \
  --i-table invertebrates_table_clean.qza \
  --m-metadata-file metadata.tsv \
  --o-feature-frequencies feature-frequencies_clean.qza \
  --o-sample-frequencies sample-frequencies_clean.qza \
  --o-summary invertebrates_table_clean.qzv
```
Let's visualize [invertebrates_table_clean.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/invertebrates_table_clean.qzv).

First we can create a new directory for storing diversity results.
```bash
mkdir -p diversity_results
```

To evaluate the impact of rarefaction on diversity estimates we can use the command `qiime diversity alpha-rarefaction` that generates interactive alpha rarefaction curves. The main parameters are:
- `--p-min-depth`, the minimum rarefaction depth to be tested
- `--p-max-depth`, the maximum rarefaction depth to be tested
- `--p-steps`, the number of rarefaction depths to include between `min-depth` and `max-depth`
- `--p-iterations`, the number of rarefied feature tables to compute at each step

```bash
qiime diversity alpha-rarefaction \
  --i-table invertebrates_table_clean.qza \
  --p-min-depth 1 \
  --p-max-depth 25000 \
  --p-steps 25 \
  --p-iterations 10 \
  --m-metadata-file metadata.tsv \
  --o-visualization diversity_results/alpha_rarefaction.qzv
```
Let's visualize [alpha_rarefaction.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/alpha_rarefaction.qzv).


Now that we have chosen a rarefaction depth we can randomly subsample the ASV table using the command `qiime feature-table rarefy`. The main parameters are:
- `--p-sampling-depth`, the total frequency that each sample should be rarefied to
- `--p-random-seed`, seed for random subsampling, for reproducibility purpose
```bash
qiime feature-table rarefy \
  --i-table invertebrates_table_clean.qza \
  --p-sampling-depth 8124 \
  --p-random-seed 1234 \
  --o-rarefied-table invertebrates_table_clean_raref.qza
```

Then we can use the rarefied table to compute alpha diversity. We will calculate three classic alpha diversity indices. 

- **Richness** ($S$, `observed_features`)

$$
S
$$

- **Shannon index** ($H'$, `shannon`)

$$
H' = -\sum_{i=1}^{S} p_i \ln p_i
$$

- **Simpson index** ($D$, `simpson`)

$$
D = 1 - \sum_{i=1}^{S} p_i^2
$$


$S =$ Number of taxa (ASVs in our case) observed in the sample

$p_i =$ Relative abundance of taxon $i$ in the sample

To perform diversity index calculation we can use `qiime diversity alpha` in a `for` loop as follow.
```bash
for div in observed_features shannon simpson; do
  qiime diversity alpha \
    --i-table invertebrates_table_clean_raref.qza \
    --p-metric $div \
    --o-alpha-diversity diversity_results/alpha_${div}_vector.qza
done
```

Now that we have some alpha diversity indices we can use them to investigate diversity variation in our data. First let's se if different groups of samples have different diversity values. We can use the command `qiime diversity alpha-group-significance` for comparing diversity between the levels of each categorical variable in our metadata.

```bash
for div in observed_features shannon simpson; do
  qiime diversity alpha-group-significance \
    --i-alpha-diversity diversity_results/alpha_${div}_vector.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization diversity_results/alpha_${div}_significance.qzv
done
```
Let's have a look at these visualizations: [alpha_observed_features_significance.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/alpha_observed_features_significance.qzv), [alpha_shannon_significance.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/alpha_shannon_significance.qzv), [alpha_simpson_significance.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/alpha_simpson_significance.qzv).


We can also compute correlations between diversity indices and numeric variables in our metadata using the command `qiime diversity alpha-group-significance`.
```bash
for div in observed_features shannon simpson; do
  qiime diversity alpha-correlation \
    --i-alpha-diversity diversity_results/alpha_${div}_vector.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization diversity_results/alpha_${div}_correlation.qzv
done
```
Let's have a look at these visualizations: [alpha_observed_features_correlation.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/alpha_observed_features_correlation.qzv), [alpha_shannon_correlation.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/alpha_shannon_correlation.qzv), [alpha_simpson_correlation.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/alpha_simpson_correlation.qzv).

### 3.2. Beta diversity

Beta diversity refers to the differences in species composition between samples or communities, typically measured as the degree of turnover or dissimilarity in taxa among them. Traditional beta diversity metrics are often not optimal for metabarcoding data because these datasets are compositional and often very sparse. As a result, distance-based approaches (e.g., Bray–Curtis) may be influenced by technical variation rather than true ecological differences. For this reason, alternative methods specifically designed to analyze compositional sequencing data in a more robust and interpretable way have been developed. The [q2-gemelli](https://github.com/biocore/gemelli) plugin contains several of this methods. First let's install the package since it is not contained in the qiime2 base distribution.
```bash
pip install gemelli
qiime dev refresh-cache
```

In this case, we use the Robust Aitchison PCA (RPCA), which can be seen as a PCA robust to the compositionality and sparsity of the data.
```bash
qiime gemelli rpca \
  --i-table invertebrates_table_clean.qza \
  --p-n-components 3 \
  --p-min-sample-count 1000 \
  --p-min-feature-count 30 \
  --p-min-feature-frequency 0 \
  --p-max-iterations 5 \
  --o-biplot diversity_results/RPCA_biplot.qza \
  --o-distance-matrix diversity_results/RPCA_distance_matrix.qza

qiime emperor biplot \
    --i-biplot diversity_results/RPCA_biplot.qza \
    --m-sample-metadata-file metadata.tsv \
    --m-feature-metadata-file taxonomy.qza \
    --o-visualization diversity_results/RPCA_biplot.qzv \
    --p-number-of-features 5

```
Let's have a look at the visualization [RPCA_biplot.qzv](https://view.qiime2.org/visualization/?src=https://raw.githubusercontent.com/MontagnaLab/InnovativeApproachesForInvertebrateBiodiversityMonitoring/main/outputs/QIIME2_visualizations/RPCA_biplot.qzv)

> [!IMPORTANT]
> RPCA is usually robust enough to be applied to not normalized ASV tables, but in real case scenario you should test it comparing the results of the RPCA performed on normalized and not normalized data. See this [tutorial](https://github.com/biocore/gemelli/blob/master/ipynb/tutorials/RPCA-moving-pictures.ipynb) for more datails.






## **exporting the results**


```bash

# export table (ASVs abundance per sample)
qiime tools export \
  --input-path table.qza \
  --output-path exp
biom convert -i exp/feature-table.biom -o table.tsv --to-tsv

# export taxonomy (ASVs taxonomic assignments)
qiime tools export \
  --input-path taxonomy.qza \
  --output-path exp
mv exp/taxonomy.tsv taxonomy.tsv

```

---
