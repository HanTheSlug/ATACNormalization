ATAC-Seq Normalization Pipeline

![Rough Normalization Schematic](https://github.com/HanTheSlug/ATACNormalization/blob/main/Norm_Schematic.png?raw=true)

Step 1: Normalize for GC content (Still working on)

Step 2: Raw TN5 Counts
Obtain Raw TN5 Counts using genome loader. Format these counts into a CSV, where each row is a sample and each column is a peak. Each entry in the CSV will be a raw TN5 count value

Step 3: TMM Normalization
Run "tmm_normalize_tn5.R" to obtain TMM normalized TN5 Counts. 
The goal of TMM normalization is to normalize for library composition and depth. 

How to run: 
Rscript path_to/tmm_normalize_tn5.R \
        -i path_to/raw_tn5_counts.csv \
        -o output_path/tmm_normalized_tn5_counts.csv

Step 3: PCA
Run "pca_normalization.R" to obtain principle components (PCs) for each sample and principle loadings for each peak. These PCs will be used to represent batch effect and subtracted from the elasticnet model. 
The optimal number of PCs to use will be determined by the BE algorithm (Buja, Andreas, and Nermin Eyuboglu. 1992. “Remarks on Parallel Analysis.” Multivariate Behavioral Research 27 (4): 509–40) or elbow plot method. 

How to run:
Rscript path_to/pca_normalization.R \
        -i path_to/tmm_normalized_tn5_counts.csv \
        -o output_path/output/
        
Step 4: ElasticNet Normalization
After selecting your PCs, run "elasticnet_normalization.py" to obtain your batch effect normalized TN5 counts. The current code takes a long time to run since I test many different alpha and l1 ratio values. Runtime can be reduced by just testing l1_ratio = 0 (L2), 0.5 (Mix), and 1 (L1). Alpha values test the regularization strength, so you can use any values you want depending on how strong you want regularization. 

How to run:
python3 path_to/elasticnet_normalization.py \
        --tn5_counts path_to/tmm_normalized_tn5_counts.csv \
        --pc path_to/output/pc_table.csv \
        --num_pc number-of-pcs \
        --out path_to/elasticnet_normalized_tn5_counts.csv

Optional Steps: Some steps I used to verify my model. 
1. Create a boxplot comparing raw vs normalized TN5 Count values. Averages should be similar after normalization, but original data behaviour and relative magnitudes should be somewhat maintained.
2. Run PCA again on the normalized TN5 counts to ensure that batches cluster and mix together well, but opposing effects (ie treatment vs. non treatment) cluster separately. 
