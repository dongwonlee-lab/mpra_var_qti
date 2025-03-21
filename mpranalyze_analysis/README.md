Here, we have Jupyter Notebooks for our MPRA analysis pipeline using MPRAnalyze.

Our pipeline consists of six steps, which are described briefly below. 
All inputs and outputs are also available for download here: 
[dataset for the Jupyter Notebooks](https://drive.google.com/file/d/1GTQ2Wd9bDJdpi5ddtUgF6-HPIQc_YG85/view?usp=sharing)

* [1_MPRAnalyze_Prepare_Annotation_Data.ipynb](1_MPRAnalyze_Prepare_Annotation_Data.ipynb): 
generate annotation tables that will be used as inputs for MPRAnalyze
    - output:
        - data/dna_annot.txt
        - data/rna_annot_pool.txt
        - data/merged_dna_annot_allelic.txt
        - data/merged_rna_annot_pool_allelic.txt


* [2_MPRAnalyze_Quantification_Analysis.ipynb](2_MPRAnalyze_Quantification_Analysis.ipynb):
quantify expression levels for all elements using MPRAnalyze
    - input:
        - data/mpra_qtigc_pgl4v1_pool1.txt
        - data/mpra_qtigc_pgl4v2_pool1.txt
        - data/mpra_qtigc_pgl4v1_pool2.txt
        - data/mpra_qtigc_pgl4v2_pool2.txt
        - data/dna_annot.txt
        - data/rna_annot_pool.txt
    - output:
        - results/pool1_mpranalyze_results.txt
        - results/pool2_mpranalyze_results.txt

* [3a_MPRAnalyze_Allelic_Cmp_Prepare_Count_Data.ipynb](3a_MPRAnalyze_Allelic_Cmp_Prepare_Count_Data.ipynb):
reformat the count data for allelic comparison analysis
    - input:
        - data/mpra_qtigc_pgl4v1_pool1.txt
        - data/mpra_qtigc_pgl4v2_pool1.txt
        - data/mpra_qtigc_pgl4v1_pool2.txt
        - data/mpra_qtigc_pgl4v2_pool2.txt
    - output:
        - data/pool1_dna_counts_allelic.txt
        - data/pool1_rna_counts_allelic.txt
        - data/pool2_dna_counts_allelic.txt
        - data/pool2_rna_counts_allelic.txt
        - data/merged_rna_annot_pool_allelic_p2.txt

* [3b_MPRAnalyze_Allelic_Comparison.r](3b_MPRAnalyze_Allelic_Comparison.r):
R script for allelic comparison analysis using the files generated in the previous step.  
This is the script for parallel running. It takes 6 arguments, and the usage is 
`Rscript 3b_MPRAnalyze_Allelic_Comparison.r <dna_annot> <rna_annot> <dna_cnt> <rna_cnt>
<pool_number> <job_number>`  Each job evaluates 10 variants.
    - input:
        - data/merged_dna_annot_allelic.txt
        - data/merged_rna_annot_pool_allelic.txt
        - data/merged_rna_annot_pool_allelic_p2.txt
        - data/pool1_dna_counts_allelic.txt
        - data/pool1_rna_counts_allelic.txt
        - data/pool2_dna_counts_allelic.txt
        - data/pool2_rna_counts_allelic.txt
    - output:
        - results/pool[1-2]/MPRAObject_allelic_[start_var_num]-[end_var_num].txt

* [4_MPRAnalyze_Merge_Allelic_Comparison_Results.ipynb](4_MPRAnalyze_Merge_Allelic_Comparison_Results.ipynb):
merge outputs of allelic comparison analysis from the previous step
    - input:
        - results/pool[1-2]/MPRAObject_allelic_[start_var_num]-[end_var_num].txt
    - output:
        - results/pool1_mpranalyze_allelic_results_merged.txt
        - results/pool2_mpranalyze_allelic_results_merged.txt

* [5_Calculate_CPM.ipynb](5_Calculate_CPM.ipynb):
calculate CPM (count per million reads) value for each element
    - input
        - data/pool1_dna_counts_allelic.txt
        - data/pool2_dna_counts_allelic.txt
    - output
        - results/pool1_dna_counts_cpm.txt
        - results/pool2_dna_counts_cpm.txt

* [6_MPRA_Enhancer_Allelic_Significance.ipynb](6_MPRA_Enhancer_Allelic_Significance.ipynb):
identify differential allelic enhancer variants 
    - input
        - results/pool1_mpranalyze_results.txt
        - results/pool1_dna_counts_cpm.txt
        - results/pool1_mpranalyze_allelic_results_merged.txt
        - results/pool2_mpranalyze_results.txt
        - results/pool2_dna_counts_cpm.txt
        - results/pool3_mpranalyze_allelic_results_merged.txt
    - output
        - pool1_final_result.txt
