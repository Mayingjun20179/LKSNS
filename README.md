# LKSNS
This is the data and code for "Prediction of Long Non-coding RNA-protein Interaction through Kernel Soft Neighborhood Similarity"
Dataset description:
five  cross  datasets (Refer to Zhang wen)  
   extracted_interaction.txt        The known lncRNA-protein interaction matrix  
   protein_ctd,mat                  Proteins' feature dataset    
   extracted_lncRNA_expression.txt  The lncRNA expression profile dataset
   extracted_lncRNA_sequence_CT.txt This is the sequence composition dataset of lncRNA           

case study dataset
  YCdata.TQ_LJ_stasebate  test dataset 1，the lncRNA-protein interaction extracted from stasebate corresponding to the benchmark dataset 
  YCdata.TQ_LJ_inter3     test dataset 2，the lncRNA-protein interaction extracted from NPInter v3.0 corresponding to the benchmark dataset
  YCdata.ZL_GY_inter2     The mapping between transcriptional ID and gene ID of lncRNA ，extracted from NONCODE v5.0
  YCdata.name             The names of lncRNA and proteins used in the Benchmark datasets extracted from NPInter v2.0 and 
  YCdata.TR_LJM           The lncRNA-protein interaction matrix in the benchmark dataset   

The code description:

Compare_different_similarity_models_five_cross.m  Performance comparison with different similarity computation models Based on five-fold crossover experiment and label propagation algorithm
compare_different_similarity_drow_figure1.m       The drawing procedure for the above program (pare_different_similarity_models_five_cross.m) results, the drawing procedure of figure 1
Comparison_new_methods_integrated_dataset.m       This program is about comparison with some relatively new methods on integrated dataset
kernel_neighbors_parameters_comparison.m          This program is about comparing the results of different model parameters, different kernel functions and different neighbor Numbers
compare_kernel_figure_4.m                         This program is based on the above program (kernel_neighbors_parameters_comparison.m ) results to draw a comparison image of different kernel functions
Neighbor_number_selection_5.m                     This program is based on the above program (kernel_neighbors_parameters_comparison.m ) results to draw a comparison image of different number of neighbors
AUPR_AUC_on_different_datasets_figure6.m          This program is based on the above program (kernel_neighbors_parameters_comparison.m ) results to draw a comparison image of different datasets
case_study.m                                      This program is a case study
