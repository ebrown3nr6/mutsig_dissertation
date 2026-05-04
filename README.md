# mutsig_dissertation
Repository of files and code for my MBiol 2026 Dissertation

All Bash commands are written in dissertation methods.

Please see below for input/output folders and files

To find code on Jupyter: /home/jovyan/shared-team/2025-masters-project/people/eleanor/CLEAN_CODE/for_github

Base folder: /home/jovyan/shared-team/2025-masters-project/people/eleanor/original_and_processed_files

    1. prepare_files.ipynb - Quality control, select singletons, remove recombination for 6-class, trinucleotide context, codon usage, checks for strand symmetry
    
        inputs: base_folder - VCF files
                            - GFF files
                            - FASTA files
                            
        outputs: base_folder/processed_c - species_strand_symmetry.csv
                                         - mutation_summary_final_c.csv
                                         - trinucleotide_context.csv

    2. tvts.ipynb - Transition versus tranversion rates; breakdown within each class
        
        inputs: base_folder/processed_c - mutation_summary_final_c.csv
        
        outputs: base_folder/processed_c/charts/tstv - tstv_composition.png
                                                     - tstv_pairwise_stats.csv
                                                     - tstv_per_class.csv
                                                     - tstv_ratio.csv
                                                     - tstv_ratio.png
    
    3. dNdS.ipynb - dN/dS rates
        
        inputs: ~/shared-team/2025-masters-project/people/eleanor/non_synon/snpEff/data/summaries - all files
                base_folder - FNA files
                            - GFF files
                
        outputs: base_folder/processed_c/charts/dnds - dnds_stats.csv
                                                     - dnds.png
                                                     - snpeff_summary_c.csv
    
    4. genic.ipynb - Proportion of genic mutations vs proportion of genic genome content
        
        inputs: base_folder - FNA files
                            - GFF files
                base_folder/processed_c - mutation_summary_final_c.csv
           
        outputs: base_folder/processed_c/charts/genic_intergenic - genic_chisq_results.csv
                                                                 - genic_descriptive_stats.csv
                                                                 - genic_difference_plot.png
                                                                 - genic_intergenic_content.png
                                                                 - genic_intergenic_proportions.png
    
    5. codon_pos.ipynb - codon position for each mutational class in each species and overall, looking at C>T rate specifically for each species and position 
        
        inputs: base_folder/processed_c - mutation_summary_final_c.csv
        
        outputs: base_folder/processed_c/charts/codon_position - codon_gof_per_species.csv
                                                               - plot1a_grouped_normalised_all_species.png
                                                               - plot1b_grouped_normalised_all_species.png
                                                               - plot2_proportional_all_species.png
                                                               
    
    6. clustering.ipynb - hierarchical and UMAP clustering
        
        inputs: base_folder/processed_c - mutation_summary_final_c.csv
        
        outputs: base_folder/processed_c/charts/clustering - cluster_assignments_6class.csv
                                                           - hierarchical_clustering_6class.csv
                                                           - pca_coordinates.csv
                                                           - pca_loadings.png
                                                           - pca_scatter.png
    
    

    7. Investigations - Investigations into SNP distribution, gap in P. aeruginosa, Finding Ts/Tv from Ruis et al., 2023
        
    General investigations:
        
        inputs: base_folder/processed_c - filtered VCF files
        
        outputs: base_folder/processed_c/charts/investigations/outlier_plots - [species]_snps_per_sample.png
                                                                             - snp_density_genome.png

    Finding Ruis Tv/Ts:
        
        inputs: base_folder/processed_c/charts/investigations - 41467_2023_42916_MOESM5_ESM.xlsx
        
        outputs: prints Ts/Tvs (does not save)
        

    8. NB_stats_model.R - Negative binomial model for 6-class *Run before mut_signatures.ipynb*
        
        inputs: base_folder/processed_c - mutation_summary_final_c.csv
        
        outputs: base_folder/processed_c/charts/mut_signatures - model1_coefficients.csv
                                                               - model1_predicted_rates.csv
                                                               - model1_species_constrasts.csv

    
    
    9. mut_signatures.ipynb - 6 Class and trinucleotide signatures *stats from R*
    
        inputs: base_folder/processed_c - mutation_summary_final_c.csv
                                        - trinucleotide_context.csv
                                        
        outputs: base_folder/processed_c/charts/mut_signatures - trinucleotide_context.png
                                                               - trinucleotide_enrichment_combined.csv
                                                               - trinuc_enrichment_96.png
                                                               - trinuc_cross_species_consistency.csv
                                                               - nb_deviation_by_class.png
                                                               - mutational_signatures.png
