# SpinalCordStratification

Main Code and Order

1) TE Quantification on ASU's AGAVE HPC using SQuIRE. EM normalized gene expression from GSE153960. Clinical data by request from NYGC. 
2) SQuIRE_PostProcessing_HiSeq.R and SQuIRE_PostProcessing_NovaSeq.R
3) SQuIRE_Meta.R
4) ClusterPrep_4Covar_Glia_RIN_Site_Tissue_HiSeq.R ; ClusterPrep_4Covar_Glia_RIN_Site_Tissue_NovaSeq.R ; ClusterPrep_4Covar_HiSeq_withControls.R ; ClusterPrep_4Covar_NovaSeq_withControls.R
5) RankEstimation_4Covar.R
6) SAKE.R
7) PostSAKE_SpinalCord.R
8) TopFeatures_vectorized_MoR.R
9A) GSEA 
9B) EnrichrPrep.R
9C) Enrichr_HGEA_pval_heatmap.R
10) CortexSpinal_Concordance_Clean.R
11) ConcordanceMeta_TISSUE_and_PLATFORM_LEVEL.R (files obtained from previous step)
12) Bootstrap_Concordance_pval.R ; Bootstrap_Concordance_pval_HiSeq.R ; Bootstrap_Concordance_pval_NovaSeq.R
13) SpinalCordStratification_survival_clinicalparameters.R
14) SpinalCord_IIDSurvival_TissueLevel.R (Cox Regression)
15) SpinalCordStratification_UnivariateAnalysis_Final.R
16) SpinalCord_SubtypeBiomarkers_GroupedViolin_CNSLevelExpression.R
17) PLSDA_SpinalCord.R
18) Python scikit learn scripts for classification

Supplemental Scripts:

-JH_MuSiC_CellDecon.R

-TruncatedStathmin2_SubtypePlot.R

-UsefulFunctions_v3.R

-VisualizeHighConcordantPatients.R

-Plot_ALSOX_markergenes_RPKM.R

-SubtypeLabelPiechart.R

-MutationFrequency.R

R and Python scripts used in the ALS Spinal Cord Stratification analysis for: Eshima J, Pennington TR, Choudrey R, Garcia JM, Fricks J, Smith BS. Elevated expression of B4GALT6, GABRA1, GAD2, GLRA3, HTR2A, PCSK1, and SLC17A6 are postmortem markers for the ALS-Ox subtype.
