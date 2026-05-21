## Pairwise Label transfer between queries to establish cross-study concordance

### Source Data: 
mouse "query" datasets: `/space/grp/rschwartz/rschwartz/get_gemma_data.nf/study_names_mouse.txt_author_true_process_samples_true/h5ad/**h5ad`
relabeling files for label harmonization: `/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/relabel_mus_musculus/*_relabel.tsv`
nextflow pipeline for original ref -> query label transfer workflow: `/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline`

### Goal

To establish cross-study concordance between the mouse query datasets, we will perform pairwise label transfer between the query datasets. Certain cell types (e.g. microglia, OPCs) are not identifiable when transferring labels from reference -> query across all query studies. We want to assess transcriptomic similarity between cell types across studies to determine if the "problem" lies with reference or query data. E.g., if queries can be used as references to transfer labels to other queries, then the issue likely lies with the reference data. If performance is poor across all query-query transfers, then the issue likely lies with the query data concordance with established cell type signatures.

### Implementation

look at the way we do label transfer in `/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline`. Specifically, look at workflows/, subworkflows, modules/, main.nf, and bin/. We want to set up the same sort of workflow, but instead of reference -> query, we will do query -> query. We will need to set up a nextflow workflow that takes in the query datasets, performs pairwise label transfer between them, and outputs the results for analysis. we want to enable the same parameters and methods explored in `/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline` (e.g. scVI vs seurat, multiple levels of granularity, variable cutoffs, etc).