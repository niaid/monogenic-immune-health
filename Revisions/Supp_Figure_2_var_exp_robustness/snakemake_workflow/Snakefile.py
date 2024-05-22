
INPUTS = {
        'cbctbnk' : "../../Snakemake/Pipeline_out/Data/TBNK/data_analysis_ready/tbnk_training_sample_level_eset.rds",
        'array' : "../../Snakemake/Pipeline_out/Data/Microarray/data_analysis_ready/eset_batch_training_sample.rds",
        'array_mod' : "../../Snakemake/Pipeline_out/Data/Microarray/analysis_output/WGCNA/array_sample_scores.rds",
        'soma' : "../../Snakemake/Pipeline_out/Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds",
        'soma_mod' : "../../Snakemake/Pipeline_out/Data/Somalogic/analysis_output/wgcna_results/scores_sample_level.rds"
        }

#DATA_TYPES = ["soma", "array", "cbctbnk"]
DATA_TYPES = INPUTS.keys()

BOOT_ITER = range(0,101)

rule all:
    input:
        expand("varpart_objects/{datatype}/{boot}_varpart.rds", datatype=DATA_TYPES, boot=BOOT_ITER)
        
rule run_varpart:
    input:
        lambda wildcards: INPUTS[wildcards.datatype]
    output:
        "varpart_objects/{datatype}/{boot}_varpart.rds",
        "metadata_subsets/{datatype}/{boot}_meta.tsv",
        "subject_boots/{datatype}/{boot}_subjects_rep.txt",
        "subject_boots/{datatype}/{boot}_subjects_norep.txt"
    singularity:
        "path_to_singularity/monogenic_0.4.sif"
    script:
        "run_varpart_snakemake.R"
