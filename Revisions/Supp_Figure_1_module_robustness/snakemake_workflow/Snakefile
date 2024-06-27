
INPUTS = {
        'array' : "../../Snakemake/Pipeline_out/Data/Microarray/data_analysis_ready/eset_batch_training_sample.rds",
        'soma' : "../../Snakemake/Pipeline_out/Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds",
        }

#DATA_TYPES = ["soma", "array", "cbctbnk"]
DATA_TYPES = INPUTS.keys()

BOOT_ITER = range(0,101)

rule all:
    input:
        expand("varpart_objects/{datatype}/{boot}_wgcna_modules.rds", datatype=DATA_TYPES, boot=BOOT_ITER)
        
rule run_varpart:
    input:
        lambda wildcards: INPUTS[wildcards.datatype]
    output:
        "varpart_objects/{datatype}/{boot}_wgcna_modules.rds",
        "varpart_objects/{datatype}/{boot}_scores_sample_level.rds",
        "varpart_objects/{datatype}/{boot}_scores_subject_level.rds"#,
        #"metadata_subsets/{datatype}/{boot}_meta.tsv",
        #"subject_boots/{datatype}/{boot}_subjects_rep.txt"
    singularity:
        "path_to_singularity/monogenic_0.4.sif"
    script:
        "wgcna.R"
