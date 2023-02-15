configfile: "config/config.yaml"

rule all:
    # All the output files
    input:
        follicular_bc='{workdir}/sce_follicular_final_bc.rds'.format(workdir=config['workdir']),
        follicular_annotated='{outdir}/sce_follicular_annotated_final_bc.rds'.format(outdir=config['outdir']),
        de_tables_reactome=expand('{outdir}/differential_expression/ReactomePA/timepoint/{{celltype_group}}/{{patient_group}}.rds'.format(outdir=config['outdir']), celltype_group=config['differential_expression_settings']['ReactomePA']['celltype_groups'], patient_group=config['differential_expression_settings']['ReactomePA']['patient_groups']),
        de_res=expand('{outdir}/differential_expression/ReactomePA/malignant_timepoint/b_cells/{{patient_group}}.rds'.format(outdir=config['outdir']), patient_group=config['differential_expression_settings']['ReactomePA']['patient_groups']),
        de_tables_fgsea=expand('{outdir}/differential_expression/fgsea/timepoint/{{celltype_group}}/{{patient_group}}.rds'.format(outdir=config['outdir']), celltype_group=config['differential_expression_settings']['fgsea']['celltype_groups'], patient_group=config['differential_expression_settings']['fgsea']['patient_groups']),
        
########### add new input files ##############

# Create SingleCellExperiment
rule create_follicular_sce:
    input:
        filtered_matrices=[config['follicular_data'][x]['filtered_matrix_path'] for x in config['follicular_data']]
    output:
        '{workdir}/sce_follicular.rds'.format(workdir=config['workdir'])
    params:
        name='create-follicular-sce',
        sample_names=[config['follicular_data'][x]['dataset'] for x in config['follicular_data']],
        timepoints=[config['follicular_data'][x]['timepoint'] for x in config['follicular_data']],
        progression_statuses=[config['follicular_data'][x]['progression_status'] for x in config['follicular_data']],
        patient_progression_statuses=[config['follicular_data'][x]['patient_progression_status'] for x in config['follicular_data']],
        patients=[config['follicular_data'][x]['patient'] for x in config['follicular_data']],
        batches=[config['follicular_data'][x]['batch'] for x in config['follicular_data']]
    log:
        '{logdir}/logs/create_follicular_sce.log'.format(logdir=config['logdir'])
    benchmark:
        '{logdir}/benchmarks/create_follicular_sce.txt'.format(logdir=config['logdir'])
    shell:
        'Rscript R/create_sce.R '
        '--sample_names {params.sample_names} '
        '--filtered_matrices {input.filtered_matrices} '
        '--timepoints {params.timepoints} '
        '--progression {params.progression_statuses} '
        '--patient_progression {params.patient_progression_statuses} '
        '--patients {params.patients} '
        '--batches {params.batches} '
        '--outfname {output} '
        '>& {log}'


# Filter and normalize SingleCellExperiment
rule preprocess_follicular_sce:
    input:
        follicular_raw='{workdir}/sce_follicular.rds'.format(workdir=config['workdir']),
    output:
        follicular_normalized='{workdir}/sce_follicular_normalized.rds'.format(workdir=config['workdir']),
    params:
        name='preprocess-follicular-sce',
        mito_thres=config['filter_settings']['mito_thres'],
        ribo_thres=config['filter_settings']['ribo_thres'],
    log:
        '{logdir}/logs/preprocess_follicular_sce.log'.format(logdir=config['logdir'])
    benchmark:
        '{logdir}/benchmarks/preprocess_follicular_sce.txt'.format(logdir=config['logdir'])
    shell:
        'Rscript R/preprocess_sce.R '
        '--sce {input.follicular_raw} '
        '--mito_thres {params.mito_thres} '
        '--ribo_thres {params.ribo_thres} '
        '--outfname {output} '
        '>& {log}'

# Cell cycle assignments
rule cyclone_follicular_sce:
    input:
        follicular_normalized='{workdir}/sce_follicular_normalized.rds'.format(workdir=config['workdir']),
    output:
        cc_result='{outdir}/cyclone/cc_result.rds'.format(outdir=config['outdir']),
    params:
        name='cyclone-follicular-sce',
    log:
        '{logdir}/logs/cyclone_follicular_sce.log'.format(logdir=config['logdir'])
    benchmark:
        '{logdir}/benchmarks/cyclone_follicular_sce.txt'.format(logdir=config['logdir'])
    threads: 15
    shell:
        'Rscript R/cyclone_sce.R '
        '--sce {input.follicular_normalized} '
        '--ncpus {threads} '
        '--outfname {output} '
        '>& {log}'

# batch correction
rule batch_correct:
    input:
        follicular_normalized='{outdir}/sce_follicular_normalized.rds'.format(outdir=config['outdir']),
    output:
        follicular_bc='{workdir}/sce_follicular_final_bc.rds'.format(workdir=config['workdir']),
    params:
        name='batch-correct',
        batch_col=config['batch_correction_settings']['batch_variable'],
        method=config['batch_correction_settings']['method'],
        conda_env='r-tensorflow',
    log:
        '{logdir}/logs/batch_correct.log'.format(logdir=config['logdir'])
    benchmark:
        '{logdir}/benchmarks/batch_correct.txt'.format(logdir=config['logdir'])
    shell:
        'Rscript R/batch_correct.R '
        '--sce {input.follicular_annotated} '
        '--conda_env {params.conda_env} '
        '--method {params.method} '
        '--batch_column {params.batch_col} '
        '--outfname {output.follicular_bc} '
        '>& {log}'

# Annotate follicular SCE and CDS with assignments
rule annotate_follicular_final:
    input:
        follicular_bc='{workdir}/sce_follicular_final_bc.rds'.format(workdir=config['workdir']),
        cc_result='{outdir}/cyclone/cc_result.rds'.format(outdir=config['outdir']),
    output:
        follicular_annotated='{outdir}/sce_follicular_annotated_final_bc.rds'.format(outdir=config['outdir']),
    params:
        name='annotate-follicular-final',
    log:
        '{logdir}/logs/annotate_follicular_final.log'.format(logdir=config['logdir'])
    benchmark:
        '{logdir}/benchmarks/annotate_follicular_final.txt'.format(logdir=config['logdir'])
    shell:
        'Rscript R/annotate_follicular_final.R '
        '--sce {input.follicular_normalized} '
        '--cyclone {input.cc_result} '
        '--outfname {output.follicular_annotated} '
        '>& {log}'

# DE between FL and DLBCL
rule differential_expression_malignant_b_by_timepoint:
    input:
        follicular_bc='{workdir}/sce_follicular_annotated_final_bc.rds'.format(workdir=config['workdir']),
    output:
        de_tables='{outdir}/differential_expression/ReactomePA/timepoint/DE_FL_DLBCL.rds'.format(outdir=config['outdir']),
    params:
        de_method=config['differential_expression_settings']['ReactomePA']['gene_method'],
        num_genes=config['differential_expression_settings']['ReactomePA']['num_genes'],
        min_gene_counts=config['differential_expression_settings']['ReactomePA']['min_gene_counts'],
        name='differential-expression-reactome-timepoint',
    log:
        '{logdir}/logs/differential_expression_timepoint/ReactomePA/DE_FL_DLBCL.log'.format(logdir=config['logdir'])
    benchmark:
        '{logdir}/benchmarks/differential_expression_timepoint/ReactomePA/DE_FL_DLBCL.txt'.format(logdir=config['logdir'])
    shell:
        'Rscript R/de_timepoint.R '
        '--sce {input.follicular_bc} '
        '--method_gene {params.de_method} '
        '--method_pathway ReactomePA '
        '--ngene {params.num_genes} '
        '--outfname {output.de_tables} '
        '>& {log}'

# DE between FL and DLBCL within patient
rule differential_expression_malignant_b_by_timepoint_by_patient:
    input:
        follicular_bc='{outdir}/sce_annotated_subset/sce_follicular_annotated_final_bc.rds'.format(outdir=config['outdir']),
    output:
        de_res='{outdir}/differential_expression/ReactomePA/malignant_b_timepoint/{{patient_group}}.rds'.format(outdir=config['outdir']),
    params:
        min_gene_counts=config['differential_expression_settings']['ReactomePA']['min_gene_counts'],
        patient_group=lambda wildcards: config['differential_expression_settings']['ReactomePA']['patient_groups'][wildcards.patient_group],
        num_genes=config['differential_expression_settings']['ReactomePA']['num_genes'],
        name='differential-expression-reactome-bcell-malignant-timepoint-{patient_group}',
        bcell_labels=config['b_cell_labels'],
    log:
        '{logdir}/logs/differential_expression_bcell_malignant_timepoint/{{patient_group}}.log'.format(logdir=config['logdir'])
    benchmark:
        '{logdir}/benchmarks/differential_expression_bcell_malignant_timepoint/{{patient_group}}.txt'.format(logdir=config['logdir'])
    shell:
        'Rscript R/de_b_cells_malignant_timepoint.R '
        '--sce {input.bcells_annotated} '
        '--min_gene_counts {params.min_gene_counts} '
        '--patients {params.patient_group} '
        '--bcell_labels {params.bcell_labels} '
        '--ngene {params.num_genes} '
        '--outfname {output.de_res} '
        '>& {log}'





