

sequencing_path_2022_07_08 = "/n/groups/kwon/data1/sequencing_run_archive_DO_NOT_EDIT/2022_07_08_MiSeq_V4/Demultiplexed/split_sample_ids_2022_07_08/"

all_samples = {
    "2022_07_08": [s for s in glob_wildcards(sequencing_path_2022_07_08 +"{sample}.fastq").sample if ((os.path.getsize(sequencing_path_2022_07_08 +"{}.fastq".format(s)) > 50000))]
}
print(all_samples)

rule target:
    input:
        # In a first run of this meta-wrapper, comment out all other inputs and only keep this one.
        # Looking at the resulting plot, adjust the `truncLen` in rule `dada2_filter_trim_pe` and then
        # rerun with all inputs uncommented.
        "results/dada2/taxonomy_silva_species.RDS",
        "results/dada2/taxonomy_decipher.RDS"

rule dada2_filter_trim_se:
    input:
        # Single-end files without primers sequences
        fwd="/n/groups/kwon/data1/sequencing_run_archive_DO_NOT_EDIT/{run}_MiSeq_V4/Demultiplexed/split_sample_ids_{run}/{sample}.fastq"
    output:
        filt="filtered-se/{run}-{sample}.fastq.gz",
        stats="reports/dada2/filter-trim-se/{run}-{sample}.tsv"
    # Even though this is an R wrapper, use named arguments in Python syntax
    # here, to specify extra parameters. Python booleans (`arg1=True`, `arg2=False`)
    # and lists (`list_arg=[]`) are automatically converted to R.
    # For a named list as an extra named argument, use a python dict
    #   (`named_list={name1=arg1}`).
    params:
        maxEE=1,
        truncLen=230,
        trimLeft=10
    log:
        "logs/dada2/filter-trim-se/{run}-{sample}.log"
    threads: 1 # set desired number of threads here
    group:
        "filter_trim_se"
    conda:
        "envs/dada2.yaml"
    wrapper:
        "0.79.0/bio/dada2/filter-trim"


def get_samples(wildcards):
    return ["filtered-se/{}-{}.fastq.gz".format(wildcards.run, s) for s in all_samples[wildcards.run]]

rule dada2_learn_errors:
    input:
    # Quality filtered and trimmed forward FASTQ files (potentially compressed)
        get_samples
    output:
        err="results/dada2/model_{run}.RDS",# save the error model
        plot="reports/dada2/errors_{run}.png",# plot observed and estimated rates
    params:
        randomize=True
    log:
        "logs/dada2/learn-errors/learn-errors_{run}.log"
    threads: 10 # set desired number of threads here
    conda:
        "envs/dada2.yaml"
    wrapper:
        "0.79.0/bio/dada2/learn-errors"

rule dada2_dereplicate_fastq:
    input:
    # Quality filtered FASTQ file
        "filtered-se/{run}-{sample}.fastq.gz"
    output:
    # Dereplicated sequences stored as `derep-class` object in a RDS file
        "uniques/{run}-{sample}.RDS"
    log:
        "logs/dada2/dereplicate-fastq/{run}-{sample}.log"
    conda:
        "envs/dada2.yaml"
    group:
        "dereplicate_fastq"
    wrapper:
        "0.79.0/bio/dada2/dereplicate-fastq"


rule dada2_sample_inference:
    input:
    # Dereplicated (aka unique) sequences of the sample
        derep = "uniques/{run}-{sample}.RDS",
        err = "results/dada2/model_{run}.RDS" # Error model
    output:
        "denoised/{run}-{sample}.RDS" # Inferred sample composition
    log:
        "logs/dada2/sample-inference/{run}-{sample}.log"
    threads: 1 # set desired number of threads here
    conda:
        "envs/dada2.yaml"
    group:
        "sample_inference"
    wrapper:
        "0.79.0/bio/dada2/sample-inference"

rule dada2_make_table_se:
    input:
    # Inferred composition
        ["denoised/{}-{}.RDS".format(r, s) for r in all_samples.keys() for s in all_samples[r]]
    output:
        "results/dada2/seqTab-se.RDS"
    params:
        names=[s for r in all_samples.keys() for s in all_samples[r]] # Sample names instead of paths
    log:
        "logs/dada2/make-table/make-table-se.log"
    threads: 10 # set desired number of threads here
    conda:
        "envs/dada2.yaml"
    wrapper:
        "0.79.0/bio/dada2/make-table"

rule dada2_remove_chimeras:
    input:
        "results/dada2/seqTab-se.RDS" # Sequence table
    output:
        "results/dada2/seqTab.nochimeras.RDS" # Chimera-free sequence table
    log:
        "logs/dada2/remove-chimeras/remove-chimeras.log"
    threads: 10 # set desired number of threads here
    conda:
        "envs/dada2.yaml"
    wrapper:
        "0.79.0/bio/dada2/remove-chimeras"

rule dada2_collapse_nomismatch:
    input:
        "results/dada2/seqTab.nochimeras.RDS" # Chimera-free sequence table
    output:
        "results/dada2/seqTab.collapsed.RDS"
    log:
        "logs/dada2/collapse-nomismatch/collapse-nomismatch.log"
    threads: 10 # set desired number of threads here
    conda:
        "envs/dada2.yaml"
    wrapper:
        "0.79.0/bio/dada2/collapse-nomismatch"

rule taxonomy_decipher:
    input:
        seqtab = "results/dada2/seqTab.collapsed.RDS",
        decipher_db = "/n/groups/kwon/joseph/SILVA_SSU_r138_2019.RData"
    output:
        taxonomy_decipher = "results/dada2/taxonomy_decipher.RDS"
    conda:
        "envs/dada2.yaml"
    threads: 10
    script:
        "scripts/assign_taxonomy_decipher.R"

rule taxonomy_silva_species:
    input:
        seqtab = "results/dada2/seqTab.collapsed.RDS",
        silva_species_db = "/n/groups/kwon/joseph/silva_species_assignment_v138.1.fa.gz"
    output:
        taxonomy_silva_species = "results/dada2/taxonomy_silva_species.RDS"
    conda:
        "envs/dada2.yaml"
    threads: 1
    script:
        "scripts/assign_taxonomy_silva_species.R"