# amplicon_for_terra

Here is a snakemake pipeline that processes raw fastq files and outputs several tables for further processing in R.

I think a path to make this run on Terra is to use the snakemake offical docker container: https://hub.docker.com/r/snakemake/snakemake

and then make a WDL that simply calls snakemake (the submit_jobs.sh file shows how snakemake gets called. https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options has all the command line options).



This workflow currently has no data or anything, but just getting to the point on Terra where snakemake gets called successfully (and gives you an error) would be awesome. Then I can add data to this repo and we can figure out how to use Terra's data storage to do that part.

