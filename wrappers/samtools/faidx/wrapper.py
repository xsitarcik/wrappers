from snakemake.shell import shell

shell("samtools faidx {snakemake.input.reference}")
