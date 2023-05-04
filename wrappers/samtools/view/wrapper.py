from snakemake.shell import shell

shell("samtools view -S -b {snakemake.input.sam} > {snakemake.output.bam}")
