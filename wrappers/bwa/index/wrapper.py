from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell("bwa index " " -p {snakemake.params.prefix}" " -a {snakemake.params.approach}" " {snakemake.input[0]}" " {log}")
