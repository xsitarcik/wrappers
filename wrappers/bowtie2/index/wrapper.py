from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "bowtie2-build "
    " --threads {snakemake.threads} "
    " {snakemake.input.reference}"
    " {snakemake.params.prefix}"
    " {log}"
)
