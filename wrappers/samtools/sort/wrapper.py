from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "samtools sort"
    " -o {snakemake.output.bam}"
    " -m 2G"
    " --threads {snakemake.threads}"
    " --output-fmt BAM"
    " --reference {snakemake.input.ref}"
    " {snakemake.input.bam}"
    " {log}"
)
