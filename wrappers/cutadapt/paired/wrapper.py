from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")


shell(
    "cutadapt"
    " --output {snakemake.output.r1}"
    " --paired-output {snakemake.output.r2}"
    " --cores {snakemake.threads}"
    " --report full"
    " --json={snakemake.output.report}"
    " {extra}"
    " {snakemake.input}"
    " {log}"
)
