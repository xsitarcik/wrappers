from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")


shell(
    "cutadapt"
    " --output {snakemake.output.r1}"
    " --paired-output {snakemake.output.r2}"
    " --cores {snakemake.threads}"
    " --action {snakemake.params.action}"
    " --overlap {snakemake.params.overlap}"
    " --times {snakemake.params.times}"
    " --error-rate {snakemake.params.error_rate}"
    " --report full"
    " --json={snakemake.output.report}"
    " {extra}"
    " {snakemake.input.r1}"
    " {snakemake.input.r2}"
    " {log}"
)
