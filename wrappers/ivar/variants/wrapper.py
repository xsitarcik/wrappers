import os

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

samtools_extra = snakemake.params.get("samtools_params", "")
ivar_extra = snakemake.params.get("ivar_params", "")

if snakemake.output.tsv.endswith(".tsv"):
    out_prefix = os.path.splitext(snakemake.output.tsv)[0]
else:
    raise ValueError("The output must be .tsv file")

shell(
    "("
    " samtools mpileup {samtools_extra} {snakemake.input.bam}"
    " |"
    " ivar variants -p {out_prefix} -r {snakemake.input.ref} {ivar_extra}"
    ") {log}"
)
