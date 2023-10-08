import os

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

samtools_extra = snakemake.params.get("samtools_params", "")
ivar_extra = snakemake.params.get("ivar_params", "")

if snakemake.output.consensus.endswith(".fa"):
    out_prefix = os.path.splitext(snakemake.output.consensus)[0]
else:
    raise ValueError("The output must be .fa file")

name_arg = ""
if name := snakemake.params.get("name", ""):
    name_arg = f"-i {name}"

shell(
    "("
    " samtools mpileup {samtools_extra} {snakemake.input.bam}"
    " |"
    " ivar consensus -p {out_prefix} {name_arg} {ivar_extra}"
    ") {log}"
)
