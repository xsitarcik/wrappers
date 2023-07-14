import os

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


if len(snakemake.input.reads) == 1:
    input_arg = f"-U {snakemake.input.reads[0]}"
else:
    input_arg = f"-1 {snakemake.input.reads[0]} -2 {snakemake.input.reads[1]}"

index = os.path.splitext(snakemake.input.index[0])[0]
extra = snakemake.params.get("extra", "")

shell(
    "(bowtie2 "
    " -x {index}"
    " {input_arg}"
    " --threads {snakemake.threads}"
    " {extra}"
    " |"
    " samtools view -bS -"
    " > {snakemake.output.bam}"
    ") {log}"
)
