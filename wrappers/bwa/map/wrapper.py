import os

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)


if len(snakemake.input.reads) == 1:
    input_reads = snakemake.input.reads[0]
elif len(snakemake.input.reads) == 2:
    input_reads = f"{snakemake.input.reads[0]} {snakemake.input.reads[1]}"
else:
    raise ValueError(f"Expected 1 or 2 reads, got {len(snakemake.input.reads)}")

index = os.path.splitext(snakemake.input.index[0])[0]

shell(
    "(bwa mem "
    " -t {snakemake.threads}"
    " -R \"$(sed 's/\\t/\\\\t/g' < {snakemake.input.read_group})\""
    " {index}"
    " {input_reads}"
    " |"
    " samtools view -bS - > {snakemake.output.bam}) {log}"
)
