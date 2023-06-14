import os
import tempfile

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)


if len(snakemake.input.reads) == 1:
    input_reads = snakemake.input.reads[0]
elif len(snakemake.input.reads) == 2:
    input_reads = f"{snakemake.input.reads[0]} {snakemake.input.reads[1]}"
else:
    raise ValueError(f"Expected 1 or 2 reads, got {len(snakemake.input.reads)}")

index = os.path.splitext(snakemake.input.index[0])[0]
filter_flags = []
if "exclude_flags" in snakemake.params:
    filter_flags.append(f"--exclude-flags {snakemake.params.exclude_flags}")
if "require_flags" in snakemake.params:
    filter_flags.append(f"--require-flags {snakemake.params.require_flags}")

filter_arg = ""
if filter_flags:
    filters = " ".join(filter_flags)
    filter_arg = f"| samtools view {filters}"

total_memory = snakemake.resources.get("mem_mb", 0)
thread_memory = int(total_memory / snakemake.threads)
memory_arg = ""
if thread_memory != 0:
    memory_arg = f"-m {thread_memory}M"

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "(bwa mem -t {snakemake.threads}"
        " -R \"$(sed 's/\\t/\\\\t/g' < {snakemake.input.read_group})\""
        " {index} {input_reads}"
        " {filter_arg}"
        " |"
        " samtools sort -o {snakemake.output.bam} {memory_arg} -@ {snakemake.threads} -T {tmpdir}"
        " ) {log}"
    )
