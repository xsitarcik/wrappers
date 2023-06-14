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
exclude_flag_param = snakemake.params.get("exclude_flag", "")
if exclude_flag_param:
    filter_flags.append(f"--exclude-flags {exclude_flag_param}")
require_flag_param = snakemake.params.get("require_flag", "")
if require_flag_param:
    filter_flags.append(f"--require-flags {require_flag_param}")

filter_arg = ""
if filter_flags:
    filters = " ".join(filter_flags)
    filter_arg = f"| sramtools view {filters}"

print(filter_arg)
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
