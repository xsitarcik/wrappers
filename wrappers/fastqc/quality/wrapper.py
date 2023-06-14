import os
import tempfile

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# fastqc allocates memory per thread, so we need to divide the allocated memory per thread
total_memory = snakemake.resources.get("mem_mb", 0)
thread_memory = int(total_memory / snakemake.threads)
memory_arg = ""
if thread_memory != 0:
    memory_arg = f"--memory {thread_memory}"

PRESUMED_SUFFIX = ".fastq.gz"
if not snakemake.input.read.endswith(PRESUMED_SUFFIX):
    raise ValueError(f"{snakemake.input.read} does not ends with {PRESUMED_SUFFIX}")


with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "fastqc --outdir {tmpdir} --extract --threads {snakemake.threads}" " {memory_arg} {snakemake.input.read} {log}"
    )

    base_name = os.path.basename(snakemake.input.read).replace(PRESUMED_SUFFIX, "")
    html_path = os.path.join(tmpdir, f"{base_name}_fastqc.html")
    zip_path = os.path.join(tmpdir, f"{base_name}_fastqc.zip")

    fastqc_datapath = os.path.join(tmpdir, f"{base_name}_fastqc", "fastqc_data.txt")
    summary_path = os.path.join(tmpdir, f"{base_name}_fastqc", "summary.txt")

    shell("mv {html_path} {snakemake.output.html}")
    shell("mv {zip_path} {snakemake.output.zip}")
    shell("mv {fastqc_datapath} {snakemake.output.qc_data}")
    shell("mv {summary_path} {snakemake.output.summary_txt}")
