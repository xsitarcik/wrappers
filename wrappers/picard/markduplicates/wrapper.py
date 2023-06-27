import tempfile

from snakemake.shell import shell

log = snakemake.log_fmt_shell()

mem_mb = snakemake.resources.get("mem_mb", 0)
if not mem_mb:
    mem_mb = snakemake.resources.get("mem_gb", 0) * 1024

java_mem_arg = ""
if mem_mb:
    java_mem_mb = round(0.75 * mem_mb)
    java_mem_arg = "-Xmx{}M".format(java_mem_mb)

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "("
        " picard MarkDuplicates"
        " {java_mem_arg}"
        " I={snakemake.input.bam}"
        " O={snakemake.output.bam}"
        " M={snakemake.output.stat}"
        " TMP_DIR={tmpdir}"
        " VALIDATION_STRINGENCY=SILENT"
        ") {log}"
    )
