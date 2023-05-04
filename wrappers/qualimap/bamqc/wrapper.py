import os

from snakemake.shell import shell

# use the same way to unset DISPLAY as in official snakemake wrappers
if os.environ.get("DISPLAY"):
    del os.environ["DISPLAY"]

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", [])
extra_str = " ".join(extra)

resources_str = ""
if "mem_mb" in snakemake.resources.keys():
    mem = snakemake.resources["mem_mb"]
    resources_str += f"--java-mem-size={mem}M"

shell(
    "qualimap bamqc"
    " -bam {snakemake.input.bam}"
    " -outdir {snakemake.output.report_dir}"
    " -nt {snakemake.threads}"
    " {extra_str}"
    " {resources_str}"
    " {log}"
)
