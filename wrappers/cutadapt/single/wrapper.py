from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("extra", "")

adapter_cmd = ""
if snakemake.params.get("anywhere_adapter", ""):
    adapter_cmd += f" --anywhere {snakemake.params.anywhere_adapter}"
if snakemake.params.get("front_adapter", ""):
    adapter_cmd += f" --front {snakemake.params.front_adapter}"
if snakemake.params.get("regular_adapter", ""):
    adapter_cmd += f" --adapter {snakemake.params.regular_adapter}"

overlap = snakemake.params.get("overlap", 3)
error_rate = snakemake.params.get("error_rate", 0.1)
times = snakemake.params.get("times", 1)
action = snakemake.params.get("action", "trim")

cut_cmd = ""
if cut_value := snakemake.params.get("head_cut", ""):
    cut_cmd += f" --cut {cut_value}"
if cut_value := snakemake.params.get("tail_cut", ""):
    cut_cmd += f" --cut -{cut_value}"

shell(
    "cutadapt"
    " {adapter_cmd}"
    " --output {snakemake.output.read}"
    " --cores {snakemake.threads}"
    " --action {action}"
    " --overlap {overlap}"
    " --times {times}"
    " --error-rate {error_rate}"
    " {cut_cmd}"
    " --minimum-length {snakemake.params.minimum_length}"
    " --quality-cutoff {snakemake.params.quality_cutoff}"
    " --report full"
    " --json={snakemake.output.report}"
    " {extra}"
    " {snakemake.input.read}"
    " {log}"
)
