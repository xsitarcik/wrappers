from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")

adapter_cmd = ""
if snakemake.params.get("anywhere_adapter", ""):
    adapter_cmd += f" --anywhere {snakemake.params.anywhere_adapter} -B {snakemake.params.anywhere_adapter}"
if snakemake.params.get("front_adapter", ""):
    adapter_cmd += f" --front {snakemake.params.front_adapter} -G {snakemake.params.front_adapter}"
if snakemake.params.get("regular_adapter", ""):
    adapter_cmd += f" --adapter {snakemake.params.regular_adapter} -A {snakemake.params.regular_adapter}"

overlap = snakemake.params.get("overlap", 3)
error_rate = snakemake.params.get("error_rate", 0.1)
times = snakemake.params.get("times", 1)
action = snakemake.params.get("action", "trim")

cut_cmd = ""
if cut_value := snakemake.params.get("head_cut", ""):
    cut_cmd += f" --cut {cut_value}"
if cut_value := snakemake.params.get("tail_cut", ""):
    cut_cmd += f" --cut -{cut_value}"

if r2_quality_cutoff := snakemake.params.get("r2_quality_cutoff", ""):
    r2_quality_cutoff_cmd = f" -Q {r2_quality_cutoff}"
else:
    r2_quality_cutoff_cmd = ""

shell(
    "cutadapt"
    " {adapter_cmd}"
    " --output {snakemake.output.r1}"
    " --paired-output {snakemake.output.r2}"
    " --cores {snakemake.threads}"
    " --action {action}"
    " --overlap {overlap}"
    " --times {times}"
    " --error-rate {error_rate}"
    " --minimum-length {snakemake.params.minimum_length}"
    " --quality-cutoff {snakemake.params.quality_cutoff}"
    " {r2_quality_cutoff_cmd}"
    " --report full"
    " --json={snakemake.output.report}"
    " {extra}"
    " {snakemake.input.r1}"
    " {snakemake.input.r2}"
    " {log}"
)
