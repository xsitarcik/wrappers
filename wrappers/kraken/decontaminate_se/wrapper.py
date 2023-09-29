from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("extra", "")

temp_read = snakemake.output.read.removesuffix(".gz")

shell(
    "("
    " extract_kraken_reads.py"
    " -k {snakemake.input.kraken_output} -r {snakemake.input.kraken_report}"
    " -s {snakemake.input.read} -o {temp_read} "
    " --exclude -t {snakemake.params.taxid} {extra} --fastq-output > {snakemake.output.std_out}"
    " &&"
    " pigz {temp_read} -9 -c -p {snakemake.threads} > {snakemake.output.read}"
    " ) {log}"
)
