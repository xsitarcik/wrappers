from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

unzipped_in_r1 = snakemake.input.r1.removesuffix(".gz")
unzipped_in_r2 = snakemake.input.r2.removesuffix(".gz")

inputs = f"{unzipped_in_r1}\n{unzipped_in_r2}"

shell(
    "("
    " pigz --decompress --force --keep --processes {snakemake.threads} {snakemake.input.r1}"
    " &&"
    " pigz --decompress --force --keep --processes {snakemake.threads} {snakemake.input.r2}"
    " &&"
    " echo {inputs:q} > {snakemake.output.pair_description}"
    " &&"
    " fastuniq -i {snakemake.output.pair_description} -o {snakemake.output.unzipped_out_r1} -p {snakemake.output.unzipped_out_r2}"
    " &&"
    " rm {unzipped_in_r1} && rm {unzipped_in_r2}"
    " &&"
    " pigz {snakemake.output.unzipped_out_r1} -9 -c -p {snakemake.threads} > {snakemake.output.r1}"
    " &&"
    " pigz {snakemake.output.unzipped_out_r2} -9 -c -p {snakemake.threads} > {snakemake.output.r2}"
    " ) {log}"
)
