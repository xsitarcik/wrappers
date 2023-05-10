from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)


shell(
    " picard CreateSequenceDictionary"
    " REFERENCE={snakemake.input.reference}"
    " OUTPUT={snakemake.output.seq_dict}"
    " {log}"
)
