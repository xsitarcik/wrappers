from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


# Number of input/output compression threads to use in addition to main thread
# from samtools index documentation: http://www.htslib.org/doc/samtools-index.html
threads = snakemake.threads - 1

shell("samtools index -@ {threads} {snakemake.input.bam} {snakemake.output.bai} {log}")
