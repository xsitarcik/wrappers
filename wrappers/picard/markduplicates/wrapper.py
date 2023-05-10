import tempfile

from snakemake.shell import shell

log = snakemake.log_fmt_shell()

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "("
        " picard MarkDuplicates"
        " I={snakemake.input.bam}"
        " O={snakemake.output.bam}"
        " M={snakemake.output.stat}"
        " TMP_DIR={tmpdir}"
        " VALIDATION_STRINGENCY=SILENT"
        ") {log}"
    )
