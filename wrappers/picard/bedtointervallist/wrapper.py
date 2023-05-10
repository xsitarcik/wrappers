import tempfile

from snakemake.shell import shell

log = snakemake.log_fmt_shell()

unique_flag = "UNIQUE=true" if snakemake.params.get("unique_only", False) else ""
sort_flag = "SORT=true" if snakemake.params.get("do_sorting", False) else ""

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "picard BedToIntervalList"
        " INPUT={snakemake.input.bed}"
        " SEQUENCE_DICTIONARY={snakemake.input.seq_dict}"
        " OUTPUT={snakemake.output.intervals}"
        " {sort_flag}"
        " {unique_flag}"
        " TMP_DIR={tmpdir}"
        " VALIDATION_STRINGENCY=SILENT"
        " {log}"
    )
