import os
import tempfile

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

sample = os.path.basename(snakemake.input.read)
sample = sample.replace(".fastq.gz", "")
sample = sample.replace(".fastq", "")

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "("
        " R_IN={snakemake.input.read};"
        " for INDEX in {snakemake.params.indices}; do"
        "  REF=`basename $INDEX`;"
        "  REF_DIR={tmpdir}/$REF;"
        "  mkdir -p $REF_DIR;"
        "  BAM_OUT=$REF_DIR/{sample}_out.bam;"
        "  R_OUT=$REF_DIR/{sample}.fastq.gz;"
        "  bwa mem "
        "   -t {snakemake.threads}"
        "   $INDEX"
        "   $R_IN"
        "  |"
        "  samtools view"
        "   -o $BAM_OUT;"
        "  samtools collate"
        "   -u -O $BAM_OUT"
        "  |"
        "  samtools fastq"
        "   -o $R_OUT"
        "   {snakemake.params.keep_param}"
        "   -0 /dev/null"
        "   -s /dev/null"
        "   -t"
        "   -n;"
        "  R_IN=$R_OUT;"
        " done;"
        " mv $R_OUT {snakemake.output.read};"
        " ) {log}"
    )
