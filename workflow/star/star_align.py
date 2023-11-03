import os
import tempfile

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


fq1 = snakemake.input.get("fq1")
assert fq1 is not None, "input-> fq1 is a required input parameter"
fq1 = (
    [snakemake.input.fq1]
    if isinstance(snakemake.input.fq1, str)
    else snakemake.input.fq1
)
fq2 = snakemake.input.get("fq2")
if fq2:
    fq2 = (
        [snakemake.input.fq2]
        if isinstance(snakemake.input.fq2, str)
        else snakemake.input.fq2
    )
    assert len(fq1) == len(
        fq2
    ), "input-> equal number of files required for fq1 and fq2"
input_str_fq1 = ",".join(fq1)
input_str_fq2 = ",".join(fq2) if fq2 is not None else ""
input_str = " ".join([input_str_fq1, input_str_fq2])


if fq1[0].endswith(".gz"):
    readcmd = "--readFilesCommand gunzip -c"
elif fq1[0].endswith(".bz2"):
    readcmd = "--readFilesCommand bunzip2 -c"
else:
    readcmd = ""

prefix = os.path.dirname(snakemake.output.get("genome_bam"))

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "STAR "
        " --runThreadN {snakemake.threads}"
        " --genomeDir {snakemake.input.idx}"
        " --sjdbGTFfile {snakemake.input.gtf}"
        " --readFilesIn {input_str}"
        " {readcmd}"
        " --outSAMunmapped Within"
        " --outSAMattributes NH HI AS NM MD"
        " --outSAMstrandField intronMotif"
        " {extra}"
        " --outTmpDir {tmpdir}/STARtmp"
        " --outFileNamePrefix {prefix}/"
        " --quantMode TranscriptomeSAM"
        " > {snakemake.output.log}"
        " {log}"
    )
