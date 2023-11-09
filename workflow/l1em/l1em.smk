rule get_l1em:
    output:
        directory("resources/L1EM"),
    shell:
        "git clone https://github.com/FenyoLab/L1EM {output}"


rule bwa_index:
    input:
        rules.make_txome.output.genome_fa,
    output:
        fa="results/{txome}/bwa_index/genome.fa",
        idx=multiext(
            "results/{txome}/bwa_index/genome.fa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "results/{txome}/bwa_index.log",
    params:
        algorithm="bwtsw",
    conda:
        "l1em.yaml"
    shell:
        """
        cp {input} {output.fa}
        bwa index -a {params.algorithm} {output.fa} > {log} 2>&1
        """


rule build_l1em_ref:
    input:
        genome=rules.bwa_index.output.fa,
        index=rules.bwa_index.output.idx,
        l1em=rules.get_l1em.output,
    output:
        "results/{txome}/build_l1em_ref.log",
    conda:
        "l1em.yaml"
    shell:
        """
        # get full paths
        genome=$(readlink -f {input.genome})
        output=$(readlink -f {output})
        cd {input.l1em}
        bash generate_L1EM_fasta_and_index.sh $genome > $output 2>&1
        """


rule l1em:
    input:
        bam=rules.samtools_sort.output,
        bai=rules.samtools_index.output,
        l1em=rules.get_l1em.output,
        fa=rules.bwa_index.output.fa,
        index=rules.bwa_index.output.idx,
        l1em_ref=rules.build_l1em_ref.output,
    output:
        full_counts="results/{txome}/{tx_sim}/{te_sim}/l1em/{sample}/full_counts.txt",
        l1hs_counts="results/{txome}/{tx_sim}/{te_sim}/l1em/{sample}/l1hs_transcript_counts.txt",
        filter_fpm="results/{txome}/{tx_sim}/{te_sim}/l1em/{sample}/filter_L1HS_FPM.txt",
    conda:
        "l1em.yaml"
    log:
        "results/{txome}/{tx_sim}/{te_sim}/l1em/{sample}/l1em.log",
    threads: 16
    shell:
        """
        outdir=$(dirname {output.filter_fpm})
        touch {log}
        l1em=$(readlink -f {input.l1em})
        bam=$(readlink -f {input.bam})
        ref=$(readlink -f {input.fa})
        log=$(readlink -f {log})
        mkdir -p $outdir && cd $outdir

        trap "rm -rf G_of_R split_fqs idL1reads L1EM" EXIT

        sed -i 's/threads=[0-9]*/threads={threads}/g' $l1em/run_L1EM.sh
        bash -e $l1em/run_L1EM.sh $bam $l1em $ref > $log 2>&1
        """
