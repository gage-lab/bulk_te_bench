rule test_sim_counts:
    input:
        txome_fa=rules.make_txome.output.txome_fa,
    output:
        counts="results/{txome}/test_sim/true_counts.tsv",
    run:
        from collections import defaultdict
        from Bio import SeqIO
        import pandas as pd

        TX = [
            "ENST00000401850.5",
            "ENST00000401959.6",
            "ENST00000401975.5",
            "ENST00000401994.5",
            "L1HS_dup1",
        ]
        counts = defaultdict(list)
        for tx in SeqIO.parse(input.txome_fa, "fasta"):
            if tx.id not in TX:
                continue
            counts["tx_id"].append(tx.id)
            for sample in range(0, 2):
                if "ENS" in tx.id:
                    counts[sample].append(20 * len(tx.seq) // 100)
                elif "L1HS" in tx.id:
                    counts[sample].append(sample + 1 * len(tx.seq) // 100)
                else:
                    counts[sample].append(0)
        pd.DataFrame(counts).to_csv(output.counts, sep="\t", index=False)


checkpoint simulate_reads:
    input:
        txome_fa=rules.make_txome.output.txome_fa,
        counts="results/{txome}/{sim}/true_counts.tsv",
    output:
        reads=directory("results/{txome}/{sim}/reads"),
    conda:
        "simulate.yaml"
    threads: 16
    params:
        strand_specific=lambda wc: config["txomes"][wc.txome]["simulations"][wc.sim][
            "strand_specific"
        ],
        readlen=lambda wc: config["txomes"][wc.txome]["simulations"][wc.sim]["readlen"],
    log:
        "results/{txome}/{sim}/simulate_reads.log",
    script:
        "simulate_reads.py"
