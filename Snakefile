split_number = 50

localrules: subsample_sequences,diff,split_nextclade_results,collect_nextclade_results,unzip_split,split_sequences,download_nextclade_binary,download_nextclade_dataset,download_sequences

wildcard_constraints:
    version="[^_]*",
    part="[^_]*",
    types="[^_]*",

rule all:
    input: "results/nextclade.tsv

rule download_sequences:
    output: "data/sequences.fasta.xz"
    params:
        url = "https://github.com/robert-koch-institut/SARS-CoV-2-Sequenzdaten_aus_Deutschland/raw/master/SARS-CoV-2-Sequenzdaten_Deutschland.fasta.xz"
    shell: "wget {params.url} {output}"

# rule subsample_sequences:
#     input: rules.download_sequences.output
#     output: "data/subsample.fasta.gz"
#     shell: "xz -dc {input} | seqkit sample -p {subsample_ratio} -o {output}"

rule xz_to_gz:
    input: "data/sequences.fasta.xz"
    output: temp("data/sequences.fasta.gz")
    shell: "xz -dc {input} | gzip -c > {output}"

rule split_sequences:
    input: "data/sequences.fasta.gz"
    output: temp(expand("split/sequences.part_{part:03d}.fasta.gz", part=range(1,split_number+1)))
    shell:
        """
        seqkit split2 {input} -p {split_number} -O split && \
        """

rule unzip_split:
    input: "split/sequences_subsample.part_{part}.fasta.gz"
    output: temp("split/sequences_subsample.part_{part}.fasta")
    shell: "gunzip -c {input} > {output}"

rule download_nextclade_dataset:        
    output: directory("data/nextclade_dataset")
    shell: "nextclade dataset get --name='sars-cov-2' --output-dir={output}"

rule run_nextclade:
    input:
        sequences = "split/all_subsample.part_{part}.fasta",
        dataset = rules.download_nextclade_dataset.output,
    output:
        output_tsv = temp("results/nextclade_results_{part}.tsv"),
    params:
        output_alignments = "data/alignments_{part}",
    threads: 8
    shell:
        """
        nextclade run \
            -j{threads} \
            --in-order \
            -i {input.sequences} \
            --input-dataset {input.dataset} \
            -t {output.output_tsv} \
            -d {params.output_alignments};
        rm -r {params.output_alignments}; \
        """

rule collect_nextclade_results:
    input: expand("results/nextclade_results_{part:03d}.tsv", part=range(1,split_number+1))
    output: "nextclade.tsv"
    shell: "keep-header {input} -- cat | dos2unix > {output}"
