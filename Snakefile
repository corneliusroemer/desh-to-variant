split_number = 50

localrules: xz_to_gz,collect_nextclade_results,unzip_split,split_sequences,download_nextclade_dataset,download_sequences

wildcard_constraints:
    version="[^_]*",
    part="[^_]*",
    types="[^_]*",

rule all:
    input: "results/nextclade.tsv"

rule download_sequences:
    output: "data/sequences.fasta.xz"
    params:
        url = "https://github.com/robert-koch-institut/SARS-CoV-2-Sequenzdaten_aus_Deutschland/raw/master/SARS-CoV-2-Sequenzdaten_Deutschland.fasta.xz"
    shell: "curl -L {params.url} >{output}"

rule split_sequences:
    input: "data/sequences.fasta.xz"
    output: temp(expand("split/stdin.part_{part:03d}", part=range(1,split_number+1)))
    shell:
        """
        xzcat -dc {input} | seqkit split2 /dev/stdin -p {split_number} -O split
        """

# rule unzip_split:
#     input: "split/sequences.part_{part}.fasta.gz"
#     output: temp("split/sequences.part_{part}.fasta")
#     shell: "gunzip -c {input} > {output}"

rule download_nextclade_dataset:        
    output: directory("data/nextclade_dataset")
    shell: "nextclade dataset get --name='sars-cov-2' --output-dir={output}"

rule run_nextclade:
    input:
        sequences = "split/stdin.part_{part}",
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
            -i {input.sequences} \
            --input-dataset {input.dataset} \
            -t {output.output_tsv} \
            --verbosity warn \
            -d {params.output_alignments};
        rm -r {params.output_alignments}; \
        """

rule collect_nextclade_results:
    input: expand("results/nextclade_results_{part:03d}.tsv", part=range(1,split_number+1))
    output: "results/nextclade.tsv"
    shell: "keep-header {input} -- cat | dos2unix > {output}"
