input_file="./samples.csv"
adapter_file="./adapters/Truseq3.PE.fa"
fasta_file="/path/to/genome.fa"
gtf_file="/path/to/genes.gtf"
bt2_index_file="/path/to/genome"
bed_file="/path/to/genes.bed"
blacklist_file="/path/to/hg38-blacklist.bed"
out_dir="./results"
work_dir="./work"

nextflow run main.nf --indir ${input_dir} \
    --adapter ${adapter_file} \
    --seq_len 42 \
    --fasta ${fasta_file} \
    --gtf ${gtf_file} \
    --bt2_index ${bt2_index_file} \
    --gene_bed ${bed_file} \
    --blacklist ${blacklist_file} \
    --macs_gsize "2.7e9" \
    --outdir ${out_dir} \
    --skip_plot_profile \
    -work-dir ${work_dir} \
    -profile docker \
    -resume