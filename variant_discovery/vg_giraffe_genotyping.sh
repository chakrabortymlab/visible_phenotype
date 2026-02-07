# vg giraffe pangenome genotyping using long reads
# Author: Alex Samano, 2026


#module load VG/1.65.0

#index graph for long read genotyping
vg autoindex --workflow lr-giraffe -g biol450_genomes.gfa -t 24 -p biol450_genomes_234X

# convert to xg
vg convert -x biol450_genomes_234X.giraffe.gbz > biol450_genomes_234X.giraffe.xg


# run given the strain ID
strain=$1
reads="/scratch/group/chakraborty_lab/alex/dmel/biol450/barcoded_pool2/"$strain"/*fq.gz"

# map reads to pangenome, output GAF
vg giraffe -b r10 -x biol450_genomes_234X.giraffe.xg -m biol450_genomes_234X.longread.withzip.min -d biol450_genomes_234X.dist  -f $reads -t 24 -o gaf > $strain.biol450_genomes_234X.gaf

# pack alignments
vg pack -x biol450_genomes_234X.giraffe.xg -a $strain.biol450_genomes_234X.gaf -Q 5 -t 24 -o $strain.biol450_genomes_234X.pack

#generate VCF with each variant in pangenome genotyped
vg call biol450_genomes_234X.giraffe.xg -k $strain.biol450_genomes_234X.pack -a -t 24 > $strain.biol450_genomes_234X.vcf