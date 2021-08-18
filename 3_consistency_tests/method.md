For this test, I compared assemblies from:
* Flye
* Raven
* Miniasm/Minipolish
* Trycycler (me)
* Trycycler (Jane)
* Trycycler (Guillaume)
* Trycycler (Kelly)
* Trycycler (Ben)
* Trycycler (Louise)
* Trycycler+Medaka+Pilon (me) - used as a proxy for the true sequence

This was to see how consistently Trycycler performs in different hands.




## Gathering assemblies

I first put the assemblies here:
```
/Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/consistency_tests/assemblies/full_assemblies
```

I didn't change them at all, except to add prefixes to the sequence names, like `flye_` or `kelly_`.




## Grouping and rotating

I then split the assemblies into per-replicon files and ensure that they were all on the same strand and with the same starting position. I did the grouping manually, but I used Trycycler cluster to help guide me:

```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/consistency_tests/assemblies/full_assemblies
mkdir reads
scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/Acinetobacter_baumannii_J9/rapid.fastq.gz reads/Acinetobacter_baumannii_J9.fastq.gz
scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/Citrobacter_koseri_MINF_9D/rapid.fastq.gz reads/Citrobacter_koseri_MINF_9D.fastq.gz
scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/Enterobacter_kobei_MSB1_1B/rapid.fastq.gz reads/Enterobacter_kobei_MSB1_1B.fastq.gz
scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/Haemophilus_M1C132_1/rapid.fastq.gz reads/Haemophilus_M1C132_1.fastq.gz
scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/Klebsiella_oxytoca_MSB1_2C/rapid.fastq.gz reads/Klebsiella_oxytoca_MSB1_2C.fastq.gz
scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/Klebsiella_variicola_INF345/rapid.fastq.gz reads/Klebsiella_variicola_INF345.fastq.gz

mkdir clustering
trycycler cluster --assemblies */Acinetobacter_baumannii_J9.fasta --reads reads/Acinetobacter_baumannii_J9.fastq.gz --out_dir clustering/Acinetobacter_baumannii_J9 2> clustering/Acinetobacter_baumannii_J9.out
trycycler cluster --assemblies */Citrobacter_koseri_MINF_9D.fasta --reads reads/Citrobacter_koseri_MINF_9D.fastq.gz --out_dir clustering/Citrobacter_koseri_MINF_9D 2> clustering/Citrobacter_koseri_MINF_9D.out
trycycler cluster --assemblies */Enterobacter_kobei_MSB1_1B.fasta --reads reads/Enterobacter_kobei_MSB1_1B.fastq.gz --out_dir clustering/Enterobacter_kobei_MSB1_1B 2> clustering/Enterobacter_kobei_MSB1_1B.out
trycycler cluster --assemblies */Haemophilus_M1C132_1.fasta --reads reads/Haemophilus_M1C132_1.fastq.gz --out_dir clustering/Haemophilus_M1C132_1 2> clustering/Haemophilus_M1C132_1.out
trycycler cluster --assemblies */Klebsiella_oxytoca_MSB1_2C.fasta --reads reads/Klebsiella_oxytoca_MSB1_2C.fastq.gz --out_dir clustering/Klebsiella_oxytoca_MSB1_2C 2> clustering/Klebsiella_oxytoca_MSB1_2C.out
trycycler cluster --assemblies */Klebsiella_variicola_INF345.fasta --reads reads/Klebsiella_variicola_INF345.fastq.gz --out_dir clustering/Klebsiella_variicola_INF345 2> clustering/Klebsiella_variicola_INF345.out
```

For each assembly, I looked for each 'true' replicon (i.e. replicons that appeared in my Trycycler assembly). Often assemblies were missing one or more contigs, or included extra contigs. When there were extra contigs, I tried to investigate their origin (e.g. cross-barcode contamination plasmid from another genome).

I then rotated all chromosomes to have a consistent starting position. I did this manually, but it was pretty easy because only non-Trycycler assemblies needed rotating (Trycycler gives consistent starting positions when it can find a starting gene).



## Distance matrix for each chromosome

I then made an MSA for each chromosome using Trycycler:
```bash
cat *.fasta > 2_all_seqs.fasta
trycycler msa -c .
rm 2_all_seqs.fasta
mv 3_msa.fasta msa.fasta
```

And used custom scripts to get variant sites and a PHYLIP distance matrix:
```bash
../../../../scripts/remove_invariant_positions_from_msa.py msa.fasta > msa_no_invariant.fasta
../../../../scripts/msa_to_distance_matrix.py msa_no_invariant.fasta > distances.phylip
```

And then made a tree from those distances in R:
```R
library(ape)
library(phangorn)
distances <- readDist('distances.phylip')
tree <- fastme.bal(distances)
tree$edge.length <- pmax(tree$edge.length, 0.0)  # set any negative branch lengths to zero
tree <- midpoint(tree)
write.tree(tree, 'distances.newick')
q(save="no")
```

I then did an unrooted layout:
```bash
~/Programs/Smoothtree/smoothtree-runner.py --random --seed 0 --angle_power 0 distances.newick tree_0.svg
~/Programs/Smoothtree/smoothtree-runner.py --random --seed 1 --angle_power 0 distances.newick tree_1.svg
~/Programs/Smoothtree/smoothtree-runner.py --random --seed 2 --angle_power 0 distances.newick tree_2.svg
~/Programs/Smoothtree/smoothtree-runner.py --random --seed 3 --angle_power 0 distances.newick tree_3.svg
~/Programs/Smoothtree/smoothtree-runner.py --random --seed 4 --angle_power 0 distances.newick tree_4.svg
~/Programs/Smoothtree/smoothtree-runner.py --random --seed 5 --angle_power 0 distances.newick tree_5.svg
```

I also normalised the distance matrix to give me per-Mbp differences:
```bash
genome_size=$(fast_count *_polished_*.fasta | cut -f3)
../../../../scripts/normalise_distance_matrix_to_mbp.py distances.phylip $genome_size > distances_per_mbp.phylip
```



## Assembly differences

I also wanted to characterise the differences between assemblies. E.g. how much substitutions vs indels, what size indels, how much indels were in homopolymers, etc.

```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/consistency_tests/assemblies/grouped_contigs
for c in Acinetobacter_baumannii_J9_chromosome Citrobacter_koseri_MINF_9D_chromosome Enterobacter_kobei_MSB1_1B_chromosome Haemophilus_M1C132_1_chromosome Klebsiella_oxytoca_MSB1_2C_chromosome Klebsiella_variicola_INF345_chromosome; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/consistency_tests/assemblies/grouped_contigs/"$c"
    printf "\n\n\n"$c"\n"
    ../../../../scripts/summarise_differences.py *_ben_*.fasta *_guillaume_*.fasta *_jane_*.fasta *_kelly_*.fasta *_louise_*.fasta *_ryan_*.fasta
done > trycycler_vs_trycycler_differences
```

```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/consistency_tests/assemblies/grouped_contigs
for c in Acinetobacter_baumannii_J9_chromosome Enterobacter_kobei_MSB1_1B_chromosome Haemophilus_M1C132_1_chromosome Klebsiella_oxytoca_MSB1_2C_chromosome Klebsiella_variicola_INF345_chromosome; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/consistency_tests/assemblies/grouped_contigs/"$c"
    printf "\n\n\n"$c"\n"
    ../../../../scripts/summarise_differences.py *_miniasm_*.fasta *_raven_*.fasta *_flye_*.fasta
done > other_vs_other_differences
```

```bash
for c in Acinetobacter_baumannii_J9_chromosome Citrobacter_koseri_MINF_9D_chromosome Enterobacter_kobei_MSB1_1B_chromosome Haemophilus_M1C132_1_chromosome Klebsiella_oxytoca_MSB1_2C_chromosome Klebsiella_variicola_INF345_chromosome; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/consistency_tests/assemblies/grouped_contigs/"$c"
    printf "\n\n\n"$c"\n"
    ../../../../scripts/summarise_differences_with_motifs.py *_ben_*.fasta.gz *_polished_*.fasta.gz
    ../../../../scripts/summarise_differences_with_motifs.py *_guillaume_*.fasta.gz *_polished_*.fasta.gz
    ../../../../scripts/summarise_differences_with_motifs.py *_jane_*.fasta.gz *_polished_*.fasta.gz
    ../../../../scripts/summarise_differences_with_motifs.py *_kelly_*.fasta.gz *_polished_*.fasta.gz
    ../../../../scripts/summarise_differences_with_motifs.py *_louise_*.fasta.gz *_polished_*.fasta.gz
    ../../../../scripts/summarise_differences_with_motifs.py *_ryan_*.fasta.gz *_polished_*.fasta.gz
done
```
