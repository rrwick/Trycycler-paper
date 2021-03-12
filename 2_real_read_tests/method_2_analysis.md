## Read analysis

I analysed the reads on Nectar, both before filtering:
```bash
cd ~/trycycler_real_read_tests
fast_count */raw_reads/ligation.fastq.gz
fast_count */raw_reads/rapid.fastq.gz
fast_count */raw_reads/illumina_[ab]_*.fastq.gz
```
and after filtering:
```bash
cd ~/trycycler_real_read_tests
fast_count */ligation.fastq.gz
fast_count */rapid.fastq.gz
fast_count */illumina_[ab]_*.fastq.gz
```

For Illumina read sets, I summed the counts/totals for all files in the set (1 and 2 for pre-filtered reads; 1, 2 and u for post-filtered reads).




## Resource usage

For Trycycler RAM usage, I took the largest value of all the steps:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /home/ubuntu/trycycler_real_read_tests/"$s"/05_trycycler_ligation
    printf $s"\t"
    for f in $(find . -name "*.time" | grep -v "medaka" | sort); do grep "Maximum resident set size" $f | grep -oP "\d+"; done | sort -n | tail -n1
done
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /home/ubuntu/trycycler_real_read_tests/"$s"/06_trycycler_rapid
    printf $s"\t"
    for f in $(find . -name "*.time" | grep -v "medaka" | sort); do grep "Maximum resident set size" $f | grep -oP "\d+"; done | sort -n | tail -n1
done
```

For Trycycler time, I took the sum of all the steps (used Excel to total them up):
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /home/ubuntu/trycycler_real_read_tests/"$s"/05_trycycler_ligation
    echo "$s"
    for f in $(find . -name "*.time" | grep -v "medaka" | sort); do grep "Elapsed (wall clock) time" $f | grep -oP "[\d:\.]{2,}"; done
    printf "\n\n"
done
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /home/ubuntu/trycycler_real_read_tests/"$s"/06_trycycler_rapid
    echo "$s"
    for f in $(find . -name "*.time" | grep -v "medaka" | sort); do grep "Elapsed (wall clock) time" $f | grep -oP "[\d:\.]{2,}"; done
    printf "\n\n"
done
```







## Transferring the assemblies to my computer

```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies
mkdir entire_assemblies
mkdir chromosome_only
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    mkdir entire_assemblies/"$s"
    mkdir chromosome_only/"$s"
done

for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/entire_assemblies/"$s"

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/01_flye_ligation/assembly.fasta flye_ligation.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/01_flye_ligation/assembly_graph.gfa flye_ligation.gfa
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/01_flye_ligation/medaka.fasta flye_medaka_ligation.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/01_flye_ligation/polished.fasta flye_medaka_pilon_ligation.fasta

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/02_flye_rapid/assembly.fasta flye_rapid.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/02_flye_rapid/assembly_graph.gfa flye_rapid.gfa
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/02_flye_rapid/medaka.fasta flye_medaka_rapid.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/02_flye_rapid/polished.fasta flye_medaka_pilon_rapid.fasta

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/03_unicycler_hybrid_ligation/assembly.fasta unicycler_hybrid_ligation.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/03_unicycler_hybrid_ligation/assembly.gfa unicycler_hybrid_ligation.gfa
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/03_unicycler_hybrid_ligation/004_bridges_applied.gfa unicycler_hybrid_ligation_bridged.gfa

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/04_unicycler_hybrid_rapid/assembly.fasta unicycler_hybrid_rapid.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/04_unicycler_hybrid_rapid/assembly.gfa unicycler_hybrid_rapid.gfa
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/04_unicycler_hybrid_rapid/004_bridges_applied.gfa unicycler_hybrid_rapid_bridged.gfa

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/05_trycycler_ligation/assembly.fasta trycycler_ligation.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/05_trycycler_ligation/medaka.fasta trycycler_medaka_ligation.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/05_trycycler_ligation/polished.fasta trycycler_medaka_pilon_ligation.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/05_trycycler_ligation/hybrid_polished.fasta trycycler_medaka_hybrid_pilon_ligation.fasta

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/06_trycycler_rapid/assembly.fasta trycycler_rapid.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/06_trycycler_rapid/medaka.fasta trycycler_medaka_rapid.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/06_trycycler_rapid/polished.fasta trycycler_medaka_pilon_rapid.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/06_trycycler_rapid/hybrid_polished.fasta trycycler_medaka_hybrid_pilon_rapid.fasta

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/07_raven_ligation/raven.fasta raven_ligation.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/07_raven_ligation/raven.gfa raven_ligation.gfa
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/07_raven_ligation/medaka.fasta raven_medaka_ligation.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/07_raven_ligation/polished.fasta raven_medaka_pilon_ligation.fasta

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/08_raven_rapid/raven.fasta raven_rapid.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/08_raven_rapid/raven.gfa raven_rapid.gfa
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/08_raven_rapid/medaka.fasta raven_medaka_rapid.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/08_raven_rapid/polished.fasta raven_medaka_pilon_rapid.fasta

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/09_miniasm_ligation/miniasm.fasta miniasm_ligation.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/09_miniasm_ligation/miniasm.gfa miniasm_ligation.gfa
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/09_miniasm_ligation/medaka.fasta miniasm_medaka_ligation.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/09_miniasm_ligation/polished.fasta miniasm_medaka_pilon_ligation.fasta

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/10_miniasm_rapid/miniasm.fasta miniasm_rapid.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/10_miniasm_rapid/miniasm.gfa miniasm_rapid.gfa
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/10_miniasm_rapid/medaka.fasta miniasm_medaka_rapid.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_real_read_tests/"$s"/10_miniasm_rapid/polished.fasta miniasm_medaka_pilon_rapid.fasta
done
```

I then manually inspected the graphs in Bandage to determine whether or not the chromosome was complete.

I then manually made chromosome-only FASTA files for each of the completed Unicycler hybrid, Flye and Trycycler assemblies. I did this by:
* Copying the FASTA and adding `_chromosome` to the end of its name.
* Removing in-sequence line breaks: `seqtk seq $f > temp; mv temp $f`
* Deleting non-chromosome contigs
* Manually adjusting the start point of each assembly to match the Unicycler hybrid start point (flipping to reverse strand if necessary with `seqtk seq -r $f > temp; mv temp $f`)



For Flye assemblies, circularisation could be an issue, so I created versions of each assembly where I manually fixed the circularisation (using the polished Trycycler assembly as a point of comparison). This had to be done before changing the starting position. Here were the adjustments I made:
* Acinetobacter_baumannii_J9 ligation:      missing 1 bp (T)
* Acinetobacter_baumannii_J9 rapid:         missing 2 bp (AA)
* Citrobacter_koseri_MINF_9D ligation:      good
* Citrobacter_koseri_MINF_9D rapid:         missing 4 bp (TAAT)
* Enterobacter_kobei_MSB1_1B ligation:      extra 1 bp (T)
* Enterobacter_kobei_MSB1_1B rapid:         missing 4 bp (AGTT)
* Haemophilus_M1C132_1 ligation:            good
* Haemophilus_M1C132_1 rapid:               extra 4 bp (TTTA)
* Klebsiella_oxytoca_MSB1_2C ligation:      missing 9 bp (ATAGGGAGT)
* Klebsiella_oxytoca_MSB1_2C rapid:         missing 1 bp (C)
* Klebsiella_variicola_INF345 ligation:     missing 4 bp (AAAT)
* Klebsiella_variicola_INF345 rapid:        missing 13 bp (CCAGGCATCAAAT)





Now that I have matched up assemblies (same replicon, same strand, same start), I can do an Edlib global alignment to do pairwise comparisons.

Flye ligation vs Flye rapid:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py flye_ligation_chromosome.fasta flye_rapid_chromosome.fasta
done
```

Flye+Medaka ligation vs Flye+Medaka rapid:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py flye_medaka_ligation_chromosome.fasta flye_medaka_rapid_chromosome.fasta
done
```

Flye+Medaka+Pilon ligation vs Flye+Medaka+Pilon rapid:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py flye_medaka_pilon_ligation_chromosome.fasta flye_medaka_pilon_rapid_chromosome.fasta
done
```

Flye ligation (fixed circ) vs Flye rapid (fixed circ):
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py flye_ligation_chromosome_fixed_circularisation.fasta flye_rapid_chromosome_fixed_circularisation.fasta
done
```

Flye+Medaka ligation (fixed circ) vs Flye+Medaka rapid (fixed circ):
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py flye_medaka_ligation_chromosome_fixed_circularisation.fasta flye_medaka_rapid_chromosome_fixed_circularisation.fasta
done
```

Flye+Medaka+Pilon ligation (fixed circ) vs Flye+Medaka+Pilon rapid (fixed circ):
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py flye_medaka_pilon_ligation_chromosome_fixed_circularisation.fasta flye_medaka_pilon_rapid_chromosome_fixed_circularisation.fasta
done
```

Trycycler ligation vs Trycycler rapid:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py trycycler_ligation_chromosome.fasta trycycler_rapid_chromosome.fasta
done
```

Trycycler+Medaka ligation vs Trycycler+Medaka rapid:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py trycycler_medaka_ligation_chromosome.fasta trycycler_medaka_rapid_chromosome.fasta
done
```

Trycycler+Medaka+Pilon ligation vs Trycycler+Medaka+Pilon rapid:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py trycycler_medaka_pilon_ligation_chromosome.fasta trycycler_medaka_pilon_rapid_chromosome.fasta
done
```

Unicycler hybrid ligation vs Unicycler hybrid rapid:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py unicycler_hybrid_ligation_chromosome.fasta unicycler_hybrid_rapid_chromosome.fasta
done
```

Miniasm ligation vs Miniasm rapid:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py miniasm_ligation_chromosome.fasta miniasm_rapid_chromosome.fasta
done
```

Miniasm+Medaka ligation vs Miniasm+Medaka rapid:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py miniasm_medaka_ligation_chromosome.fasta miniasm_medaka_rapid_chromosome.fasta
done
```

Miniasm+Medaka+Pilon ligation vs Miniasm+Medaka+Pilon rapid:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py miniasm_medaka_pilon_ligation_chromosome.fasta miniasm_medaka_pilon_rapid_chromosome.fasta
done
```

Raven ligation vs Raven rapid:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py raven_ligation_chromosome.fasta raven_rapid_chromosome.fasta
done
```

Raven+Medaka ligation vs Raven+Medaka rapid:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py raven_medaka_ligation_chromosome.fasta raven_medaka_rapid_chromosome.fasta
done
```

Raven+Medaka+Pilon ligation vs Raven+Medaka+Pilon rapid:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py raven_medaka_pilon_ligation_chromosome.fasta raven_medaka_pilon_rapid_chromosome.fasta
done
```




## Repeat analysis

I was also curious what the longest repeat was in each chromosome. I did this by downloading the Unicycler hybrid Nanopore bridged graph (all of which had a complete chromosome) and looked for the longest bridge segment that was in the chromosomal loop.
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/entire_assemblies/"$s"
    printf $s"\t"
    ../../../../scripts/find_longest_chromosomal_repeat.py unicycler_hybrid_ligation_bridged.gfa
    printf $s"\t"
    ../../../../scripts/find_longest_chromosomal_repeat.py unicycler_hybrid_rapid_bridged.gfa
done
```
If ligation and rapid gave two different values, I took the smaller of the two.


I also did the same thing using MUMmer:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    nucmer --maxmatch --nosimplify --prefix=seq_seq trycycler_medaka_ligation_chromosome.fasta trycycler_medaka_ligation_chromosome.fasta 2> /dev/null
    show-coords -H -T -r seq_seq.delta | awk '{if ($2-$1 < 1000000 && $2-$1 > 1000 && $1 != $3 && $2 != $4) print $0"\t"($2-$1+1);}' | sort -nk10,10 | tail -n1
    rm seq_seq.delta
        printf $s"\t"
    nucmer --maxmatch --nosimplify --prefix=seq_seq trycycler_medaka_rapid_chromosome.fasta trycycler_medaka_rapid_chromosome.fasta 2> /dev/null
    show-coords -H -T -r seq_seq.delta | awk '{if ($2-$1 < 1000000 && $2-$1 > 1000 && $1 != $3 && $2 != $4) print $0"\t"($2-$1+1);}' | sort -nk10,10 | tail -n1
    rm seq_seq.delta
done
```





## ALE

For ALE, I used the unfiltered (before fastp) Illumina reads and I aligned them using different Bowtie2 parameters which were consistent with the ALE docs. I also used my chromosome-only FASTA files, because the presence/absence of other contigs really changes the ALE score.

```bash
# Since all real genomes had similar insert size distributions, I used these values for all of them:
insert_min=158
insert_max=723

for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /home/ubuntu/trycycler_real_read_tests/chromosome_only/"$s"
    i1=/home/ubuntu/trycycler_real_read_tests/"$s"/raw_reads/illumina_1.fastq.gz
    i2=/home/ubuntu/trycycler_real_read_tests/"$s"/raw_reads/illumina_2.fastq.gz
    for a in *.fasta; do
        ale="${a/fasta/ale}"
        bowtie2-build "$a" "$a" > /dev/null
        bowtie2 -1 "$i1" -2 "$i2" -x "$a" --threads 32 -I "$insert_min" -X "$insert_max" -a --fr --end-to-end --very-sensitive | samtools sort > illumina_alignments.bam
        ALE illumina_alignments.bam "$a" "$ale"
        rm illumina_alignments.bam
        rm "$ale".param
        find . -name "*.bt2" -delete
    done
done
```



## Reapr

For REAPR, I followed the 'Any size genome, where only one library of reads is available' instructions from their manual.

```bash
export PATH=/home/ubuntu/Reapr_1.0.18:/home/ubuntu/Reapr_1.0.18/src:"$PATH"

for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /home/ubuntu/trycycler_real_read_tests/chromosome_only/"$s"
    
    i1=/home/ubuntu/trycycler_real_read_tests/"$s"/raw_reads/illumina_1.fastq.gz
    i2=/home/ubuntu/trycycler_real_read_tests/"$s"/raw_reads/illumina_2.fastq.gz

    # REAPR needs all reads to be the same length. So we chop the reads down to 145 bp and then keep pairs where both reads are 145 bp.
    seqtk trimfq -L 145 "$i1" > reads_1.fastq
    seqtk trimfq -L 145 "$i2" > reads_2.fastq
    seqtk mergepe reads_1.fastq reads_2.fastq | paste - - - - - - - - | awk 'BEGIN {FS="\t"}; {if (length($2)==145 && length($6)==145) print $0;}' > temp
    cat temp | cut -f 1-4 | tr "\t" "\n" > reads_1.fastq
    cat temp | cut -f 5-8 | tr "\t" "\n" > reads_2.fastq
    rm temp

    for a in *.fasta; do
        reapr smaltmap -n 32 "$a" reads_1.fastq reads_2.fastq long_mapped.bam
        reapr perfectmap "$a" reads_1.fastq reads_2.fastq 450 perfect  # I used a mean insert size of 450 for all isolates as they were pretty consistent.
        reapr pipeline "$a" long_mapped.bam "$a"_reapr perfect
        rm long_mapped.bam long_mapped.bam.bai
        rm perfect.hist perfect.perfect_cov.gz perfect.perfect_cov.gz.tbi
        rm "$a".fai "$a"_reapr.run-pipeline.sh
    done
    rm reads_1.fastq reads_2.fastq
done
```

I found REAPR's results to be quite confusing. It almost always gave errors, even for genomes that I'm pretty sure are in good shape (e.g. the Medaka+Pilon _A. baumannii_ genomes). It also gave tons of warnings for each genome. And the 'error free bases' metric which sounds nice gives low numbers (e.g. 97%) so I don't know how to interpret them.

Another issue was REAPR was inconsistent - it gave different answers for identical sequences. I saw this for the Flye sequences which didn't need circularisation correction (Citrobacter_koseri_MINF_9D ligation and Haemophilus_M1C132_1 ligation) - these have the same sequence for the non-corrected and corrected sequences, but REAPR was giving slightly different results. I thought it might be due to the random seed in SMALT: it was set to use `smalt map -r 0` which seeds based on the current time. I hacked the `task_smaltmap.pl` file of REAPR to instead use `smalt map -r 1` (a constant seed which is actually what they did in the REAPR paper), but that didn't fix the problem.

For these reasons, I decided against including REAPR results in my paper, but I'll leave the methods here for completeness.



## Ideel

I ran [ideel](https://github.com/mw55309/ideel) on each of the genomes by:
* Downloading UniProt/TrEMBL release 2020_05: https://www.uniprot.org/statistics/TrEMBL
* Building a Diamond database (v2.0.4)
* Cloning the ideel repo
* Making a directory called 'genomes' with the assemblies in there with .fa file extension
* Editing the `Snakefile` to point to the Diamond database
* Running it: `snakemake --cores 16`

This made a `.data` file for each genome which contained one row for each protein annotation and two columns:
* The amino acid length in the genome
* The amino acid length in the best UniProt/TrEMBL match

I then got the fraction of these where the genome's aa sequence was at least 95% the length of the best UniProt/TrEMBL match:
```bash
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    cd /home/ubuntu/trycycler_real_read_tests/chromosome_only/"$s"
    printf "\n\n"$s"\n"
    for i in *.ideel; do
        full_count=$(cat $i | wc -l)
        long_count=$(cat $i | awk '{if ($1/$2 >= 0.95) print $0;}' | wc -l)
        long_fraction=$(echo "scale=5; "$long_count"/"$full_count | bc)
        print $i"\t"$long_fraction
    done
done
```



## Error positions

First I ran a script to get the Flye starting position, as that's needed for the next command:
```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/Acinetobacter_baumannii_J9
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Acinetobacter_baumannii_J9/flye_ligation.fasta flye_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Acinetobacter_baumannii_J9/flye_rapid.fasta flye_rapid_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Acinetobacter_baumannii_J9/flye_medaka_ligation.fasta flye_medaka_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Acinetobacter_baumannii_J9/flye_medaka_rapid.fasta flye_medaka_rapid_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Acinetobacter_baumannii_J9/flye_medaka_pilon_ligation.fasta flye_medaka_pilon_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Acinetobacter_baumannii_J9/flye_medaka_pilon_rapid.fasta flye_medaka_pilon_rapid_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/Citrobacter_koseri_MINF_9D
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Acinetobacter_baumannii_J9/flye_ligation.fasta flye_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Acinetobacter_baumannii_J9/flye_rapid.fasta flye_rapid_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Acinetobacter_baumannii_J9/flye_medaka_ligation.fasta flye_medaka_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Acinetobacter_baumannii_J9/flye_medaka_rapid.fasta flye_medaka_rapid_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Acinetobacter_baumannii_J9/flye_medaka_pilon_ligation.fasta flye_medaka_pilon_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Acinetobacter_baumannii_J9/flye_medaka_pilon_rapid.fasta flye_medaka_pilon_rapid_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/Enterobacter_kobei_MSB1_1B
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Enterobacter_kobei_MSB1_1B/flye_ligation.fasta flye_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Enterobacter_kobei_MSB1_1B/flye_rapid.fasta flye_rapid_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Enterobacter_kobei_MSB1_1B/flye_medaka_ligation.fasta flye_medaka_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Enterobacter_kobei_MSB1_1B/flye_medaka_rapid.fasta flye_medaka_rapid_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Enterobacter_kobei_MSB1_1B/flye_medaka_pilon_ligation.fasta flye_medaka_pilon_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Enterobacter_kobei_MSB1_1B/flye_medaka_pilon_rapid.fasta flye_medaka_pilon_rapid_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/Haemophilus_M1C132_1
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Haemophilus_M1C132_1/flye_ligation.fasta flye_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Haemophilus_M1C132_1/flye_rapid.fasta flye_rapid_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Haemophilus_M1C132_1/flye_medaka_ligation.fasta flye_medaka_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Haemophilus_M1C132_1/flye_medaka_rapid.fasta flye_medaka_rapid_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Haemophilus_M1C132_1/flye_medaka_pilon_ligation.fasta flye_medaka_pilon_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Haemophilus_M1C132_1/flye_medaka_pilon_rapid.fasta flye_medaka_pilon_rapid_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/Klebsiella_oxytoca_MSB1_2C
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_oxytoca_MSB1_2C/flye_ligation.fasta flye_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_oxytoca_MSB1_2C/flye_rapid.fasta flye_rapid_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_oxytoca_MSB1_2C/flye_medaka_ligation.fasta flye_medaka_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_oxytoca_MSB1_2C/flye_medaka_rapid.fasta flye_medaka_rapid_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_oxytoca_MSB1_2C/flye_medaka_pilon_ligation.fasta flye_medaka_pilon_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_oxytoca_MSB1_2C/flye_medaka_pilon_rapid.fasta flye_medaka_pilon_rapid_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/Klebsiella_variicola_INF345
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_variicola_INF345/flye_ligation.fasta flye_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_variicola_INF345/flye_rapid.fasta flye_rapid_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_variicola_INF345/flye_medaka_ligation.fasta flye_medaka_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_variicola_INF345/flye_medaka_rapid.fasta flye_medaka_rapid_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_variicola_INF345/flye_medaka_pilon_ligation.fasta flye_medaka_pilon_ligation_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_variicola_INF345/flye_medaka_pilon_rapid.fasta flye_medaka_pilon_rapid_chromosome.fasta
```

This will generate the data that my R script can then visualise:
```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/Acinetobacter_baumannii_J9
../../../../scripts/error_positions.py flye_ligation_chromosome.fasta flye_rapid_chromosome.fasta 493819 3586520 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_ligation_chromosome_fixed_circularisation.fasta flye_rapid_chromosome_fixed_circularisation.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_medaka_ligation_chromosome.fasta flye_medaka_rapid_chromosome.fasta 493891 3587557 > flye_medaka_chromosome.errors
../../../../scripts/error_positions.py flye_medaka_ligation_chromosome_fixed_circularisation.fasta flye_medaka_rapid_chromosome_fixed_circularisation.fasta > flye_medaka_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_medaka_pilon_ligation_chromosome.fasta flye_medaka_pilon_rapid_chromosome.fasta 493883 3587582 > flye_medaka_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_medaka_pilon_ligation_chromosome_fixed_circularisation.fasta flye_medaka_pilon_rapid_chromosome_fixed_circularisation.fasta > flye_medaka_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_ligation_chromosome.fasta trycycler_rapid_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_medaka_ligation_chromosome.fasta trycycler_medaka_rapid_chromosome.fasta > trycycler_medaka_chromosome.errors
../../../../scripts/error_positions.py trycycler_medaka_pilon_ligation_chromosome.fasta trycycler_medaka_pilon_rapid_chromosome.fasta > trycycler_medaka_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_ligation_chromosome.fasta unicycler_hybrid_rapid_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/Citrobacter_koseri_MINF_9D
../../../../scripts/error_positions.py flye_ligation_chromosome.fasta flye_rapid_chromosome.fasta 4622071 1893123 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_ligation_chromosome_fixed_circularisation.fasta flye_rapid_chromosome_fixed_circularisation.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_medaka_ligation_chromosome.fasta flye_medaka_rapid_chromosome.fasta 4623420 1893801 > flye_medaka_chromosome.errors
../../../../scripts/error_positions.py flye_medaka_ligation_chromosome_fixed_circularisation.fasta flye_medaka_rapid_chromosome_fixed_circularisation.fasta > flye_medaka_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_medaka_pilon_ligation_chromosome.fasta flye_medaka_pilon_rapid_chromosome.fasta 4623342 1893782 > flye_medaka_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_medaka_pilon_ligation_chromosome_fixed_circularisation.fasta flye_medaka_pilon_rapid_chromosome_fixed_circularisation.fasta > flye_medaka_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_ligation_chromosome.fasta trycycler_rapid_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_medaka_ligation_chromosome.fasta trycycler_medaka_rapid_chromosome.fasta > trycycler_medaka_chromosome.errors
../../../../scripts/error_positions.py trycycler_medaka_pilon_ligation_chromosome.fasta trycycler_medaka_pilon_rapid_chromosome.fasta > trycycler_medaka_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_ligation_chromosome.fasta unicycler_hybrid_rapid_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/Enterobacter_kobei_MSB1_1B
../../../../scripts/error_positions.py flye_ligation_chromosome.fasta flye_rapid_chromosome.fasta 871102 4662345 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_ligation_chromosome_fixed_circularisation.fasta flye_rapid_chromosome_fixed_circularisation.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_medaka_ligation_chromosome.fasta flye_medaka_rapid_chromosome.fasta 871355 4664175 > flye_medaka_chromosome.errors
../../../../scripts/error_positions.py flye_medaka_ligation_chromosome_fixed_circularisation.fasta flye_medaka_rapid_chromosome_fixed_circularisation.fasta > flye_medaka_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_medaka_pilon_ligation_chromosome.fasta flye_medaka_pilon_rapid_chromosome.fasta 871329 4664141 > flye_medaka_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_medaka_pilon_ligation_chromosome_fixed_circularisation.fasta flye_medaka_pilon_rapid_chromosome_fixed_circularisation.fasta > flye_medaka_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_ligation_chromosome.fasta trycycler_rapid_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_medaka_ligation_chromosome.fasta trycycler_medaka_rapid_chromosome.fasta > trycycler_medaka_chromosome.errors
../../../../scripts/error_positions.py trycycler_medaka_pilon_ligation_chromosome.fasta trycycler_medaka_pilon_rapid_chromosome.fasta > trycycler_medaka_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_ligation_chromosome.fasta unicycler_hybrid_rapid_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/Haemophilus_M1C132_1
../../../../scripts/error_positions.py flye_ligation_chromosome.fasta flye_rapid_chromosome.fasta 515457 515413 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_ligation_chromosome_fixed_circularisation.fasta flye_rapid_chromosome_fixed_circularisation.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_medaka_ligation_chromosome.fasta flye_medaka_rapid_chromosome.fasta 515593 515595 > flye_medaka_chromosome.errors
../../../../scripts/error_positions.py flye_medaka_ligation_chromosome_fixed_circularisation.fasta flye_medaka_rapid_chromosome_fixed_circularisation.fasta > flye_medaka_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_medaka_pilon_ligation_chromosome.fasta flye_medaka_pilon_rapid_chromosome.fasta 515565 515563 > flye_medaka_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_medaka_pilon_ligation_chromosome_fixed_circularisation.fasta flye_medaka_pilon_rapid_chromosome_fixed_circularisation.fasta > flye_medaka_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_ligation_chromosome.fasta trycycler_rapid_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_medaka_ligation_chromosome.fasta trycycler_medaka_rapid_chromosome.fasta > trycycler_medaka_chromosome.errors
../../../../scripts/error_positions.py trycycler_medaka_pilon_ligation_chromosome.fasta trycycler_medaka_pilon_rapid_chromosome.fasta > trycycler_medaka_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_ligation_chromosome.fasta unicycler_hybrid_rapid_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/Klebsiella_oxytoca_MSB1_2C
../../../../scripts/error_positions.py flye_ligation_chromosome.fasta flye_rapid_chromosome.fasta 5743486 1325984 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_ligation_chromosome_fixed_circularisation.fasta flye_rapid_chromosome_fixed_circularisation.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_medaka_ligation_chromosome.fasta flye_medaka_rapid_chromosome.fasta 5745312 1326505 > flye_medaka_chromosome.errors
../../../../scripts/error_positions.py flye_medaka_ligation_chromosome_fixed_circularisation.fasta flye_medaka_rapid_chromosome_fixed_circularisation.fasta > flye_medaka_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_medaka_pilon_ligation_chromosome.fasta flye_medaka_pilon_rapid_chromosome.fasta 5745404 1326543 > flye_medaka_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_medaka_pilon_ligation_chromosome_fixed_circularisation.fasta flye_medaka_pilon_rapid_chromosome_fixed_circularisation.fasta > flye_medaka_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_ligation_chromosome.fasta trycycler_rapid_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_medaka_ligation_chromosome.fasta trycycler_medaka_rapid_chromosome.fasta > trycycler_medaka_chromosome.errors
../../../../scripts/error_positions.py trycycler_medaka_pilon_ligation_chromosome.fasta trycycler_medaka_pilon_rapid_chromosome.fasta > trycycler_medaka_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_ligation_chromosome.fasta unicycler_hybrid_rapid_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/real_read_tests/assemblies/chromosome_only/Klebsiella_variicola_INF345
../../../../scripts/error_positions.py flye_ligation_chromosome.fasta flye_rapid_chromosome.fasta 451098 451022 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_ligation_chromosome_fixed_circularisation.fasta flye_rapid_chromosome_fixed_circularisation.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_medaka_ligation_chromosome.fasta flye_medaka_rapid_chromosome.fasta 451272 451236 > flye_medaka_chromosome.errors
../../../../scripts/error_positions.py flye_medaka_ligation_chromosome_fixed_circularisation.fasta flye_medaka_rapid_chromosome_fixed_circularisation.fasta > flye_medaka_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_medaka_pilon_ligation_chromosome.fasta flye_medaka_pilon_rapid_chromosome.fasta 451256 451247 > flye_medaka_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_medaka_pilon_ligation_chromosome_fixed_circularisation.fasta flye_medaka_pilon_rapid_chromosome_fixed_circularisation.fasta > flye_medaka_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_ligation_chromosome.fasta trycycler_rapid_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_medaka_ligation_chromosome.fasta trycycler_medaka_rapid_chromosome.fasta > trycycler_medaka_chromosome.errors
../../../../scripts/error_positions.py trycycler_medaka_pilon_ligation_chromosome.fasta trycycler_medaka_pilon_rapid_chromosome.fasta > trycycler_medaka_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_ligation_chromosome.fasta unicycler_hybrid_rapid_chromosome.fasta > unicycler_hybrid_chromosome.errors
```
