## Read analysis


I analysed the reads on Nectar, both before filtering:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    fast_count ~/trycycler_simulated_read_tests/"$s"/raw_reads/nanopore.fastq.gz
done
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    fast_count ~/trycycler_simulated_read_tests/"$s"/raw_reads/illumina_1.fastq.gz
    fast_count ~/trycycler_simulated_read_tests/"$s"/raw_reads/illumina_2.fastq.gz
done
```
and after filtering:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    fast_count ~/trycycler_simulated_read_tests/"$s"/nanopore.fastq.gz
done
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    fast_count ~/trycycler_simulated_read_tests/"$s"/illumina_1.fastq.gz
    fast_count ~/trycycler_simulated_read_tests/"$s"/illumina_2.fastq.gz
    fast_count ~/trycycler_simulated_read_tests/"$s"/illumina_u.fastq.gz
done
```




## Resource usage

For Trycycler RAM usage, I took the largest value of all the steps:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /home/ubuntu/trycycler_simulated_read_tests/"$s"/04_trycycler
    printf $s"\t"
    for f in $(find . -name "*.time" | grep -v "medaka" | sort); do grep "Maximum resident set size" $f | grep -oP "\d+"; done | sort -n | tail -n1
done
```

For Trycycler time, I took the largest sum of all the steps (used Excel to total them up):
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /home/ubuntu/trycycler_simulated_read_tests/"$s"/04_trycycler
    echo "$s"
    for f in $(find . -name "*.time" | grep -v "medaka" | sort); do grep "Elapsed (wall clock) time" $f | grep -oP "[\d:\.]{2,}"; done
    printf "\n\n"
done
```

Subsequent analyses were done on my laptop.






Transferring the assemblies to my computer:
```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/entire_assemblies
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    mkdir $s
done

for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/entire_assemblies/"$s"

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/01_unicycler_illumina_only/002_overlaps_removed.gfa unicycler_illumina_only_no_overlap.gfa

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/02_flye/assembly.fasta flye.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/02_flye/assembly_graph.gfa flye.gfa
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/02_flye/polished.fasta flye_pilon.fasta

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/03_unicycler_hybrid/assembly.fasta unicycler_hybrid.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/03_unicycler_hybrid/assembly.gfa unicycler_hybrid.gfa
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/03_unicycler_hybrid/004_bridges_applied.gfa unicycler_hybrid_bridged.gfa

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/04_trycycler/assembly.fasta trycycler.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/04_trycycler/polished.fasta trycycler_pilon.fasta

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/05_raven/raven.fasta raven.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/05_raven/raven.gfa raven.gfa
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/05_raven/polished.fasta raven_pilon.fasta

    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/06_miniasm/miniasm.fasta miniasm.fasta
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/06_miniasm/miniasm.gfa miniasm.gfa
    scp assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/06_miniasm/polished.fasta miniasm_pilon.fasta
done
```

I then manually inspected the graphs in Bandage to determine whether or not the chromosome was complete.

I then manually made chromosome-only FASTA files for each of the completed Unicycler hybrid, Flye and Trycycler assemblies. I did this by:
* Copying the FASTA and adding `_chromosome` to the end of its name.
* Removing in-sequence line breaks: `seqtk seq $f > temp; mv temp $f`
* Deleting non-chromosome contigs
* Manually adjusting the start point of each assembly to match the Unicycler hybrid start point (flipping to reverse strand if necessary with `seqtk seq -r $f > temp; mv temp $f`)




For Flye assemblies, circularisation could be an issue, so I created versions of each assembly where I manually fixed the circularisation (using the reference genome as a point of comparison). This had to be done before changing the starting position. Here were the adjustments I made:
* Campylobacter_jejuni:       extra 3 bp (ATG)
* Escherichia_coli:           good
* Klebsiella_pneumoniae:      missing 7 bp (TTTGATG)
* Listeria_monocytogenes:     extra 1 bp (A)
* Mycobacterium_tuberculosis: good
* Neisseria_meningitidis:     missing 9 bp (AGACGGCAT)
* Pseudomonas_aeruginosa:     missing 1 bp (C)
* Salmonella_enterica:        missing 6 bp (GGAGTT)
* Staphylococcus_aureus:      missing 9 bp (ATGAACATT)
* Streptococcus_pneumoniae:   extra 2 bp (AA)






Now that I have matched up assemblies (same replicon, same strand, same start), I can do an Edlib global alignment to do pairwise comparisons.

Unicycler vs reference:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py unicycler_hybrid_chromosome.fasta reference_chromosome.fasta
done
```

Flye vs reference:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py flye_chromosome.fasta reference_chromosome.fasta
done
```

Flye+Pilon vs reference:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py flye_pilon_chromosome.fasta reference_chromosome.fasta
done
```

Flye (fixed circ) vs reference:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py flye_chromosome_fixed_circularisation.fasta reference_chromosome.fasta
done
```

Flye+Pilon (fixed circ) vs reference:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py flye_pilon_chromosome_fixed_circularisation.fasta reference_chromosome.fasta
done
```

Trycycler vs reference:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py trycycler_chromosome.fasta reference_chromosome.fasta
done
```

Trycycler+Pilon vs reference:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py trycycler_pilon_chromosome.fasta reference_chromosome.fasta
done
```

Miniasm vs reference:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py miniasm_chromosome.fasta reference_chromosome.fasta
done
```

Miniasm+Pilon vs reference:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py miniasm_pilon_chromosome.fasta reference_chromosome.fasta
done
```

Raven vs reference:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py raven_chromosome.fasta reference_chromosome.fasta
done
```

Raven+Pilon vs reference:
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/"$s"
    printf $s"\t"
    ../../../../scripts/pairwise_align.py raven_pilon_chromosome.fasta reference_chromosome.fasta
done
```





## Repeat analysis

I was also curious what the longest repeat was in each chromosome. I did this by downloading the Unicycler hybrid Nanopore bridged graph (all of which had a complete chromosome) and looked for the longest bridge segment that was in the chromosomal loop.
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/entire_assemblies/"$s"
    printf $s"\t"
    ../../../../scripts/find_longest_chromosomal_repeat.py unicycler_hybrid_bridged.gfa
done
```







## Error positions

First I ran a script to get the Flye starting position, as that's needed for the next command:
```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Campylobacter_jejuni
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Campylobacter_jejuni/flye.fasta flye_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Campylobacter_jejuni/flye_pilon.fasta flye_pilon_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Escherichia_coli
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Escherichia_coli/flye.fasta flye_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Escherichia_coli/flye_pilon.fasta flye_pilon_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Klebsiella_pneumoniae
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_pneumoniae/flye.fasta flye_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Klebsiella_pneumoniae/flye_pilon.fasta flye_pilon_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Listeria_monocytogenes
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Listeria_monocytogenes/flye.fasta flye_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Listeria_monocytogenes/flye_pilon.fasta flye_pilon_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Mycobacterium_tuberculosis
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Mycobacterium_tuberculosis/flye.fasta flye_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Mycobacterium_tuberculosis/flye_pilon.fasta flye_pilon_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Neisseria_meningitidis
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Neisseria_meningitidis/flye.fasta flye_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Neisseria_meningitidis/flye_pilon.fasta flye_pilon_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Pseudomonas_aeruginosa
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Pseudomonas_aeruginosa/flye.fasta flye_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Pseudomonas_aeruginosa/flye_pilon.fasta flye_pilon_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Salmonella_enterica
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Salmonella_enterica/flye.fasta flye_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Salmonella_enterica/flye_pilon.fasta flye_pilon_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Staphylococcus_aureus
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Staphylococcus_aureus/flye.fasta flye_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Staphylococcus_aureus/flye_pilon.fasta flye_pilon_chromosome.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Streptococcus_pneumoniae
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Streptococcus_pneumoniae/flye.fasta flye_chromosome.fasta
../../../../scripts/get_flye_starting_pos.py ../../entire_assemblies/Streptococcus_pneumoniae/flye_pilon.fasta flye_pilon_chromosome.fasta
```

This will generate the data that my R script can then visualise:
```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Campylobacter_jejuni
../../../../scripts/error_positions.py flye_chromosome.fasta reference_chromosome.fasta 399732 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_pilon_chromosome.fasta reference_chromosome.fasta 399748 > flye_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_pilon_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_chromosome.fasta reference_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_pilon_chromosome.fasta reference_chromosome.fasta > trycycler_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_chromosome.fasta reference_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Escherichia_coli
../../../../scripts/error_positions.py flye_chromosome.fasta reference_chromosome.fasta 4035146 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_pilon_chromosome.fasta reference_chromosome.fasta 4035217 > flye_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_pilon_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_chromosome.fasta reference_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_pilon_chromosome.fasta reference_chromosome.fasta > trycycler_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_chromosome.fasta reference_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Klebsiella_pneumoniae
../../../../scripts/error_positions.py flye_chromosome.fasta reference_chromosome.fasta 262580 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_pilon_chromosome.fasta reference_chromosome.fasta 262615 > flye_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_pilon_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_chromosome.fasta reference_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_pilon_chromosome.fasta reference_chromosome.fasta > trycycler_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_chromosome.fasta reference_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Listeria_monocytogenes
../../../../scripts/error_positions.py flye_chromosome.fasta reference_chromosome.fasta 2443199 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_pilon_chromosome.fasta reference_chromosome.fasta 2443199 > flye_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_pilon_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_chromosome.fasta reference_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_pilon_chromosome.fasta reference_chromosome.fasta > trycycler_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_chromosome.fasta reference_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Mycobacterium_tuberculosis
../../../../scripts/error_positions.py flye_chromosome.fasta reference_chromosome.fasta 1055732 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_pilon_chromosome.fasta reference_chromosome.fasta 1055743 > flye_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_pilon_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_chromosome.fasta reference_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_pilon_chromosome.fasta reference_chromosome.fasta > trycycler_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_chromosome.fasta reference_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Neisseria_meningitidis
../../../../scripts/error_positions.py flye_chromosome.fasta reference_chromosome.fasta 1731009 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_pilon_chromosome.fasta reference_chromosome.fasta 1731009 > flye_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_pilon_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_chromosome.fasta reference_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_pilon_chromosome.fasta reference_chromosome.fasta > trycycler_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_chromosome.fasta reference_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Pseudomonas_aeruginosa
../../../../scripts/error_positions.py flye_chromosome.fasta reference_chromosome.fasta 4788514 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_pilon_chromosome.fasta reference_chromosome.fasta 4788526 > flye_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_pilon_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_chromosome.fasta reference_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_pilon_chromosome.fasta reference_chromosome.fasta > trycycler_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_chromosome.fasta reference_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Salmonella_enterica
../../../../scripts/error_positions.py flye_chromosome.fasta reference_chromosome.fasta 2800759 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_pilon_chromosome.fasta reference_chromosome.fasta 2801957 > flye_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_pilon_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_chromosome.fasta reference_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_pilon_chromosome.fasta reference_chromosome.fasta > trycycler_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_chromosome.fasta reference_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Staphylococcus_aureus
../../../../scripts/error_positions.py flye_chromosome.fasta reference_chromosome.fasta 492985 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_pilon_chromosome.fasta reference_chromosome.fasta 492992 > flye_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_pilon_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_chromosome.fasta reference_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_pilon_chromosome.fasta reference_chromosome.fasta > trycycler_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_chromosome.fasta reference_chromosome.fasta > unicycler_hybrid_chromosome.errors

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/assemblies/chromosome_only/Streptococcus_pneumoniae
../../../../scripts/error_positions.py flye_chromosome.fasta reference_chromosome.fasta 1795679 > flye_chromosome.errors
../../../../scripts/error_positions.py flye_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py flye_pilon_chromosome.fasta reference_chromosome.fasta 1795751 > flye_pilon_chromosome.errors
../../../../scripts/error_positions.py flye_pilon_chromosome_fixed_circularisation.fasta reference_chromosome.fasta > flye_pilon_chromosome_fixed_circularisation.errors
../../../../scripts/error_positions.py trycycler_chromosome.fasta reference_chromosome.fasta > trycycler_chromosome.errors
../../../../scripts/error_positions.py trycycler_pilon_chromosome.fasta reference_chromosome.fasta > trycycler_pilon_chromosome.errors
../../../../scripts/error_positions.py unicycler_hybrid_chromosome.fasta reference_chromosome.fasta > unicycler_hybrid_chromosome.errors
```
