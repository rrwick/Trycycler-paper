I did these assemblies on my biggest Nectar instance (assembler_benchmarking_1).

Make the directories (ran on Nectar):
```bash
mkdir ~/trycycler_simulated_read_tests
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    mkdir ~/trycycler_simulated_read_tests/"$s"
    mkdir ~/trycycler_simulated_read_tests/"$s"/raw_reads
done
```

Transfer the reads to Nectar (ran on my laptop):
```bash
for s in Campylobacter_jejuni Escherichia_coli Klebsiella_pneumoniae Listeria_monocytogenes Mycobacterium_tuberculosis Neisseria_meningitidis Pseudomonas_aeruginosa Salmonella_enterica Staphylococcus_aureus Streptococcus_pneumoniae; do
    cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/"$s"
    scp *.fastq.gz assembler_benchmarking_1:/home/ubuntu/trycycler_simulated_read_tests/"$s"/raw_reads
done
```


From this point onward, the commands are per-sample:
```bash
s=Campylobacter_jejuni
cd ~/trycycler_simulated_read_tests/"$s"

# These genome sizes are from the reference sequences.
declare -A genome_sizes
genome_sizes[Campylobacter_jejuni]=1641481
genome_sizes[Escherichia_coli]=4641652
genome_sizes[Klebsiella_pneumoniae]=5682322
genome_sizes[Listeria_monocytogenes]=2944528
genome_sizes[Mycobacterium_tuberculosis]=4411532
genome_sizes[Neisseria_meningitidis]=2272360
genome_sizes[Pseudomonas_aeruginosa]=6264404
genome_sizes[Salmonella_enterica]=4951383
genome_sizes[Staphylococcus_aureus]=2821361
genome_sizes[Streptococcus_pneumoniae]=2038615
genome_size=$genome_sizes[$s]
```


Prep the reads:
```bash
fastp --in1 raw_reads/illumina_1.fastq.gz --in2 raw_reads/illumina_2.fastq.gz --out1 illumina_1.fastq.gz --out2 illumina_2.fastq.gz --unpaired1 illumina_u.fastq.gz --unpaired2 illumina_u.fastq.gz
mv fastp.* raw_reads
filtlong --min_length 1000 --keep_percent 95 raw_reads/nanopore.fastq.gz | gzip > nanopore.fastq.gz
```


Do an Illumina-only Unicycler assembly:
```bash
/usr/bin/time -v -o unicycler.time unicycler -1 illumina_1.fastq.gz -2 illumina_2.fastq.gz -s illumina_u.fastq.gz -o 01_unicycler_illumina_only --threads 32 --no_correct --no_pilon --no_rotate
mv unicycler.time 01_unicycler_illumina_only
```
If necessary, fix up the multiplicity in the 001_best_spades_graph.gfa file.


Do a Unicycler hybrid assembly:
```bash
mkdir 03_unicycler_hybrid
cp 01_unicycler_illumina_only/001_best_spades_graph.gfa 03_unicycler_hybrid
/usr/bin/time -v -o unicycler.time unicycler -1 illumina_1.fastq.gz -2 illumina_2.fastq.gz -s illumina_u.fastq.gz -l nanopore.fastq.gz -o 03_unicycler_hybrid --threads 32 --no_correct
mv unicycler.time 03_unicycler_hybrid
```


Do Flye v2.7.1 assembly:
```bash
/usr/bin/time -v -o flye.time flye --nano-raw nanopore.fastq.gz --genome-size "$genome_size" --threads 32 --plasmids --out-dir 02_flye &> flye.out
mv flye.time flye.out 02_flye
```


Polish Flye assemblies with Pilon:
```bash
cd 02_flye

# Take these values from the Unicycler hybrid assembly log:
insert_min=197
insert_max=670

before=assembly
after=round_1
bowtie2-build "$before".fasta "$before".fasta
bowtie2 -1 ../illumina_1.fastq.gz -2 ../illumina_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
samtools index illumina_alignments.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=round_1
after=round_2
bowtie2-build "$before".fasta "$before".fasta
bowtie2 -1 ../illumina_1.fastq.gz -2 ../illumina_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
samtools index illumina_alignments.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls round_*.fasta | tail -n1)
mv $last_round polished.fasta  # change to last round
cd ..
```


Do a Trycycler assembly:
```bash
mkdir 04_trycycler
cp nanopore.fastq.gz 04_trycycler/reads.fastq.gz
mkdir 04_trycycler/assemblies
target_depth=50
mean_length=$(seqtk comp nanopore.fastq.gz | awk '{count++; bases += $2} END{print bases/count}')
read_count=$(echo $target_depth"*"$genome_size"/"$mean_length | bc)

date > 04_trycycler/sampling.time
for i in {00..11}; do
    seqtk sample -s "$i" nanopore.fastq.gz "$read_count" | paste - - - - | shuf | tr '\t' '\n' > sample_"$i".fastq
done
date >> 04_trycycler/sampling.time
for i in {00..02}; do
    /usr/bin/time -v -o 04_trycycler/assemblies/assembly_"$i".time flye --nano-raw sample_"$i".fastq --genome-size "$genome_size" --threads 32 --plasmids --out-dir assembly_"$i"; cp assembly_"$i"/assembly.fasta 04_trycycler/assemblies/assembly_"$i".fasta; rm -r assembly_"$i"
done
for i in {03..05}; do
    /usr/bin/time -v -o 04_trycycler/assemblies/assembly_"$i".time miniasm_and_minipolish.sh sample_"$i".fastq 32 > assembly_"$i".gfa; any2fasta assembly_"$i".gfa > 04_trycycler/assemblies/assembly_"$i".fasta; rm assembly_"$i".gfa
done
for i in {06..08}; do
    /usr/bin/time -v -o 04_trycycler/assemblies/assembly_"$i".time raven --threads 32 sample_"$i".fastq > 04_trycycler/assemblies/assembly_"$i".fasta; rm raven.cereal
done
for i in {09..11}; do
    /usr/bin/time -v -o 04_trycycler/assemblies/assembly_"$i".time wtdbg2.pl -o assembly_"$i" -g "$genome_size" -t 32 -x ont sample_"$i".fastq; mv assembly_"$i".cns.fa 04_trycycler/assemblies/assembly_"$i".fasta; rm assembly_"$i".*
done
rm sample_*.fastq

/usr/bin/time -v -o 04_trycycler/cluster.time ~/Trycycler/trycycler-runner.py cluster --reads 04_trycycler/reads.fastq.gz --assemblies 04_trycycler/assemblies/*.fasta --out_dir 04_trycycler/trycycler --threads 32 &> 04_trycycler/cluster.out
# Choose clusters

/usr/bin/time -v -o 04_trycycler/trycycler/cluster_001/reconcile.time ~/Trycycler/trycycler-runner.py reconcile --reads 04_trycycler/reads.fastq.gz --cluster_dir 04_trycycler/trycycler/cluster_001 --threads 32
/usr/bin/time -v -o 04_trycycler/trycycler/cluster_002/reconcile.time ~/Trycycler/trycycler-runner.py reconcile --reads 04_trycycler/reads.fastq.gz --cluster_dir 04_trycycler/trycycler/cluster_002 --threads 32
/usr/bin/time -v -o 04_trycycler/trycycler/cluster_003/reconcile.time ~/Trycycler/trycycler-runner.py reconcile --reads 04_trycycler/reads.fastq.gz --cluster_dir 04_trycycler/trycycler/cluster_003 --threads 32
# Fix reconciliation issues

for c in 04_trycycler/trycycler/cluster_*; do
    /usr/bin/time -v -o "$c"/msa.time ~/Trycycler/trycycler-runner.py msa --cluster_dir "$c" --threads 32
done

/usr/bin/time -v -o 04_trycycler/partition.time ~/Trycycler/trycycler-runner.py partition --reads 04_trycycler/reads.fastq.gz --cluster_dirs 04_trycycler/trycycler/cluster_* --threads 32

for c in 04_trycycler/trycycler/cluster_*; do
    /usr/bin/time -v -o "$c"/consensus.time ~/Trycycler/trycycler-runner.py consensus --cluster_dir "$c" --threads 32 --verbose
done

cat 04_trycycler/trycycler/cluster_*/7_final_consensus.fasta > 04_trycycler/assembly.fasta
```


Polish the Trycycler assembly with Pilon:
```bash
cd 04_trycycler

# Take these values from the Unicycler hybrid assembly log:
insert_min=197
insert_max=670

before=assembly
after=round_1
bowtie2-build "$before".fasta "$before".fasta
bowtie2 -1 ../illumina_1.fastq.gz -2 ../illumina_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
samtools index illumina_alignments.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=round_1
after=round_2
bowtie2-build "$before".fasta "$before".fasta
bowtie2 -1 ../illumina_1.fastq.gz -2 ../illumina_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
samtools index illumina_alignments.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls round_*.fasta | tail -n1)
mv $last_round polished.fasta  # change to last round
cd ..
```









## Other assemblers

I later decided to throw a couple other assemblers into the mix, so here are their assembly commands:


### Raven v1.2.2

```bash
/usr/bin/time -v -o raven.time raven --disable-checkpoints --graphical-fragment-assembly raven.gfa -t 32 nanopore.fastq.gz 1> raven.fasta 2> raven.out
mkdir 05_raven; mv raven.time raven.fasta raven.gfa raven.out 05_raven
```

Polish the Raven assembly with Pilon:
```bash
cd 05_raven

# Take these values from the Unicycler hybrid assembly log:
insert_min=371
insert_max=469

before=raven
after=round_1
bowtie2-build "$before".fasta "$before".fasta
bowtie2 -1 ../illumina_1.fastq.gz -2 ../illumina_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
samtools index illumina_alignments.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=round_1
after=round_2
bowtie2-build "$before".fasta "$before".fasta
bowtie2 -1 ../illumina_1.fastq.gz -2 ../illumina_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
samtools index illumina_alignments.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls round_*.fasta | tail -n1)
mv $last_round polished.fasta  # change to last round
cd ..
```


### Miniasm/Minipolish v0.3/v0.1.3

```bash
/usr/bin/time -v -o miniasm.time miniasm_and_minipolish.sh nanopore.fastq.gz 32 1> miniasm.gfa 2> miniasm.out
any2fasta miniasm.gfa > miniasm.fasta
mkdir 06_miniasm; mv miniasm.time miniasm.fasta miniasm.gfa miniasm.out 06_miniasm
```

Polish the Miniasm/Minipolish assembly with Pilon:
```bash
cd 06_miniasm

# Take these values from the Unicycler hybrid assembly log:
insert_min=371
insert_max=469

before=miniasm
after=round_1
bowtie2-build "$before".fasta "$before".fasta
bowtie2 -1 ../illumina_1.fastq.gz -2 ../illumina_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
samtools index illumina_alignments.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=round_1
after=round_2
bowtie2-build "$before".fasta "$before".fasta
bowtie2 -1 ../illumina_1.fastq.gz -2 ../illumina_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
samtools index illumina_alignments.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls round_*.fasta | tail -n1)
mv $last_round polished.fasta  # change to last round
cd ..
```
