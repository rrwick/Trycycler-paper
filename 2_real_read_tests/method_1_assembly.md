## Prep the reads

I did these assemblies on my biggest Nectar instance (assembler_benchmarking_1).

Make the directories:
```bash
mkdir ~/trycycler_real_read_tests
for s in Acinetobacter_baumannii_J9 Citrobacter_koseri_MINF_9D Enterobacter_kobei_MSB1_1B Haemophilus_M1C132_1 Klebsiella_oxytoca_MSB1_2C Klebsiella_variicola_INF345 Serratia_marcescens_17-147-1671; do
    mkdir ~/trycycler_real_read_tests/"$s"
    mkdir ~/trycycler_real_read_tests/"$s"/raw_reads
done
```

Transfer the reads from MASSIVE:
```bash
cd ~/trycycler_real_read_tests/Acinetobacter_baumannii_J9/raw_reads
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-10_small_plasmid_ligation_02/fastq/barcode01.fastq.gz ligation.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-09_small_plasmid_rapid_02/fastq/barcode01.fastq.gz rapid.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-J9_1.fastq.gz illumina_1.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-J9_2.fastq.gz illumina_2.fastq.gz

cd ~/trycycler_real_read_tests/Citrobacter_koseri_MINF_9D/raw_reads
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-10_small_plasmid_ligation_02/fastq/barcode02.fastq.gz ligation.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-09_small_plasmid_rapid_02/fastq/barcode02.fastq.gz rapid.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-MINF-9D_1.fastq.gz illumina_1.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-MINF-9D_2.fastq.gz illumina_2.fastq.gz

cd ~/trycycler_real_read_tests/Enterobacter_kobei_MSB1_1B/raw_reads
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-10_small_plasmid_ligation_02/fastq/barcode03.fastq.gz ligation.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-09_small_plasmid_rapid_02/fastq/barcode03.fastq.gz rapid.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-MSB1-1B_1.fastq.gz illumina_1.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-MSB1-1B_2.fastq.gz illumina_2.fastq.gz

cd ~/trycycler_real_read_tests/Haemophilus_M1C132_1/raw_reads
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-10_small_plasmid_ligation_02/fastq/barcode04.fastq.gz ligation.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-09_small_plasmid_rapid_02/fastq/barcode04.fastq.gz rapid.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-MIC132-1_1.fastq.gz illumina_1.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-MIC132-1_2.fastq.gz illumina_2.fastq.gz

cd ~/trycycler_real_read_tests/Klebsiella_oxytoca_MSB1_2C/raw_reads
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-10_small_plasmid_ligation_02/fastq/barcode05.fastq.gz ligation.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-09_small_plasmid_rapid_02/fastq/barcode05.fastq.gz rapid.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-MSB1-2C_1.fastq.gz illumina_1.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-MSB1-2C_2.fastq.gz illumina_2.fastq.gz

cd ~/trycycler_real_read_tests/Klebsiella_variicola_INF345/raw_reads
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-10_small_plasmid_ligation_02/fastq/barcode07.fastq.gz ligation.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-09_small_plasmid_rapid_02/fastq/barcode07.fastq.gz rapid.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-INF345_1.fastq.gz illumina_1.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-INF345_2.fastq.gz illumina_2.fastq.gz

cd ~/trycycler_real_read_tests/Serratia_marcescens_17-147-1671/raw_reads
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-10_small_plasmid_ligation_02/fastq/barcode08.fastq.gz ligation.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/minion/2020-06-09_small_plasmid_rapid_02/fastq/barcode08.fastq.gz rapid.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-17-147-1671_1.fastq.gz illumina_1.fastq.gz
scp rwic0002@m3-dtn1.massive.org.au:/mnt/vault2_holt/vault/instruments/other_illumina/NextSeq_ACBD_200217_Klebs_colistin_small_plasmids/FASTQ/plasmids-17-147-1671_2.fastq.gz illumina_2.fastq.gz
```


From this point onward, the commands are per-sample:
```bash
s=Acinetobacter_baumannii_J9
cd ~/trycycler_real_read_tests/"$s"

# These genome size estimates are based on my previous analysis (for the small plasmid paper):
declare -A genome_sizes
genome_sizes[Acinetobacter_baumannii_J9]=3950000
genome_sizes[Citrobacter_koseri_MINF_9D]=4833000
genome_sizes[Enterobacter_kobei_MSB1_1B]=5094000
genome_sizes[Haemophilus_M1C132_1]=2112000
genome_sizes[Klebsiella_oxytoca_MSB1_2C]=5986000
genome_sizes[Klebsiella_variicola_INF345]=5953000
genome_sizes[Serratia_marcescens_17-147-1671]=5884000
genome_size=$genome_sizes[$s]
```


Filter the Nanopore reads:
```bash
filtlong --min_length 1000 --keep_percent 95 raw_reads/ligation.fastq.gz | gzip > ligation.fastq.gz
filtlong --min_length 1000 --keep_percent 95 raw_reads/rapid.fastq.gz | gzip > rapid.fastq.gz
```


There were a couple rapid reads of highly repetitive junk that made Flye hang for days, so I had to track them down and remove them:
```bash
zcat rapid.fastq.gz | paste - - - - | grep -v -P "5ef0247e-3106-4203-852f-79d78afdbdde|06ae7b30-16ea-4238-8523-ad85a1e47b0f|2621b126-195b-4a75-a7d1-3d286113aec0" | tr "\t" "\n" | gzip > rapid_filter.fastq.gz
mv rapid_filter.fastq.gz rapid.fastq.gz
```


Split the Illumina reads randomly into two sets. Here are the a and b counts for each isolate:
* Acinetobacter_baumannii_J9: 1557613,1557613
* Citrobacter_koseri_MINF_9D: 1142150,1142150
* Enterobacter_kobei_MSB1_1B: 1023434,1023433
* Haemophilus_M1C132_1: 991749,991749
* Klebsiella_oxytoca_MSB1_2C: 1818842,1818841
* Klebsiella_variicola_INF345: 833289,833289
* Serratia_marcescens_17-147-1671: 792947,792947
```bash
cd raw_reads
a_count=1557613
b_count=1557613
seqtk mergepe illumina_1.fastq.gz illumina_2.fastq.gz | paste - - - - - - - - | shuf > temp_merged.fastq
head -n "$a_count" temp_merged.fastq | cut -f 1-4 | tr "\t" "\n" | gzip > illumina_a_1.fastq.gz
head -n "$a_count" temp_merged.fastq | cut -f 5-8 | tr "\t" "\n" | gzip > illumina_a_2.fastq.gz
tail -n "$b_count" temp_merged.fastq | cut -f 1-4 | tr "\t" "\n" | gzip > illumina_b_1.fastq.gz
tail -n "$b_count" temp_merged.fastq | cut -f 5-8 | tr "\t" "\n" | gzip > illumina_b_2.fastq.gz
rm temp_merged.fastq
cd ..
```


And then QC the Illumina reads:
```bash
fastp --in1 raw_reads/illumina_a_1.fastq.gz --in2 raw_reads/illumina_a_2.fastq.gz --out1 illumina_a_1.fastq.gz --out2 illumina_a_2.fastq.gz --unpaired1 illumina_a_u.fastq.gz --unpaired2 illumina_a_u.fastq.gz
mv fastp.html raw_reads/fastp_a.html
mv fastp.json raw_reads/fastp_a.json
fastp --in1 raw_reads/illumina_b_1.fastq.gz --in2 raw_reads/illumina_b_2.fastq.gz --out1 illumina_b_1.fastq.gz --out2 illumina_b_2.fastq.gz --unpaired1 illumina_b_u.fastq.gz --unpaired2 illumina_b_u.fastq.gz
mv fastp.html raw_reads/fastp_b.html
mv fastp.json raw_reads/fastp_b.json
```




## Unicycler assemblies

Do an Illumina-ligation Unicycler assembly:
```bash
/usr/bin/time -v -o unicycler.time unicycler -1 illumina_a_1.fastq.gz -2 illumina_a_2.fastq.gz -s illumina_a_u.fastq.gz -o 03_unicycler_hybrid_ligation --threads 32 --no_correct --no_pilon --no_rotate
mv unicycler.time 03_unicycler_hybrid_ligation/unicycler.time1
rm 03_unicycler_hybrid_ligation/002_overlaps_removed.gfa 03_unicycler_hybrid_ligation/003_bridges_applied.gfa 03_unicycler_hybrid_ligation/004_final_clean.gfa 03_unicycler_hybrid_ligation/assembly.fasta 03_unicycler_hybrid_ligation/assembly.gfa
# If necessary, fix up the multiplicity in the 001_best_spades_graph.gfa file.
/usr/bin/time -v -o unicycler.time unicycler -1 illumina_a_1.fastq.gz -2 illumina_a_2.fastq.gz -s illumina_a_u.fastq.gz -l ligation.fastq.gz -o 03_unicycler_hybrid_ligation --threads 32 --no_correct
mv unicycler.time 03_unicycler_hybrid_ligation/unicycler.time2
```

For the `Klebsiella_variicola_INF345` assembly, Unicycler made a strange k-mer choice, going with a dead-end and misassembly filled k127 assembly. So I re-ran that one with `--kmers 27,47,63,77,89,99,107,115,121` to avoid the k127 graph, and it chose k77 instead.

For Pilon polishing steps later, I've noted the Illumina read insert range (as determined by Unicycler) for each of the isolates:
* Acinetobacter_baumannii_J9: 157-712
* Citrobacter_koseri_MINF_9D: 159-721
* Enterobacter_kobei_MSB1_1B: 160-719
* Haemophilus_M1C132_1: 160-758
* Klebsiella_oxytoca_MSB1_2C: 157-716
* Klebsiella_variicola_INF345: 156-704
* Serratia_marcescens_17-147-1671: 158-734


Do an Illumina-rapid Unicycler assembly:
```bash
/usr/bin/time -v -o unicycler.time unicycler -1 illumina_b_1.fastq.gz -2 illumina_b_2.fastq.gz -s illumina_b_u.fastq.gz -o 04_unicycler_hybrid_rapid --threads 32 --no_correct --no_pilon --no_rotate
mv unicycler.time 04_unicycler_hybrid_rapid/unicycler.time1
rm 04_unicycler_hybrid_rapid/002_overlaps_removed.gfa 04_unicycler_hybrid_rapid/003_bridges_applied.gfa 04_unicycler_hybrid_rapid/004_final_clean.gfa 04_unicycler_hybrid_rapid/assembly.fasta 04_unicycler_hybrid_rapid/assembly.gfa
# If necessary, fix up the multiplicity in the 001_best_spades_graph.gfa file.
/usr/bin/time -v -o unicycler.time unicycler -1 illumina_b_1.fastq.gz -2 illumina_b_2.fastq.gz -s illumina_b_u.fastq.gz -l rapid.fastq.gz -o 04_unicycler_hybrid_rapid --threads 32 --no_correct
mv unicycler.time 04_unicycler_hybrid_rapid/unicycler.time2
```

For Pilon polishing steps later, I've noted the Illumina read insert range (as determined by Unicycler) for each of the isolates:
* Acinetobacter_baumannii_J9: 157-713
* Citrobacter_koseri_MINF_9D: 159-721
* Enterobacter_kobei_MSB1_1B: 160-717
* Haemophilus_M1C132_1: 160-757
* Klebsiella_oxytoca_MSB1_2C: 157-716
* Klebsiella_variicola_INF345: 156-705
* Serratia_marcescens_17-147-1671: 158-755






## Flye assemblies

Do Flye assemblies:
```bash
/usr/bin/time -v -o flye.time flye --nano-raw ligation.fastq.gz --genome-size "$genome_size" --threads 32 --plasmids --out-dir 01_flye_ligation &> flye.out
mv flye.time flye.out 01_flye_ligation
/usr/bin/time -v -o flye.time flye --nano-raw rapid.fastq.gz --genome-size "$genome_size" --threads 32 --plasmids --out-dir 02_flye_rapid &> flye.out
mv flye.time flye.out 02_flye_rapid
```


Polish Flye assemblies with Medaka:
```bash
cd ~
. medaka/bin/activate

cd ~/trycycler_real_read_tests/"$s"/01_flye_ligation
/usr/bin/time -v -o medaka.time medaka_consensus -i ../ligation.fastq.gz -d assembly.fasta -o medaka -m r941_min_high_g360 -t 32
mv medaka/consensus.fasta medaka.fasta
rm -r medaka *.fai *.mmi

cd ~/trycycler_real_read_tests/"$s"/02_flye_rapid
/usr/bin/time -v -o medaka.time medaka_consensus -i ../rapid.fastq.gz -d assembly.fasta -o medaka -m r941_min_high_g360 -t 32
mv medaka/consensus.fasta medaka.fasta
rm -r medaka *.fai *.mmi

deactivate
```


Polish Flye assemblies with Pilon:
```bash
# Since all real genomes had similar insert size distributions, I used these values for all of them:
insert_min=158
insert_max=723


cd 01_flye_ligation

before=medaka
after=round_1
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_a_1.fastq.gz -2 ../illumina_a_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_a_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=round_1
after=round_2
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_a_1.fastq.gz -2 ../illumina_a_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_a_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls round_*.fasta | tail -n1)
mv "$last_round" polished.fasta
rm round_*.fasta
cd ..


cd 02_flye_rapid

before=medaka
after=round_1
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_b_1.fastq.gz -2 ../illumina_b_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_b_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=round_1
after=round_2
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_b_1.fastq.gz -2 ../illumina_b_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_b_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls round_*.fasta | tail -n1)
mv "$last_round" polished.fasta
rm round_*.fasta
cd ..
```




## Trycycler assemblies

Do a ligation Trycycler assembly:
```bash
mkdir 05_trycycler_ligation
cp ligation.fastq.gz 05_trycycler_ligation/reads.fastq.gz
mkdir 05_trycycler_ligation/assemblies
target_depth=50
mean_length=$(seqtk comp ligation.fastq.gz | awk '{count++; bases += $2} END{print bases/count}')
read_count=$(echo $target_depth"*"$genome_size"/"$mean_length | bc)

for i in {00..11}; do
    seqtk sample -s "$i" ligation.fastq.gz "$read_count" | paste - - - - | shuf | tr '\t' '\n' > sample_"$i".fastq
done
for i in {00..02}; do
    /usr/bin/time -v -o 05_trycycler_ligation/assemblies/assembly_"$i".time flye --nano-raw sample_"$i".fastq --genome-size "$genome_size" --threads 32 --plasmids --out-dir assembly_"$i"; cp assembly_"$i"/assembly.fasta 05_trycycler_ligation/assemblies/assembly_"$i".fasta; rm -r assembly_"$i"
done
for i in {03..05}; do
    /usr/bin/time -v -o 05_trycycler_ligation/assemblies/assembly_"$i".time miniasm_and_minipolish.sh sample_"$i".fastq 32 > assembly_"$i".gfa; any2fasta assembly_"$i".gfa > 05_trycycler_ligation/assemblies/assembly_"$i".fasta; rm assembly_"$i".gfa
done
for i in {06..08}; do
    /usr/bin/time -v -o 05_trycycler_ligation/assemblies/assembly_"$i".time raven --threads 32 sample_"$i".fastq > 05_trycycler_ligation/assemblies/assembly_"$i".fasta; rm raven.cereal
done
for i in {09..11}; do
    /usr/bin/time -v -o 05_trycycler_ligation/assemblies/assembly_"$i".time wtdbg2.pl -o assembly_"$i" -g "$genome_size" -t 32 -x ont sample_"$i".fastq; mv assembly_"$i".cns.fa 05_trycycler_ligation/assemblies/assembly_"$i".fasta; rm assembly_"$i".*
done
rm sample_*.fastq

/usr/bin/time -v -o 05_trycycler_ligation/cluster.time ~/Trycycler/trycycler-runner.py cluster --reads 05_trycycler_ligation/reads.fastq.gz --assemblies 05_trycycler_ligation/assemblies/*.fasta --out_dir 05_trycycler_ligation/trycycler --threads 32 &> 05_trycycler_ligation/cluster.out
# Choose clusters

/usr/bin/time -v -o 05_trycycler_ligation/trycycler/cluster_001/reconcile.time ~/Trycycler/trycycler-runner.py reconcile --reads 05_trycycler_ligation/reads.fastq.gz --cluster_dir 05_trycycler_ligation/trycycler/cluster_001 --threads 32
/usr/bin/time -v -o 05_trycycler_ligation/trycycler/cluster_002/reconcile.time ~/Trycycler/trycycler-runner.py reconcile --reads 05_trycycler_ligation/reads.fastq.gz --cluster_dir 05_trycycler_ligation/trycycler/cluster_002 --threads 32
/usr/bin/time -v -o 05_trycycler_ligation/trycycler/cluster_003/reconcile.time ~/Trycycler/trycycler-runner.py reconcile --reads 05_trycycler_ligation/reads.fastq.gz --cluster_dir 05_trycycler_ligation/trycycler/cluster_003 --threads 32
# Fix reconciliation issues

for c in 05_trycycler_ligation/trycycler/cluster_*; do
    /usr/bin/time -v -o "$c"/msa.time ~/Trycycler/trycycler-runner.py msa --cluster_dir "$c" --threads 32
done

/usr/bin/time -v -o 05_trycycler_ligation/partition.time ~/Trycycler/trycycler-runner.py partition --reads 05_trycycler_ligation/reads.fastq.gz --cluster_dirs 05_trycycler_ligation/trycycler/cluster_* --threads 32

for c in 05_trycycler_ligation/trycycler/cluster_*; do
    /usr/bin/time -v -o "$c"/consensus.time ~/Trycycler/trycycler-runner.py consensus --cluster_dir "$c" --threads 32
done

cat 05_trycycler_ligation/trycycler/cluster_*/7_final_consensus.fasta > 05_trycycler_ligation/assembly.fasta
```


Polish the ligation Trycycler assembly with Medaka:
```bash
cd ~
. medaka/bin/activate
cd ~/trycycler_real_read_tests/"$s"/05_trycycler_ligation/trycycler

for c in cluster_*; do
    /usr/bin/time -v -o "$c"/medaka.time medaka_consensus -i "$c"/4_reads.fastq -d "$c"/7_final_consensus.fasta -o "$c"/medaka -m r941_min_high_g360 -t 32
    mv "$c"/medaka/consensus.fasta "$c"/8_medaka.fasta
    rm -r "$c"/medaka "$c"/*.fai "$c"/*.mmi
done
deactivate

cd ~/trycycler_real_read_tests/"$s"
cat 05_trycycler_ligation/trycycler/cluster_*/8_medaka.fasta > 05_trycycler_ligation/medaka.fasta
```


Polish the ligation Trycycler assembly with Pilon:
```bash
cd 05_trycycler_ligation

# Since all real genomes had similar insert size distributions, I used these values for all of them:
insert_min=158
insert_max=723

before=medaka
after=round_1
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_a_1.fastq.gz -2 ../illumina_a_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_a_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=round_1
after=round_2
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_a_1.fastq.gz -2 ../illumina_a_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_a_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls round_*.fasta | tail -n1)
mv "$last_round" polished.fasta
rm round_*.fasta
cd ..
```


Do a rapid Trycycler assembly:
```bash
mkdir 06_trycycler_rapid
cp rapid.fastq.gz 06_trycycler_rapid/reads.fastq.gz
mkdir 06_trycycler_rapid/assemblies
target_depth=50
mean_length=$(seqtk comp rapid.fastq.gz | awk '{count++; bases += $2} END{print bases/count}')
read_count=$(echo $target_depth"*"$genome_size"/"$mean_length | bc)

for i in {00..11}; do
    seqtk sample -s "$i" rapid.fastq.gz "$read_count" | paste - - - - | shuf | tr '\t' '\n' > sample_"$i".fastq
done
for i in {00..02}; do
    /usr/bin/time -v -o 06_trycycler_rapid/assemblies/assembly_"$i".time flye --nano-raw sample_"$i".fastq --genome-size "$genome_size" --threads 32 --plasmids --out-dir assembly_"$i"; cp assembly_"$i"/assembly.fasta 06_trycycler_rapid/assemblies/assembly_"$i".fasta; rm -r assembly_"$i"
done
for i in {03..05}; do
    /usr/bin/time -v -o 06_trycycler_rapid/assemblies/assembly_"$i".time miniasm_and_minipolish.sh sample_"$i".fastq 32 > assembly_"$i".gfa; any2fasta assembly_"$i".gfa > 06_trycycler_rapid/assemblies/assembly_"$i".fasta; rm assembly_"$i".gfa
done
for i in {06..08}; do
    /usr/bin/time -v -o 06_trycycler_rapid/assemblies/assembly_"$i".time raven --threads 32 sample_"$i".fastq > 06_trycycler_rapid/assemblies/assembly_"$i".fasta; rm raven.cereal
done
for i in {09..11}; do
    /usr/bin/time -v -o 06_trycycler_rapid/assemblies/assembly_"$i".time wtdbg2.pl -o assembly_"$i" -g "$genome_size" -t 32 -x ont sample_"$i".fastq; mv assembly_"$i".cns.fa 06_trycycler_rapid/assemblies/assembly_"$i".fasta; rm assembly_"$i".*
done
rm sample_*.fastq

/usr/bin/time -v -o 06_trycycler_rapid/cluster.time ~/Trycycler/trycycler-runner.py cluster --reads 06_trycycler_rapid/reads.fastq.gz --assemblies 06_trycycler_rapid/assemblies/*.fasta --out_dir 06_trycycler_rapid/trycycler --threads 32 &> 06_trycycler_rapid/cluster.out
# Choose clusters

/usr/bin/time -v -o 06_trycycler_rapid/trycycler/cluster_001/reconcile.time ~/Trycycler/trycycler-runner.py reconcile --reads 06_trycycler_rapid/reads.fastq.gz --cluster_dir 06_trycycler_rapid/trycycler/cluster_001 --threads 32
/usr/bin/time -v -o 06_trycycler_rapid/trycycler/cluster_002/reconcile.time ~/Trycycler/trycycler-runner.py reconcile --reads 06_trycycler_rapid/reads.fastq.gz --cluster_dir 06_trycycler_rapid/trycycler/cluster_002 --threads 32
/usr/bin/time -v -o 06_trycycler_rapid/trycycler/cluster_003/reconcile.time ~/Trycycler/trycycler-runner.py reconcile --reads 06_trycycler_rapid/reads.fastq.gz --cluster_dir 06_trycycler_rapid/trycycler/cluster_003 --threads 32
# Fix reconciliation issues

for c in 06_trycycler_rapid/trycycler/cluster_*; do
    /usr/bin/time -v -o "$c"/msa.time ~/Trycycler/trycycler-runner.py msa --cluster_dir "$c" --threads 32
done

/usr/bin/time -v -o 06_trycycler_rapid/partition.time ~/Trycycler/trycycler-runner.py partition --reads 06_trycycler_rapid/reads.fastq.gz --cluster_dirs 06_trycycler_rapid/trycycler/cluster_* --threads 32

for c in 06_trycycler_rapid/trycycler/cluster_*; do
    /usr/bin/time -v -o "$c"/consensus.time ~/Trycycler/trycycler-runner.py consensus --cluster_dir "$c" --threads 32
done

cat 06_trycycler_rapid/trycycler/cluster_*/7_final_consensus.fasta > 06_trycycler_rapid/assembly.fasta
```


Polish the rapid Trycycler assembly with Medaka:
```bash
cd ~
. medaka/bin/activate
cd ~/trycycler_real_read_tests/"$s"/06_trycycler_rapid/trycycler

for c in cluster_*; do
    /usr/bin/time -v -o "$c"/medaka.time medaka_consensus -i "$c"/4_reads.fastq -d "$c"/7_final_consensus.fasta -o "$c"/medaka -m r941_min_high_g360 -t 32
    mv "$c"/medaka/consensus.fasta "$c"/8_medaka.fasta
    rm -r "$c"/medaka "$c"/*.fai "$c"/*.mmi
done
deactivate

cd ~/trycycler_real_read_tests/"$s"
cat 06_trycycler_rapid/trycycler/cluster_*/8_medaka.fasta > 06_trycycler_rapid/medaka.fasta
```


Polish the rapid Trycycler assembly with Pilon:
```bash
cd 06_trycycler_rapid

# Since all real genomes had similar insert size distributions, I used these values for all of them:
insert_min=158
insert_max=723

before=medaka
after=round_1
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_b_1.fastq.gz -2 ../illumina_b_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_b_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=round_1
after=round_2
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_b_1.fastq.gz -2 ../illumina_b_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_b_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls round_*.fasta | tail -n1)
mv "$last_round" polished.fasta
rm round_*.fasta
cd ..
```











## Pilon experiment

As an experiment, I also tried Pilon polishing with its `--nanopore` option, to see if that did any better:

Polish the ligation Trycycler assembly with hybrid Pilon:
```bash
cd 05_trycycler_ligation

# Since all real genomes had similar insert size distributions, I used these values for all of them:
insert_min=158
insert_max=723

before=medaka
after=hybrid_round_1
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_a_1.fastq.gz -2 ../illumina_a_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_a_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
minimap2 -a -x map-ont -t 32 "$before".fasta ../ligation.fastq.gz | samtools sort > nanopore_alignments.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
samtools index nanopore_alignments.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --nanopore nanopore_alignments.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=hybrid_round_1
after=hybrid_round_2
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_a_1.fastq.gz -2 ../illumina_a_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_a_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
minimap2 -a -x map-ont -t 32 "$before".fasta ../ligation.fastq.gz | samtools sort > nanopore_alignments.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
samtools index nanopore_alignments.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --nanopore nanopore_alignments.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls hybrid_round_*.fasta | tail -n1)
mv "$last_round" hybrid_polished.fasta
rm hybrid_round_*.fasta
cd ..
```

Polish the rapid Trycycler assembly with hybrid Pilon:
```bash
cd 06_trycycler_rapid

# Since all real genomes had similar insert size distributions, I used these values for all of them:
insert_min=158
insert_max=723

before=medaka
after=hybrid_round_1
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_b_1.fastq.gz -2 ../illumina_b_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_b_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
minimap2 -a -x map-ont -t 32 "$before".fasta ../rapid.fastq.gz | samtools sort > nanopore_alignments.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
samtools index nanopore_alignments.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --nanopore nanopore_alignments.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=hybrid_round_1
after=hybrid_round_2
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_b_1.fastq.gz -2 ../illumina_b_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_b_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
minimap2 -a -x map-ont -t 32 "$before".fasta ../rapid.fastq.gz | samtools sort > nanopore_alignments.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
samtools index nanopore_alignments.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --nanopore nanopore_alignments.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls hybrid_round_*.fasta | tail -n1)
mv "$last_round" hybrid_polished.fasta
rm hybrid_round_*.fasta
cd ..
```

However, I found this to be worse than Illumina-only Pilon polishing, so I didn't use these assemblies.





## Other assemblers

I later decided to throw a couple other assemblers into the mix, so here are their assembly commands:


### Raven

Assemblies:
```bash
/usr/bin/time -v -o raven.time raven --disable-checkpoints --graphical-fragment-assembly raven.gfa -t 32 ligation.fastq.gz 1> raven.fasta 2> raven.out
mkdir 07_raven_ligation; mv raven.time raven.fasta raven.gfa raven.out 07_raven_ligation

/usr/bin/time -v -o raven.time raven --disable-checkpoints --graphical-fragment-assembly raven.gfa -t 32 rapid.fastq.gz 1> raven.fasta 2> raven.out
mkdir 08_raven_rapid; mv raven.time raven.fasta raven.gfa raven.out 08_raven_rapid
```

Medaka:
```bash
cd ~
. medaka/bin/activate

cd ~/trycycler_real_read_tests/"$s"/07_raven_ligation
/usr/bin/time -v -o medaka.time medaka_consensus -i ../ligation.fastq.gz -d raven.fasta -o medaka -m r941_min_high_g360 -t 32
mv medaka/consensus.fasta medaka.fasta
rm -r medaka *.fai *.mmi

cd ~/trycycler_real_read_tests/"$s"/08_raven_rapid
/usr/bin/time -v -o medaka.time medaka_consensus -i ../rapid.fastq.gz -d raven.fasta -o medaka -m r941_min_high_g360 -t 32
mv medaka/consensus.fasta medaka.fasta
rm -r medaka *.fai *.mmi

deactivate
```

Pilon (ligation):
```bash
cd ~/trycycler_real_read_tests/"$s"/07_raven_ligation

# Since all real genomes had similar insert size distributions, I used these values for all of them:
insert_min=158
insert_max=723

before=medaka
after=round_1
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_a_1.fastq.gz -2 ../illumina_a_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_a_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=round_1
after=round_2
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_a_1.fastq.gz -2 ../illumina_a_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_a_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls round_*.fasta | tail -n1)
mv "$last_round" polished.fasta
rm round_*.fasta
cd ..
```

Pilon (rapid):
```bash
cd ~/trycycler_real_read_tests/"$s"/08_raven_rapid

# Since all real genomes had similar insert size distributions, I used these values for all of them:
insert_min=158
insert_max=723

before=medaka
after=round_1
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_b_1.fastq.gz -2 ../illumina_b_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_b_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=round_1
after=round_2
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_b_1.fastq.gz -2 ../illumina_b_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_b_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls round_*.fasta | tail -n1)
mv "$last_round" polished.fasta
rm round_*.fasta
cd ..
```


### Miniasm/Minipolish

Assemblies:
```bash
/usr/bin/time -v -o miniasm.time miniasm_and_minipolish.sh ligation.fastq.gz 32 1> miniasm.gfa 2> miniasm.out
any2fasta miniasm.gfa > miniasm.fasta
mkdir 09_miniasm_ligation; mv miniasm.time miniasm.fasta miniasm.gfa miniasm.out 09_miniasm_ligation

/usr/bin/time -v -o miniasm.time miniasm_and_minipolish.sh rapid.fastq.gz 32 1> miniasm.gfa 2> miniasm.out
any2fasta miniasm.gfa > miniasm.fasta
mkdir 10_miniasm_rapid; mv miniasm.time miniasm.fasta miniasm.gfa miniasm.out 10_miniasm_rapid
```

Medaka:
```bash
cd ~
. medaka/bin/activate

cd ~/trycycler_real_read_tests/"$s"/09_miniasm_ligation
/usr/bin/time -v -o medaka.time medaka_consensus -i ../ligation.fastq.gz -d miniasm.fasta -o medaka -m r941_min_high_g360 -t 32
mv medaka/consensus.fasta medaka.fasta
rm -r medaka *.fai *.mmi

cd ~/trycycler_real_read_tests/"$s"/10_miniasm_rapid
/usr/bin/time -v -o medaka.time medaka_consensus -i ../rapid.fastq.gz -d miniasm.fasta -o medaka -m r941_min_high_g360 -t 32
mv medaka/consensus.fasta medaka.fasta
rm -r medaka *.fai *.mmi

deactivate
```

Pilon (ligation):
```bash
cd ~/trycycler_real_read_tests/"$s"/09_miniasm_ligation

# Since all real genomes had similar insert size distributions, I used these values for all of them:
insert_min=158
insert_max=723

before=medaka
after=round_1
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_a_1.fastq.gz -2 ../illumina_a_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_a_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=round_1
after=round_2
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_a_1.fastq.gz -2 ../illumina_a_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_a_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls round_*.fasta | tail -n1)
mv "$last_round" polished.fasta
rm round_*.fasta
cd ..
```

Pilon (rapid):
```bash
cd ~/trycycler_real_read_tests/"$s"/10_miniasm_rapid

# Since all real genomes had similar insert size distributions, I used these values for all of them:
insert_min=158
insert_max=723

before=medaka
after=round_1
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_b_1.fastq.gz -2 ../illumina_b_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_b_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

before=round_1
after=round_2
bowtie2-build "$before".fasta "$before".fasta > /dev/null
bowtie2 -1 ../illumina_b_1.fastq.gz -2 ../illumina_b_2.fastq.gz -x "$before".fasta --threads 32 -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U ../illumina_b_u.fastq.gz -x "$before".fasta --threads 32 --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
seqtk seq "$after".fasta > temp; mv temp "$after".fasta

# Repeat until no more changes (up to a maximum of 5 rounds).

last_round=$(ls round_*.fasta | tail -n1)
mv "$last_round" polished.fasta
rm round_*.fasta
cd ..
```





## Medaka before Trycycler

I later decided to test an alternative approach: running Medaka on Trycycler's input assemblies.

```bash
cd ~/trycycler_real_read_tests/"$s"
mkdir 11_medaka_trycycler_ligation
cp ligation.fastq.gz 11_medaka_trycycler_ligation/reads.fastq.gz
cp -r 05_trycycler_ligation/assemblies 11_medaka_trycycler_ligation/assemblies
target_depth=50
mean_length=$(seqtk comp ligation.fastq.gz | awk '{count++; bases += $2} END{print bases/count}')
read_count=$(echo $target_depth"*"$genome_size"/"$mean_length | bc)

for i in {00..11}; do
    seqtk sample -s "$i" ligation.fastq.gz "$read_count" | paste - - - - | shuf | tr '\t' '\n' > sample_"$i".fastq
done

mkdir 11_medaka_trycycler_ligation/medaka_assemblies
for i in {00..11}; do
    medaka_consensus -i sample_"$i".fastq -d 11_medaka_trycycler_ligation/assemblies/assembly_"$i".fasta -o medaka_temp -m r941_min_high_g360 -t 32
    mv medaka_temp/consensus.fasta 11_medaka_trycycler_ligation/medaka_assemblies/assembly_"$i".fasta
    rm -r medaka_temp 11_medaka_trycycler_ligation/assemblies/*.fai 11_medaka_trycycler_ligation/assemblies/*.mmi
done
rm sample_*.fastq

~/programs/Trycycler/trycycler-runner.py cluster --reads 11_medaka_trycycler_ligation/reads.fastq.gz --assemblies 11_medaka_trycycler_ligation/medaka_assemblies/*.fasta --out_dir 11_medaka_trycycler_ligation/trycycler --threads 32 &> 11_medaka_trycycler_ligation/cluster.out
# Choose clusters

~/programs/Trycycler/trycycler-runner.py reconcile  --threads 32 --reads 11_medaka_trycycler_ligation/reads.fastq.gz --cluster_dir 11_medaka_trycycler_ligation/trycycler/cluster_001
~/programs/Trycycler/trycycler-runner.py reconcile  --threads 32 --reads 11_medaka_trycycler_ligation/reads.fastq.gz --cluster_dir 11_medaka_trycycler_ligation/trycycler/cluster_002
~/programs/Trycycler/trycycler-runner.py reconcile  --threads 32 --reads 11_medaka_trycycler_ligation/reads.fastq.gz --cluster_dir 11_medaka_trycycler_ligation/trycycler/cluster_003
# Fix reconciliation issues

for c in 11_medaka_trycycler_ligation/trycycler/cluster_*; do
    ~/programs/Trycycler/trycycler-runner.py msa --cluster_dir "$c" --threads 32
done

~/programs/Trycycler/trycycler-runner.py partition --reads 11_medaka_trycycler_ligation/reads.fastq.gz --cluster_dirs 11_medaka_trycycler_ligation/trycycler/cluster_* --threads 32

for c in 11_medaka_trycycler_ligation/trycycler/cluster_*; do
    ~/programs/Trycycler/trycycler-runner.py consensus --cluster_dir "$c" --threads 32
done
cat 11_medaka_trycycler_ligation/trycycler/cluster_*/7_final_consensus.fasta > 11_medaka_trycycler_ligation/assembly.fasta

for c in 11_medaka_trycycler_ligation/trycycler/cluster_*; do
    medaka_consensus -i "$c"/4_reads.fastq -d "$c"/7_final_consensus.fasta -o "$c"/medaka -m r941_min_high_g360 -t 32
    mv "$c"/medaka/consensus.fasta "$c"/8_medaka.fasta
    rm -r "$c"/medaka "$c"/*.fai "$c"/*.mmi
done
cat 11_medaka_trycycler_ligation/trycycler/cluster_*/8_medaka.fasta > 11_medaka_trycycler_ligation/medaka.fasta
```

```bash
cd ~/trycycler_real_read_tests/"$s"
mkdir 12_medaka_trycycler_rapid
cp rapid.fastq.gz 12_medaka_trycycler_rapid/reads.fastq.gz
cp -r 06_trycycler_rapid/assemblies 12_medaka_trycycler_rapid/assemblies
target_depth=50
mean_length=$(seqtk comp rapid.fastq.gz | awk '{count++; bases += $2} END{print bases/count}')
read_count=$(echo $target_depth"*"$genome_size"/"$mean_length | bc)

for i in {00..11}; do
    seqtk sample -s "$i" rapid.fastq.gz "$read_count" | paste - - - - | shuf | tr '\t' '\n' > sample_"$i".fastq
done

mkdir 12_medaka_trycycler_rapid/medaka_assemblies
for i in {00..11}; do
    medaka_consensus -i sample_"$i".fastq -d 12_medaka_trycycler_rapid/assemblies/assembly_"$i".fasta -o medaka_temp -m r941_min_high_g360 -t 32
    mv medaka_temp/consensus.fasta 12_medaka_trycycler_rapid/medaka_assemblies/assembly_"$i".fasta
    rm -r medaka_temp 12_medaka_trycycler_rapid/assemblies/*.fai 12_medaka_trycycler_rapid/assemblies/*.mmi
done
rm sample_*.fastq

~/programs/Trycycler/trycycler-runner.py cluster --reads 12_medaka_trycycler_rapid/reads.fastq.gz --assemblies 12_medaka_trycycler_rapid/medaka_assemblies/*.fasta --out_dir 12_medaka_trycycler_rapid/trycycler --threads 32 &> 12_medaka_trycycler_rapid/cluster.out
# Choose clusters

~/programs/Trycycler/trycycler-runner.py reconcile --threads 32 --reads 12_medaka_trycycler_rapid/reads.fastq.gz --cluster_dir 12_medaka_trycycler_rapid/trycycler/cluster_001
~/programs/Trycycler/trycycler-runner.py reconcile --threads 32 --reads 12_medaka_trycycler_rapid/reads.fastq.gz --cluster_dir 12_medaka_trycycler_rapid/trycycler/cluster_002
~/programs/Trycycler/trycycler-runner.py reconcile --threads 32 --reads 12_medaka_trycycler_rapid/reads.fastq.gz --cluster_dir 12_medaka_trycycler_rapid/trycycler/cluster_003
# Fix reconciliation issues

for c in 12_medaka_trycycler_rapid/trycycler/cluster_*; do
    ~/programs/Trycycler/trycycler-runner.py msa --cluster_dir "$c" --threads 32
done

~/programs/Trycycler/trycycler-runner.py partition --reads 12_medaka_trycycler_rapid/reads.fastq.gz --cluster_dirs 12_medaka_trycycler_rapid/trycycler/cluster_* --threads 32

for c in 12_medaka_trycycler_rapid/trycycler/cluster_*; do
    ~/programs/Trycycler/trycycler-runner.py consensus --cluster_dir "$c" --threads 32
done
cat 12_medaka_trycycler_rapid/trycycler/cluster_*/7_final_consensus.fasta > 12_medaka_trycycler_rapid/assembly.fasta

for c in 12_medaka_trycycler_rapid/trycycler/cluster_*; do
    medaka_consensus -i "$c"/4_reads.fastq -d "$c"/7_final_consensus.fasta -o "$c"/medaka -m r941_min_high_g360 -t 32
    mv "$c"/medaka/consensus.fasta "$c"/8_medaka.fasta
    rm -r "$c"/medaka "$c"/*.fai "$c"/*.mmi
done
cat 12_medaka_trycycler_rapid/trycycler/cluster_*/8_medaka.fasta > 12_medaka_trycycler_rapid/medaka.fasta
```
