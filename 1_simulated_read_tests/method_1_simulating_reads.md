I downloaded a report of NCBI reference genomes (got the link from https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/#reference_genomes):
```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests
wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prok_reference_genomes.txt
```
I then loaded it in Excel and sorted by genome count to get the most common species and took the top 10. Sometimes there were more than one representative for a species (e.g. E. coli) and in these cases I took the one with the lowest chromosome RefSeq accession (i.e. the oldest).

I then went to the NCBI Assembly site (https://www.ncbi.nlm.nih.gov/assembly) and searched for each of the chromosome accessions to get an assembly RefSeq accession and download the assembly. I put these in `/Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes` in a subdirectory named for the species.


Make them one-line per sequence:
```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes
for f in */*.fna; do seqtk seq $f > temp; mv temp $f; done
```

I then manually added 'circular=true' and 'depth=' to each FASTA header, where the depth was 1 for chromosomes and something length-dependent for plasmids - I just made up the values, with small plasmids getting higher depths.









## Illumina read simulation

To get ART to simulate reads over circular genomes and to adjust for per-replicon read depth, I ran the references through a script (`prep_for_illumina_read_simulation.py`) which multiplied the sequences.

For the ART commands, I used the 7 built-in read profiles which had a read length of at least 100 bp (excluding the HiSeq 1000 profile which gave bad results). Since I had 10 genomes, 4 of these profiles were used 1 time and 3 were used twice. I took depths from 40x to 130x in 10x increments (a total of 10 depths), randomly shuffled between the genomes. Mean fragment lengths were from 300 bp to 660 bp in 50 bp increments, randomly shuffled. The fragment length stdevs were randomly sampled (from a uniform distribution) from the range of 20 to 1/4 the mean fragment length.

```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Campylobacter_jejuni
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys MSv1 --len 250 --fcov 09 --mflen 340 --sdev 053 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Escherichia_coli
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HS25 --len 150 --fcov 11 --mflen 420 --sdev 021 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Klebsiella_pneumoniae
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HSXt --len 150 --fcov 04 --mflen 500 --sdev 122 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Listeria_monocytogenes
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys MSv3 --len 250 --fcov 13 --mflen 300 --sdev 030 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Mycobacterium_tuberculosis
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HS25 --len 125 --fcov 07 --mflen 460 --sdev 031 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Neisseria_meningitidis
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HS20 --len 100 --fcov 12 --mflen 540 --sdev 045 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Pseudomonas_aeruginosa
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HS25 --len 125 --fcov 08 --mflen 660 --sdev 051 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Salmonella_enterica
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HS20 --len 100 --fcov 05 --mflen 580 --sdev 061 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Staphylococcus_aureus
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HS25 --len 150 --fcov 06 --mflen 380 --sdev 081 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Streptococcus_pneumoniae
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HSXn --len 150 --fcov 10 --mflen 620 --sdev 117 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta
```



## Nanopore read simulation

I used MUMmer to get the longest repeat in each genome (so I could ensure that the long reads were longer):
```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes
for f in */*.fna; do
    printf $f"\t"
    nucmer --maxmatch --nosimplify --prefix=seq_seq $f $f 2> /dev/null
    show-coords -H -T -r seq_seq.delta | awk '{if ($2-$1 < 1000000 && $2-$1 > 1000 && $1 != $3 && $2 != $4) print $0"\t"($2-$1+1);}' | sort -nk10,10 | tail -n1
    rm seq_seq.delta
done
```


For Badread parameters:
* mean read length was sampled from 5k over the max repeat length to 40k
* read length stdev was sampled from 1/4 of the mean length to the mean length
* max identity was sampled from 95 to 100
* mean identity was sampled from 90 to the max identity minus one
* identity stdev was sampled from 1 to 5
* junk and random percentages were set to 0 (so Flye didn't run too slowly)
* glitch rate was sampled from 1k to 10k
* glitch size and skip were sampled (separately) from 0 to 50
* start and end adapter sequence lengths were sampled (separately) from 0 to 50
* depth ranged from 100x to 290x in steps of 10x (shuffled)



```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Campylobacter_jejuni
badread simulate --reference *.fna --quantity 200x --length 20788,06361 --identity 94,097,2 --start_adapter_seq 26 --end_adapter_seq 15 --junk 0 --random 0 --chimeras 1.2 --glitches 7083,40,25 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Escherichia_coli
badread simulate --reference *.fna --quantity 280x --length 26103,08384 --identity 93,095,5 --start_adapter_seq 02 --end_adapter_seq 40 --junk 0 --random 0 --chimeras 1.2 --glitches 5766,28,08 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Klebsiella_pneumoniae
badread simulate --reference *.fna --quantity 210x --length 14688,12402 --identity 91,095,5 --start_adapter_seq 02 --end_adapter_seq 44 --junk 0 --random 0 --chimeras 2.1 --glitches 3422,44,04 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Listeria_monocytogenes
badread simulate --reference *.fna --quantity 150x --length 30814,24000 --identity 96,100,3 --start_adapter_seq 09 --end_adapter_seq 16 --junk 0 --random 0 --chimeras 1.1 --glitches 8367,28,14 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Mycobacterium_tuberculosis
badread simulate --reference *.fna --quantity 270x --length 25285,19637 --identity 93,096,4 --start_adapter_seq 02 --end_adapter_seq 19 --junk 0 --random 0 --chimeras 2.5 --glitches 7377,10,28 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Neisseria_meningitidis
badread simulate --reference *.fna --quantity 250x --length 37314,27176 --identity 96,098,5 --start_adapter_seq 27 --end_adapter_seq 14 --junk 0 --random 0 --chimeras 0.6 --glitches 6710,44,32 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Pseudomonas_aeruginosa
badread simulate --reference *.fna --quantity 260x --length 35001,18994 --identity 93,100,5 --start_adapter_seq 21 --end_adapter_seq 02 --junk 0 --random 0 --chimeras 1.9 --glitches 6726,25,27 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Salmonella_enterica
badread simulate --reference *.fna --quantity 120x --length 35097,26894 --identity 91,097,2 --start_adapter_seq 43 --end_adapter_seq 12 --junk 0 --random 0 --chimeras 2.3 --glitches 8267,31,28 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Staphylococcus_aureus
badread simulate --reference *.fna --quantity 190x --length 36588,29376 --identity 93,095,4 --start_adapter_seq 23 --end_adapter_seq 14 --junk 0 --random 0 --chimeras 0.4 --glitches 4138,03,07 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Streptococcus_pneumoniae
badread simulate --reference *.fna --quantity 230x --length 13506,05927 --identity 93,096,3 --start_adapter_seq 43 --end_adapter_seq 02 --junk 0 --random 0 --chimeras 2.2 --glitches 7045,14,11 2> badread.out | gzip > nanopore.fastq.gz
```







































Commands for the top 11-20 (didn't use):
```bash
cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Bacillus_cereus
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HS10 --len 100 --fcov 13.5 --mflen 540 --sdev 114 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Burkholderia_pseudomallei
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HS20 --len 100 --fcov 08.0 --mflen 300 --sdev 022 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Clostridioides_difficile
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HS25 --len 150 --fcov 07.0 --mflen 400 --sdev 074 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Enterococcus_faecalis
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HSXn --len 150 --fcov 10.0 --mflen 200 --sdev 038 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Enterococcus_faecium
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HSXt --len 150 --fcov 06.0 --mflen 480 --sdev 086 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Helicobacter_pylori
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys MSv3 --len 250 --fcov 11.5 --mflen 320 --sdev 037 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Mycobacteroides_abscessus
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HS25 --len 150 --fcov 09.0 --mflen 240 --sdev 042 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Streptococcus_agalactiae
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HS10 --len 100 --fcov 11.0 --mflen 440 --sdev 043 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Streptococcus_suis
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HS25 --len 125 --fcov 08.5 --mflen 340 --sdev 060 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Vibrio_cholerae
../../../scripts/prep_for_illumina_read_simulation.py *.fna > temp.fasta
~/Applications/art/art_illumina --paired --seqSys HS25 --len 150 --fcov 07.5 --mflen 280 --sdev 048 --noALN --in temp.fasta --out illumina &> art.out
mv illumina1.fq illumina_1.fastq; mv illumina2.fq illumina_2.fastq; gzip *.fastq; rm temp.fasta

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Bacillus_cereus
badread simulate --reference *.fna --quantity 180x --length 29638,28127 --identity 90,097,2 --start_adapter_seq 24 --end_adapter_seq 32 --junk 0 --random 0 --chimeras 0.1 --glitches 8962,07,48 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Burkholderia_pseudomallei
badread simulate --reference *.fna --quantity 160x --length 27261,23051 --identity 94,099,3 --start_adapter_seq 39 --end_adapter_seq 17 --junk 0 --random 0 --chimeras 0.2 --glitches 7113,43,31 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Clostridioides_difficile
badread simulate --reference *.fna --quantity 110x --length 39136,36652 --identity 94,097,4 --start_adapter_seq 49 --end_adapter_seq 42 --junk 0 --random 0 --chimeras 0.9 --glitches 1332,30,44 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Enterococcus_faecalis
badread simulate --reference *.fna --quantity 290x --length 34232,22624 --identity 90,100,5 --start_adapter_seq 25 --end_adapter_seq 22 --junk 0 --random 0 --chimeras 0.7 --glitches 5410,14,29 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Enterococcus_faecium
badread simulate --reference *.fna --quantity 140x --length 21634,19123 --identity 93,098,2 --start_adapter_seq 01 --end_adapter_seq 14 --junk 0 --random 0 --chimeras 0.0 --glitches 4757,15,40 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Helicobacter_pylori
badread simulate --reference *.fna --quantity 100x --length 31867,11989 --identity 93,099,3 --start_adapter_seq 39 --end_adapter_seq 39 --junk 0 --random 0 --chimeras 1.2 --glitches 9720,50,45 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Mycobacteroides_abscessus
badread simulate --reference *.fna --quantity 240x --length 39269,38384 --identity 95,096,5 --start_adapter_seq 13 --end_adapter_seq 13 --junk 0 --random 0 --chimeras 1.9 --glitches 5747,47,05 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Streptococcus_agalactiae
badread simulate --reference *.fna --quantity 130x --length 17298,07349 --identity 90,098,4 --start_adapter_seq 34 --end_adapter_seq 42 --junk 0 --random 0 --chimeras 1.1 --glitches 3538,50,03 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Streptococcus_suis
badread simulate --reference *.fna --quantity 170x --length 21926,12410 --identity 92,096,3 --start_adapter_seq 31 --end_adapter_seq 37 --junk 0 --random 0 --chimeras 2.0 --glitches 3726,27,00 2> badread.out | gzip > nanopore.fastq.gz

cd /Users/ryan/Dropbox/Uni_research/Projects/Trycycler/PAPER/GitHub_repo/simulated_read_tests/genomes/Vibrio_cholerae
badread simulate --reference *.fna --quantity 220x --length 17132,06145 --identity 91,095,3 --start_adapter_seq 20 --end_adapter_seq 21 --junk 0 --random 0 --chimeras 0.3 --glitches 5316,00,29 2> badread.out | gzip > nanopore.fastq.gz
```
