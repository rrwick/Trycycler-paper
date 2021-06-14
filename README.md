<p align="center"><img src="https://github.com/rrwick/Trycycler/blob/main/images/logo.png" alt="Trycycler" width="70%"></p>

This repo contains supplementary data for our paper describing Trycycler: [Wick RR, Judd LM, Cerdeira LT, Hawkey J, Méric G, Vezina B, Wyres KL, Wyres KL, Holt KE. Trycycler: consensus long-read assemblies for prokaryote genomes. bioRxiv. 2021. doi:TBA.](https://www.biorxiv.org/content/TBA)

If you're interested in Trycycler itself (i.e. not the paper), then head over to [its GitHub repo](https://github.com/rrwick/Trycycler).

In this repo, you'll find the following subdirectories:
* [`1_simulated_read_tests`](1_simulated_read_tests): contains the methods, data and results for the simulated read tests. These are the tests where we assembled simulated read sets and assessed accuracy by comparing the assemblies back to the original reference genome.
* [`2_real_read_tests`](2_real_read_tests): contains the methods, data and results for the real read tests. These are the tests where we assembled two independent real read sets for each genome and assessed accuracy by comparing the two assemblies to each other.
* [`3_consistency_tests`](3_consistency_tests): contains the methods, data and results for the consistency tests. These are the tests where we had multiple different people assemble the same read sets using Trycycler to assess how consistently Trycycler performs in different hands.
* [`figures`](figures): contains all figures (both main text and supplementary) for the paper.
* [`scripts`](scripts): contains custom Python scripts we wrote to conduct the analyses in the paper.
