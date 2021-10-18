<p align="center"><img src="https://github.com/rrwick/Trycycler/blob/main/images/logo_transparent.png" alt="Trycycler" width="70%"></p>

This repo contains supplementary data for our paper describing Trycycler: [Wick RR, Judd LM, Cerdeira LT, Hawkey J, Méric G, Vezina B, Wyres KL, Holt KE. Trycycler: consensus long-read assemblies for bacterial genomes. Genome Biology. 2021. doi:10.1186/s13059-021-02483-z.](https://doi.org/10.1186/s13059-021-02483-z)

If you're interested in Trycycler itself (i.e. not the paper), then head over to [its GitHub repo](https://github.com/rrwick/Trycycler).

This repo has the following subdirectories:
* [`1_simulated_read_tests`](1_simulated_read_tests): the methods, data and results for the simulated read tests. These are the tests where we assembled simulated read sets and assessed accuracy by comparing the assemblies back to the original reference genome.
* [`2_real_read_tests`](2_real_read_tests): the methods, data and results for the real read tests. These are the tests where we assembled two independent real read sets for each genome and assessed accuracy by comparing the two assemblies to each other.
* [`3_consistency_tests`](3_consistency_tests): the methods, data and results for the consistency tests. These are the tests where we had multiple different people assemble the same read sets using Trycycler to assess how consistently Trycycler performs in different hands.
* [`figures`](figures): contains all figures (both main text and supplementary) for the paper.
* [`scripts`](scripts): contains custom Python scripts we wrote to conduct the analyses in the paper.

Summarised results for each of the sections are available in the following supplementary tables, each an Excel file with multiple worksheets:
* [table_s1_simulated_read_results.xlsx](table_s1_simulated_read_results.xlsx): results for the simulated read tests
* [table_s2_real_read_results.xlsx](table_s2_real_read_results.xlsx): results for the real read tests
* [table_s3_consistency_test_results.xlsx](table_s3_consistency_test_results.xlsx): results for the multi-user consistency tests

The datasets used in this paper (reads, assemblies and references) are too large for a GitHub repo, so you can download them from here: [bridges.monash.edu/articles/dataset/Trycycler_paper_dataset/14890734](https://bridges.monash.edu/articles/dataset/Trycycler_paper_dataset/14890734)
