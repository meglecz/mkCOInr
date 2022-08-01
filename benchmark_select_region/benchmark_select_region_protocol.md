# Benchmark select_region.pl

## Download from NCBI

### Download complete mitogenomes
~~~
conda activate mkcoinr
cd ~/mkCOInr

nsdpy -r "complete[Title] AND genome[Title] AND Mitochondrion[Filter]" -T -v --cds

mv NSDPY_results/yyy-mm-dd_hh-mm-ss benchmark_select_region/download_mito_nsdpy 
mv report.tsv benchmark_select_region/download_mito_nsdpy 
rmdir NSDPY_results
~~~

### Download complete chloroplast genomes
~~~
conda activate mkcoinr
cd ~/mkCOInr

nsdpy -r "complete[Title] AND genome[Title] AND Chloroplast[Filter]" -T -v --cds

mv NSDPY_results/yyy-mm-dd_hh-mm-ss benchmark_select_region/download_chloro_nsdpy 
mv report.tsv benchmark_select_region/download_chloro_nsdpy 
rmdir NSDPY_results
~~~


## Make test datasets

### Positive mitochondrial dataset (COI sequences)

#### Select COI sequences from complete genome CDS file and clean them 
Eliminate identical sequences of the same taxID.
Clean tax names and taxids.
~~~
perl scripts/format_ncbi.pl -cds benchmark_select_region/download_mito_nsdpy/sequences.fasta -taxids benchmark_select_region/download_mito_nsdpy/TaxIDs.txt -taxonomy COInr/taxonomy.tsv -outdir benchmark_select_region/test_dataset/mitogenomes_COI/1_format_ncbi
~~~

#### Check_sequences

- add lineage
- select sequence of length between 1100 and 2000 bp

According to the length distribution of the COI genes, there are sequences that clearly uncomplete and complete genome is an error in the title (lenght <1000).
From 1100 there are COI genes with no evidence that it is incomplete
Longest sequences are around 2000.

~~~
perl benchmark_select_region/scripts/check_sequences.pl -tsv benchmark_select_region/test_dataset/mitogenomes_COI/1_format_ncbi/ncbi_sequences.tsv -outdir benchmark_select_region/test_dataset/mitogenomes_COI/2_check_sequences -taxonomy COInr/taxonomy.tsv -taxrank phylum -min_length 1100 -max_length 2000
~~~

#### Move the cleaned file that will be used as a **positive test dataset**
~~~
mv benchmark_select_region/test_dataset/mitogenomes_COI/2_check_sequences/sequences_checked.tsv benchmark_select_region/test_dataset/positive_mito.tsv
~~~

### Negative mitochondrial  (non-COI sequences)

#### Select all other than COI genes from mitogenomes
- Start from mitogenomes downloaded previously
- Eliminate COI sequences and clean the rest (@regexp_gene = ('cox*.*[1i]'), @regexp_protein = ('cytochrome', 'cox*.*[1i]'))
- Eliminate identical sequences of the same taxID.
- Eliminate coding sequences with introns
- Select sequence length between 700 and 2000 nt
- Clean tax names and taxids
- Keep the original cds seqid, so it is possible to get back to the original file 
- Keep only protein-coding sequences

~~~
perl benchmark_select_region/scripts/format_ncbi_negative.pl -cds benchmark_select_region/download_mito_nsdpy/sequences.fasta -taxids benchmark_select_region/download_mito_nsdpy/TaxIDs.txt -taxonomy COInr/taxonomy.tsv -outdir benchmark_select_region/test_dataset/mitogenomes_nonCOI/1_format_ncbi -min_length 700 -max_length 2000
~~~

#### Check_sequences 
- Add lineage
- Select sequences of length between 700 and 2000 bp
- Random sample to have as many sequences as in the positive dataset (30062)

~~~
perl benchmark_select_region/scripts/check_sequences.pl -tsv benchmark_select_region/test_dataset/mitogenomes_nonCOI/1_format_ncbi/ncbi_sequences.tsv -outdir benchmark_select_region/test_dataset/mitogenomes_nonCOI/2_check_sequences -taxonomy COInr/taxonomy.tsv -taxrank  phylum -min_length 700 -max_length 2000 -seq_n 30062
~~~

#### Move the cleaned file that will be used as a **negative-mito test dataset**
~~~
mv benchmark_select_region/test_dataset/mitogenomes_nonCOI/2_check_sequences/random_sample.tsv benchmark_select_region/test_dataset/negative_mito.tsv
~~~

### Negative chloroplast (non-COI sequences)

#### Select other then COI-like annotated genes from chloroplasts

~~~
perl benchmark_select_region/scripts/format_ncbi_negative.pl -cds benchmark_select_region/download_chloro_nsdpy/sequences.fasta -taxids benchmark_select_region/download_chloro_nsdpy/TaxIDs.txt -taxonomy COInr/taxonomy.tsv -outdir benchmark_select_region/test_dataset/chlorogenomes_nonCOI/1_format_ncbi -min_length 700 -max_length 2000
~~~

#### Check_sequences
- Add lineage
- Select sequence of length between 700 and 2000 bp
- Random sample to have as many sequences as in the positive dataset (30062)

~~~
perl benchmark_select_region/scripts/check_sequences.pl -tsv benchmark_select_region/test_dataset/chlorogenomes_nonCOI/1_format_ncbi/ncbi_sequences.tsv -outdir benchmark_select_region/test_dataset/chlorogenomes_nonCOI/2_check_sequences -taxonomy COInr/taxonomy.tsv -taxrank phylum -min_length 700 -max_length 2000 -seq_n 30062
~~~

#### Move the cleaned file that will be used as a **negative-chloro test dataset**
~~~
mv benchmark_select_region/test_dataset/chlorogenomes_nonCOI/2_check_sequences/random_sample.tsv benchmark_select_region/test_dataset/negative_chloro.tsv
~~~


## Trim sequences with select_region.pl with e_pcr option

Run select_region.pl with the following combination of parameters for each test datasets (positive-mito, negative-mito, negative-chloro).

The following parameters were constant:
e_pcr: 1
fw: TCHACNAAYCAYAARGAYATHGG
rv: TANACYTCNGGRTGNCCRAARAAYCA
min_amplicon_length: 600
max_amplicon_length: 700
bait_fas: 0
tcov: 1

The primers amplify the ca. 658 bp barcoding fragment of the COI gene.


| trim_error | min_overlap | identity |
| ---------- | ----------- | -------- |
| 0.3        | 10          | 0.7      |
| 0.3        | 10          | 0.6      |
| 0.3        | 20          | 0.7      |
| 0.3        | 20          | 0.6      |
| 0.2        | 10          | 0.7      |
| 0.2        | 10          | 0.6      |
| 0.2        | 20          | 0.7      |
| 0.2        | 20          | 0.6      |

Run select_region in a loop
~~~
cd ~/mkCOInr

perl benchmark_select_region/scripts/run_scripts.pl -param_file benchmark_select_region/run_select_region_E_pcr.tsv -script_name scripts/select_region.pl

~~~
Count trimmed and untrimmed sequences
~~~
perl benchmark_select_region/scripts/count_trimmed.pl  -dir /home/meglecz/mkCOInr/benchmark_select_region/test_dataset -motif _epcr_ -e_pcr 1
~~~



## Make bait files

Use the sequences of the positive test dataset trimmed by the following parameters of select_region.pl 

e_pcr: 1
fw: TCHACNAAYCAYAARGAYATHGG
rv: TANACYTCNGGRTGNCCRAARAAYCA
trim_error:0.3
min_amplicon_length: 600
max_amplicon_length: 700
min_overlap: 20
bait_fas: 0
tcov: 1
identity: 0.7

From the trimmed sequences, select 1 or 5 random sequences per
- order
- class
- phylum

Repeat each bait selection 10 times. That makes 3 x 2 x 10 = 60 bait files

~~~
perl benchmark_select_region/scripts/run_scripts.pl -param_file benchmark_select_region/run_make_bait.tsv -script_name benchmark_select_region/scripts/make_bait.pl
~~~


## Trim sequences with select_region.pl with bait_file option

- Use each of the 60 bait files to select the target region using the **bait_file option** instead of the E_pcr option of the **select_region.pl** script.
- Run select_region.pl on the **3 test datasets** (positive_mito, negative_mito, negative chloro) using each of the 60 bait files.
- Use 0.7 and 0.6 as identity threshold for usearch_global

That makes 60 x 3 x 2 = 360 runs:
~~~
cd ~/mkCOInr
perl benchmark_select_region/scripts/run_scripts.pl -param_file benchmark_select_region/run_select_region_bait.tsv -script_name scripts/select_region.pl
~~~

Count trimmed and untrimmed sequences
```
perl benchmark_select_region/scripts/count_trimmed.pl  -dir benchmark_select_region/positive_mito_bait_select_region/ -motif positive_mito_ -e_pcr 0

perl benchmark_select_region/scripts/count_trimmed.pl  -dir benchmark_select_region/negative_mito_bait_select_region/ -motif negative_mito_ -e_pcr 0

perl benchmark_select_region/scripts/count_trimmed.pl  -dir benchmark_select_region/negative_chloro_bait_select_region/ -motif negative_chloro_ -e_pcr 0

```

The number of False Positives (trimmed sequences in the negative dataset) is constant, irrespective of the bait file and sequence identity threshold.

The number of False Negatives (untrimmed sequences in the positive dataset) varies in function of the diversity in the bait files (1 or 5 sequences from each order, class or phylum), and in function of the identity threshold.

### Make bxplt of Sensitivity using different bait sets and identity thresholds

The proportion of trimmed sequences of the positive-mito dataset is the sensitivity of the select_region.pl, since we know that all sequences are complete COI sequences in this test dataset.

~~~
Rscript benchmark_select_region/scripts/bxplt.R
~~~

![Sensitivity](positive_mito_bait_select_region/Sensitivity_per_bait_file.png)



