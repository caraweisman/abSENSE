# abSENSE: a method to interpret undetected homologs

- [abSENSE: a method to interpret undetected homologs](#absense--a-method-to-interpret-undetected-homologs)
  * [__INTRODUCTION__](#--introduction--)
  * [GUI WEBSITE AVAILABLE FOR SIMPLE ANALYSES](#gui-website-available-for-simple-analyses)
  * [1. RUNNING abSENSE: THE BASICS](#1-running-absense--the-basics)
    + [Quickstart: main analysis](#quickstart--main-analysis)
    + [Quickstart: visualization](#quickstart--visualization)
  * [2. ADVANCED OPTIONS](#2-advanced-options)
  * [3: OUTPUT FILES](#3--output-files)
    + [a) `Detection_failure_probabilities`](#a---detection-failure-probabilities-)
    + [b) `Predicted_bitscores`](#b---predicted-bitscores-)
    + [c) `Bitscore_99PI_upperbound_predictions`](#c---bitscore-99pi-upperbound-predictions-)
    + [d) `Bitscore_99PI_lowerbound_predictions`](#d---bitscore-99pi-lowerbound-predictions-)
    + [e) `Parameter_values`](#e---parameter-values-)
    + [f) `Run_info`](#f---run-info-)
  * [4: INPUT FILE FORMATS](#4--input-file-formats)
    + [a) The bitscore file](#a--the-bitscore-file)
    + [b) The distance file](#b--the-distance-file)
    + [c) The gene length file](#c--the-gene-length-file)
    + [d) The database length file](#d--the-database-length-file)
  * [5: DATA FOR FUNGAL AND INSECT GENES IN Weisman, C. M., Murray, A. W., & Eddy, S. R. (2020). Many but not all lineage-specific genes can be explained by homology detection failure, 12, 263–29. http://doi.org/10.1101/2020.02.27.968420](#5--data-for-fungal-and-insect-genes-in-weisman--c-m--murray--a-w----eddy--s-r--2020--many-but-not-all-lineage-specific-genes-can-be-explained-by-homology-detection-failure--12--263-29-http---doiorg-101101-20200227968420)
  * [6: QUESTIONS?](#6--questions-)

## __INTRODUCTION__ 

Welcome! This is the public Git repository for abSENSE, a method that calculates the probability that a homolog of a given gene would fail to be detected by a homology search (using BLAST or a similar method) in a given species, even if the homolog were present and evolving normally. 

The result of this calculation informs how one interprets the result of a homology search failing to find homologs of a gene in some species. One possibility to explain such a result is that the gene is _actually absent_ from the genome in that species: a biological, and potentially interesting (e.g. if due to a gene loss or the birth of a new gene), result. 

A second explanation, often ignored, is that the homolog _is_ present in the genome of that species, but that the homology search merely lacks statistical power to detect it. Here, the apparent absense of the homolog is a technical/statistical limitation, and does not reflect underlying biology.

By calculating the probability that your homology search would fail to detect a homolog _even if one were present_ and _even if it were evolving normally_ (e.g. no rate accelerations on a specific branch, potentially suggestive of biologically interesting changes), abSENSE informs the interpretation of a negative homology search result. If abSENSE finds that there is a high probability of a homolog being undetected even if present, you may not be as inclined to invoke a biological explanation for the result: the null model of a failure of the homology search is sufficient to explain what you observe.

The method is explained in further detail in the paper (https://www.biorxiv.org/content/10.1101/2020.02.27.968420v1). There, it is applied to the specific case of lineage-specific genes, for which homologs appear absent in all species outside of a narrow lineage. The method itself is applicable to any case in which a homolog appears absent (e.g. a single species missing a homolog that one might interpret as a gene loss), and likewise, this code is applicable to all such cases. 

In this repo, you will find:

__a)__ Code that can perform this analysis using the required input data (see below)

__b)__ All input data used in the fungal and insect lineages discussed in (PAPER CITATION) that can be used as input to the analysis code. 

## GUI WEBSITE AVAILABLE FOR SIMPLE ANALYSES

We have a GUI website that can perform basic single-gene analyses that don't require the advanced command line options available here (see below). If this is your case, you may find it easier to use than downloading and running this code on the command line.

The website is pre-loaded with the input data required to perform these analyses for all of the fungal and insect genes analyzed in (PAPER CITATION). It can also perform analyses on user-provided input data. 

The website is available at: http://eddylab.org/abSENSE/

## 1. RUNNING abSENSE: THE BASICS

### Quickstart: main analysis

The main analysis script, `Run_abSENSE.py`, calculates the probabilities of homologs of genes from some "focal" species being undetected in a set of N other species. It can perform this analysis for an arbitrary number of genes in the taxon at one time. It requires a minimum of two input files:

i) A file containing the bitscores of homologs of each gene to be analyzed in at least three of the species (including the focal species itself, so two others). 

ii) A file containing the N evolutionary distances, in substitutions/site, between the focal species and each other species. The distance between the focal species and itself should be 0. (If you don't already have such distances, a description of how to calculate them relatively painlessly can be found in PAPER CITATION.)

Examples of both of these files for a subset of genes from S. cerevisiae and their orthologs 11 other fungal species (the same species analyzed in CITATION) can be found in the folder Quickstart\_Examples: the bitscore file is Fungi\_Example\_Bitscores, and the distance file is Fungi\_Distances. They exemplify the formatting required for abSENSE to run (explained in more detail below).

To run abSENSE on a given bitscore and distance file, type:

```bash
python Run_abSENSE.py --distfile (NAME OF DISTANCE FILE) --scorefile (NAME OF SCORE FILE)
``` 

For example, to run abSENSE on the example fungal genes, type: 

```bash
python Run_abSENSE.py --distfile Quickstart_Examples/Fungi_Distances --scorefile Quickstart_Examples/Fungi_Example_Bitscores
```

For each gene in the input bitscore file, the following will be computed:

a) The probabilities of a homolog being undetected in each species (in the file `Detection_Failure_Probabilities`);

b) The expected bitscores of homologs in each species (in the file `Predicted_bitscores`);

c) The 99\% confidence interval around this bitscore in each species (low and high bounds listed in separate files: `Bitscore_99PI_lowerbound_predictions` and `Bitscore_99PI_higherbound_predictions`)

These results will be output to a set of tab-delimited files in a separate directory, by default named with the start time of the analysis. (You can specify the name with a command line option, see below). Additional information on output files is below.


### Quickstart: visualization

A supplemental visualization script, `Plot_abSENSE.py`, performs the same analysis as above, but for one gene at a time, and also produces a visualization of the results (see PAPER CITATION).
It is run in the same way, except that it also requires specifying which single gene in the bitscore input file you wish to analyze.

To run abSENSE on gene GENEID contained in a given bitscore with a given distance file, type:

```bash
python Plot_abSENSE.py --distfile (NAME OF DISTANCE FILE) --scorefile (NAME OF SCORE FILE) --gene GENEID
```

For example, to analyze the S. cerevisiae gene Uli1, listed in the bitscore file unde its RefSeq ID (NP_116682.3), type:

```bash
python Plot_abSENSE.py --distfile Quickstart\_Examples/Fungi_Distances --scorefile Quickstart\_Examples/Fungi\_Example\_Bitscores --gene NP\_116682.3
```

The same results as above will be computed, but now they will be output to the terminal, and then the visualiation will be shown.


## 2. ADVANCED OPTIONS

You can specify advanced options with the additional command line options. You can view them all by typing

```bash
python Run_abSENSE.py --help
```

They are:

- `--out`: The prefix of the directory to which your results will be output. Default is the time at which you ran the analysis (to avoid overwriting of results). 

- `--Eval`: The E-value threshold to be used (above this value, homologs will be considered undetected). Default is 0.001 (fairly permissive).

- `--genelenfile`: Allows you to specify a file containing the lengths (in aa) of all genes in the bitscore file to be analyzed. Default is 400 amino acids (~average protein size in many species) for all proteins.

abSENSE predicts a bitscore, which is then converted to an E-value to determine detectability; this conversation technically requires knowledge of both the size of the database in which the search occurs (see below) and the length of the gene being searched. Because the conversion between these values and E-value is logarithmic, though, only fairly large changes in these values substantially affect results. 

Examples of such files containing the lengths of all S. cerevisiae and D. melanogaster genes can be found in `Fungi_Data/S_cer_Protein_Lengths` and `Insect_Data/D_mel_Protein_Lengths`. The format required is described more below.

- `--dblenfile`: Allows you to specify a file containing the sizes (in aa) of the database in which the homology search for each of your N species is performed. Default is 400 amino acids * 20,000 amino acids / gene = 8,000,000 amimo acids (~average protein and proteome size in many species) for all species.

abSENSE predicts a bitscore, which is then converted to an E-value to determine detectability; this conversation technically requires knowledge of both the size of the database in which the search occurs and the length of the gene being searched (see above). Because the conversion between these values and E-value is logarithmic, though, only fairly large changes in these values substantially affect results. 

Examples containing the lengths of all S. cerevisiae and D. melanogaster genes can be found in `Fungi_Data/Fungi_Database_Lengths` and `Insect_Data/Insect_Database\_Lengths`. The format required is described more below.

- `--predall`: Default is False. When True, Causes abSENSE to calculate the probability of homologs being undetected, the expected bitscores, and 99\% confidence intervals not only in species in which homologs were actually undetected, but also for species in which homologs have been found. This is obviously not the main use case, and is especially uninformative when those homologs and their bitscores have been used in the prediction itself (see below). May potentially be useful to see if a homolog in one species, although detected, seems to be behaving anomalously compared to those in other species (eg due to rate acceleration). 

- `--includeonly`: Allows you to restrict the species whose bitscores are used in predicting bitscores in other species. Mainly useful to do control-type analyses, such as Figure 5 in (PAPER CITATION), to show that predictions made from only a subset of the data are nonetheless reliable. If not specified, abSENSE uses all available bitscores in the prediction. 


For example, to run an analysis on all S. cerevisiae proteins in the selected fungal species in which the lengths of each S. cerevisiae protein and the sizes of each species' search database (their annotated proteomes as indicated in the supplement of https://www.biorxiv.org/content/10.1101/2020.02.27.968420v1) are specified:

```
python Run_abSENSE.py --distfile Fungi_Data/Fungi_Distances --scorefile Fungi_Data/Fungi_Bitscores --genelenfile Fungi_Data/S_cer_Protein_Lengths --dblenfile Fungi_Data/Fungi_Database_Lengths
```

The visualization script `Plot_abSENSE.py` takes all of the same options, with the exception of again requiring that the gene to be analyzed is specified by the `--gene` option, and also that the `--genelenfile` option is instead `--genelen`, after which should be entered an integer corresponding to the length of the gene. (With a single gene, it's hardly worth requiring a whole file: just give the number.)


## 3: OUTPUT FILES

`Run_abSENSE.py` outputs six output files. Examples resulting from running abSENSE on the provided insect and fungal data can be found in `Fungi_Data/Fungi_abSENSE_Results/` and `Insect_Data/Insect_abSENSE_Results/` respectively.

### a) `Detection_failure_probabilities`

The central output of the program. For each gene in the analysis, this contains the predicted probability that a homolog in each species would be undetected at the specified E-value by a homology search, even if the homolog were present.

By default, this is only calculated in species in which the gene was not detected. Results for species in which homologs were detected are therefore listed as "detected". The setting --predall will calculate this value for all species, even those in which a homolog was in fact detected.

If not enough data for a gene was provided to generate a bitscore prediction (bitscores of homologs from at least three species are needed), the results will read `not_enough_data`.

### b) `Predicted_bitscores`

For each gene in the analysis, this contains the predicted (maximum likliehood) bitscore of a homolog in each species.

By default, bitscores are only predicted in species in which the gene was not detected. Results for species in which homologs were detected are therefore listed as "detected". The setting --predall will calculates this value for all species, even those for which a homolog was in fact detected. Here, the known bitscore (often used in the prediction process; see the option --includeonly) will be shown alongside the prediction. If the known bitscore was used in the fitting process, of course, these will usually be quite similar! 

If not enough data for a gene was provided to generate a bitscore prediction (bitscores of homologs from at least three species are needed), the results will read `not_enough_data`.

### c) `Bitscore_99PI_upperbound_predictions`

For each gene in the analysis, this contains the upper bound of the 99\% confidence interval for the bitscore of a homolog in each species.

By default, this is only calculated in species in which the gene was not detected. Results for species in which homologs were detected are therefore listed as "detected". The setting --predall will calculates this value for all species, even those for which a homolog was in fact detected. Here, the known bitscore (often used in the prediction process; see the option --includeonly) will be shown alongside the prediction. If the known bitscore was used in the fitting process, of course, these will usually be quite similar! 

If not enough data for a gene was provided to generate a bitscore prediction (bitscores of homologs from at least three species are needed), the results will read `not_enough_data`.

### d) `Bitscore_99PI_lowerbound_predictions`

For each gene in the analysis, this contains the lower bound of the 99\% confidence interval for the bitscore of a homolog in each species.

By default, this is only calculated in species in which the gene was not detected. Results for species in which homologs were detected are therefore listed as "detected". The setting `--predall` will calculates this value for all species, even those for which a homolog was in fact detected. Here, the known bitscore (often used in the prediction process; see the option `--includeonly`) will be shown alongside the prediction. If the known bitscore was used in the fitting process, of course, these will usually be quite similar! 

If not enough data for a gene was provided to generate a bitscore prediction (bitscores of homologs from at least three species are needed), the results will read `not_enough_data`.

### e) `Parameter_values`

For each gene in the analysis, this contains the best-fit (maximum likelihood) values of the a and b parameters. (See PAPER CITATION for explanation.) 

These a and b parameters are calculated from bitscores of homologs in species included in the prediction process. If the command line option `--includeonly` is used, this will be only the species specified by that option. By default, all provided bitscores are used.

If not enough data for a gene was provided to generate a bitscore prediction (bitscores of homologs from at least three species are needed), the results will read `not_enough_data`.

### f) `Run_info`

Contains information about the analysis, including names of input files, options/settings used, and the analysis time.

## 4: INPUT FILE FORMATS

Required files:

### a) The bitscore file

For an analysis of M genes in N species (including the focal species), the bitscore file should be a tab-delimited file of N+1 columns by M+1 rows.
The first row should begin with a blank entry, and should be followed by N entries containing the names of the N species in your analysis. These names should match those in the distance file (below) exactly.
The remaining M rows should each begin with the name/identifier of the gene from the focal species to be analyzed, followed by the bitscore of that gene against its homolog in the species indicated at the top of the given column. For species in which homologs are undetected, this value should be __0__. For species in which a homolog is detected, but the orthology is unclear and so you wish to exclude it from being used in the fit (see PAPER CITATION), this value should be __'N/A'__.
Two examples are provided: Fungi\_Data/Fungi_\Bitscores and Insect_\Data/Insect\_Bitscores.

### b) The distance file

For an analysis with N species (including the focal species), the distance file should be a tab-delimited file of 2 columns by N rows.
Entries in the first column should contain the name of each species in your analysis. These names should match those in the bitscore file (above) exactly.
Entries in the second column should contain the evolutionary distance between each species in the indicated column and the focal species. (The distance between the focal species and itself should always be 0.)
Two examples are provided: Fungi\_Data/Fungi\_Distances and Insect\_Data/Insect\_Distances.

Optional files:

### c) The gene length file

For an analysis of M genes, the gene length file should be a tab-delimited file of 2 columns by M rows.
The first column should contain the names/identifiers of the gene from the focal species to be analyzed. These should match exactly the names in the bitscore file (above).
The second column should contain the length in amino acids of that gene.

If you don't already have such a file, here is a command to make one (named OUTPUTFILENAME), from a FASTA file (SEQFILENAME) containing all of the sequences: (It requires the easel package, which comes with the HMMER software, available at http://hmmer.org/documentation.html.)

```
esl-seqstat -a (SEQFILENAME) | awk '{print $2 "\t" $3}' | tac | sed -e '1,7d' | tac > (OUTFILENAME)
```

### d) The database length file

For an analysis of N species, the database length file should be a tab-delimited file of 2 columns by N rows.
The first column should contain the names of the N species. These should match exactly the names in the bitscore and distance files (above).
The second column should contain the sizes, in amino acids, of the database in which the homology search for each species is performed. For example, if you are searching for a gene in each of the species' annotated proteomes, it should be the size of that species' proteome in amino acids. If instead you are searching against a pan-specific database, like for example NR, for all species, it should be the size of that database in amino acids.

If you don't already have such a file, you can make one easily if you have the databases themselves in eg FASTA format (again requiring the easel package, downloadable with HMMER as in c) above:
just run the command
__esl-seqstat (FASTA FILE)__
on each database; this will report the total length in aa of each database file. You can then put these into a tab-delimited file manually.

## 5: DATA FOR FUNGAL AND INSECT GENES IN Weisman, C. M., Murray, A. W., & Eddy, S. R. (2020). Many but not all lineage-specific genes can be explained by homology detection failure, 12, 263–29. http://doi.org/10.1101/2020.02.27.968420

Also in this folder are all input data required to perform abSENSE analyses on the fungal and insect lineages analyzed in (paper citation). The fungal files are Fungi_Data/Fungi_Bitscores and Fungi_Data/Fungi_Distances; the insect files are Insect_Data/Insect_Bitscores and Insect_Data/Insect_Distances. 

As discussed in (CITATION), the bitscore files contain the bitscores of orthologs of all S. cerevisiae genes in the current RefSeq annotation in 11 other fungal species. The insect bitscore files contain the bitscores of orthologs of all D. melanogaster genes in the current RefSeq annotation in 21 other fungal species. 


## 6: QUESTIONS?

If you have questions, concerns, or suspect a bug, please contact me at weisman@g.harvard.edu. 

