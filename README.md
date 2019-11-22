# abSENSE: a method to interpret undetected homologs

__0. INTRODUCTION__ 

Welcome! This is the public Git repository for abSENSE, a method that calculates the probability that a homolog of a given gene would fail to be detected by a BLAST (or similar methods) homology search in a given species if the gene is evolving normally. 

The result of this calculation informs how one interprets homologs apparently being missing in some species. One possibility to explain such a result is that the gene is actually absent in that species, which can be biologically interesting (e.g. if due to a gene loss or the birth of a new gene). 

But another explanation, often not considered, is that the homolog /is/ present in that species, but that the given homology search method just lacks statistical power to detect it: that is, that the apparent absense of the homolog is due to a technical/statistical limitation, rather than representing true biology.

By calculating the probability that a homology search would fail to detect a homolog /even if one is present/ and /even if it is evolving normally/ (e.g. no rate accelerations on a specific branch, which could be suggestive of biologically interesting changes), abSENSE informs how seriously you take such an apparently absent homolog. If abSENSE calculates a high probability of being undetected, you may not be as inclined to invoke a biological explanation for the result: the null model of homology search failure is sufficient to explain what you observe.

The method is explained in further detail in the paper (citation). There, it is applied to the specific case of lineage-specific genes, for which homologs appear absent in all species outside of a narrow lineage. The method itself, however, is applicable to any case in which a homolog appears absent (e.g. a single species missing a homolog that one might interpret as a gene loss). This code is applicable to all such cases. 

__1. RUNNING abSENSE: THE BASICS__

The main analysis script, which will generate bitscore and detectability predictions for an arbitrary number of genes, is __Run\_abSENSE.py__. To predict the bitscores/E-values of a gene from some particular "focal" species in N other species, it needs at minimum two input files: 

i) A file containing the bitscores of homologs of the genes to be analyzed in at least three of the  species (including the focal species itself, so two others). 

ii) A file containing the N evolutionary distances, in substitutions/site, between the focal species and each other species. (The distance between the focal species and itself is 0.)

Examples of both of these files can be found for the fungal and insect lineages analyzed in (paper citation). The fungal files are Fungi_Data/Fungi_Bitscores and Fungi_Data/Fungi_Distances; the insect files are Insect_Data/Insect_Bitscores and Insect_Data/Insect_Distances. They exemplify the formatting required for abSENSE to run (explained in more detail below).

To run abSENSE on a given bitscore and distance file, simply type:

__python abSENSE_run.py --distfile (NAME OF DISTANCE FILE) --scorefile (NAME OF SCORE FILE)__
    
For example, to run abSENSE on the example fungal data (containing all S. cerevisiae genes in the current RefSeq annotation and their orthologs in 11 other fungal species), type: 

__python abSENSE_run.py --distfile Fungi_Data/Fungi_Distances --scorefile Fungi_Data/Fungi_Bitscores__

For each gene in the input bitscore file, the following will be computed:

a) The expected bitscore in each species in which no homolog was detected;

b) The 99\% confidence interval around this bitscore in each of these species;

c) The overall probability of a homolog being undetected in each of these species. 

These results will be output to a set of tab-delimited files in a separate directory (by default named with the start time of the analysis; you can specify the name with a command line option, see below).



There is also a supplemental visualization script, __Plot_abSENSE.py__, which performs the same analysis as above, but for a single gene at a time, and also produces a visualization of the results (see PAPER CITATION).
It is run in the same way as above, except that it also requires specifying which gene in the bitscore input file you wish to analyze. To run abSENSE on a given bitscore and distance file, on GENEID contained in the bitscore file, simply type:

__python abSENSE_run.py --distfile (NAME OF DISTANCE FILE) --scorefile (NAME OF SCORE FILE) --gene GENEID__

For example, to analyze the S. cerevisiae gene Uli1, listed in the bitscore file unde its RefSeq ID (NP_116682.3), type:

__python abSENSE_run.py --distfile Fungi_Data/Fungi_Distances --scorefile Fungi_Data/Fungi_Bitscores --gene NP_116682.3__

The same results as above will be computed, but here they will be output to the terminal, after which the visualiation will be shown.



__2.ADVANCED OPTIONS__

You can specify advanced options with additional command line options. You can view them all by typing

__python abSENSE_run.py --help__

They are:

__--out__: The prefix of the directory to which your results will be output. If not specified, defaults to the time at which you ran the analysis (to avoid overwriting of results). 

__--Eval__: The E-value threshold to be used (above this value, homologs will be considered undetected). If not specified, the default value of 0.001 (fairly permissive) will be used.

__--genelenfile__: Allows you to provide a file containing the lengths (in aa) of all genes in the bitscore file to be analyzed. abSENSE predicts a bitscore, which is then converted to an E-value to determine detectability; this conversation technically requires knowledge of both the size of the database in which the search occurs (see below) and the length of the gene being searched. Because the conversion between these values and E-value is logarithmic, though, only fairly large changes in these values substantially affect results. The default value is 400 amino acids (~average protein size in many species).

__--dblenfile__: Allows you to provide a file containing the sizes (in aa) of the database in which the homology search for each of your N species is performed. abSENSE predicts a bitscore, which is then converted to an E-value to determine detectability; this conversation technically requires knowledge of both the size of the database in which the search occurs and the length of the gene being searched (see above). Because the conversion between these values and E-value is logarithmic, though, only fairly large changes in these values substantially affect results. The default value is 400 amino acids * 20,000 amino acids / gene = 8,000,000 amimo acids (~average protein and proteome size in many species)..

__--predall__: Default is False. When True, Causes abSENSE to predict bitscores, confidence intervals, and probability of being undetected not only in species in which homologs are actually undetected, but also for species in which homologs are detected. This is obviously not the main use case, and is especially uninformative when those homologs and their bitscores have been used in the prediction itself (see below). May potentially be useful to see if a homolog in one species, although detected, seems to be behaving anomalously compared to those in other species (eg due to rate acceleration). 

__--includeonly__: Allows you to restrict the species whose bitscores are used in predicting bitscores in other species. Mainly useful to do control-type analyses, such as Figure 5 in (PAPER CITATION), to show that predictions made from only a subset of the data are nonetheless reliable. If not specified, abSENSE uses all available bitscores in the prediction. 



__N: FILE FORMAT__
