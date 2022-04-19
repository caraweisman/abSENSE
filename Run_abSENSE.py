from __future__ import print_function
import numpy as np
import glob
import sys
import random
from scipy.optimize import curve_fit
from scipy.stats import chi2
import inspect
import math
import re
import sys
import glob
import os
import warnings
import argparse
from datetime import datetime
from scipy import stats

#!/usr/bin/env python

###### Define user inputs #####

now = datetime.now()

starttime = now.strftime("%m.%d.%Y_%H.%M")

parser = argparse.ArgumentParser(description='abSENSE arguments:')

parser.add_argument("--distfile", type=str, required=True, help="Required. Name of file containing pairwise evolutionary distances between focal species and each of the other species")
parser.add_argument("--scorefile", type=str, required=True, help="Required. Name of file containing bitscores between focal species gene and orthologs in other species")
parser.add_argument("--Eval", default=0.001, type=float, help="Optional. E-value threshold. Scientific notation (e.g. 10E-5) accepted. Default 0.001.")
parser.add_argument("--includeonly", type=str, help="Optional. Species whose orthologs' bitscores will be included in fit; all others will be omitted. Default is all species. Format as species names, exactly as in input files, separated by commas (no spaces).")
parser.add_argument("--genelenfile", type=str, help="Optional. File containing lengths (aa) of all genes to be analyzed. Used to accurately calculate E-value threshold. Default is 400aa for all genes. Only large deviations will qualitatively affect results.")
parser.add_argument("--dblenfile", type=str, help="Optional. File containing size (aa) of databases on which the anticipated homology searches will be performed. Species-specific. Used to accurately calculate E-value threshold. Default is 400aa/gene * 20,000 genes for each species, intended to be the size of an average proteome. Only large deviations will significantly affect results.")
parser.add_argument("--predall", type=bool, default=False, help="Optional. True: Predicts bitscores and P(detectable) of homologs in all species, including those in which homologs were actually detected. Default is False: only make predictions for homologs that seem to be absent.")
now = datetime.now()
starttime = now.strftime("%m.%d.%Y_%H.%M")
parser.add_argument("--out", type=str, default='abSENSE_results_' + starttime, help="Optional. Name of directory for output data. Default is date and time when analysis was run.")
args = parser.parse_args()

###### User input processing, general preprocessing, and checks #####

distancefilecheck = glob.glob(args.distfile)
if len(distancefilecheck) == 0:
        sys.exit('Distance file with that name not found! Is it in the current directory? If not, specify directory information. Quitting. \n')
else:
        distancefile = np.transpose(np.genfromtxt(args.distfile, dtype=str, delimiter='\t'))

scorefilecheck = glob.glob(args.scorefile)
if len(scorefilecheck) == 0:
        sys.exit('Bitscore file with that name not found! Is it in the current directory? If not, specify directory information. Quitting. \n')
else:
        bitscores = np.genfromtxt(args.scorefile, dtype=str, delimiter='\t')
        
        
speciesorder = distancefile[0]
rawdistances = distancefile[1].astype(float)
genelist = bitscores[1:,0] # skip header

ethresh = args.Eval


if args.genelenfile == None:
        genelengthfilefound = False
        defgenelen = float(400)
else:
        genelengthfilefound = True
        
if args.dblenfile == None:
        speciesdblengthfilefound = False
        defdbsize = float(8000000)
        speciesdblengths = np.transpose(np.vstack((speciesorder, [float(defdbsize)]*len(speciesorder))))
else:
        speciesdblengthfilefound = True


if genelengthfilefound ==  True:
        genelengthfilecheck = glob.glob(args.genelenfile)
        genelengthfilename = args.genelenfile
        if len(genelengthfilecheck) != 1:
                sys.exit('Gene length file with that name not found! Is it in the current directory? If not, specify directory information. Quitting. \n')
        else:
                genelengths = np.genfromtxt(args.genelenfile, dtype = str, delimiter='\t')
                

if speciesdblengthfilefound == True:
        speciesdblengthfilecheck = glob.glob(args.dblenfile)
        speciesdblengthfilename = args.dblenfile
        if len(speciesdblengthfilecheck) != 1:
                sys.exit('Species database size file with that name not found! Is it in the current directory? If not, specify directory information. Quitting. \n')
        else:
                speciesdblengths = np.genfromtxt(args.dblenfile, dtype = str, delimiter='\t')
        
## Determine order of species in the bitscore file, in case they're not the same as in the distance file, from which they're extracted
## Also determine the locations in the ordering to be used of the species to omit from curve fit, if any
ordervec = [] # first index is location of first species in speciesorder in the file; and so on
if args.includeonly != None:
        pred_specs = re.split(',', args.includeonly)
        pred_spec_locs = []
else:
        pred_specs = []
        pred_spec_locs = []
for i in range(0, len(speciesorder)):
    found = False
    for j in range(0, len(bitscores[0])):
        if speciesorder[i] in bitscores[0][j]:
            ordervec.append(j)
            found = True
            if speciesorder[i] in pred_specs:
                pred_spec_locs.append(i)
    if found == False:
            sys.exit('One or more species names in header of bitscore file do not match species names in header of distance file! The first I encountered was ' + speciesorder[i] + '. Quitting. \n')

invordervec = [] # first index is location of first species in file in speciesorder; and so on
for i in range(0, len(bitscores[0])):
        for j in range (0, len(speciesorder)):
                if bitscores[0][i] in speciesorder[j]:
                        invordervec.append(j)


## Find species db sizes in the right order, from either file or manual input (done above)
speciestotallengths = []
for i in range(0, len(speciesorder)):
        found = False
        for j in range(0, len(speciesdblengths)): 
                if speciesorder[i] in speciesdblengths[j][0]:
                        speciestotallengths.append(float(speciesdblengths[j][1]))
                        found = True
        if found == False:
                if speciesdblengthfilefound == True:
                        sys.exit('One or more species names in your database size file do not match species names in distance file! The first I encountered was ' + speciesorder[i] + '. Quitting. \n')



###### Define functions ######

## curve to fit
def func(x, a, b):
        return a*np.exp(-b*x)

def isfloat(s):
        try:
                float(s)
                return True
        except ValueError:
                return False


## function to, where possible, use maximum likelihood estimates of a and b parameter plus estimated covariance matrix to directly sample from the probability distribution of a and b (assume Gaussian with mean of max likelihood estimates and given covariance structure)
def parameter_CI_find(mla, mlb, covar):
        
        testavals = []
        testbvals = []

        if True not in np.isinf(covar):
        
                for i in range(0, 200):
                        a = np.random.multivariate_normal([mla, mlb], covar)[0]
                        b = np.random.multivariate_normal([mla, mlb], covar)[1]
                        testavals.append(a)
                        testbvals.append(b)

                if len(testavals) > 0:
                        return testavals, testbvals 
                else:
                        return 'failed'
        else:
                return 'failed'

## function to take each of the sampled a, b values and use them to sample directly from the distribution of scores taking into account the Gaussian noise (a function of distance, a, b) 
## this gives an empirical estimate of the prediction interval 
def PI_find(testavals, testbvals, currx):

        # sample from score distribution: Gaussian with mean a, b and noise determined by distance (currx), a, b
        PIsamples = []
        for i in range(0, len(testavals)):
                detval = func(currx, testavals[i], testbvals[i])
                estnoise = np.sqrt(testavals[i]*(1-math.exp(-1*testbvals[i]*currx))*(math.exp(-1*testbvals[i]*currx)))
                if estnoise > 0:
                        parpairvals = []
                        for j in range(0, 200):
                                PIsamples.append(detval + np.random.normal(0, estnoise))
                                parpairvals.append(detval + np.random.normal(0, estnoise))
                else:
                        PIsamples.append(detval)
        # compute mean of sample
        mean = np.mean(PIsamples)
        # compute std dev of sample
        std = np.std(PIsamples)

        # empirically determine, from sampled scores, how many are below detectability threshold 
        undetcount = 0
        for i in range(0, len(PIsamples)):
                if PIsamples[i] < bitthresh: 
                        undetcount = undetcount + 1

        # compute fraction of sampled scores below threshold = P(undetected) = empriical "p-value"
        emppval = float(undetcount)/float(len(PIsamples))

        # calculate this analytically from std estimate
        pval = stats.norm.cdf(bitthresh, mean, std)
        
        # calculate 99% CI 
        (lowint, highint) = stats.norm.interval(0.99, mean, std)

        return lowint, highint, pval
        
###### Output preparation ######

# make directory named for control file into which all output files will be put
outdirectory = args.out
os.system('mkdir ' + outdirectory)
                
mloutputfile = open(outdirectory + '/' + 'Predicted_bitscores', 'w')
lowboundoutputfile = open(outdirectory + '/' + 'Bitscore_99PI_lowerbound_predictions', 'w')
highboundoutputfile = open(outdirectory + '/' + 'Bitscore_99PI_higherbound_predictions', 'w')
pvaloutputfile = open(outdirectory + '/' + 'Detection_failure_probabilities', 'w')
outputfileparams = open(outdirectory + '/' + 'Parameter_values', 'w')

runinfofile = open(outdirectory + '/' + 'Run_info', 'w')

# Write brief summary header to output files

mloutputfile.write('# This file contains maximum likelihood bitscore predictions for each tested gene in each species' + '\n')
lowboundoutputfile.write('# This file contains the lower bound of the 99% bitscore prediction interval for each tested gene in each species' + '\n')
highboundoutputfile.write('# This file contains the upper bound of the 99% bitscore prediction interval for each tested gene in each species' + '\n')
pvaloutputfile.write('# This file contains the probability of a homolog being undetected at the specified significance threshold (see run info file) in each tested gene in each species' + '\n')
outputfileparams.write('# This file contains the best-fit parameters (performed using only bitscores from species not omitted from the fit; see run info file) for a and b for each gene' + '\n')

# Write run information to run info output file

runinfofile.write('abSENSE analysis run on ' + starttime + '\n')
runinfofile.write('Input bitscore file: ' + args.scorefile + '\n')
runinfofile.write('Input distance file: ' + args.distfile + '\n')
if genelengthfilefound ==  True:
        runinfofile.write('Gene length file: ' + args.genelenfile + '\n')
elif genelengthfilefound ==  False:
        runinfofile.write('Gene length (for all genes): ' + str(defgenelen) + ' (default)' + '\n')
if speciesdblengthfilefound == True:
        runinfofile.write('Database length file: ' + args.dblenfile + '\n')
elif speciesdblengthfilefound == False:
        runinfofile.write('Database length (for all species): ' + str(defdbsize) + ' (default)' + '\n')
runinfofile.write('Species used in fit: ')
for i in range(0, len(pred_specs)):
        runinfofile.write(pred_specs[i] + ' ')
if len(pred_specs) == 0:
        runinfofile.write('All (default)')
runinfofile.write('\n')
runinfofile.write('E-value threshold: ' + str(ethresh) + ' (default)' + '\n')
runinfofile.close()

# Ignore warning that sometimes happen as a result of stochastic sampling but that doesn't affect overall computation
warnings.filterwarnings("ignore", message="invalid value encountered in sqrt")


###### Main execution and output ######

# Write headers to file (first column will have the gene name; subsequent columns have prediction info)
mloutputfile.write('Gene')
lowboundoutputfile.write('Gene')
highboundoutputfile.write('Gene')
pvaloutputfile.write('Gene')
outputfileparams.write('Gene')
for i in range(0, len(speciesorder)):
        mloutputfile.write('\t' + speciesorder[invordervec[i]])
        lowboundoutputfile.write('\t' + speciesorder[invordervec[i]])
        highboundoutputfile.write('\t' + speciesorder[invordervec[i]])
        pvaloutputfile.write('\t' + speciesorder[invordervec[i]])
mloutputfile.write('\n')
lowboundoutputfile.write('\n')
highboundoutputfile.write('\n')
pvaloutputfile.write('\n')
outputfileparams.write('\t' + 'a' + '\t' + 'b' + '\n')


print('Running!')

for i in range(0, len(genelist)):
        # report current position and gene name
        print('gene', i, 'out of', str(len(genelist)), ':', genelist[i])
        
        # print current gene to output file
        mloutputfile.write(genelist[i])
        lowboundoutputfile.write(genelist[i])
        highboundoutputfile.write(genelist[i])
        outputfileparams.write(genelist[i])
        pvaloutputfile.write(genelist[i])

        # make new arrays to put truncated (not using omitted species) scores, distances
        genebitscores = []
        truncdistances = []

        # make new arrays to indicate which species/distances are ambiguous orthology but have a homolog; don't predict these later
        ambigdists =  [] 
        
        # if gene length input file given, look for length of current gene
        # if not given, assume default value (specified above)
        lengthfound = False
        if genelengthfilefound ==  True:
                for z in range(1, len(genelengths)): # skip header
                        if genelist[i] in genelengths[z]:
                                seqlen = float(genelengths[z][1])
                                lengthfound = True
                                break
                if lengthfound == False:
                        sys.exit('Gene ' + genelist[i] + ' not found in specified gene length file! Quitting \n')
        elif genelengthfilefound ==  False:
                seqlen = float(defgenelen)

        # put scores for current gene in bitscore file in right order
        orderedscores = []
        for k in range(0, len(ordervec)): # ordervec starts at 1 
                orderedscores.append(bitscores[i+1][ordervec[k]]) ## # i + 1 because header skipped in gene list formation, so one behind now
        # append score of species and corresponding distance of species to gene-specific distance, score vectors if:
        # score isn't 0 (can't distinguish absence from loss from bad data etc)
        # score isn't 'N/A' or some other string (ie indicator of unclear orthology or otherwise absent, or generally unclear what it is)
        # score isn't from species that is to be excluded from fit
        for k in range(0, len(orderedscores)):
                if len(pred_spec_locs) > 0:
                        if isfloat(orderedscores[k]) == True and orderedscores[k] != '0' and k in pred_spec_locs:
                                genebitscores.append(float(orderedscores[k]))
                                truncdistances.append(rawdistances[k])
                        elif orderedscores[k] == 'N/A':
                                ambigdists.append(rawdistances[k])
                else:
                        if isfloat(orderedscores[k]) == True and orderedscores[k] != '0':
                                genebitscores.append(float(orderedscores[k]))
                                truncdistances.append(rawdistances[k])
                        elif orderedscores[k] == 'N/A':
                                ambigdists.append(rawdistances[k])
        
        if len(truncdistances) > 2: 
                try: 
                        (a, b), covar = curve_fit(func, truncdistances, genebitscores, bounds=((-np.inf, 0), (np.inf, np.inf)))
                except RuntimeError:
                        for j in range(0, len(rawdistances)):
                                mloutputfile.write('\t' + 'analysis_error')
                                highboundoutputfile.write('\t' + 'analysis_error')
                                lowboundoutputfile.write('\t' + 'analysis_error')
                                pvaloutputfile.write('\t' + 'analysis_error')
                        mloutputfile.write('\n')
                        highboundoutputfile.write('\n')
                        lowboundoutputfile.write('\n')
                        outputfileparams.write('\t' + 'analysis_error' + '\t' + 'analysis_error' + '\n')
                        pvaloutputfile.write('\n')
                        continue
                parout = parameter_CI_find(a, b, covar) 
                if parout != 'failed':
                        testavals, testbvals = parout
                        for j in range(0, len(rawdistances)):
                                bitthresh = -1*math.log(ethresh/(seqlen*speciestotallengths[invordervec[j]]), 2)
                                prediction = round(func(rawdistances[invordervec[j]], a,b),2)
                                lowprediction, highprediction, pval = PI_find(testavals, testbvals, rawdistances[invordervec[j]])
                                if rawdistances[invordervec[j]] not in truncdistances and rawdistances[invordervec[j]] not in ambigdists:
                                        mloutputfile.write('\t' + str(prediction))
                                        highboundoutputfile.write('\t' + str(round(highprediction,2)))
                                        lowboundoutputfile.write('\t' + str(round(lowprediction,2)))
                                        pvaloutputfile.write('\t' + str(round(pval,2)))
                                elif rawdistances[invordervec[j]] in truncdistances: #  could make new vector, ambigdists, and test  whether dist is in it herre; if so, don't predict
                                        if args.predall == True:
                                                realscore = genebitscores[truncdistances.index(rawdistances[invordervec[j]])]
                                                mloutputfile.write('\t' + str(prediction) + '(Ortholog_detected:' + str(realscore) + ')')
                                                highboundoutputfile.write('\t' + str(round(highprediction,2)) + '(Ortholog_detected)')
                                                lowboundoutputfile.write('\t' + str(round(lowprediction,2)) + '(Ortholog_detected)')
                                                pvaloutputfile.write('\t' + str(round(pval,2)) + '(Ortholog_detected)')
                                        elif args.predall == False:
                                                mloutputfile.write('\t' + 'Ortholog_detected')
                                                highboundoutputfile.write('\t' + 'Ortholog_detected')
                                                lowboundoutputfile.write('\t' + 'Ortholog_detected')
                                                pvaloutputfile.write('\t' + 'Ortholog_detected')
                                elif rawdistances[invordervec[j]] in ambigdists:
                                        if args.predall == True:
                                                mloutputfile.write('\t' + str(prediction) + 'Homolog_detected(orthology_ambiguous)')
                                                highboundoutputfile.write('\t' + str(round(highprediction,2)) + 'Homolog_detected(orthology_ambiguous)')
                                                lowboundoutputfile.write('\t' + str(round(lowprediction,2)) + 'Homolog_detected(orthology_ambiguous)')
                                                pvaloutputfile.write('\t' + str(round(pval,2)) + 'Homolog_detected(orthology_ambiguous)')
                                        elif args.predall == False:
                                                mloutputfile.write('\t' + 'Homolog_detected(orthology_ambiguous)')
                                                highboundoutputfile.write('\t' + 'Homolog_detected(orthology_ambiguous)')
                                                lowboundoutputfile.write('\t' + 'Homolog_detected(orthology_ambiguous)')
                                                pvaloutputfile.write('\t' + 'Homolog_detected(orthology_ambiguous)')
                                
                        mloutputfile.write('\n')
                        highboundoutputfile.write('\n')
                        lowboundoutputfile.write('\n')
                        pvaloutputfile.write('\n')
                        outputfileparams.write('\t' + str(a) + '\t' +  str(b))
                        outputfileparams.write('\n')
                else:
                        for j in range(0, len(rawdistances)):
                                prediction = round(func(rawdistances[j], a,b),2) 
                                mloutputfile.write('\t' + str(prediction))
                                highboundoutputfile.write('\t' + 'analysis_error')
                                lowboundoutputfile.write('\t' + 'analysis_error')
                                pvaloutputfile.write('\t' + 'analysis_error')
                        mloutputfile.write('\n')
                        highboundoutputfile.write('\n')
                        outputfileparams.write('\t' + 'analysis_error' + '\t' + 'analysis_error' + '\n')
                        lowboundoutputfile.write('\n')
                        pvaloutputfile.write('\n')
        else:
                for j in range(0, len(rawdistances)):
                        mloutputfile.write('\t' + 'not_enough_data')
                        highboundoutputfile.write('\t' + 'not_enough_data')
                        lowboundoutputfile.write('\t' + 'not_enough_data')
                        pvaloutputfile.write('\t' + 'not_enough_data')
                mloutputfile.write('\n')
                highboundoutputfile.write('\n')
                lowboundoutputfile.write('\n')
                pvaloutputfile.write('\n')
                outputfileparams.write('\t' + 'not_enough_data' + '\t' + 'not_enough_data' + '\n')
mloutputfile.close()
highboundoutputfile.close()
lowboundoutputfile.close()
outputfileparams.close()
pvaloutputfile.close()

