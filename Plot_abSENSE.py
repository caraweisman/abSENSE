import numpy as np
import glob
import sys
import matplotlib.pyplot as plt
import random
from scipy.optimize import curve_fit
from scipy.stats import chi2
import inspect
import math
import re
import sys
import glob
import datetime
import os
import warnings
from cycler import cycler
from dill.source import getsource
from scipy import stats
import argparse
from datetime import datetime



###### Define user inputs #####

now = datetime.now()

starttime = now.strftime("%m.%d.%Y_%H.%M")

parser = argparse.ArgumentParser(description='abSENSE arguments:')

parser.add_argument("--distfile", type=str, required=True, help="Required. Name of file containing pairwise evolutionary distances between focal species and each of the other species")
parser.add_argument("--scorefile", type=str, required=True, help="Required. Name of file containing bitscores between focal species gene and orthologs in other species")
parser.add_argument("--gene", type=str, required=True, help="Required. Name of gene to be analyzed. Must match exactly the name of a gene in SCOREFILE (in first column).")
parser.add_argument("--Eval", default=0.001, type=float, help="Optional. E-value threshold. Scientific notation (e.g. 10E-5) accepted. Default 0.001.")
parser.add_argument("--includeonly", type=str, help="Optional. Species whose orthologs' bitscores will be included in fit; all others will be omitted. Default is all species. Format as species names, exactly as in input files, separated by commas (no spaces).")
parser.add_argument("--genelen", default=400, type=float, help="Optional. Length (aa) of the gene to be analyzed. Used to accurately calculate E-value threshold. Default is 400aa. Only large deviations will qualitatively affect results.")
parser.add_argument("--dblenfile", type=str, help="Optional. File containing size (aa) of databases on which the anticipated homology searches will be performed. Species-specific. Used to accurately calculate E-value threshold. Default is 400aa/gene * 20,000 genes for each species, intended to be the size of an average proteome. Only large deviations will significantly affect results.")
parser.add_argument("--predall", type=bool, default=False, help="Optional. True: Predicts bitscores and P(detectable) of homologs in all species, including those in which homologs were actually detected. Default is False: only make predictions for homologs that seem to be absent.")
now = datetime.now()
starttime = now.strftime("%m.%d.%Y_%H.%M")
args = parser.parse_args()

###### User input processing, general preprocessing, and checks #####

distancefilecheck = glob.glob(args.distfile)
if len(distancefilecheck) == 0:
        sys.exit('Distance file with that name not found! Is it in the current directory? If not, specify directory information. Quitting. \n')
else:
        distancefile = np.genfromtxt(args.distfile, dtype=str, delimiter='\t')

scorefilecheck = glob.glob(args.scorefile)
if len(scorefilecheck) == 0:
        sys.exit('Bitscore file with that name not found! Is it in the current directory? If not, specify directory information. Quitting. \n')
else:
        bitscores = np.genfromtxt(args.scorefile, dtype=str, delimiter='\t')
        
        
speciesorder = distancefile[0]
rawdistances = distancefile[1].astype(float)

genelist = bitscores[1:,0] # skip header

ethresh = args.Eval

seqlen = args.genelen

gene = args.gene

if args.dblenfile == None:
        speciesdblengthfilefound = False
        defdbsize = float(8000000)
        speciesdblengths = np.transpose(np.vstack((speciesorder, [float(defdbsize)]*len(speciesorder))))
else:
        speciesdblengthfilefound = True

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
            sys.exit('One or more species names in distance file do not match species names in distance file! The first I encountered was' + speciesorder[i] + '. Quitting. \n')

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


# function to, where possible, use maximum likelihood estimates of a and b parameter plus estimated covariance matrix to directly sample from the probability distribution of a and b (assume Gaussian with mean of max likelihood estimates and given covariance structure)
def parameter_CI_find(mla, mlb, covar):
        
        testavals = []
        testbvals = []

        if True not in np.isinf(covar):
        
                for i in range(0, 20): # 200
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

# function to take each of the sampled a, b values and use them to sample directly from the distribution of scores taking into account the Gaussian noise (a function of distance, a, b) 
# this gives an empirical estimate of the prediction interval 
def PI_find(testavals, testbvals, currx):

        # sample from score distribution: Gaussian with mean a, b and noise determined by distance (currx), a, b
        PIsamples = []
        for i in range(0, len(testavals)):
                detval = func(currx, testavals[i], testbvals[i])
                estnoise = np.sqrt(testavals[i]*(1-math.exp(-1*testbvals[i]*currx))*(math.exp(-1*testbvals[i]*currx)))
                if estnoise > 0:
                        parpairvals = []
                        for j in range(0, 20): # 200
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

        # compute fraction of sampled scores below threshold = P(undetected) = "p-value"
        pval = float(undetcount)/float(len(PIsamples))

        return mean - 3*std, mean + 3*std, pval



###### Main execution ######

# Nice terminal output formatting
print '\n'
print '\n'
print '\n'

print 'Running! (May take a moment...)'

# Ignore warning that sometimes happen as a result of stochastic sampling but that doesn't affect overall computation
warnings.filterwarnings("ignore", message="invalid value encountered in sqrt")

# make new arrays to put truncated (not using omitted species) scores, distances
genebitscores = []
truncdistances = []       

##print rawdistances
##speciesorderbydist = [x for _,x in sorted(zip(rawdistances,speciesorder))] # defined from bitscore fi
##rawdistancesbydist = sorted(rawdistances)
##print rawdistances

predbitscores = []
preddistances = []
allbitscores = []
alldistances = []
logscores = []

tocalcdists = []
tocalcspecs = []

ambigspecs = []
ambigdists = []

absentspecs = []
absentdists = []

for j in range(0, len(bitscores)): 
        if gene in bitscores[j][0]:
                scores = bitscores[j][1:]
                #orderedscores = [x for _,x in sorted(zip(rawdistances,scores))]
        	# put scores for current gene in bitscore file in right order
                orderedscores = []
                for k in range(0, len(ordervec)): # ordervec starts at 1 
                        orderedscores.append(bitscores[j][ordervec[k]]) ## # i + 1 because header skipped in gene list formation, so one behind now
                for k in range(0, len(orderedscores)):
                        tocalc = True
                        if isfloat(orderedscores[k]) == True and orderedscores[k] != '0': 
                                allbitscores.append(float(orderedscores[k]))
                                alldistances.append(rawdistances[k])
                                logscores.append(math.log(float(orderedscores[k])))
                                if len(pred_spec_locs) > 0:
                                        if k in pred_spec_locs:
                                                tocalc = False
                                                predbitscores.append(float(orderedscores[k]))
                                                preddistances.append(rawdistances[k])
                                else:
                                        predbitscores.append(float(orderedscores[k]))
                                        tocalc = False
                                        preddistances.append(rawdistances[k])
                        elif orderedscores[k] == 'N/A':
                                ambigspecs.append(speciesorder[k])
                                ambigdists.append(rawdistances[k])
                                tocalc = False
                        elif orderedscores[k] == '0':
                                absentspecs.append(speciesorder[k])
                                absentdists.append(rawdistances[k])                                       
                        if tocalc == True:
                                tocalcdists.append(rawdistances[k])
                                tocalcspecs.append(speciesorder[k])
if len(predbitscores) < 3:
        sys.exit('Too few species with detected homologs! There must be at least 3. Quitting. \n') 
                                
smoothx = np.linspace(min(rawdistances), max(rawdistances), 1000)
predictions = []
highs = []
lows = []

notinfitpvals = []

# define average database size for general line
avgdbsize = np.mean(speciesdblengths[:,1].astype(float))
globalbitthresh = -1*math.log(ethresh/(avgdbsize*seqlen), 2)
bitthresh = globalbitthresh

try: 
        (a, b), covar = curve_fit(func, preddistances, predbitscores)
        slope, intercept, r_value, p_value, std_err = stats.linregress(alldistances,logscores)
except RuntimeError:
        print 'Runtime Error'
parout = parameter_CI_find(a, b, covar) 
if parout != 'failed':
        testavals, testbvals = parout
        for j in range(0, len(smoothx)):
                predictions.append(round(func(smoothx[j], a,b),2))
                lowprediction, highprediction, pval = PI_find(testavals, testbvals, smoothx[j])
                highs.append(highprediction)
                lows.append(lowprediction)
        print '\n'
        print '\n'
        print '\n'
        print 'Results:'

        print '\n'
        print 'Best-fit a parameter from bitscores included in fit (black points): ', round(a,2)
        print 'Best-fit b parameter from bitscores included in fit (black points): ', round(b,2)
        print 'r squared from all bitscores (black points and orange points, if any):', round(r_value**2,2)
        print '\n'
        for j in range(0, len(tocalcdists)):
                for k in range(1, len(speciesdblengths)):
                        # now use species-specific db length
                        if tocalcspecs[j] in speciesdblengths[k][0]:
                                dblen = float(speciesdblengths[k][1])
                                bitthresh = -1*math.log(ethresh/(dblen*seqlen), 2)
                lowprediction, highprediction, pval = PI_find(testavals, testbvals, tocalcdists[j])
                print tocalcspecs[j], ':'
                print 'Maximum likelihood bitscore prediction: ', str(round(func(tocalcdists[j], a,b),2))
                print '99th percentile (high) prediction: ', str(round(highprediction, 2))
                print '1st percentile (low) prediction: ', str(round(lowprediction, 2))
                print 'Probability of homolog being undetected: ', pval
                print '\n'



startheight = 0.88

rawdistancesbydist = sorted(rawdistances)
speciesorderbydist = [x for _,x in sorted(zip(rawdistances,speciesorder))]

labels = speciesorderbydist
afont = {'fontname':'Arial'}
plt.title('Gene: ' + gene + '\n' + 'a = ' + str(round(a,1)) + ', b = ' + str(round(b,2)) + '\n' + '$r^2$ = ' + str(round(r_value**2, 2)), color='black', fontsize=13, fontweight='bold', **afont)
plt.scatter(preddistances, predbitscores, s=40, c='black', label='Bitscores of detected orthologs used in fit')
plt.plot(smoothx, predictions, c='red', label='Predicted bitscore')
plt.plot(smoothx, highs, c='black')
plt.plot(smoothx, lows, c='black')
plt.fill_between(smoothx, highs, lows, facecolor='blue', alpha=0.2, label='99% confidence interval')
plt.ylabel('Bitscore', fontsize=13, labelpad=10, **afont)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().tick_params(axis='x', width=2, length=7, direction='inout')
plt.xticks(rawdistancesbydist, labels, fontsize=10, rotation=90) #labels
totrange = max(smoothx) - min(smoothx)
inc = float(totrange)/100
for i in range(0, len(absentspecs)):
        if i == 0:
                plt.axvspan(absentdists[i] - inc, absentdists[i] + inc, facecolor='#fc8123', alpha=0.3, label='No homolog detected in species', capstyle='round')
                plt.gca().get_xticklabels()[labels.index(absentspecs[i])].set_color('#fc8123')
                plt.gca().get_xticklabels()[labels.index(absentspecs[i])].set_weight('bold')
                
        else:
                plt.axvspan(absentdists[i] - inc, absentdists[i] + inc, facecolor='#fc8123', alpha=0.3, capstyle='round')
                plt.gca().get_xticklabels()[labels.index(absentspecs[i])].set_color('#fc8123')
                plt.gca().get_xticklabels()[labels.index(absentspecs[i])].set_weight('bold')
for i in range(0, len(ambigspecs)):
        if i == 0:
                plt.gca().get_xticklabels()[labels.index(ambigspecs[i])].set_color('#a3a29b')
                
        else:
                plt.gca().get_xticklabels()[labels.index(ambigspecs[i])].set_color('#a3a29b')
plt.yticks(fontsize=10)
plt.axhline(y=globalbitthresh, linestyle='dashed', c='black', label='Detectability threshold')
plt.xlim([-inc, max(rawdistances)+inc])
plt.ylim([0, max(predictions)+max(predictions)/10])
handles, labels = plt.gca().get_legend_handles_labels()
if len(absentspecs) > 0: 
        order = [3, 0, 4, 1, 2]
else:
        order = [2,0, 3,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize=9)
plt.tight_layout()
plt.show()










