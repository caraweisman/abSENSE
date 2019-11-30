import numpy as np

inputfile = np.transpose(np.genfromtxt('Insect_Distances', dtype=str, delimiter='\t'))

outfile = open('newdists', 'w')
for i in range(0, len(inputfile)):
    outfile.write(inputfile[i][0])
    for j in range(1, len(inputfile[i])):
        outfile.write('\t')
        outfile.write(inputfile[i][j])
        outfile.write('\n')
outfile.close()
