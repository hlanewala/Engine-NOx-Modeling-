import subprocess
import shutil
import csv
import os
import string
from numpy import numarray
import numpy
 
## open the path       
path = os.getcwd()
os.chdir(path)

## search for names of the files of interest
dirList=os.listdir(path)
i = 0
while 1>0:
    badguy = string.split(dirList[i], 'JT15D_Idle')
    if (len(badguy) == 1):
        dirList.pop(i)
        i = i-1
    i += 1
    if i == len(dirList):
        break 

dirList.sort(key=lambda x: float(x.split('_')[2][0:]))

NOx = numpy.empty([len(dirList)+1,2],dtype=object)
NOx[0][0] = 'SMD'
NOx[0][1] = 'Stoich Mass'
LoopCounter = 1
for fname in dirList:
    p = string.rsplit(fname, sep='_')    
    if (len(p) == 4):
        with open(fname, 'r+') as f:
            data = f.readlines()
            
        data = string.split(data[4], ',')
        NOx[LoopCounter][0] = float(p[2])
        NOx[LoopCounter][1] = data[7]        
        LoopCounter += 1

## Writing NOx vs. SMD data to a csv file
s = 'Stoich Mass vs. SMD.csv'
csv_file = open(str(s),'wb')
csv_writer = csv.writer(csv_file)
csv_writer.writerows(NOx)
csv_file.close()


