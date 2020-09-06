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
    badguy = string.split(dirList[i], 'PW305')
    if (len(badguy) == 1):
        dirList.pop(i)
        i = i-1
    i += 1
    if i == len(dirList):
        break 

dirList.sort(key=lambda x: float(x.split('_')[1][0:]))

Emis = numpy.empty([len(dirList)+1,3],dtype=object)
Emis[0][0] = 'Condition'
Emis[0][1] = 'NOx'
Emis[0][2] = 'CO'
LoopCounter = 1
for fname in dirList:    
    p = string.rsplit(fname, sep='_')
    if (len(p) == 5):
        with open(fname, 'r+') as f:
            data = f.readlines()
            
        data1 = string.split(data[len(data)-2], ',')
        Emis[LoopCounter][0] = float(p[1])
        Emis[LoopCounter][1] = data1[2]
        Emis[LoopCounter][2] = data1[3] 
        LoopCounter += 1

## Writing Emis vs. SMD data to a csv file
s = 'Emissions.csv'
csv_file = open(str(s),'wb')
csv_writer = csv.writer(csv_file)
csv_writer.writerows(Emis)
csv_file.close()


