'''
Date: 1/21/22
Author: Andrew Martin
Title: Sequencing Exercise - Base calling
Python: 3.6.5
Description: parses two biochem methods signal intensity data to determine basecalls
for each cycle. Records relevant quality analytics and outputs to a text file.
'''
import csv
import math
from decimal import *

#Results class with four methods
class Results():
    def __init__(self):
        self.basecount = {1:{0:0,1:0,2:0,3:0,4:0,5:0},2:{0:0,1:0,2:0,3:0,4:0,5:0}}
        self.call = {1:[],2:[]}
        self.intensitysum = {1:{1:0,2:0,3:0,4:0},2:{1:0,2:0,3:0,4:0}}
        self.probecount = 0
        self.quality = {1:{'perror':0,'phred':0},2:{'perror':0,'phred':0}}
        self.avgcall = {}
        
    #convert to base or no-read
    def Basecall(self, base):
        basecall = {0:'N',1:'A',2:'C',3:'G',4:'T'}
        return basecall[base]
    
    #extract base with highest signal intensity and assign the corresponding base/no-read
    def getResults(self, cycle):       
        
        maxrow = {1:max(cycle[1][1:]),2:max(cycle[2][1:])}
        base = {1:'',2:''}
        for i in base:
            base[i] = cycle[i].index(maxrow[i])
        for i in cycle:
                
            if cycle[i][1:] == ['0.0','0.0','0.0','0.0']:
                self.basecount[i][0] += 1
                base[i] = 0
            else:
                self.intensitysum[i][base[i]] += Decimal(maxrow[i])
                self.basecount[i][base[i]] += 1
                if cycle[i][0] != self.Basecall(base[i]):
                    self.basecount[i][5] += 1
            self.call[i].append(self.Basecall(base[i]))
            

        self.probecount += 1
    
    #calculate percent error and quality scores along with average signal intensities for each channel
    def calculate(self):
        for i in self.quality:
            self.quality[i]['perror'] = ((self.basecount[i][0]+self.basecount[i][5])/self.probecount)*100
            self.quality[i]['phred'] = -10*(math.log10(self.quality[i]['perror']/100))
        
        #take sums of each called base (j) per channel (i) per cycle (k) and divide by the number of times that channel was called
        self.avgcall = {k:[j/self.basecount[k][i] for (i,j) in v.items()] for (k,v) in self.intensitysum.items()}
    
    #write results to a text file
    def printResults(self, resultsFile):
        with open(resultsFile,'w') as w:
            for i in self.call:
                w.write('CYCLE: '+str(i)+'\n')
                w.writelines(self.call[i])
                w.write('\n')
                w.write('number of reads: '+str(self.probecount-self.basecount[i][0])+'\n')
                w.write('number of no-reads:'+str(self.basecount[i][0])+'\n')
                w.write('number of incorrect reads:'+str(self.basecount[i][5])+'\n')
                w.write('percent error:'+str(self.quality[i]['perror'])+ '%\n')
                w.write('quality score:'+str(self.quality[i]['phred'])+'\n')
                w.write('avg call intensity:'+str((self.avgcall[i][0]+self.avgcall[i][1]+self.avgcall[i][2]+self.avgcall[i][3])/4)+'\n')
                w.write('gc content:'+str((self.basecount[i][2]+self.basecount[i][3])/self.probecount*100)+ '%\n\n')

#main function
def main():
    with open('sequencing_data_biochem1.csv', newline='') as f:
        results = Results()
        cycle = dict()
        reader = csv.reader(f)
        header = next(reader)
        
        #split into cycles one and two as [ref,A,C,G,T]
        for row in reader:
            cycle[1] = row[0:5]
            cycle[2] = row[6:11]
            results.getResults(cycle)
        results.calculate()  
        results.printResults('results1.txt')

    with open('sequencing_data_biochem2.csv', newline='') as a:
        results = Results()
        cycle = dict()
        reader = csv.reader(a)
        header = next(reader)
        
        #split and apply correction of rotating one channel to the right
        for row in reader:
            cycle[1] = row[0:5]
            cycle[2] = row[6:11]
            for i in cycle:
                cycle[i][1],cycle[i][2],cycle[i][3],cycle[i][4]=cycle[i][4],cycle[i][1],cycle[i][2],cycle[i][3]
            results.getResults(cycle)
        results.calculate()
        results.printResults('results2.txt')

if __name__ == "__main__":
    main()
    
