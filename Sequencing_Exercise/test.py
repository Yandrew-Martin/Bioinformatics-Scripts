import re

header = ['ref_1','A_1','C_1','G_1','T_1','call_1','ref_2','A_2','C_2','G_2','T_2']
regex = r'[ATGC]_'

test = re.compile(regex)
newlist = list(filter(test.match, header))
print(newlist)

with open('sequencing_data_biochem2.csv', newline='') as a:
    results = Results()
    reader = csv.reader(a)
    header = next(reader)
    correction = {1:2,2:3,3:4,4:1}
    
    for row in reader:
        cycle1 = row[1:5]
        cycle2 = row[7:11]
        maxrow1 = max(cycle1)
        maxrow2 = max(cycle2)
        
        if cycle1 == ['0.0','0.0','0.0','0.0']:
            results.basecount1[0] += 1
            base1 = 0
        else:
            base1 = correction[cycle1.index(maxrow1)+1]
            results.intensitysum1[base1] += Decimal(maxrow1)
            results.basecount1[base1] += 1
        if cycle2 == ['0.0','0.0','0.0','0.0']:
            results.basecount2[0] += 1
            base2 = 0
        else:
            base2 = correction[cycle2.index(maxrow2)+1]
            results.intensitysum2[base2] += Decimal(maxrow2)
            results.basecount2[base1] += 1
            
        results.call1.append(basecall[base1])
        results.call2.append(basecall[base2])
        
        if row[0] != basecall[base1]:
            results.basecount1[5] += 1
        if row[6] != basecall[base2]:
            results.basecount2[5] += 1
        results.probecount += 1
        
    perror1 = (results.basecount1[5]/results.probecount)*100
    phred1 = -10*(math.log10(perror1/100))    
    perror2 = (results.basecount2[5]/results.probecount)*100
    phred2 = -10*(math.log10(perror2/100))
    
    print(results.call1)
    print(results.basecount1)
    print('percent error:', perror1, '%')
    print('phred quality score:', phred1)
    print(results.call2)
    print(results.basecount2)
    print('percent error:', perror2, '%')
    print('phred quality score:', phred2)
    
    for i in results.intensityavg1:
        results.intensityavg1[i] = results.intensitysum1[i]/results.basecount1[i]
        results.intensityavg2[i] = results.intensitysum2[i]/results.basecount2[i]
        
    print((results.intensityavg1[1]+results.intensityavg1[2]+results.intensityavg1[3]+results.intensityavg1[3])/4)
    print((results.intensityavg2[1]+results.intensityavg2[2]+results.intensityavg2[3]+results.intensityavg2[3])/4)
    print('cycle 1 gc content:',((results.basecount1[2]+results.basecount1[3])/results.probecount), '%')
    print('cycle 2 gc content:',((results.basecount2[2]+results.basecount2[3])/results.probecount), '%')
    
    
with open('sequencing_data_biochem2.csv', newline='') as a:
    count1 = {0:0,1:0,2:0,3:0,4:0,5:0}
    count2 = {0:0,1:0,2:0,3:0,4:0,5:0}
    sum1 = {1:0,2:0,3:0,4:0}
    sum2 = {1:0,2:0,3:0,4:0}
    avg1 = {1:0,2:0,3:0,4:0}
    avg2 = {1:0,2:0,3:0,4:0}
    call1 = ''
    call2 = ''
    probecount = 0
    reader2 = csv.reader(a)
    header2 = next(reader2)
    
            
    for row in reader2:
        cycle1 = row[1:5]
        cycle2 = row[7:11]
        maxrow1 = max(cycle1)
        maxrow2 = max(cycle2)
        base1 = correction[cycle1.index(maxrow1)+1]
        base2 = correction[cycle2.index(maxrow2)+1]
        sum1[base1] += Decimal(maxrow1)
        sum2[base2] += Decimal(maxrow2)
        if cycle1 == ['0.0','0.0','0.0','0.0']:
            count1[0] += 1
            base1 = 0
        else:
            count1[base1] += 1
        if cycle2 == ['0.0','0.0','0.0','0.0']:
            count2[0] += 1
            base2 = 0
        else:
            count2[base1] += 1
        call1 += basecall[base1]
        call2 += basecall[base2]
        if row[0] != basecall[base1]:
            count1[5] += 1
        if row[6] != basecall[base2]:
            count2[5] += 1
        probecount += 1
        
    perror1 = (count1[5]/probecount)*100
    phred1 = -10*(math.log10(perror1/100))    
    perror2 = (count2[5]/probecount)*100
    phred2 = -10*(math.log10(perror2/100))
    
    print(call1)
    print(count1)
    print('percent error:', perror1, '%')
    print('phred quality score:', phred1)
    print(call2)
    print(count2)
    print('percent error:', perror2, '%')
    print('phred quality score:', phred2)
    print('cycle 1 gc content:',((count1[2]+count1[3])/probecount), '%')
    print('cycle 2 gc content:',((count2[2]+count2[3])/probecount), '%')

    for i in avg1:
        avg1[i] = sum1[i]/count1[i]
        avg2[i] = sum2[i]/count2[i]

    print((avg1[1]+avg1[2]+avg1[3]+avg1[3])/4)
    print((avg2[1]+avg2[2]+avg2[3]+avg2[3])/4)