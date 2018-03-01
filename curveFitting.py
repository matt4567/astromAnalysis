from __future__ import division
import numpy
import pickle
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import operator
#from astropy.time import time

dataSets = []
ps = []

superX1 = []
color = None

cal1Mag = None
cal2Mag = None


# V - band: 31, 29, 06
# B - band: 9, 13, 15
def inputData():
    '''Inport data and set colour band'''
   #V - band start
    global color
     
  #  color = "V"
    color = "B"

    if (color == "B"):
        global cal2Mag
        cal2Mag = 14.277
        global cal1Mag
        cal1Mag = 14.61

        file_in = open('../18_02_15/summary.obs')
        data1 = open('../mystack/summary.obs')
        data2 = open('../18_02_13/summary.obs')
        data3 = open('../18_02_09/summary.obs')
        datain = [data1, data3, data2]

    if (color == "V"):
        global cal1Mag
        cal1Mag = 13.658
        global cal2Mag
        cal2Mag = 13.323
        file_in = open('../18_02_06/summary.obs')
        data2 = open('../18_01_31_new/summary.obs')
        data3 = open('../18_01_29_new/summary.obs')
        datain = [data3, data2]
        

#Set period of orbits
    period = float(raw_input('What prop are you looking at? '))

    print period
    return file_in, datain, period

def extractData(data):
    '''Pull out data from summary.obs files and calculate apparent magnitudes'''
    time = []
    magnitude = []
    trigger = False
    trigger2 = False
    for line in data:
        if (trigger):
            print line
            
            time.append(float(line.split()[0]))
          #  magA = float(line.split()[1]) - float(line.split()[3])
          #  magB = float(line.split()[1]) - float(line.split()[5])
        #    magTot = (magA + magB) / 2
            delta1 = cal1Mag - float(line.split()[3])
          
            delta2 = cal2Mag - float(line.split()[5])
            avDelta = (delta1 + delta2) / 2
            magTot = float(line.split()[1]) + avDelta
      
            magnitude.append(magTot)
        else:
            #stops running first line
            trigger = True
    data.close()

    return time, magnitude

def gaus(x,a,x0,sigma):
    '''define gaussian distribution'''
    return a*exp(-(x-x0)**2/(2*sigma**2))

def fitData(time, magnitude, minimum):
    '''Fit data to Guassian distribution'''
 #   if (upper == None):
 #       time, magnitude = extractData(file_in)
 #       index, value = max(enumerate(time), key=operator.itemgetter(1))
 #       lower = index - 30
 #       upper = index + 40
 #      
 #       
 #   else: time, magnitude = extractData(file_in)
##    time, magnitude = genBinnedLines([time, magnitude])
    print len(time)
 #   print lower, upper
 #   lower = int(lower * (len(time)))
 #   upper = int(upper * (len(time)))

    print minimum, "mini"
#set bounds for guassian fitting
    lower = int(minimum - len(time) * 0.5)
    upper = int(minimum + len(time) * 0.2)


 #   print lower, upper

  
    x = ar(time[lower:upper])

    y = ar(magnitude[lower:upper])
    plt.plot(x, y)
    plt.show()
    total = 0
    for i, val in enumerate(x):
        total += val * y[i]

    mean = total / sum (y)


    sigma = sum(y*(x-mean)**2)/sum(y) 
                  
    popt,pcov = curve_fit(gaus,x,y,p0=[10,mean,sigma])



    index, value = max(enumerate(gaus(x, *popt)), key=operator.itemgetter(1))

                       

    
    plt.plot(x,y,'b+:',label='data')
    plt.plot(x,gaus(x,*popt),'ro:',label='fit')
   
    plt.legend()
    plt.title('Gaussian fit')
    plt.xlabel('Times')
    plt.ylabel('Magnitudes')
    plt.show()

    start = time.index(x[index])
    origin = start + index
    return origin
    return x[index]
    

def findRelMagnitudesAndPlot(file_in, datain, period, minimum):
    '''Create (and save) the lightcurve both B and V'''

    
    time, magnitude = extractData(file_in)
 #   timeBinned, magnitudeBinned = genBinnedLines([time, magnitude])
    fig, ax = plt.subplots()
   # ax.plot(timeBinned, magnitudeBinned)
    lightCurveMags = []
    lightCurveTimes = []

    lightCurveMags.extend(magnitude)
    lightCurveTimes.extend(time)

    print len(lightCurveMags)
  
    ##ax.plot(time, magnitude, 'o', ms=2)

   #recentMin = fitData(file_in, .3, .6)
  #  print recentMin

    
    for i, data in enumerate(datain):

        time2, magnitude2 = extractData(data)

        numPeriods = abs((time2[0] - time[0])) / period
        time2 = [x + ((math.ceil(numPeriods)) * period) for x in time2]
        lightCurveMags, lightCurveTimes = combArrays(lightCurveTimes, lightCurveMags, time2, magnitude2)
        print len(lightCurveMags)

    cutoff = lightCurveTimes[0] + period
    print cutoff, 'cutoff'
   
    lightCurveTimesCut = [x for x in lightCurveTimes if x < cutoff]
    lightCurveTimesRollover = lightCurveTimes[len(lightCurveTimesCut):]
    lightCurveMagsRollover = lightCurveMags[len(lightCurveTimesCut):]
    lightCurveTimesRollover = [x - period for x in lightCurveTimesRollover]
    lightCurveMagsCut = lightCurveMags[:len(lightCurveTimesCut)]
    lightCurveCompleteMags, lightCurveCompleteTimes = combArrays(lightCurveTimesCut, lightCurveMagsCut,\
                                                                 lightCurveTimesRollover, lightCurveMagsRollover)
    print lightCurveCompleteTimes[-1] - lightCurveCompleteTimes[0]
    lightCurveCompleteTimes = findMinPos(lightCurveCompleteTimes, lightCurveCompleteMags, period)
    lightCurveCompleteTimes, lightCurveCompleteMags = genBinnedLines([lightCurveCompleteTimes, lightCurveCompleteMags])

    ax.plot(lightCurveCompleteTimes, lightCurveCompleteMags, 'o', ms = 2)
  # 
    if (color == "B"):
        plt.title("Complete light curve B - band")
        numpy.save("lCTimesB", lightCurveCompleteTimes)
        numpy.save("lCMagB", lightCurveCompleteMags)
    if (color == "V"):
        plt.title("Complete light curve V - band")
   #     numpy.save("lCTimesV", lightCurveCompleteTimes)
    #    numpy.save("lCMagV", lightCurveCompleteMags)
    plt.gca().invert_yaxis()
    plt.show()

    return lightCurveCompleteTimes, lightCurveCompleteMags
    

def findMinPos(lCTimes, lCMags, period):
    '''Find minimum/start point for light curves'''
    minVal = max(lCMags)
    print minVal
    index = lCMags.index(minVal)
    print len(lCTimes)
    timeCut = lCTimes[index]
    
    lCTimesShifted = [x + period if x < timeCut else x for x in lCTimes]
    print len(lCTimesShifted)
   # origin = fitData(lCTimes, lCMags, index)
    return lCTimesShifted
    
    


def genBinnedLines(data):
    '''Bin the data in groups to sharpen and reduce spread'''
  
    timeBinned = []
    magnitudeBinned = []
    n = 3
    curTotal = 0
   
    for i in data[0]:
      
        if n > 0:
            curTotal += i
            n -= 1
         
        else:
            n = 3
          
            timeBinned.append(curTotal / n)
           
            curTotal = 0

    n = 3
    curTotal = 0
    for m in data[1]:
   
        if n > 0:
            curTotal += m
            n -= 1
        else:
            n = 3
            magnitudeBinned.append(curTotal / n)
  
        
            curTotal = 0

  

    return timeBinned, magnitudeBinned
            
    

def onclick(event):
    '''define superX1 on pyplot.plot click'''
    global superX1
    superX1.append(event.xdata)

   

def findPeriod(tGuess, alpha, times, magnitudes):
    '''Iterative method for finding minimum of lightcurve'''
    offsetMean = 100

    while (offsetMean > .010):
        value = min(times, key=lambda x:abs(x-tGuess))
       
        #print times
        index = times.index(value)
      
        offset = []
       

        if (tGuess > times[index]):
            index += 1
        for i in range(20):
            offset.append(magnitudes[index + i] - magnitudes[index - i])
            
        offsetMean = numpy.mean(offset)
        tGuess -= offsetMean * alpha

    return tGuess


    


def combArrays(times1, vals1, times2, vals2):
    '''Combine 2 arrays at the correct times'''
    combVals = []
    combTimes = []
    while (len(vals1) != 0 or len(vals2) != 0):
    #    print len(vals2), len(times2)
        if (len(times1) == 0):
            combVals.append(vals2[0])
            combTimes.append(times2[0])
            vals2.remove(vals2[0])
            times2.remove(times2[0])

        elif(len(times2) == 0):
            combVals.append(vals1[0])
            combTimes.append(times1[0])
            vals1.remove(vals1[0])
            times1.remove(times1[0])

        elif (times1[0] <= times2[0]):
            combVals.append(vals1[0])
            combTimes.append(times1[0])
            vals1.remove(vals1[0])
            times1.remove(times1[0])
            
      
            
        
        else:
            combVals.append(vals2[0])
            combTimes.append(times2[0])
            vals2.remove(vals2[0])
            times2.remove(times2[0])
    return combVals, combTimes

def getLightCurveAtTime(time, lightCurveTimes, lightCurveMags):
    '''Find the position in the lightCurve at the given time'''
    mjdtime = time.mjd[0]
    deltaTime = int(mjdtime - lightCurveTimes[0])

    newTimes = [x + deltaTime for x in lightCurveTimes]

    plt.plot(newTimes, lightCurveMags)
    plt.plot(mjdtime, label = "Time")
    plt.show()
    

    


file_in, datain, ps          = inputData()
lcTimes, lcMags = findRelMagnitudesAndPlot(file_in, datain, ps, 0)

#time = ['2018-02-25T22:00:00.123456789']

#getLightCurveAtTime(time, lcTimes, lcMags)
