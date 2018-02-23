from __future__ import division
import numpy
import pickle
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import operator

dataSets = []
ps = []

superX1 = []
color = None

# V - band: 31, 29, 06
# B - band: 9, 13, 15
def inputData():
   #V - band start
  #  color = "V"
  #  file_in = open('../18_02_06/summary.obs')
  #  data2 = open('../18_01_31_new/summary.obs')
  #  data3 = open('../18_01_29_new/summary.obs')
   #B - band start
    color = "B"
    file_in = open('../18_02_15/summary.obs')
    data2 = open('../18_02_13/summary.obs')
    data3 = open('../18_02_09/summary.obs')

    #   minVal1 = fitData(file_in, 30/484, 120/484)
    
#    print "time 1", minVal1
#    file_in = open('../18_01_31_new/summary.obs')
  #  minVal2 = fitData(file_in, 310/484, 450/484)
  #  print "time 2", minVal2
   # time1 = optimise(x, guassFit)




    #trial_in = open('../18_01_31_new/summary.obs')
    #trial_in = open('../18_01_29_new/results.diff')
  #  data1 = open('../18_02_02/summary.obs')

    datain = [data3, data2]

   # period = minVal2 - minVal1
    period = float(raw_input('What prop are you looking at? '))
  # period = period / prop
    print period
 ## file_in = open('../18_02_06/summary.obs')
    
   
    return file_in, datain, period

def extractData(data):
    time = []
    magnitude = []
    trigger = False
    trigger2 = False
    for line in data:
        if (trigger):
            
            time.append(float(line.split()[0]))
            magA = float(line.split()[1]) - float(line.split()[3])
            magB = float(line.split()[1]) - float(line.split()[5])
            magTot = (magA + magB) / 2
                     
            magnitude.append(magTot)
        else:
            trigger = True
    data.close()
   
    return time, magnitude

def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def fitData(file_in, lower, upper):
    if (upper == None):
        time, magnitude = extractData(file_in)
        index, value = max(enumerate(time), key=operator.itemgetter(1))
        lower = index - 30
        upper = index + 40
       
        
    else: time, magnitude = extractData(file_in)
#    time, magnitude = genBinnedLines([time, magnitude])
    print len(time)
    print lower, upper
    lower = int(lower * (len(time)))
    upper = int(upper * (len(time))) 

    print lower, upper

  
    x = ar(time[lower:upper])

    y = ar(magnitude[lower:upper])
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
    return x[index]
    

def findRelMagnitudesAndPlot(file_in, datain, period, minimum):

    
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
   #     timeBinned, magnitudeBinned = genBinnedLines([time2, magnitude2])

        #todo: add a boolean which allows you to check the data in advance to set distribution around max point for more automation

     #  timeVal = fitData(data, None, None)
     #  print timeVal

     #  noPeriods = (recentMin - timeVal) / period
       #print noPeriods, '<<--num periods'
        
     
        
        
     #   pos = int(len(time) * period)
        numPeriods = abs((time2[0] - time[0])) / period
        time2 = [x + ((math.ceil(numPeriods)) * period) for x in time2]
      # timeBinned = [x + ((math.ceil(numPeriods)) * period) for x in timeBinned]
      # ax.plot(timeBinned, magnitudeBinned)
        #
        #ax.plot(time2, magnitude2, 'o', ms=2)
       

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
    lightCurveCompleteTimes, lightCurveCompleteMags = genBinnedLines([lightCurveCompleteTimes, lightCurveCompleteMags])
    ax.plot(lightCurveCompleteTimes, lightCurveCompleteMags, 'o', ms = 2)
    if (color == "B"):
        
        ax.title("Complete light curve B - band")
    if (color == "V"):
        ax.title("Complete light curve V - band")
    plt.gca().invert_yaxis()
    plt.show()
    


def genBinnedLines(data):
  
  
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

    global superX1
    superX1.append(event.xdata)

   

def findPeriod(tGuess, alpha, times, magnitudes):
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


    
def optimise(x, gaussFit):

  #  x = [-i for i in x]

   # time, magnitude = extractData(file_in)
  #  timeBinned, magnitudeBinned = genBinnedLines([time, magnitude])
    fig, ax = plt.subplots()
    ax.plot(x, gaussFit)
  
#    ax.plot(time, magnitude, 'o', ms=2)
    
    plt.gca().invert_yaxis()

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    fig.canvas.mpl_disconnect(cid)
    
    return superX1[0]

    print 'First time input: ', superX1[0]
    
    t1 = findPeriod(superX1[0], 0.1, timeBinned, magnitudeBinned)

    print 'First time optimsed: ', t1
    

    print 'Second time input: ', superX1[1]
    
    t2 = findPeriod(superX1[1], 0.1, timeBinned, magnitudeBinned)

    print 'Second time optimised: ', t2

    proportion = float(raw_input("What prop of period do you have? "))
    
    period = t2 - t1
    period = period / proportion
    print "the period is", period
    return period

def combArrays(times1, vals1, times2, vals2):
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



file_in, datain, ps          = inputData()
findRelMagnitudesAndPlot(file_in, datain, ps, 0)
