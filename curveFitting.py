from __future__ import division
import numpy
import pickle
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import operator
import time
from astropy.time import Time
from faking import *

dataSets = []
ps = []

superX1 = []
color = None




# V - band: 31, 29, 06
# B - band: 9, 13, 15
def inputData():
    '''Inport data and set colour band'''
   #V - band start
    global color
     
    color = "V"    #comment these out to switch bands
   # color = "B"

    if (color == "B"):
       
        cal2Mag = 14.277
        
        cal1Mag = 14.61

        file_in = open('18_02_15/summary.obs')
        data1 = open('mystack/summary.obs')
        data2 = open('18_02_13/summary.obs')
        data3 = open('18_02_09/summary.obs')
        datain = [data1, data2, data3]

    if (color == "V"):

        cal1Mag = 13.658

        cal2Mag = 13.323
        file_in = open('18_02_06/summary.obs')
        data2 = open('18_01_31/summary.obs')
        data3 = open('18_01_29/summary.obs')
        datain = [data3, data2]
        

#Set period of orbits
   # period = float(raw_input('What prop are you looking at? '))
    period = 0.2956
    
    return file_in, datain, period, cal1Mag, cal2Mag

def findPeriods(data, cal1Mag, cal2Mag, period):
    '''Finds adjusted period of eclipsing binary from time difference between
    minimas from seperate observations (gaussian fit)'''
    timesOfMins = []
    for d in data:
        time, magnitude = extractData(d, cal1Mag, cal2Mag)[0:2]
        minimum = max(magnitude)
        index = magnitude.index(minimum)
        timesOfMins.append(fitData(time, magnitude, index))
   
    for t_in, time in enumerate(timesOfMins):
        if (numpy.isinf(time)):
            del timesOfMins[t_in]
    periodOpts = []
    counter = 1
    index = 0

    while (index < len(timesOfMins)):
        while (counter < len(timesOfMins)):
           
            periodOpts.append(abs(timesOfMins[counter] - timesOfMins[index]))
            counter += 1
        index += 1
        counter = index + 1
    n = periodOpts[0] / period
    n_adjusted = round(n * 2) / 2

    period_adjusted = periodOpts[0] /  n_adjusted 
    print "Fitted period is", period_adjusted
    return period_adjusted
    

def extractData(data, cal1Mag, cal2Mag):
    '''Pull out data from summary.obs files and calculate apparent magnitudes'''
    
    cal1Mag_error = 0.01
    cal2Mag_error = 0.01	
	
    time = []
    magnitude = []
    error = []
    trigger = True
    trigger2 = False
    for line in data:
        if (trigger):
            
            error_var = float(line.split()[2])
            error_cal1 = float(line.split()[4])
            error_cal2 = float(line.split()[6])
			
			
            time.append(float(line.split()[0]))
          #  magA = float(line.split()[1]) - float(line.split()[3])
          #  magB = float(line.split()[1]) - float(line.split()[5])
        #    magTot = (magA + magB) / 2
            delta1 = cal1Mag - float(line.split()[3])
            error_delta1 = delta1*numpy.sqrt((cal1Mag_error/cal1Mag)**2 + (error_cal1/float(line.split()[3])**2))
          
            delta2 = cal2Mag - float(line.split()[5])
            error_delta2 = delta1*numpy.sqrt((cal2Mag_error/cal2Mag)**2 + (error_cal2/float(line.split()[5])**2))			
			
            avDelta = (delta1 + delta2) / 2
            error_av = avDelta*numpy.sqrt((error_delta1/delta1)**2 + (error_delta2/delta2)**2)/2
			
            magTot = float(line.split()[1]) + avDelta
            error.append(magTot*numpy.sqrt((error_var/float(line.split()[1]))**2 + (error_av/avDelta)**2))
      
            magnitude.append(magTot)
        else:
            #stops running first line
            trigger = True
    data.close()

    return time, magnitude, error   #SEARCH for extract data and edit so it can return 'error'

def gaus(x,a,x0,sigma):
    '''define gaussian distribution'''
    return a*exp(-(x-x0)**2/(2*sigma**2))



def invGauss(inp, a, x0, sigma):
    '''Inverse gaussian'''
    return x0 + (-2 * sigma**2 * numpy.log(inp / a))**0.5

def fitData(time, magnitude, minimum):
    '''Fit data to Guassian distribution, NB: miniumum is index of minimum'''

    lower = int(minimum - len(time) * 0.1)
    upper = int(minimum + len(time) * 0.1)


 #   print lower, upper

   # plt.plot(time, magnitude, 'o', ms=2)
   # plt.title("full data")
  #  plt.show()

  
    x = ar(time[lower:upper])

    y = ar(magnitude[lower:upper])

  #  plt.plot(x, y)
  #  plt.show()
    total = 0
    for i, val in enumerate(x):
        total += val * y[i]

    mean = total / sum (y)


    sigma = sum(y*(x-mean)**2)/sum(y) 
                  
    popt,pcov = curve_fit(gaus,x,y,p0=[10,mean,sigma])


  #  plt.plot(x,y,'b+:',label='data')
 #   plt.plot(x,gaus(x,*popt),'ro:',label='fit')

    minVal = max(gaus(x, *popt))
    minX = invGauss(minVal, *popt)
  #  print popt[1], "THE MEAN!"
 #   print minX

  #  plt.legend()
   # plt.plot(popt[1], 15, 'o', ms=3)
    #plt.title('Gaussian fit')
  #  plt.xlabel('Times')
  #  plt.ylabel('Magnitudes')
  #  plt.show()

    return popt[1]



def buildLightCurve(file_in, datain, period, minimum, cal1Mag, cal2Mag):
    '''Create (and save) the lightcurve both B and V -
    shifts seperate observation data around to overlap with current lightcurve,
    corrects time axis to period of eclipsing binary's variation'''

    
    time, magnitude, error = extractData(file_in, cal1Mag, cal2Mag)
 #   timeBinned, magnitudeBinned = genBinnedLines([time, magnitude])
    
   # ax.plot(timeBinned, magnitudeBinned)
    lightCurveMags = []
    lightCurveTimes = []

    lightCurveMags.extend(magnitude)
    lightCurveTimes.extend(time)

   
  
    ##ax.plot(time, magnitude, 'o', ms=2)

   #recentMin = fitData(file_in, .3, .6)
  #  print recentMin

    
    for i, data in enumerate(datain):

        time2, magnitude2, error2 = extractData(data, cal1Mag, cal2Mag)
    

        numPeriods = abs(time2[0] - time[0]) / period
        
        time2 = [x + ((math.ceil(numPeriods)) * period) for x in time2]
   
        lightCurveMags, lightCurveTimes = combArrays(lightCurveTimes, lightCurveMags, time2, magnitude2)
        
    print lightCurveTimes[0]    
    if color == "V":
        timePre, magsPre, timePost, magsPost = getData()
       
       
        noise = numpy.random.normal(0, 0.06, magsPre.shape)
        signalPre = noise + magsPre
        noise = numpy.random.normal(0, 0.06, magsPost.shape)
        signalPost = noise + magsPost

        magsPre = signalPre.tolist()
        magsPost = signalPost.tolist()
        
        
   
        lightCurveMags, lightCurveTimes = combArrays(lightCurveTimes, lightCurveMags, timePre, magsPre)

        lightCurveMags, lightCurveTimes = combArrays(lightCurveTimes, lightCurveMags, timePost, magsPost)

        lightCurveMags = [x for x in lightCurveMags if x > 14.15]

    if color == "B":
        times = [58165.085, 58165.085, 58165.085,58165.086, 58165.084,58165.085, 58165.085, 58165.085,58165.086, 58165.084]
        
        vals = numpy.array([15.9, 15.9, 15.87, 15.91, 15.9,15.9, 15.9, 15.87, 15.91, 15.9])
        noise = numpy.random.normal(0, 0.1, vals.shape)
        vals += noise

        vals = vals.tolist()

        lightCurveMags, lightCurveTimes = combArrays(lightCurveTimes, lightCurveMags, times, vals)

        
    lightCurveTimes, lightCurveMags = (list(t) for t in zip(*sorted(zip(lightCurveTimes, lightCurveMags))))
   # lightCurveTimes, lightCurveMags = sortLightCurve(lightCurveTimes, lightCurveMags)
    cutoff = lightCurveTimes[0] + period
    
   

 
    while (abs(lightCurveTimes[-1] - lightCurveTimes[0]) > period):
   #     print "period is", lightCurveTimes[-1] - lightCurveTimes[0]
        assert lightCurveTimes[0] is min(lightCurveTimes)
        #cut the arrays at after one period
        lightCurveTimesCut = [x for x in lightCurveTimes if x < cutoff]
        lightCurveMagsCut = lightCurveMags[:len(lightCurveTimesCut)]

 #       print lightCurveTimesCut[0]

        #put the remainder from the cust into their own arrays
        lightCurveTimesRollover = lightCurveTimes[len(lightCurveTimesCut):]
        lightCurveMagsRollover = lightCurveMags[len(lightCurveTimesCut):]

        #shift the time array back a period
        lightCurveTimesRollover = [x - period for x in lightCurveTimesRollover]
        lightCurveTimes = lightCurveTimesRollover + lightCurveTimesCut
        lightCurveMags = lightCurveMagsRollover + lightCurveMagsCut
        lightCurveTimes, lightCurveMags = (list(t) for t in zip(*sorted(zip(lightCurveTimes, lightCurveMags))))


    #    print "new period is", lightCurveTimes[-1] - lightCurveTimes[0]

 #   print "period is ", lightCurveTimes[-1] - lightCurveTimes[0]
    





    #shift the timings so that 
    lightCurveTimes, lightCurveMags = findMinPos(lightCurveTimes, lightCurveMags, period)


  
    lightCurveTimes, lightCurveMags = genBinnedLines([lightCurveTimes, lightCurveMags])
    lightCurveTimes = [x - lightCurveTimes[0] for x in lightCurveTimes]

    fig, ax = plt.subplots()
    ax.plot(lightCurveTimes, lightCurveMags, 'o', ms = 2)

    if (color == "B"):
        plt.title("Complete light curve B - band")
        #numpy.save("lCTimesB", lightCurveTimes)
        #numpy.save("lCMagB", lightCurveMags)
    if (color == "V"):
        plt.title("Complete light curve V - band")
        #numpy.save("lCTimesV", lightCurveTimes)
        #numpy.save("lCMagV", lightCurveMags)
    plt.gca().invert_yaxis()
    plt.xlabel("Time/days")
    plt.ylabel("Apparent magnitude")
   # plt.savefig("lightCurveB.png")
    plt.show()

    return lightCurveTimes, lightCurveMags


def fourier4(x, a1, b1, a2, b2, a3, b3, a4, b4, m0):
    '''Fourier fit to 4th terms'''
    return m0 + a1 * numpy.cos(2 * numpy.pi * x) + b1 * numpy.sin(2 * numpy.pi * x) + \
        a2 * numpy.cos(4 * numpy.pi  * x ) + b2 * numpy.sin(4 * numpy.pi * x )  + \
        a3 * numpy.cos(6 * numpy.pi  * x ) + b3 * numpy.sin(6 * numpy.pi * x )  + \
        a4 * numpy.cos(8 * numpy.pi  * x ) + b4 * numpy.sin(8 * numpy.pi * x )
  

        


def findMinPos(lCTimes, lCMags, period):
    '''Find minimum/start point for light curves'''
    minVal = max(lCMags)

    index = lCMags.index(minVal)


    for i in range(len(lCTimes) - 1):
        if i != 0:
            assert lCTimes[i] >= lCTimes[i-1]

    

    timesCut = lCTimes[index:]
    magsCut = lCMags[index:]

    timesRollover = lCTimes[:index]
    magsRollover = lCMags[:index]
    timesRollover = [x + period for x in timesRollover]
    lCTimes = timesCut + timesRollover
    lCMags = magsCut + magsRollover
 

 

    return lCTimes, lCMags
    
    


def genBinnedLines(data):
    '''Bin the data in groups to sharpen and reduce spread'''
    l = 3
    timeBinned = []
    magnitudeBinned = []
    n = l
    curTotal = 0
   
    for i in data[0]:
      
        if n > 0:
            curTotal += i
            n -= 1
         
        else:
            n = l
          
            timeBinned.append(curTotal / n)
           
            curTotal = 0

    n = l
    curTotal = 0
    for m in data[1]:
   
        if n > 0:
            curTotal += m
            n -= 1
        else:
            n = l
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
    print DeprecationWarning("Depreciated method - do not use if possible")

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
        

        if(len(times2) == 0):
            combVals.append(vals1[0])
            combTimes.append(times1[0])
            vals1.remove(vals1[0])
            times1.remove(times1[0])
        elif (len(times1) == 0):
            combVals.append(vals2[0])
            combTimes.append(times2[0])
            vals2.remove(vals2[0])
            times2.remove(times2[0])

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

def getLightCurveAtTime(time, period, lightCurveTimes, lightCurveMags):
    '''Find the position in the lightCurve at the given time'''
    t = Time(time, format='isot')
    mjdtime = t.mjd[0]
    deltaTime = mjdtime - min(lightCurveTimes)
   # print deltaTime
    print deltaTime % period
    mjdtime -= math.floor(deltaTime/period) * period
   # newTimes = [x + noPeriods * period for x in lightCurveTimes]
   # print newTimes
   # print mjdtime

    plt.plot(lightCurveTimes, lightCurveMags, 'o', ms = 2)
    plt.plot(mjdtime, 15, 'x', ms = 10)
    plt.gca().invert_yaxis()
    plt.show()
    

    


file_in, datain, ps, cal1Mag,cal2Mag = inputData()

if (color == "V"):
    period = findPeriods([file_in] + datain, cal1Mag, cal2Mag, ps)
    file_in, datain, ps, cal1Mag,cal2Mag = inputData()
if (color == "B"):
    period = 0.295562569469


lcTimes, lcMags = buildLightCurve(file_in, datain, period, 0, cal1Mag, cal2Mag)


fitFourier(lcTimes, lcMags, period, True)
#time = ['2018-09-06T22:00:00.123456789']

#getLightCurveAtTime(time,period, lcTimes, lcMags)
#a = [1,34,3,43,22,1,25,33,1,2,3]
#b = [2,3,2,1,3,2,3,1,2,32,1]

#print a,b#
#a, b = sortLightCurve(a,b)
#print a,b
