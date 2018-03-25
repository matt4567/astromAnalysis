import matplotlib.pyplot as plt
import math
import numpy
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
input_file = open("summary.obs")


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

    return time, magnitude, error

def fourier4(x, a1, b1, a2, b2, a3, b3, a4, b4, m0):
    
    return m0 + a1 * numpy.cos(2 * numpy.pi * x) + b1 * numpy.sin(2 * numpy.pi * x) + \
        a2 * numpy.cos(4 * numpy.pi  * x ) + b2 * numpy.sin(4 * numpy.pi * x )  + \
        a3 * numpy.cos(6 * numpy.pi  * x ) + b3 * numpy.sin(6 * numpy.pi * x )  + \
        a4 * numpy.cos(8 * numpy.pi  * x ) + b4 * numpy.sin(8 * numpy.pi * x )
  

def fitFourier(times, mags, period, plot):
    # fits
    x = ar(times)
    
  #  print x
    x = x / period
    
    y = ar(mags)
   # popt,pcov = curve_fit(gaus,x,y,p0=[10,mean,sigma])
    popt, pcov = curve_fit(fourier4, x, y, p0=[1, 1, 1, 1, 1, 1, 1,1, -14.5])
 #   print pcov
    if plot:
        print "a1", popt[0], "+/-", numpy.sqrt(pcov[0,0])
        print "a2",  popt[2],"+/-", numpy.sqrt(pcov[2,2])
        print "a3",  popt[4],"+/-", numpy.sqrt(pcov[4,4])
        print "a4",  popt[6],"+/-", numpy.sqrt(pcov[6,6])
        print "b1", popt[1],"+/-", numpy.sqrt(pcov[1,1])
        print "b2",  popt[3],"+/-", numpy.sqrt(pcov[3,3])
        print "b3",  popt[5],"+/-", numpy.sqrt(pcov[5,5])
        print "b4",  popt[7],"+/-", numpy.sqrt(pcov[7,7])
        print "m0", popt[8], "+/-", numpy.sqrt(pcov[8,8])
        print popt
  #  print popt

    # further plots
    newMags = fourier4(x,*popt)
    deltaT = x[2] - x[1]
    extTimes = numpy.arange(x[-1], x[-1] + deltaT * 50, deltaT)
    extTimes = ar(extTimes.tolist())
    extMags = fourier4(extTimes, *popt)
    extTimes = [a * period for a in extTimes]

    beforeTimes =  numpy.arange(x[0] - deltaT*20, x[0], deltaT)
    beforeTimes = ar(beforeTimes.tolist())
    beforeMags = fourier4(beforeTimes, *popt)
    beforeTimes = [a * period for a in beforeTimes]

    if plot == True:
        p2, = plt.plot(x, newMags, 'o', ms = 3, label="Fourier fit")
        p3, = plt.plot(x, y, 'ro:', ms = 2, label = "Our data")

        plt.legend()
        plt.gca().invert_yaxis()
        #plt.title("Fourier fitting of light curve")
        plt.xlabel("Phase")
        plt.ylabel("Apparent magnitude")
		

        plt.show()

    return beforeTimes, beforeMags, extTimes, extMags


def getData():
    
    time, mags = extractData(input_file, 13.658, 13.323)
   

    return fitFourier(time, mags, 0.295609, False)

    
