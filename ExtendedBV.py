from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

def gaus(x, a, x0, sigma):
    '''define gaussian distribution'''
    return a*exp(-(x-x0)**2/(2*sigma**2))

def invGauss(inp, a, x0, sigma):
    '''define inverse gaussian distribution'''
    return x0 + (-2 * sigma**2 * np.log(inp / a))**0.5
    

def fitData(time, magnitude, minimum, scalefactor, colour):
    '''Fit data to Guassian distribution'''

    #print minimum, len(time)
    lower = int(minimum - len(time) * scalefactor)
    upper = int(minimum + len(time) * scalefactor)

    
    #print lower, upper, minimum
    # plt.figure()
    # plt.plot(time, magnitude, 'o', ms=2)
    # plt.title("full data")
    # plt.show()

  
    x = ar(time[lower:upper])

    y = ar(magnitude[lower:upper])
    # plt.figure()
    # plt.plot(x, y, 'o')
    # plt.show()
    total = 0
    for i, val in enumerate(x):
        total += val * y[i]

    # print time
    # print magnitude
    # print time[lower:upper]
    # print magnitude[lower:upper]
    # print x
    # print y
    
    mean = total / sum (y)


    sigma = sum(y*(x-mean)**2)/sum(y) 
                  
    popt,pcov = curve_fit(gaus,x,y,p0=[10,mean,sigma])
    plt.figure()
    plt.plot(x,y,'b+:',label='data', color = "red")
    plt.plot(x,gaus(x,*popt),'ro:',label='fit', color = colour)

    #minVal = max(gaus(x, *popt))
    #minX = invGauss(minVal, *popt)

    plt.legend()
    #plt.title('Gaussian fit')
    plt.xlabel('Times')
    plt.ylabel('Magnitudes')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()

    print popt 
    print pcov
    
    #return popt,pcov
    
    
VBand_mag = np.load("lCMagV.npy").tolist()
VBand_time = np.load("lCTimesV.npy").tolist()
VBand_mag.extend(np.load("lCMagV.npy").tolist())
VBand_time.extend((VBand_time[-1]+np.load("lCTimesV.npy")).tolist())




BBand_mag = np.load("lCMagB.npy").tolist()
BBand_time = np.load("lCTimesB.npy").tolist()
BBand_mag.extend(np.load("lCMagB.npy").tolist())
BBand_time.extend((BBand_time[-1]+np.load("lCTimesB.npy")).tolist())





VBand_mag = VBand_mag[50:-50]
VBand_time = VBand_time[50:-50]
BBand_mag = BBand_mag[50:-50]
BBand_time = BBand_time[50:-50]

V_minimum = VBand_mag.index(max(VBand_mag))
B_minimum = BBand_mag.index(max(BBand_mag))

fitV = fitData(VBand_time, VBand_mag, V_minimum, 0.1, "green")
fitB = fitData(BBand_time, BBand_mag, B_minimum, 0.1, "blue")


#V_min_time = fitV[1]
#B_min_time = fitB[1]

#V_min_mag_1 = fitV[0]*exp(-(0)**2/(2*fitV[2]**2))
#B_min_mag_1 = fitB[0]*exp(-(0)**2/(2*fitB[2]**2))

#print V_min_mag_1
#print B_min_mag_1

# plt.figure()
# #plt.plot(VBand_time,VBand_mag, 'o', ms = 2)
# plt.plot(BBand_time,BBand_mag, 'o', ms = 2)
# plt.gca().invert_yaxis()
# plt.tight_layout()
# plt.show()