import numpy
import matplotlib.pyplot as plt
import heapq
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
#import curveFitting

def inputData(file):
    file_in = open(file)
    return file_in

def extractDataAs(data):
    time = []
    magnitude = []
    trigger = False
    trigger2 = False
    for line in data:
        if (trigger):

            time.append(float(line.split()[0]))
            magnitude.append(float(line.split()[2]))

        else:

            trigger = True


    return time,magnitude

def extractData(data):

    cal2Mag = 14.277
        
    cal1Mag = 14.61
    
    time = []
    magnitude = []
    trigger = False
    trigger2 = False
    for line in data:
        if (trigger):
            
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
            
           
            trigger = True
    data.close()

    return time, magnitude



def savefiles(file_in):
    time, magnitude = extractData(file_in)
    #numpy.save("asassnTimes", time)
    #numpy.save("asassnMagV", magnitude)



def sort_asassn(file_in):
    P = 0.2956
    time, magnitude = extractDataAs(file_in)
    #time_ignore, orig_magnitude = extractData(file_in)
    orig_times = []
    minMag = min(magnitude)
    min_index = magnitude.index(minMag)
    t_0 = time[min_index]
    for i in range(len(time)):
        orig_times.append(time[i])
        n = (time[i]-t_0)/P
        time[i] = time[i] - (int(n)*P)

    time, magnitude = (list(t) for t in zip(*sorted(zip(time, magnitude))))
        
   # plt.figure()
   # plt.plot(time, magnitude, 'o', ms = 3)
   # plt.gca().invert_yaxis()

    #plt.figure()
    #plt.plot(orig_times, magnitude,'o', ms=3)
    #plt.gca().invert_yaxis()
   # plt.show()

    return time, magnitude

def find_min_n(data, n):
    time, magnitude = data
    time = [x - 2400000 for x in time]  #only for asassn data
    mins = heapq.nlargest(n,magnitude)
    min_times = []
    for i in range(len(mins)):
        index = magnitude.index(mins[i])
        min_times.append(time[index])
        
    return min_times


def gaus(x, a, x0, sigma):
    '''define gaussian distribution'''
    return a*exp(-(x-x0)**2/(2*sigma**2))

def invGauss(inp, a, x0, sigma):
    '''define inverse gaussian distribution'''
    return x0 + (-2 * sigma**2 * numpy.log(inp / a))**0.5
    

def fitData(time, magnitude, minimum, scalefactor):
    '''Fit data to Guassian distribution'''

    #print minimum, len(time)
    lower = int(minimum - len(time) * scalefactor)
    upper = int(minimum + len(time) * scalefactor)


   # print lower, upper, minimum

    #plt.plot(time, magnitude, 'o', ms=2)
    #plt.title("full data")
    #plt.show()

  
    x = ar(time[lower:upper])

    y = ar(magnitude[lower:upper])

    #plt.plot(x, y, 'o')
    #plt.show()
    total = 0
    for i, val in enumerate(x):
        total += val * y[i]

    #print time
    #print magnitude
    #print time[lower:upper]
    #print magnitude[lower:upper]
    #print x
    #print y
    
    mean = total / sum (y)


    sigma = sum(y*(x-mean)**2)/sum(y) 
                  
    popt,pcov = curve_fit(gaus,x,y,p0=[10,mean,sigma])

   # plt.plot(x,y,'b+:',label='data')
   # plt.plot(x,gaus(x,*popt),'ro:',label='fit')

    minVal = max(gaus(x, *popt))
    minX = invGauss(minVal, *popt)

   # plt.legend()
   # plt.title('Gaussian fit')
   # plt.xlabel('Times')
   # plt.ylabel('Magnitudes')
   # plt.show()

    return popt[1]

def O_C():
    P_ours = 0.296097392087


    #Set up the two lots of data from our v band and the corrected asassn data
    time_ours, magnitude_ours = extractData(inputData("18_01_31/summary.obs"))
    our_minimum = magnitude_ours.index(max(magnitude_ours))
    our_curve_min = fitData(time_ours, magnitude_ours, our_minimum, 0.1)

    time_as, magnitude_as = sort_asassn(inputData("lc507905.dat"))
    time_as = [x - 2400000 for x in time_as] 
    as_minimum = magnitude_as.index(max(magnitude_as))
    as_curve_min = fitData(time_as, magnitude_as, as_minimum, 0.218)

    #find the raw data 4 minima and correct them using the difference found with the Gaussian fit, thus find observed minima
    min_interest_raw = find_min_n(extractDataAs(inputData("lc507905.dat")),3)
    min_interest_corrected = find_min_n(sort_asassn(inputData("lc507905.dat")),3)
    #print min_interest_raw
    #print min_interest_corrected
    observed_mins = []
    time_diff = []
    for mini in min_interest_corrected:
        #print mini
        time_diff.append(as_curve_min - mini)

    #print time_diff

    for i in range(len(min_interest_raw)):
        if time_diff[i] < 0:
            observed_mins.append(min_interest_raw[i] + time_diff[i])
        else:
            observed_mins.append(min_interest_raw[i] - time_diff[i])

    #print observed_mins
    
    #Find the number of periods (as a float) between minima, extract the int from this to find the calculated value for the minima
    calculated_mins = []
    period_ints = []
    for mini in observed_mins:
        period_number = (our_curve_min - mini)/P_ours
        period_ints.append(round(period_number*2)/2)
        calculated_mins.append((our_curve_min - (P_ours * round(period_number * 2)/2)))

    #print calculated_mins    
    calc_mins = numpy.array(calculated_mins)
    obs_mins = numpy.array(observed_mins)
    #print calc_mins
   # print obs_mins


    #O_C = obs_mins - calc_mins

   
    O_C = []
    for i in range(len(period_ints)):
        O_C.append((obs_mins[i] - calc_mins[i])/period_ints[i])
    #print O_C

    obs_mins, O_C = (list(t) for t in zip(*sorted(zip(obs_mins, O_C))))

    
    some_data_x = []
    some_data_y = []
    some_data_x.extend(numpy.ndarray.tolist(numpy.linspace(obs_mins[0],obs_mins[1])))
    some_data_x.extend(numpy.ndarray.tolist(numpy.linspace(obs_mins[1],obs_mins[2])))
    some_data_y.extend(numpy.ndarray.tolist(numpy.linspace(O_C[0],O_C[1])))
    some_data_y.extend(numpy.ndarray.tolist(numpy.linspace(O_C[1],O_C[2])))

    print len(some_data_x)
    print len(some_data_y)
    
    plt.figure()
    plt.plot(some_data_x, some_data_y,'o', color = 'blue')
    
    plt.figure()
    #popt,pcov = curve_fit(expo,obs_mins,O_C,p0=[P_ours,5.00,60000])
    #print popt, pcov
    #plt.plot(obs_mins,O_C,'b+:',label='data')
    #plt.plot(obs_mins,expo(obs_mins,*popt),'ro:',label='fit')

    popt,pcov = curve_fit(expo,some_data_x,some_data_y,p0=[P_ours,2.00,100000])
    print popt, pcov
    plt.plot(some_data_x,some_data_y,'b+:',label='data')
    plt.plot(some_data_x,expo(some_data_x,*popt),'ro:',label='fit')
    
    plt.figure()
    plt.plot(observed_mins, O_C, 'o')
    plt.xlabel("HJD")
    plt.ylabel("Observed - Calculated")
    plt.show()          
        
        
def expo(t,P_0,tau,t_0):
    return P_0*numpy.exp(tau/(t-t_0))




    
time_ours, magnitude_ours = extractData(inputData("18_01_31/summary.obs"))
our_minimum = magnitude_ours.index(max(magnitude_ours))
#plt.plot(time_as, magnitude_as,'o')
#plt.show()

#sort_asassn(inputData("lc507905.dat"))
#print find_min_n(inputData("lc507905.dat"),4)
#print as_minimum
#print fitData(time_ours, magnitude_ours, our_minimum)
O_C()
