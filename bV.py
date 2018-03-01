from __future__ import division
import matplotlib.pyplot as plt
import numpy
v_Mags = numpy.load("lCMagV.npy")
v_Times = numpy.load("lCTimesV.npy")
b_Mags = numpy.load("lCMagB.npy")
b_Times = numpy.load("lCTimesB.npy")



while(len(b_Mags) != len(v_Mags)):
    
    diff = -len(v_Mags) + len(b_Mags)
    if (diff > 0):
        print diff
        dropper = int(len(b_Mags) / diff)
        print len(b_Mags)
        b_Mags = [x for i, x in enumerate(b_Mags) if i % dropper != 0]
        b_Times = [x for i, x in enumerate(b_Times) if i % dropper != 0]
        print len(v_Mags)
        print len(b_Mags)
    else:
        diff = diff * -1
        dropper = int(len(v_Mags) / diff)
        print len(b_Mags)
        v_Mags = [x for i, x in enumerate(v_Mags) if i % dropper != 0]
        v_Times = [x for i, x in enumerate(v_Times) if i % dropper != 0]
        print len(v_Mags)
        print len(b_Mags)
              
plt.subplot(221)

plt.plot(v_Times, v_Mags, 'o', ms = 2)
plt.gca().invert_yaxis()
plt.title("V - colour magnitude")

plt.subplot(222)
plt.plot(b_Times, b_Mags, 'o', ms = 2)
plt.gca().invert_yaxis()
plt.title("B - colour magnitude")

plt.subplot(223)
bMinusV = [i - j for i in b_Mags for j in v_Mags]
plt.plot(bMinusV)
plt.gca().invert_yaxis()
plt.title("B - V")

plt.subplot(224)
temperatures = []
for bv in bMinusV:
    temperatures.append(4600 * (1 / (0.92 * (bv) + 1.7) + (1 / 0.92 * bv + 0.62)))
plt.plot(temperatures)
plt.title("Temperature Variation")
plt.show()
