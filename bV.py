from __future__ import division
import matplotlib.pyplot as plt
import numpy
v_Mags = numpy.load("lCMagV.npy")
v_Times = numpy.load("lCTimesV.npy")
b_Mags = numpy.load("lCMagB.npy")
b_Times = numpy.load("lCTimesB.npy")

#len(b_Mags) != len(v_Mags)
find = False

while(False):
    
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
bminusV = []
aVals = []
scale = []
diff = len(v_Mags) - len(b_Mags)
print diff
if diff > 0:
    scale = b_Times
    for i,t in enumerate(b_Times):
        a = min(v_Times, key=lambda x:abs(x-t))
        index = numpy.where(v_Times == a)
      #  print b_Mags[i], v_Mags[index[0]][0], b_Mags[i] - v_Mags[index[0]][0]
        
        bminusV.append(b_Mags[i] - v_Mags[index[0]][0])

if diff < 0:
    scale = v_Times
    for i,t in enumerate(v_Times):
        a = min(b_Times, key=lambda x:abs(x-t))
       
        
        if (a > 0.2769 and a < 0.27694):
            bminusV.append(0)
            continue

        
       
        index = numpy.where(b_Times == a)
     #   print t, b_Times[index[0]]
      #  print index
        print b_Mags[index[0]][0], v_Mags[i], b_Mags[index[0]][0] - v_Mags[i]
        bminusV.append(b_Mags[index[0]][0] - v_Mags[i])
        #
        #print bminusV[-1]

if diff == 0:
    scale = v_Times
    for i, t in enumerate(b_Mags):
        bminusV.append(t - v_Mags[i])


print len(v_Times), len(bminusV)
print "B-V color ", numpy.mean(bminusV)

diff = 50
diffB = (len(b_Times) // diff) - 1
diffV = (len(v_Times) // diff) - 1
bminusV = []
for i in range(diff):

    bminusV.append(numpy.mean(b_Mags[(i * diffB): (i + 1) * diffB]) - numpy.mean(v_Mags[i * diffV: (i + 1) * diffV]))


print bminusV


plt.suptitle("B-V analysis")
plt.subplot(221)

plt.plot(v_Times, v_Mags, 'o', ms = 2)
plt.gca().invert_yaxis()
plt.title("V - colour magnitude")

plt.subplot(222)
plt.plot(b_Times, b_Mags, 'o', ms = 2)
plt.gca().invert_yaxis()
plt.title("B - colour magnitude")

plt.subplot(223)


plt.plot(bminusV, 'o', ms = 2)

plt.title("B - V")

plt.subplot(224)
temperatures = []
for bv in bminusV:
    temperatures.append(4600 * (1 / (0.92 * (bv) + 1.7) + (1 / 0.92 * bv + 0.62)))
plt.plot(temperatures, 'o', ms =2)
plt.title("Temperature Variation")
plt.show()
#plt.savefig("bVplot.png")
