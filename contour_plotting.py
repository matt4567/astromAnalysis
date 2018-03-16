import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


def inputData(file):
    file_in = open(file)
    return file_in


def extractData(data):
    for line in data:
        
        x = line.split()
        break 
    
    y = []
    triggerY = False

    Z = np.zeros((16,16))
    
    for i, line in enumerate(data):
            
            y.append(float(line.split()[0]))
            Z[i,:] = line.split()[1:]
 
        
    x = np.asarray(x)
    y = np.asarray(y)

    #print x
    #print y
    
    X, Y = np.meshgrid(x, y)
    data.close()

    #print np.shape(X)
    #print np.shape(Y)
    #print np.shape(Z)
    return X, Y, Z

#print extractData(inputData("NFCSD.txt"))
plt.figure()
x, y, z = extractData(inputData("NFCSD.txt"))
print np.shape(x)

CS = plt.contour(x, y, z)
plt.clabel(CS,inline=1, fontsize=10)
plt.show()
