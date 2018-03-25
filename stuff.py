def binarySearch(val, data):
    '''Totally pointless binary search method I wrote - literally ignore it completely it does nothing'''
    first = 0
    last = len(data) - 1

    while first <= last:
        print "value", val

        midpoint = (first + last) // 2
        print "midpoint", midpoint
        if data[midpoint] == val:

            return midpoint

        else:
            if val < data[midpoint]:

                last = midpoint - 1
                midpoint -= 1
            else:

                first = midpoint + 1
                midpoint += 1

    return midpoint


def sortLightCurve(time, mags):
    '''Another pointless sorting method'''
    sortTime = []
    sortMags = []

    for i, t in enumerate(time):

        if sortTime != []:

            index = binarySearch(t, sortTime)

            print "index returned", index

            sortTime.insert(index, t)
            sortMags.insert(index, mags[i])
            print sortTime

        else:
            sortTime.append(t)
            sortMags.append(mags[i])

    return sortTime, sortMags