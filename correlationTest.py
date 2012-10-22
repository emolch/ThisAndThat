from matplotlib.pyplot import *
from pyrocko import io, trace
import numpy as np

infileDir = "/Users/mariuskriegerowski/hiwi/MineDemo/Demodataset_new.mseed"

allTracesInput = io.load(infileDir)

for trac in allTracesInput:
    print trac.channel

def lineFromAnglePoint(x,y,angle):
    xvals = [ x-1.*cos(angle), x+1.*cos(angle) ]
    yvals = [ y-1.*sin(angla), y+1.*sin(angle) ]
    return [xvals, yvals]

def plotStations():

    openedInFile = open('/Users/mariuskriegerowski/hiwi/MineDemo/Stations.dat', 'r')

    inl = openedInFile.readlines()

    stats=[]
    lat=[]
    lon=[]
    depth=[]
    relLat=[]
    relLon=[]

    for line in inl:
        tmp=line.split()
        stats.append(tmp[0])
        lat.append(tmp[1])
        lon.append(tmp[2])
        depth.append(tmp[3])
        relLat.append(tmp[4])
        relLon.append(tmp[5])

        statPlot = plot(relLat, relLon, 'o')    
        annotate(tmp[0],xy=(float(tmp[4]),float(tmp[5])+0.05))
    return statPlot

    


def chopAllTraces(traces, tmin, dt):
    ''' Chop all traces and return chopped list '''
    outTraces=[]
    [outTraces.append(trac.chop(tmin, tmin+dt, inplace=False)) for trac in traces]
    return outTraces


def getPairStationTool(InputTraces):
    '''Sort out station pairs
    :param list InputTraces: Traces, each instance of pyrocko.trace
    :param return: List contains Sets contains stationpairs.
    '''
    outTraces = []
    doit = True
    while doit:
        try: 
            trac1 = InputTraces.pop()
            for val in InputTraces.__reversed__():
                if val.station == trac1.station:
                    trac2 = InputTraces.pop(InputTraces.index(val))
            outTraces.append((trac1,trac2)) 
        except IndexError:
            return outTraces

def getTracesCorrelation(trac1, trac2):
    '''
    :param trac1: a pyrocko.trace
    :param trac2: a pyrocko.trace
    :param return: zero lag cross correlation value
    '''
    return np.correlate(trac1.get_ydata(), trac2.get_ydata())

def RotateAllTraces(degrange, Traces):
    '''
    :param degrange:    list of degree values to rotate
    :param Traces:      list of Traces to rotate
    :return:            dict containing degrees as keys and 
                                        rotatedTraces
    '''
    rotTrac={}
    for deg in degrange:
        rotatedTraces = trace.rotate(TestTraces, deg, ['n', 'e'], ['p1rot', 'p2rot']) 
        rotTrac[deg]=rotatedTraces
    return rotTrac


'''
pyrocko.trace documentation bei rotate ist nicht klar, dass inchannels und 
outchannels Listen sein muessen
'''
TestTraces = chopAllTraces(allTracesInput, allTracesInput[0].tmin+4.7,1)
degrees = range(-180,180,5)
AllRotatedTraces = RotateAllTraces(degrees, TestTraces)
correlationStatVal={}

for deg in degrees:
    for stationPair in getPairStationTool(AllRotatedTraces[deg]):
        try:
            correlationStatVal[stationPair[0].station].append([deg, getTracesCorrelation(stationPair[0], stationPair[1])])
        except KeyError:
            correlationStatVal[stationPair[0].station]=[]

for key in correlationStatVal.keys():

    tmpDeg=[]; tmpVal=[]
    oldmaxval = 0
    for line in correlationStatVal[key]:
        tmpDeg.append(line[0]); tmpVal.append(line[1])
        maxval = max(tmpVal) 

        if maxval > oldmaxval:
            merke = line[0]
            oldmaxval = maxval
    print '{0} {1} {2}'.format(key, merke+45, merke+180+45)
    #figure()
    #plot(tmpDeg, tmpVal)
    #title(key)
    #show()
statPlot = plotStations()
show()

sys.exit(0)
