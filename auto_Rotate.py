from pyrocko import pile, trace, util, io, pile_viewer
from pyrocko.snuffling import Param, Snuffling, Switch
import numpy as np

infileDir = "/Users/mariuskriegerowski/hiwi/MineDemo/Demodataset_new.mseed"

allTracesInput = io.load(infileDir)

for trac in allTracesInput:
    print trac.channel

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
        rotatedTraces = trace.rotate(Traces, deg, ['n', 'e'], ['p1rot', 'p2rot']) 
        rotTrac[deg]=rotatedTraces
    return rotTrac

def printFormattedOutput(rotdict):
    for key in rotdict.iterkeys():
        print '\n', key
        for value in rotdict[key]:
            print value  



'''
pyrocko.trace documentation bei rotate ist nicht klar, dass inchannels und 
outchannels Listen sein muessen
'''

class autoRotate(Snuffling):

    def setup(self):
        self.set_name('autoRotate')
        self.set_live_update(False)
        self.add_parameter(Param('T_0', 't_0', 0., 0., 100.))

    def call(self):
        
        self.cleanup()
        currentPile = self.get_pile()
        currentPileViewer = self.get_viewer()
        markers = currentPileViewer.selected_markers()

        tminpick=[]; tmaxpick=[]
        for phase_marker in markers:
            tminpick.append(phase_marker.get_tmin())
            tmaxpick.append(phase_marker.get_tmax())
       
        choppedTraces=currentPile.chopper(tmin=min(tminpick), tmax=min(tmaxpick) )
        
        degrees = range(-180,180,5)
        for traces in choppedTraces:
            AllRotatedTraces = RotateAllTraces(degrees, traces)
            correlationStatVal={}

            for deg in degrees:
                for stationPair in getPairStationTool(AllRotatedTraces[deg]):
                    try:
                        correlationStatVal[stationPair[0].station].append([deg, getTracesCorrelation(stationPair[0], stationPair[1])])
                    except KeyError:
                        correlationStatVal[stationPair[0].station]=[]
            
            
            for trace in traces:
                trace.set_codes(network='rot')

            printFormattedOutput(correlationStatVal)        
        self.add_traces(traces)     
        


def __snufflings__():
    return [ autoRotate() ]
