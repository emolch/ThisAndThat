from pyrocko import trace, util, io
import numpy
import time
import glob

traces=([])
t_start = util.str_to_time('2012-01-01 00:00:00')

infiles = glob.glob("HM*")

for i in infiles:
    Y_data=[]
    X_data=[]
    print i[5]
    
    file = open(i,"r")#.readline()#split()
    for line in file:
        inf=line.split()
        X_data.append(float(inf[0]))
        Y_data.append(float(inf[1]))
    
    traces.append(trace.Trace(station=i[0:4], channel=i[5], tmin=util.str_to_time('2012-02-06 12:00:00.'), deltat=0.005, ydata=numpy.array(Y_data))) 

io.save(traces,'Demodataset.mseed')



