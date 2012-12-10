#!/usr/bin/env python
import glob
from obspy.core import read
from obspy.core import UTCDateTime
from obspy.signal import cornFreq2Paz
from obspy.signal.array_analysis import sonic
import pickle
import urllib

st = pickle.load(urllib.urlopen("http://examples.obspy.org/agfa.dump"))

for file in glob.glob('*.z'):
    st = read(file)
    tr = st[0]
    msg = "%s %s %f %f" % (tr.stats.station, str(tr.stats.starttime),
                           tr.data.mean(), tr.data.std())
    print msg
