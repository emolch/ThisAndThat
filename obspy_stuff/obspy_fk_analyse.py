from obspy.core import UTCDateTime, read
from obspy.signal import cornFreq2Paz
from obspy.signal.array_analysis import *
import pickle
import urllib

# Load data
#st = pickle.load(urllib.urlopen("http://examples.obspy.org/agfa.dump"))

st = read('/home/zmaw/u254061/hiwi/MineDemo/Demodataset.mseed')
import pdb 
pdb.set_trace()
for sti in st:
    sti.stats.coordinates={'longitude':3.4, 'latitude':33.4, 'elevation':0}
st[2].stats.coordinates={'longitude':3.4, 'latitude':35.4, 'elevation':0}
st[1].stats.coordinates={'longitude':2.4, 'latitude':34.4, 'elevation':0}
st[0].stats.coordinates={'longitude':3.4, 'latitude':34.4, 'elevation':0}

# Instrument correction to 1Hz corner frequency
paz1hz = cornFreq2Paz(1.0, damp=0.707)
#st.simulate(paz_remove='self', paz_simulate=paz1hz)

# Execute sonic
kwargs = dict(
    # slowness grid: X min, X max, Y min, Y max, Slow Step
    sll_x=-3.0, slm_x=3.0, sll_y=-3.0, slm_y=3.0, sl_s=0.03,
    # sliding window properties
    win_len=1.0, win_frac=0.05,
    # frequency properties
    frqlow=1.0, frqhigh=8.0, prewhiten=0,
    # restrict output
    semb_thres=-1e9, vel_thres=-1e9, verbose=True, timestamp='mlabday',
    stime=UTCDateTime("20120206120000"), etime=UTCDateTime("20120206120006")
)
out = sonic(st, **kwargs)

# Plot
import matplotlib.pyplot as plt
labels = 'rel.power abs.power baz slow'.split()

fig = plt.figure()
for i, lab in enumerate(labels):
    ax = fig.add_subplot(4, 1, i + 1)
    ax.scatter(out[:, 0], out[:, i + 1], c=out[:, 1], alpha=0.6,
               edgecolors='none')
    ax.set_ylabel(lab)
    ax.xaxis_date()

fig.autofmt_xdate()
fig.subplots_adjust(top=0.95, right=0.95, bottom=0.2, hspace=0)
plt.show()
