from pyrocko.snuffling import Param, Snuffling, Switch, pile, Choice
from pyrocko import io, util, model
import pickle
import urllib
try:
    from obspy.core import UTCDateTime, read, stream, trace
    from obspy.signal import cornFreq2Paz
    from obspy.signal.array_analysis import sonic
    try:
        from obspy.signal import array_analysis
    except:
        print 'using deprecated obspy.signal.sonic'
        pass
except ImportError:
    print 'please note: cannot import each required obspy module needed for fk_analysis snuffling'
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

def prepare_time_string(t):
    gm_time = util.time.gmtime(t)
    return '%s%s%s%s%s%s'%(str(gm_time.tm_year), 
            str(gm_time.tm_mon).zfill(2), 
            str(gm_time.tm_mday).zfill(2),
            str(gm_time.tm_hour).zfill(2), 
            str(gm_time.tm_min).zfill(2), 
            str(gm_time.tm_sec).zfill(2)) 

class fk(Snuffling):
    
    def setup(self):
        self.set_name('fk Analysis')
        self.add_parameter(Param('Number of Steps','numberOfFraction',4,4,40))
        self.set_live_update(False)

    def call(self):
        self.cleanup()

        # get time range visible in viewer
        self.viewer = self.get_viewer()
        self.tmin, self.tmax = self.viewer.get_time_range()
        print 'tmin: %s, tmax: %s'%(self.tmin, self.tmax)
        pile = self.get_pile()
        self.outTracs = [] 
        for trac in pile.chopper(tmin=self.tmin, tmax=self.tmax):
            self.outTracs.extend(trac)
        
        self.traces2analize = []
        for trac in self.outTracs:
            if trac.channel=='BHZ':
                self.traces2analize.append(trac)
                
        io.save(self.traces2analize, 'tmp_fkIn.mseed')

        self.start_fk_anlysis()


    def start_fk_anlysis(self):

        # todo ... zeitlimits
        tmax = prepare_time_string(self.tmax-1)
        tmin = prepare_time_string(self.tmin+1)

        st= read('tmp_fkIn.mseed')

        # tried to avoid file io
        #allTraces = []
        #for trac in self.outTracs:
        #    allTraces.append(trace.Trace(trac))
        #st = stream.Stream(traces=allTraces)

        # Set station coordinates
        for n,sti  in enumerate(st):
            lat, lon = self.viewer.station_latlon(self.traces2analize[n])
            sti.stats.coordinates={'longitude': lat, 'latitude': lon, 'elevation':0}

        # Instrument correction to 1Hz corner frequency
        paz1hz = cornFreq2Paz(1.0, damp=0.707)
        st.simulate(paz_remove='self', paz_simulate=paz1hz)

        # Execute sonic
        kwargs = dict(
            # slowness grid: X min, X max, Y min, Y max, Slow Step
            sll_x=-3.0, slm_x=3.0, sll_y=-3.0, slm_y=3.0, sl_s=0.03,
            # sliding window properties
            win_len=0.8, win_frac=0.05,
            # frequency properties
            frqlow=1.0, frqhigh=8.0, prewhiten=0,
            # restrict output
            semb_thres=-1e9, vel_thres=-1e9, verbose=True, timestamp='mlabday',
            stime=UTCDateTime(tmin), etime=UTCDateTime(tmax)
        )
        #out = sonic(st, **kwargs)
        out = array_analysis.array_processing(st, **kwargs)



        cmap = cm.hot_r
        pi = np.pi

        #
        # make output human readable, adjust backazimuth to values between 0 and 360
        t, rel_power, abs_power, baz, slow = out.T
        baz[baz < 0.0] += 360

        # choose number of fractions in plot (desirably 360 degree/N is an integer!)
        N = int(self.numberOfFraction)
        abins = np.arange(N + 1) * 360. / N
        sbins = np.linspace(0, 3, N + 1)

        # sum rel power in bins given by abins and sbins
        hist, baz_edges, sl_edges = np.histogram2d(baz, slow,
                bins=[abins, sbins], weights=rel_power)

        # transform to gradient
        baz_edges = baz_edges / 180 * np.pi

        # add polar and colorbar axes
        fig = plt.figure(figsize=(8, 8))
        cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
        ax = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)

        dh = abs(sl_edges[1] - sl_edges[0])
        dw = abs(baz_edges[1] - baz_edges[0])

        # circle through backazimuth
        for i, row in enumerate(hist):
            print hist.max()
            bars = ax.bar(left=(pi / 2 - (i + 1) * dw) * np.ones(N),
                          height=dh * np.ones(N),
                          width=dw, bottom=dh * np.arange(N),
                          color=cmap(row / hist.max()))
                            

        ax.set_xticks([pi / 2, 0, 3. / 2 * pi, pi])
        ax.set_xticklabels(['N', 'E', 'S', 'W'])

        # set slowness limits
        ax.set_ylim(0, 3)
        ColorbarBase(cax, cmap=cmap,
                     norm=Normalize(vmin=hist.min(), vmax=hist.max()))

        plt.show()

def __snufflings__():
    return [ fk() ]
        
