from pyrocko.snuffling import Param, Snuffling, Switch, pile, Choice
from pyrocko import io, util, model, trace
import pickle
import urllib
try:
    from obspy.core import UTCDateTime, stream
    from obspy.core import trace as obspy_trace

    from obspy.signal import cornFreq2Paz
    from obspy.signal.array_analysis import sonic
    try:
        from obspy.signal import array_analysis
    except:
        print 'using deprecated obspy.signal.sonic'
        pass

except ImportError:
    print 'please note: cannot import each required obspy module needed for fk_analysis snuffling'
    raise

from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np


def p2o_trace(ptrace, station):
    '''Convert Pyrocko trace to ObsPy trace.'''

    otr = obspy_trace.Trace(
            data = ptrace.get_ydata(),
            header=dict(
                network = ptrace.network,
                station = ptrace.station,
                location = ptrace.location,
                channel = ptrace.channel,
                delta = ptrace.deltat,
                starttime = UTCDateTime(ptrace.tmin),
                coordinates = dict(
                    latitude = station.lat,
                    longitude = station.lon,
                    elevation = station.elevation/1000. )))

    return otr

class FK(Snuffling):
    
    def setup(self):
        self.set_name('FK Analysis')
        self.add_parameter(Param('Number of Steps','numberOfFraction',32,4,40))
        self.set_live_update(False)

    def call(self):
        self.cleanup()
        viewer = self.get_viewer()

        if viewer.lowpass is None or viewer.highpass is None:
            self.fail('highpass and lowpass in viewer must be set!')

        traces = []
        for trs in self.chopper_selected_traces(fallback=True):
            for tr in trs:
                tr.downsample_to(1/20.)
                tr.lowpass(4, viewer.lowpass)
                tr.highpass(4, viewer.highpass)

            traces.extend(trs)

        if not traces:
            self.fail('no traces selected')

        tmin, tmax = trace.minmaxtime(traces, key=lambda x: None)[None]

        try:
            obspy_traces = [ p2o_trace(tr, viewer.get_station(viewer.station_key(tr)) ) for tr in traces ]
        except KeyError:
            self.fail('station information missing')

        st = stream.Stream(traces=obspy_traces)

        smax = 0.5

        # Execute sonic
        kwargs = dict(
            # slowness grid: X min, X max, Y min, Y max, Slow Step
            sll_x=-smax, slm_x=smax, sll_y=-smax, slm_y=smax, sl_s=smax/20.,
            # sliding window properties
            win_len=5.0, win_frac=0.1,
            # frequency properties
            frqlow=viewer.highpass, frqhigh=viewer.lowpass, prewhiten=0,
            # restrict output
            semb_thres=-1.0e9, vel_thres=-1.0e9, verbose=True, timestamp='mlabday',
            stime=UTCDateTime(tmin), etime=UTCDateTime(tmax)
        )
        #out = sonic(st, **kwargs)
        out = array_analysis.array_processing(st, **kwargs)

        cmap = cm.hot_r
        pi = np.pi

        #
        # make output human readable, adjust backazimuth to values between 0 and 360
        t, rel_power, abs_power, baz, slow = out.T
        baz[baz < 0.0] += 360.

        # choose number of fractions in plot (desirably 360 degree/N is an integer!)
        N = int(self.numberOfFraction)
        abins = np.arange(N + 1) * 360. / N
        sbins = np.linspace(0., smax, N + 1)

        # sum rel power in bins given by abins and sbins
        hist, baz_edges, sl_edges = np.histogram2d(baz, slow,
                bins=[abins, sbins], weights=rel_power)

        # transform to gradient
        baz_edges = baz_edges / 180 * np.pi

        # add polar and colorbar axes
        fig = self.pylab(get='figure')
        cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
        ax = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)

        dh = abs(sl_edges[1] - sl_edges[0])
        dw = abs(baz_edges[1] - baz_edges[0])

        # circle through backazimuth
        for i, row in enumerate(hist):
            bars = ax.bar(left=(pi / 2 - (i + 1) * dw) * np.ones(N),
                          height=dh * np.ones(N),
                          width=dw, bottom=dh * np.arange(N),
                          color=cmap(row / hist.max()))
                            

        ax.set_xticks([pi / 2, 0, 3. / 2 * pi, pi])
        ax.set_xticklabels(['N', 'E', 'S', 'W'])

        # set slowness limits
        ax.set_ylim(0., smax)
        ColorbarBase(cax, cmap=cmap,
                     norm=Normalize(vmin=hist.min(), vmax=hist.max()))

def __snufflings__():
    return [ FK() ]
        
