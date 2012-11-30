from pyrocko.snuffling import Param, Snuffling, Switch, pile, Choice
from pyrocko import util, pile_viewer
from tunguska import gfdb, receiver, seismosizer, source, orthodrome 
from tunguska.phase import Taper, Timing
from numpy import array, pi, savetxt, complex, savetxt, prod, ones

class STS2:
    ''' Apply the STS2's transfer function which is deduced from the 
poles, zeros and gain of the transfer tunction. The Green's function database (gdfb) which is required for synthetic seismograms and the rake of the focal mechanism can be chosen and changed within snuffler.
Two gfdbs are needed.
Three synthetic seismograms of an STS2 seismometer will be the result
'''

    def evaluate(self,freqs):
        
        # transform the frequency to angular frequency.
        w = 2j*pi*freqs

        Poles = array([-3.7e-2+3.7e-2j, -3.7e-2-3.7e-2j, 
                       -2.51e2, -1.31e2+4.67e2j, -1.31e2-4.67e2])
        Zeros = array([0, 0, 0])
        K = 6.16817e7

        # Multiply factored polynomials of the transfer function's numerator
        # and denominator.
        a = ones(freqs.size,dtype=complex)*K
        for i_z in Zeros:
            a *= w-i_z
        for i_p in Poles:
            a /= w-i_p
        return a

class extendedSnuffling(Snuffling):
    def __init__(self):
        Snuffling.__init__(self)

    def pre_destroy(self):
        '''overloaded here'''
        self.cleanup()
        if self._tempdir is not None:
            import shutil
            shutil.rmtree(self._tempdir)  
        try:
            self.mydel()
        except AttributeError:
            pass

class ParaEditCp_TF_GTTG(extendedSnuffling):

    def mydel(self):
        ''' Terminates seismosizer'''
        if self.seis is not None:
            self.seis.close()

    def setup(self):
        
        # Give the snuffling a name:
        self.set_name('EQmodelling')
        
        # Add scrollbars of the parameters that you desire to adjust.
        # 1st argument: Description that appears within the snuffling.
        # 2nd argument: Name of parameter as used in the following code.
        # 3rd-5th argument: default, start, stop.
        self.add_parameter(Param('Source Depth [km]','source_depth', 10., 1., 100.))
        self.add_parameter(Param('Distance North [km]','d_north', 10., 100., 500.))
        self.add_parameter(Param('Distance East[km]','d_east', 10., 100., 500.))

        self.set_live_update(False)


    def call(self):
         
        self.cleanup()
        
        # Set up receiver configuration.
        tab = '''
        HH  3. 3. 0
        '''.strip()

        receivers = []
        station, lat, lon, depth = tab.split()

        d_north=self.d_north*1000
        d_east=self.d_east*1000
        origin_lat, origin_lon = orthodrome.ne_to_latlon_alternative_method(float(lat), float(lon), d_north, d_east)
        r = receiver.Receiver(lat,lon, components='neu', name='.%s.' % station)
        receivers.append(r)

        # Composition of the source
        otime = util.str_to_time('2000-01-1 00:00:00')
        
        db = gfdb.Gfdb('/media/exupery2/gemini-iasp91-20000km/db')

        seis = seismosizer.Seismosizer(hosts=['localhost'])
        seis.set_database(db)
        seis.set_effective_dt(db.dt)
        seis.set_local_interpolation('bilinear')
        seis.set_receivers(receivers)
        seis.set_source_location( origin_lat, origin_lon, otime)
        seis.set_source_constraints (0, 0, 0, 0 ,0 ,-1)
        self.seis = seis        
        seis = None
        # Change strike within snuffler with the added scroll bar.
        #strike = self.strike

        # Other focal mechism parameters are constants
        #dip = 90; strike=0; rake = 0; moment = 7.00e20; source_depth = self.source_depth*1000; 

        risetime=3
        mxx=0.
        myy=0.
        mzz=0.
        mxy=0.
        mxz=-1.
        myz=0.
        Sourceparams = dict(zip(['mxx', 'myy', 'mzz', 'mxy', 'mxz', 'myz', 'depth'],[mxx, myy, mzz, mxy, mxz, myz,  self.source_depth*1000]))
        s = source.Source(sourcetype='moment_tensor',
        sourceparams=Sourceparams)
        #print source.source_infos('moment_tensor')
        self.seis.set_source(s)
        recs = self.seis.get_receivers_snapshot( which_seismograms = ('syn',), which_spectra=(), which_processing='tapered')
        
        trs = []
        for rec in recs:
            rec.save_traces_mseed(filename_tmpl='%(whichset)s_%(network)s_%(station)s_%(location)s_%(channel)s.mseed' )
            trs.extend(rec.get_traces())
        
        # Define fade in and out, band pass filter and cut off fader for the TF.
        tfade = 0
        freqlimit = (0.005,.006,1,1.2)
        cut_off_fading = 50
        ntraces = []
        
        for tr in trs:
            TF = STS2()
            
            # Save synthetic trace after transfer function was applied.
            trace_filtered = tr.transfer(tfade, freqlimit, TF, cut_off_fading)            
            # Set new codes to the filtered trace to make it identifiable.
            rename={'e':'BHE','n':'BHN','u':'BHZ'}
            trace_filtered.set_codes(channel=rename[trace_filtered.channel], network='', station='HHHA', location='syn')
            ntraces.append(trace_filtered)            
            
        self.add_traces(ntraces)        

def __snufflings__():
    return [ ParaEditCp_TF_GTTG() ]


