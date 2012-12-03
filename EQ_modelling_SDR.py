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
        Zeros = array([0, 0])
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
        self.set_name('EQmodelling_SDR')
       
        self.db = gfdb.Gfdb('/scratch/local1/auto_Rotate_Test/gfdb_sediment/db')
        
        # Add scrollbars of the parameters that you desire to adjust.
        # 1st argument: Description that appears within the snuffling.
        # 2nd argument: Name of parameter as used in the following code.
        # 3rd-5th argument: default, start, stop.
        self.add_parameter(Param('Source Depth [m]','source_depth', self.db.firstz, self.db.firstz, self.db.nz*self.db.dz-self.db.dz))
        self.add_parameter(Param('Distance North [m]','d_north', self.db.firstx, self.db.firstx, self.db.nx*self.db.dx-self.db.dx))
        self.add_parameter(Param('Distance East [m]','d_east', self.db.firstx, self.db.firstx, self.db.nx*self.db.dx-self.db.dx))
        self.add_parameter(Param('Strike [deg]','strike', 0., 0.,360.))
        self.add_parameter(Param('Dip [deg]','dip', 0., 0., 90.))
        self.add_parameter(Param('Rake [deg]','rake', 0., -180.,180.))
        self.add_parameter(Param('Time fade [s]','tfade', 5.,2.,15))
        self.add_parameter(Switch('Save miniseed traces in .snufflings?','save_mseed', False))
        
        self.set_live_update(False)


    def call(self):
         
        self.cleanup()
        
        # Set up receiver configuration.
        tab = '''
        HH  3. 3. 0
        '''.strip()

        receivers = []
        station, lat, lon, depth = tab.split()

        d_north=self.d_north
        d_east=self.d_east
        origin_lat, origin_lon = orthodrome.ne_to_latlon_alternative_method(float(lat), float(lon), d_north, d_east)
        r = receiver.Receiver(lat,lon, components='neu', name='.%s.' % station)
        receivers.append(r)

        # Composition of the source
        otime = util.str_to_time('2000-01-1 00:00:00')
        db = self.db
        #db = gfdb.Gfdb('/media/exupery2/gemini-iasp91-20000km/db')

        seis = seismosizer.Seismosizer(hosts=['localhost'])
        seis.set_database(db)
        seis.set_effective_dt(db.dt)
        seis.set_local_interpolation('bilinear')
        seis.set_receivers(receivers)
        seis.set_source_location( origin_lat, origin_lon, otime)
        seis.set_source_constraints (0, 0, 0, 0 ,0 ,-1)
        self.seis = seis        
        seis = None

        risetime=3; moment=1.
        s = source.Source('bilateral',
        sourceparams_str='0 0 0 %g %g %g %g %g 0 0 0 0 1 %g' % (self.source_depth, moment, self.strike, self.dip, self.rake, risetime))
        self.seis.set_source(s)
        recs = self.seis.get_receivers_snapshot( which_seismograms = ('syn',), which_spectra=(), which_processing='tapered')
        
        trs = []
        for rec in recs:
            if self.save_mseed is True:
                rec.save_traces_mseed(filename_tmpl='%(whichset)s_%(network)s_%(station)s_%(location)s_%(channel)s.mseed' )
            trs.extend(rec.get_traces())
        self.add_traces(trs)
        # Define fade in and out, band pass filter and cut off fader for the TF.
        tfade = self.tfade
        freqlimit = (0.005,.006,1,1.2)
        cut_off_fading = 50
        ntraces = []
        
        for tr in trs:
            TF = STS2()
            
            # Save synthetic trace after transfer function was applied.
            trace_filtered = tr.transfer(tfade, freqlimit, TF, cut_off_fading)            
            # Set new codes to the filtered trace to make it identifiable.
            rename={'e':'BHE','n':'BHN','u':'BHZ'}
            trace_filtered.set_codes(channel=rename[trace_filtered.channel], network='STS2', station='HH', location='syn')
            ntraces.append(trace_filtered)            
        #self.add_traces(ntraces)        

def __snufflings__():
    return [ ParaEditCp_TF_GTTG() ]


