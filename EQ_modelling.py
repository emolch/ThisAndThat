from pyrocko.snuffling import Param, Snuffling, Switch, pile, Choice
from pyrocko import util, pile_viewer
from tunguska import gfdb, receiver, seismosizer, source, orthodrome 
from tunguska.phase import Taper, Timing
from numpy import array, pi, savetxt, complex, savetxt, prod, ones
from math import acos, cos, asin, sin, pi, radians, degrees

def lat_lon_from_dist_azi(olat, olon, dist, azim):
    '''
    :param olat: Latitude of origin, type: float
    :param olon: Longitude of otigin, type: float
    :param dist: Distance from origin, type: float
    :param azim: Azimuth from origin, type: float
    :return lat, lon: new latitude and longitude 
    '''
    rE=6371000.785
    olat=radians(olat); olon=radians(olon); azi=radians(azim)
    b=dist/rE
    a=acos(cos(b)*cos(pi/2-olat)+sin(pi/2-olat)*sin(b)*cos(azi))
    B=asin(sin(b)*sin(azi)/sin(a))
    return 90-degrees(a), degrees(B)+degrees(olon)


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
        self.set_name('EQmodelling')
       
        #self.db = gfdb.Gfdb('/scratch/local1/auto_Rotate_Test/gfdb_sediment/db')
        self.db = gfdb.Gfdb('/scratch/local1/gfdb_sediment_cut/db')
        
        # Add scrollbars of the parameters that you desire to adjust.
        # 1st argument: Description that appears within the snuffling.
        # 2nd argument: Name of parameter as used in the following code.
        # 3rd-5th argument: default, start, stop.
        gfdb_maxrange=self.db.nx*self.db.dx-self.db.dx
        gfdb_minrange=self.db.firstx

        self.add_parameter(Param('Source Depth [m]','source_depth', self.db.firstz, self.db.firstz, self.db.nz*self.db.dz-self.db.dz))
        self.add_parameter(Param('Distance [m]','dist', gfdb_maxrange/2, gfdb_minrange, gfdb_maxrange))
        self.add_parameter(Param('Azimuth [deg]','azi', 0.,-180., 180.))
        self.add_parameter(Param('Mxx=Myy=Mzz [Nm]','mxx', 1.,-1.,1.))
        self.add_parameter(Param('Mxy [Nm]','mxy', 0.,-1.,1.))
        self.add_parameter(Param('Myz [Nm]','myz', 0.,-1.,1.))
        self.add_parameter(Param('Mxz [Nm]','mxz', 0.,-1.,1.))
        self.add_parameter(Param('Time fade [s]','tfade', 5.,0.,15))
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

        origin_lat, origin_lon = lat_lon_from_dist_azi(float(lat), float(lon), self.dist, self.azi)
        r = receiver.Receiver(lat,lon, components='neu', name='.%s.' % station)
        receivers.append(r)

        # Composition of the source
        otime = util.str_to_time('2000-01-1 00:00:00')
        db = self.db

        seis = seismosizer.Seismosizer(hosts=['localhost'])
        seis.set_database(db)
        seis.set_effective_dt(db.dt)
        seis.set_local_interpolation('bilinear')
        seis.set_receivers(receivers)
        seis.set_source_location( origin_lat, origin_lon, otime)
        seis.set_source_constraints (0, 0, 0, 0 ,0 ,-1)
        self.seis = seis        
        seis = None

        risetime=3
        mxx=self.mxx
        myy=self.mxx
        mzz=self.mxx
        mxy=self.mxy
        mxz=self.mxz
        myz=self.myz
        Sourceparams = dict(zip(['mxx', 'myy', 'mzz', 'mxy', 'mxz', 'myz', 'depth', 'rise-time'],
            [mxx, myy, mzz, mxy, mxz, myz,  self.source_depth, risetime]))
        s = source.Source(sourcetype='moment_tensor', sourceparams=Sourceparams)
        self.seis.set_source(s)
        recs = self.seis.get_receivers_snapshot( which_seismograms = ('syn',), which_spectra=(), which_processing='tapered')
        
        trs = []
        for rec in recs:
            if self.save_mseed is True:
                rec.save_traces_mseed(filename_tmpl='$HOME/.snufflings/%(whichset)s_%(network)s_%(station)s_%(location)s_%(channel)s.mseed' )
            trs.extend(rec.get_traces())
        self.add_traces(trs)
        # Define fade in and out, band pass filter and cut off fader for the TF.
        tfade = self.tfade
        freqlimit = (0.005,.006,1,1.2)
        cut_off_fading = 500
        ntraces = []
        
        for tr in trs:
            TF = STS2()
            
            # Save synthetic trace after transfer function was applied.
            trace_filtered = tr.transfer(tfade, freqlimit, TF, cut_off_fading)            
            # Set new codes to the filtered trace to make it identifiable.
            rename={'e':'BHE','n':'BHN','u':'BHZ'}
            trace_filtered.set_codes(channel=rename[trace_filtered.channel], network='STS2', station='HH', location='syn')
            ntraces.append(trace_filtered)            
        self.add_traces(ntraces)        

def __snufflings__():
    return [ ParaEditCp_TF_GTTG() ]


