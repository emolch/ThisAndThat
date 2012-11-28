#!/usr/bin/env python

from pyrocko import io
import glob
import sys, os
import progressbar
import getopt

short_options = 'hnslci:'
long_options = ['help', 'network=', 'station=', 'location=', 'channel=', 'inplace=']

def usage():
    sys.exit(('''USAGE: {0} [--network=] [--station=] [--location=] [--channel=] [--inplace=(True/False)] FILE1 FILE2 .....

The new code needs to follow each individual argument.

WARNING:
    Changing the NSLC-code may have unwanted consequences as for example
    undistinguishable traces of one station after components have falsely been 
    set the same.

    Hence, a first tryout without --inplace=True is highly recommended!


EXAMPLE: 
    ./reset_NSLC_code --network=IV --station=RAS Demodataset.mseed

(marius.kriegerowski(at)zmaw.de) ''').format(sys.argv[0]))

opts=[]
args=[]

try:
    opts, args = getopt.getopt(sys.argv[1:], short_options, long_options)
except getopt.GetoptError:
    print "ERR: at least one of these options is not available"
    usage()
    sys.exit()

if opts==[]:
    usage()

new_network = None
new_station = None
new_location= None
new_channel = None
inplace = False

for o, a in opts:
    if o in ("--help", "-h"):
        print "HELP"
        usage()
    elif o in ("--network", "-n"):
        print "set network to:", a
        new_network = a
    elif o in ("--station", "-s"):
        print "set station to:", a
        new_station = a
    elif o in ("--location", "-l"):
        print "set location to:", a
        new_location = a
    elif o in ("--channel", "-c"):
        print "set channel to:", a
        new_channel = a
    elif o in ("--inplace", "-i"):
        print "inplace: ", a 
        inplace = bool(a)

for directory in args:
    all_files = glob.glob(directory)

    for file_name in all_files:
        try:
            print '(re-)setting NSLC for: {0}'.format(file_name)
            loaded_file = io.load(file_name, getdata=True)
            for trac in loaded_file:
                if new_channel is not None:
                    trac.set_channel(new_channel)
                if new_network is not None:
                   trac.set_network(new_network)
                if new_station is not None:
                   trac.set_station(new_station)
                if new_location is not None:
                   trac.set_location(new_location)
            if inplace is False:
                file_name = file_name+'.renamed'
            elif inplace is True:
                file_name = file_name

            io.save(loaded_file, file_name)
            os.chmod(file_name, 0775)
        except:
            print 'Error while reading/writing {0}'.format(file_name)
            pass
        


