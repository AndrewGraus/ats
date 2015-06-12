import sys,os
import h5py
import numpy as np
import argparse

def get_keys(dat, time_range=None):
    """Collect all time index keys, optionally only those in a given 2-tuple of start and end time."""
    keys = dat[dat.keys()[0]].keys()
    keys.sort(lambda a,b: int.__cmp__(int(a),int(b)))

    if time_range is not None:
        all_times = get_times(dat, keys)
        keys = [key for key,time in zip(keys,all_times) if time_range[0] <= time <= time_range[1]]
    return keys

def get_times(dat, keys=None):
    """Get the times in the file, optionally at a given list of time index keys."""
    a_field = dat.keys()[0]

    if keys is None:
        keys = dat[dat.keys()[0]].keys()
        keys.sort(lambda a,b: int.__cmp__(int(a),int(b)))

    times = [dat[a_field][key].attrs['Time'] for key in keys]
    return times

def get_keys_and_times(dat, time_range=None):
    """Get the time index keys and times in a file, optionally only within a range of (start,end) times."""
    if time_range is not None:
        all_keys, all_times = get_keys_and_times(dat, None)
        keys = [key for key,time in zip(all_keys,all_times) if time_range[0] <= time <= time_range[1]]
        times = [time for time in all_times if time_range[0] <= time <= time_range[1]]        
    else:
        keys = dat[dat.keys()[0]].keys()
        keys.sort(lambda a,b: int.__cmp__(int(a),int(b)))
        times = get_times(dat, keys)

    return keys, times

    

def readATS(directory='.', base="visdump_data.h5", inds=None, time_range=None):
    """Main access point to data."""
    dat = h5py.File(os.path.join(directory,base),'r')

    if inds is None:
        keys, times = get_keys_and_times(dat, time_range)
        return keys, times, dat
    else:
        keys = get_keys(dat, time_range)
        keys = [keys[ind] for ind in inds]
        times = get_times(dat, keys)
        return keys, times, dat

def getSurfaceData(keys, dat, name):
    if not name.endswith(".cell.0"):
        name = name + ".cell.0"
    return np.array([dat[name][key][0] for key in keys])

def subsetFile(directory=".", base="visdump_data.h5", outfile="my_visdump_data.h5", inds=None, interval=1, time_range=None, names=None):
    """Read one file, write another"""
    keys, times, dat = readATS(directory,base,inds,time_range)
    if interval > 1:
        keys = keys[::interval]
        times = times[::interval]
    
    print "Transfering %d times to %s"%(len(times),outfile)
    
    if names is None:
        names = dat.keys()

    
    out = h5py.File(outfile)
    for name in names:
        grp = out.create_group(name)
        for key in keys:
            grp.create_dataset(key, data=dat[name][key][:])

    out.create_dataset('Times', data=np.array(times))
    out.close()


    
