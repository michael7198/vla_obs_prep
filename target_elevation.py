###make plots of elevation as function of LST

import numpy as np, matplotlib.pyplot as plt, argparse
from distutils.util import strtobool
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from astropy import units as u
from astropy.utils.iers import conf
conf.auto_max_age = None

plt.interactive(1)

##########################################################
##########################################################
####parameters



##########################################################
##########################################################
####functions

def extract_targetposdict(data, namecol='Name', acol='RA', dcol='DEC'):
    'produce a dictionary of target names and positions to use with elevation code'
    targets = {}
    for row in data:
        name = row[namecol]
        ra = row[acol]
        dec = row[dcol]
        targets[name] = SkyCoord(ra=ra, dec=dec, unit='deg')
        
    return targets
    

def make_time_list(location,
                   date=None,
                   start='00:00:00',
                   end='00:00:00',
                   step=15*u.min):
    'make list of time objects to use'
    if date is None:
        date = Time.now().to_value('iso').split(' ')[0]
    t_start = Time(f'{date} {start}', format='iso',
                   scale='utc', location=location)
    if start != end:
        t_end = Time(f'{date} {end}', format='iso',
                     scale='utc', location=location)
    else:
        t_end = t_start + 24*u.hour
    
    tlist = [t_start]
    t1 = t_start
    while t1 <= t_end:
        t1 = t1 + step
        tlist.append(t1)
    
    return tlist


def elevation_from_time(coord, obstime):
    'obtain the elevation of target based on observing time and sky coord'
    pos_altaz = coord.transform_to(AltAz(obstime=obstime.to_value('iso'),
                                         location=obstime.location))
    
    return pos_altaz.alt


def lst_from_utc(obstime, kind='apparent'):
    'get local sidereal time from utc for observatory'
    
    lst = obstime.sidereal_time(kind=kind)
    
    return lst


def plot_alt_v_time(time, elevation, label=None,
                    figsize=(8, 6), fontsize=14,
                    xtix=np.arange(0, 25, 2),
                    ytix=np.arange(0, 100, 15),
                    lincol='C0', lw=2, ls='-',
                    warncolors=['C1', 'C1', 'C3', 'C3'],
                    warnx=[[0,24], [0, 24], [0,24], [0,24]],
                    warny=[[8,15], [80, 85], [0,8], [85,90]],
                    warnalpha=0.3, plot_grid=True,
                    xlab='LST [hours]',
                    ylab='Alt [deg]'):
    'plot time vs elevation plot'
    
    ###link time and elevation in table and sort so no wrapping on plot

    
    plt.figure(figsize=figsize)
    for i in range(len(warnx)):
        plt.fill_between(warnx[i], warny[i][0], warny[i][1],
                         color=warncolors[i], alpha=warnalpha)
    if args.multiplot:
        for i in range(len(time)):
            ttab = Table({'time': time[i], 'alt': elevation[i]})
            ttab.sort('time')
            time2 = np.array(ttab['time'])
            elevation2 = np.array(ttab['alt'])
            plt.plot(time2, elevation2, lw=lw, ls=ls, label=label[i].split('on')[0].strip())
        plt.legend()
    else:
        ttab = Table({'time': time, 'alt': elevation})
        ttab.sort('time')
        time = np.array(ttab['time'])
        elevation = np.array(ttab['alt'])
        plt.plot(time, elevation, c=lincol, lw=lw, ls=ls)
    plt.xlim(0, 24)
    plt.ylim(0, 90)
    plt.xlabel(xlab, fontsize=fontsize)
    plt.ylabel(ylab, fontsize=fontsize)
    plt.xticks(xtix)
    plt.yticks(ytix)
    
    if plot_grid==True:
        plt.grid(ls=':')
        
    if not args.multiplot:
        plt.title(label, fontsize=fontsize+2)
    
    return



def parse_args():
    'parse command line arguments'
    parser = argparse.ArgumentParser(description="make plots of elevation vs LST at the VLA on a given data for a list of targets")
    parser.add_argument("targets",
                        help="table file containing sky coords of targets")
    parser.add_argument("obsdate",
                        help="date of observations (yyyy-mm-dd)")
    parser.add_argument("--obsloc", action='store',
                        type=str, default='vla',
                        help="location of observations (yyyy-mm-dd)")
    parser.add_argument("--namecol", action='store',
                        type=str, default='Name',
                        help="Object name column name in targets")
    parser.add_argument("--RAcol", action='store',
                        type=str, default='RA',
                        help="RA column name in targets")
    parser.add_argument("--DEcol", action='store',
                        type=str, default='DEC',
                        help="Decl. column name in targets")
    parser.add_argument("--outdir", action='store',
                        type=str, default='.',
                        help="directory to write files to")
    parser.add_argument("--use_old_IERS", action='store',
                        type=str, default='False',
                        help="use old IERS predictions if IERS server down")
    parser.add_argument("--multiplot", action='store', type=str, default='False', help='plot all targets on one graph?')
    
    args = parser.parse_args()
    
    ###change datatypes where appropriate
    args.obsloc = EarthLocation.of_site(args.obsloc)
    args.use_old_IERS = bool(strtobool(args.use_old_IERS))
    args.outdir = args.outdir.removesuffix('/')
    args.multiplot = bool(strtobool(args.multiplot))
    return args




##########################################################
##########################################################
####main

    
if __name__ == '__main__':
    args = parse_args()
    data = Table.read(args.targets)
    targets = extract_targetposdict(data=data, namecol=args.namecol,
                                    acol=args.RAcol, dcol=args.DEcol)
    if args.use_old_IERS == True:
        conf.auto_max_age = None
    ###iterate through targets
    lst_list = []
    elevs_list = []
    plot_label = []
    for tname in list(targets.keys()):
        pos = targets[tname]
        times = make_time_list(location=args.obsloc,
                               date=args.obsdate)
        lst_list.append([lst_from_utc(obstime=t).value for t in times])
        elevs_list.append([elevation_from_time(coord=pos, obstime=t).value for t in times])
        plot_label.append(f'{tname}  on  {args.obsdate}')
    if args.multiplot:
        plot_alt_v_time(time=lst_list, elevation=elevs_list,
                        label=plot_label)
        plt.savefig(f'{args.outdir}/lst_v_alt_all.png', dpi=100)
        plt.close()
    else:
        for ind, tname in enumerate(list(targets.keys())):
            plot_alt_v_time(time=lst_list[ind], elevation=elevs_list[ind],
                        label=plot_label[ind])
            plt.savefig(f'{args.outdir}/lst_v_alt_{tname}.png', dpi=100)
            plt.close()


