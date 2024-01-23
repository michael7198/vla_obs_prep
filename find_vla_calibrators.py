###parse VLA calibrator list to find suitable calibrators within given radius

import numpy as np, argparse
from distutils.util import strtobool
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack, vstack, join, unique
from astropy import units as u
from collections import Counter


#########################################################
#########################################################
###parameters

calfile = '../../VLA_collaborator_list.txt'
caldelim = '\n \n'


#target_pos = SkyCoord(ra=[64.1972], dec=[74.4827], unit='deg')
#target_file = ('../../../observing_proposals/VLA24B_lensed_radio_sources/'
#               + 'data/targets_v1_22Jan.csv')
target_file = ('../../../observing_proposals/VLA24B_lensed_radio_sources/'
               + 'data/target-scheduling.csv')

#########################################################
#########################################################
###functions

def parse_nrao_list(filename, delim='\n \n'):
    'parse NRAO list of VLA calisbrators and output list of sources + info'
    
    clist = open(filename).read().split(delim)
    clist = [clist[i].split('\n') for i in range(len(clist))]
    clist = clist[1:]
    
    return clist


def parse_source(source):
    'parse an individual source from calibrator list'
    
    ###data starts at row 5
    head = source[:5]
    data = source[5:]
    dname = f'J{head[0].split()[0]}'
    ra = head[0].split()[3]
    dec = head[0].split()[4]
    if dec[0]!='-':
        dec = '+' + dec
    co = SkyCoord(ra+dec)
    
    caldata = {}
    bands = []
    for d in data:
        row = d.split()
        bands.append(row[1])
        if len(row) > 6:
            caldata[row[1]] = {'A': row[2], 'B': row[3], 'C': row[4], 'D': row[5], 'Flux': row[6]}
            if len(row) > 7:
                caldata[row[1]]['UVmin max'] = ' ; '.join(row[7])


    cdict = {dname: caldata}
    sdict = {'Name': dname, 'RAJ2000': co.ra.value, 'DEJ2000': co.dec.value}
    
    ###add info on whether source is a calibrator for desired band
    calband = {'P': False,
               'L': False,
               'C': False,
               'X': False,
               'K': False,
               'Q': False,
               'U': False}
    
    for key in list(calband.keys()):
        if key in bands:
            calband[key] = True
        sdict[key] = calband[key]

    return cdict, sdict


def create_datasets_for_querying(callist):
    'output skycoord table and dict of multi level data'
    sources, caldict = [], []
    for cal in callist:
        cdata = parse_source(cal)
        caldict.append(cdata[0])
        sources.append(cdata[1])
    
    sources = Table(sources)
    caldict = {k: v for d in caldict for k, v in d.items()}
    
    return sources, caldict


def find_calibrators(targets, cal_sources, cal_infodict,
                     tname='Name', tra='ra', tdec='dec',
                     cname='Name', cra='RAJ2000', cdec='DEJ2000',
                     cband='L', worst_calib_level='S', limitoncalq=True,
                     configs=['A'], max_sep=10*u.deg,
                     posunits_t='deg', posunits_c='deg',
                     selected_configs_only=True):
    'input target list to identify all potential calibrators for targets'
    
    ###subset only those calibrators that are usable for the selected band
    suitable_calibrators = cal_sources[cal_sources[cband]==True]
    
    ###setup skycoords for matching
    tcat = SkyCoord(ra=targets[tra], dec=targets[tdec], unit=posunits_t)
    ccat = SkyCoord(ra=suitable_calibrators[cra],
                    dec=suitable_calibrators[cdec],
                    unit=posunits_c)
    
    ###search within radius
    searchres = tcat.search_around_sky(ccat, seplimit=max_sep)
    
    
    tlite = targets[[tname]]
    clite = suitable_calibrators[[cname]]
    tlite.rename_column(name=tname, new_name='target')
    clite.rename_column(name=cname, new_name='calibrator')
    
    targetcals = hstack([tlite[searchres[1]], clite[searchres[0]]])
    targetcals['sep'] = np.round(searchres[2], 4)
    
    ###add in calibration info on config of choice
    clist, acon, bcon, ccon, dcon, flux = [], [], [], [], [], []
    for c in list(np.unique(targetcals['calibrator'])):
        clist.append(c)
        cinfo = cal_infodict[c][cband]
        cikeys = list(cinfo.keys())
        if 'A' in cikeys:
            acon.append(cinfo['A'])
        else:
            acon.append('-')
        if 'B' in cikeys:
            bcon.append(cinfo['B'])
        else:
            bcon.append('-')
        if 'C' in cikeys:
            ccon.append(cinfo['C'])
        else:
            ccon.append('-')
        if 'D' in cikeys:
            dcon.append(cinfo['D'])
        else:
            dcon.append('-')
        if 'Flux' in cikeys:
            flux.append(np.round(np.float64(cinfo['Flux']), 2))
        else:
            flux.append(-99)
    
    calinfo = Table({'calibrator': clist, f'Aconfig_{cband}band': acon,
                     f'Bconfig_{cband}band': bcon,
                     f'Cconfig_{cband}band': ccon,
                     f'Dconfig_{cband}band': dcon,
                     'Flux': flux*u.Jy})
    
    ###lighten table
    if selected_configs_only==True:
        if 'A' not in configs:
            calinfo.remove_column(f'Aconfig_{cband}band')
        if 'B' not in configs:
            calinfo.remove_column(f'Bconfig_{cband}band')
        if 'C' not in configs:
            calinfo.remove_column(f'Cconfig_{cband}band')
        if 'D' not in configs:
            calinfo.remove_column(f'Dconfig_{cband}band')
        
    
    targetcals = join(targetcals, calinfo,
                      keys='calibrator',
                      join_type='left')
                      
    ###select calibrators of cetain quality, use highest reslotion config
    if limitoncalq==True:
        ##define best config for calibration
        bestconfig = 'A'
        if 'A' not in configs:
            bestconfig = 'B'
            if 'B' not in configs:
                bestconfig = 'C'
                if 'C' not in configs:
                    bestconfig = 'D'
        
        bccol = f'{bestconfig}config_{cband}band'
        if worst_calib_level == 'P':
            cqfilt = (targetcals[bccol]=='P')
        elif worst_calib_level == 'S':
            cqfilt = ((targetcals[bccol]=='P') | (targetcals[bccol]=='S'))
        elif worst_calib_level == 'W':
            cqfilt = ((targetcals[bccol]=='P') | (targetcals[bccol]=='S')
                      | (targetcals[bccol]=='W'))
        elif worst_calib_level == 'C':
            cqfilt = ((targetcals[bccol]=='P') | (targetcals[bccol]=='S')
                      | (targetcals[bccol]=='W') | (targetcals[bccol]=='C'))
        else:
            cqfilt = np.ones(len(targetcals)).astype(bool)
        targetcals = targetcals[cqfilt]
        
    return targetcals
 

def optimise_calibrators(tcdata, tcol='target', ccol='calibrator',
                         verbose=False, scol11='sep'):
    'minimise number of calibrators used in a target list'
    target_list = np.unique(tcdata[tcol])
    n_targets = len(target_list)
    caldata = tcdata.copy() ###copy to be edited as iterating through
    data_tc, targets_with_calib = [], [] ##set up empty lists to append to
    
    print(f'optimising calibrators for {n_targets} target')
    print()
    
    ###iterate through calibrators until all targets have a calibrator
    while n_targets > 0:
        ccount = Counter(caldata[ccol])
        ccount = Table({ccol: list(ccount.keys()),
                        'n_targets': list(ccount.values())})
        caldata = join(caldata, ccount, keys=ccol, join_type='left')
        caldata.sort('n_targets', reverse=True) ##prioritise calibrators for many targets
        if np.max(caldata['n_targets']) == 1:
            caldata.sort(scol11)
        c0 = caldata[ccol][0]
        t0 = caldata[caldata[ccol]==c0]
        data_tc.append(t0)
        targets_used = list(t0[tcol])
        targets_with_calib = targets_with_calib + targets_used
        targets_left = Table({tcol: [t for t in target_list if t not in targets_with_calib]})
        if verbose == True:
            print(f'{len(targets_used)} targets for calibrator')
            print(f'{len(targets_left)} targets remaining to find calibrators for')
            print()
        if len(targets_left)>0:
            caldata = join(caldata, targets_left, keys=tcol, join_type='right')
            caldata.remove_column('n_targets')
        n_targets = len(targets_left)
    
    data_tc = vstack(data_tc)
    
    return data_tc


def parse_args():
    'parse command line arguments'
    parser = argparse.ArgumentParser(description="find VLA complex gain calibrators for target list. Also see instructions on VLA calibration at https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/calibration")
    parser.add_argument("targets",
                        help="table file containing sky coords of targets")
    parser.add_argument("--namecol", action='store',
                        type=str, default='Name',
                        help="Object name column name in targets")
    parser.add_argument("--RAcol", action='store',
                        type=str, default='RA',
                        help="RA column name in targets")
    parser.add_argument("--DEcol", action='store',
                        type=str, default='DEC',
                        help="Decl. column name in targets")
    parser.add_argument("--calibrators", action='store',
                        type=str, default='VLA_collaborator_list.txt',
                        help="text file of VLA calibrators from https://science.nrao.edu/facilities/vla/observing/callist")
    parser.add_argument("--band", action='store',
                        type=str, default='L',
                        help="VLA band to find calibrator for")
    parser.add_argument("--quality", action='store',
                        type=str, default='S',
                        help="worst acceptable calibrator quality")
    parser.add_argument("--limit_quality", action='store',
                        type=str, default='True',
                        help="only find potential calibrators with worst acceptable quality or better")
    parser.add_argument("--configs", action='store',
                        type=str, default='A,B,C,D',
                        help="comma separated list of VLA configs to find calibrators for")
    parser.add_argument("--maxsep", action='store',
                        type=str, default='10 deg',
                        help="maximum angular separation between target and calibrator")
    parser.add_argument("--outdir", action='store',
                        type=str, default='.',
                        help="directory to write files to")
    parser.add_argument("--fname_opt", action='store',
                        type=str, default="optimised_calibrator_list.fits",
                        help="filename for optimised calibrator list")
    parser.add_argument("--fname_pot", action='store',
                        type=str, default="potential_calibrators.fits",
                        help="filename for potential calibrators list")
    parser.add_argument("--output_optimised", action='store',
                        type=str, default='True',
                        help="output a minimal list of calibrators for target list")
    
    args = parser.parse_args()
    
    ###change datatypes where appropriate
    args.configs = args.configs.split(',')
    args.maxsep = u.Quantity(args.maxsep)
    args.output_optimised = bool(strtobool(args.output_optimised))
    args.limit_quality = bool(strtobool(args.limit_quality))
    
    return args


#########################################################
#########################################################
###main


if __name__ == '__main__':
    args = parse_args()
    targets = Table.read(args.targets)
    calibrators = parse_nrao_list(args.calibrators)
    clist, cdict = create_datasets_for_querying(callist=calibrators)
    target_calibrators = find_calibrators(targets=targets, cal_sources=clist,
                                          cal_infodict=cdict, tname=args.namecol,
                                          tra=args.RAcol, tdec=args.DEcol,
                                          cband=args.band,
                                          worst_calib_level=args.quality,
                                          limitoncalq=args.limit_quality,
                                          configs=args.configs, max_sep=args.maxsep)
    ###write to file
    target_calibrators.write('/'.join([args.outdir, args.fname_pot]))
    if args.output_optimised == True:
        optlist = optimise_calibrators(tcdata=target_calibrators)
        optlist.write('/'.join([args.outdir, args.fname_opt]))
