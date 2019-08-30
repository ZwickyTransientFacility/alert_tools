import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time

def open_avro(fname):
    with open(fname,'rb') as f:
        freader = fastavro.reader(f)
        schema = freader.writer_schema
        for packet in freader:
            return packet

def make_dataframe(packet):
    dfc = pd.DataFrame(packet['candidate'], index=[0])
    df_prv = pd.DataFrame(packet['prv_candidates'])
    dflc = pd.concat([dfc,df_prv], ignore_index=True,sort=True)	
    dflc.objectId = packet['objectId']	
    dflc.candid = packet['candid']
    return dflc

def dcmag(dflc, match_radius_arcsec=1.5, star_galaxy_threshold = 0.4,band=2):
    if (dflc.loc[0,'distpsnr1'] > match_radius_arcsec) & (dflc.loc[0,'sgscore1'] < star_galaxy_threshold):
        print('Object is not a variable star.')
        return dflc
                
    else:
        dflc=dflc.fillna(np.nan)
        def robust_median(x):
            if len(x) == 0:
                return np.nan
            else:
                return np.median(x[np.isfinite(x)])
        grp = dflc.groupby(['fid','field','rcid'])
        impute_magnr = grp['magnr'].agg(robust_median)
        #print(impute_magnr)
        impute_sigmagnr = grp['sigmagnr'].agg(robust_median)
        #print(impute_sigmagnr)
        for idx, grpi in grp:
            w = np.isnan(grpi['magnr'])
            w2 = grpi[w].index
            dflc.loc[w2,'magnr'] = impute_magnr[idx]
            dflc.loc[w2,'sigmagnr'] = impute_sigmagnr[idx]
        dflc['sign'] = 2* (dflc['isdiffpos'] == 't') - 1
        
        dflc['dc_mag'] = -2.5 * np.log10(10**(-0.4*dflc['magnr']) + dflc['sign'] * 10**(-0.4*dflc['magpsf']))   #u
        dflc['dc_sigmag'] = np.sqrt(
            (10**(-0.4*dflc['magnr'])* dflc['sigmagnr']) **2. + 
            (10**(-0.4*dflc['magpsf']) * dflc['sigmapsf'])**2.) / 10**(-0.4*dflc['magnr']) + dflc['sign'] * 10**(-0.4*dflc['magpsf']) #u
        dflc['dc_mag_ulim'] = -2.5 * np.log10(10**(-0.4*dflc['magnr']) + 10**(-0.4*dflc['diffmaglim']))   #v
        dflc['dc_mag_llim'] = -2.5 * np.log10(10**(-0.4*dflc['magnr']) - 10**(-0.4*dflc['diffmaglim']))   #v2
	
	
        return dflc
        
def band_amplitude(dflc, band=2):
    if 'dc_mag' in dflc.columns.values:
       mag_key= 'dc_mag'
    else:
       mag_key= 'magpsf'

    z = dflc[dflc.fid==band]
    ampli=z[mag_key].max()-z[mag_key].min()
    print('Max:',z[mag_key].max())
    print('Min:',z[mag_key].min())
    print('Amplitude:',ampli)
    print('Is amplitude > 1.0 mag?',ampli>=1)
    
    return ampli
        
        




def plot_dc_lightcurve(dflc, days_ago=True):
    
    plt.rcParams["figure.figsize"] = (10,7)
    
    filter_color = {1:'green', 2:'red', 3:'pink'}
    if days_ago:
        now = Time.now().jd
        t = dflc.jd - now
        xlabel = 'Days Ago'
    else:
        t = dflc.jd
        xlabel = 'Time (JD)'
    
    fig=plt.figure()
    for fid, color in filter_color.items():
        # plot detections in this filter:
        w = (dflc.fid == fid) & ~dflc.magpsf.isnull()
        if np.sum(w):
            plt.errorbar(t[w],dflc.loc[w,'dc_mag'], dflc.loc[w,'dc_sigmag'],fmt='.',color=color)
        wnodet = (dflc.fid == fid) & dflc.magpsf.isnull()
        if np.sum(wnodet):
            plt.scatter(t[wnodet],dflc.loc[wnodet,'dc_mag_ulim'], marker='v',color=color,alpha=0.25)
            plt.scatter(t[wnodet],dflc.loc[wnodet,'dc_mag_llim'], marker='^',color=color,alpha=0.25)

    
    plt.gca().invert_yaxis()
    plt.xlabel(xlabel)
    plt.ylabel('Magnitude')
	
	return fig

	
def get_dc_mag(dflc, band=2, days_ago=True):
    if 'dc_mag' not in dflc.columns.values:
        dflc = dcmag(dflc)
    else:
        dflc=dflc

    amplitude=band_amplitude(dflc, band=band)
	
    fig=plot_dc_lightcurve(dflc, days_ago=days_ago)
	
	return dflc, amplitude, fig
	
	

