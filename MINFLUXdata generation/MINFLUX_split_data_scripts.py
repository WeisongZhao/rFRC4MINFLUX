# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import json
from scipy.ndimage import gaussian_filter, generic_filter
import csv
# from skimage.io import imsave, imshow


def x2idx(x, olr):
    o,l,r = olr
    return (x-o)*r/l+(r-1)/2


def olr2extent(olr):
    o, l, r = olr
    return (o-l/2, o+l/2)


def olrFromValues(array, res=None, ps=None):

    mx = np.nanmax(array, axis=-1, keepdims=True)
    mn = np.nanmin(array, axis=-1, keepdims=True)
    
    l = mx - mn
    o = (mx + mn) / 2
    
    if ps is not None:
        r = np.ceil(l/ps)
        l = r * ps
    else:
        r = np.ones_like(l) * np.ceil(res)

    return np.hstack((o,l,r))


def tophat_filter(img, size):        
    gg = np.meshgrid((*(np.arange(size)-(size-1)/2,)*img.ndim), indexing='ij')
    footprint = np.linalg.norm(gg, axis=0) <= (size)/2
    n = footprint.sum()
        
    return generic_filter(
            img,
            lambda x: x.sum()/n,
            footprint=footprint,
            mode='reflect',
            cval=0.0,
            origin=0,
            ), n

    
def nonlin_filter(img, sigma_px):
    return np.sqrt(gaussian_filter(img, sigma=sigma_px))


def density(X, Y, olrX, olrY):
    den, _, _ = np.histogram2d(
                    X, Y,
                    bins=(olrX[2],
                          olrY[2]),
                    range=(olr2extent(olrX),
                           olr2extent(olrY)),
                    weights=None,
                    )
    return den


def loc2png(loc, fn, pixelsize=0.5e-9, olrs=None, method='tophat', filtersize=2e-9, oglow=None):
    
    if olrs is None:
        olrs = olrFromValues(loc.T, ps=pixelsize)
        olrs[:,1] *= 1
        olrs[:,2] *= 1
    T = loc.T
    img = density(*loc.T, *olrs)

    if method == 'tophat':    
        den, n = tophat_filter(img, size=2*filtersize/pixelsize)
        if oglow is not None:
            den = np.minimum(den, oglow/n)
    elif method == 'nonlin':
        den = nonlin_filter(img, sigma_px=filtersize/pixelsize)
    elif method == 'density':
        den = img
        if oglow is not None:
            den = np.minimum(den, oglow)
    else:
        raise NotImplementedError()

    plt.imshow(den, cmap='gray')
    plt.imsave('{}.png'.format(fn), den, cmap='gray')
    return den



class mfxdta:

    _dtype = np.dtype([
        ('tim', float),
        ('tid', int),
        ('loc', float, 2),
        ('eco', int),
        ('efo', float),
        ('cfr', float),
        ])

    def __init__(self, size):
        self._data = np.zeros(size, dtype=__class__._dtype)

    def __setitem__(self, field, data): self._data[field] = data
    def __getitem__(self, field): return self._data[field]
        
    @classmethod
    def _mfx2dict(cls, a):
        if len(np.shape(a)) > 0:
            return [cls._mfx2dict(k) for k in list(a)]
        if a.dtype.names is not None:
            return {n: cls._mfx2dict(a[n]) for n in a.dtype.names}
        return a.item()
    
    def to_json(self, fn):
        with open(fn, 'w') as f: json.dump(self._mfx2dict(self._data), f, indent='\t')
    
    @classmethod
    def from_json(cls, fn):
        with open(fn, 'r') as f: data = json.load(f)
        
        new = cls(len(data))
        
        for n in cls._dtype.names:
            new[n] = [X[n] for X in data]
                    
        return new
    @classmethod
    def apply_cutoffs(self, efo=np.inf, cfr=np.inf):
        j = (mfx['cfr'] < cfr) * (mfx['efo'] < efo)
        mfx_ = self(j.sum())
        mfx_[:] = mfx[j]
        return mfx_
    

    @classmethod
    def aggregate(self, nph, strict=False, stats=False):
        
        mfx = self._data
        crit = mfx['eco']
        rn = -1
        rix = np.zeros(len(crit), dtype=int)
        tid = None
        
        for i,c in enumerate(crit):
        
            if mfx['tid'][i] != tid:
                s = 0
                tid = mfx['tid'][i]
                rn += 1
                
            rix[i] = rn
            s += c
        
            if s >= nph: tid = None
                
        mfx_agg = mfxdta(rn)
        mfx_agg['cfr'] = np.nan
        
        for i in np.arange(rn):            
            l = mfx[np.where(rix==i)]
            mfx_agg['tid'][i] = l['tid'][0]
            mfx_agg['tim'][i] = l['tim'][0]
            mfx_agg['eco'][i] = l['eco'].sum()
            mfx_agg['efo'][i] = mfx_agg['eco'][i] / (l['eco']/l['efo']).sum()        
            mfx_agg['loc'][i] = (l['loc'] * l['eco'][:,np.newaxis]).sum(0) / l['eco'].sum()
        
        return mfx_agg
    
    
    
        
   
if __name__ == '__main__':

    # ==============================nup96===================================
    mfx = mfxdta.from_json('nup96.json')
    mfx_odd = mfxdta(1130)
    mfx_even = mfxdta(1130)
    flag_odd = 0
    flag_even = 0
    for j in range(2260):
        if j%2 == 0:
            mfx_even[flag_even] = (mfx[j])   
            flag_even = flag_even + 1
        else:
            mfx_odd[flag_odd] = (mfx[j])
            flag_odd = flag_odd +1 
            
   
    mfx_even.to_json('nup_even.json')    
    mfx_odd.to_json('nup_odd.json')  

# ====================================nup96_odd&even_filter=============================
    mfx = mfxdta.from_json('nup_even.json')
    loc2png(mfx['loc'], pixelsize=1e-9, method='nonlin', fn='nup_even_wo_filter', filtersize = 2.11932e-9)
    mfx = mfxdta.from_json('nup_odd.json')
    loc2png(mfx['loc'], pixelsize=1e-9, method='nonlin', fn='nup_odd_wo_filter', filtersize = 2.11932e-9)
    mfx = mfxdta.from_json('nup_even.json')
    mfx = mfx.apply_cutoffs(efo=50e3)
    loc2png(mfx['loc'], pixelsize=1e-9, method='nonlin', fn='nup_even_w_first_filter', filtersize = 2.11932e-9)
    mfx = mfxdta.from_json('nup_odd.json')
    mfx = mfx.apply_cutoffs(efo=50e3)
    loc2png(mfx['loc'], pixelsize=1e-9, method='nonlin', fn='nup_odd_w_first_filter', filtersize = 2.11932e-9)
    mfx = mfxdta.from_json('nup_even.json')
    mfx = mfx.apply_cutoffs(efo=50e3)
    mfx = mfx.aggregate(2100)
    loc2png(mfx['loc'], pixelsize=1e-9, method='nonlin', fn='nup_even_w_both_filter', filtersize = 2.11932e-9)
    mfx = mfxdta.from_json('nup_odd.json')
    mfx = mfx.apply_cutoffs(efo=50e3)
    mfx = mfx.aggregate(2100)
    loc2png(mfx['loc'], pixelsize=1e-9, method='nonlin', fn='nup_odd_w_both_filter', filtersize = 2.11932e-9)
    # ==============================spectrin===========================
    mfx = mfxdta.from_json('spectrin.json')
    mfx_odd = mfxdta(5411)
    mfx_even = mfxdta(5412)
    flag_odd = 0
    flag_even = 0
    for i in range(10823):
        if i%2 == 0:
            mfx_even[flag_even] = (mfx[i])   
            flag_even = flag_even + 1
        else:
            mfx_odd[flag_odd] = (mfx[i])
            flag_odd = flag_odd +1 
            
    mfx_even.to_json('spectrin_even.json')    
    mfx_odd.to_json('spectrin_odd.json')  
# =====================================spectrin_filter===========================    
    mfx = mfxdta.from_json('spectrin_even.json')
    loc2png(mfx['loc'], pixelsize=2e-9, method='nonlin', fn='spectrin_even_wo_filter', filtersize = 4.47412e-9)
    mfx = mfxdta.from_json('spectrin_odd.json')
    loc2png(mfx['loc'], pixelsize=2e-9, method='nonlin', fn='spectrin_odd_wo_filter', filtersize = 4.47412e-9)
    mfx = mfxdta.from_json('spectrin_even.json')
    mfx = mfx.apply_cutoffs(efo=90e3, cfr=0.6)
    loc2png(mfx['loc'], pixelsize=2e-9, method='nonlin', fn='spectrin_even_w_first_filter', filtersize = 4.47412e-9)
    mfx = mfxdta.from_json('spectrin_odd.json')
    mfx = mfx.apply_cutoffs(efo=90e3, cfr=0.6)
    loc2png(mfx['loc'], pixelsize=2e-9, method='nonlin', fn='spectrin_odd_w_first_filter', filtersize = 4.47412e-9)
    mfx = mfxdta.from_json('spectrin_even.json')
    loc2png(mfx['loc'], pixelsize=2e-9, method='nonlin', fn='spectrin_even_wo_filter_V2', filtersize = 2.11932e-9)
    mfx = mfxdta.from_json('spectrin_odd.json')
    loc2png(mfx['loc'], pixelsize=2e-9, method='nonlin', fn='spectrin_odd_wo_filter_V2', filtersize = 2.11932e-9)
    mfx = mfxdta.from_json('spectrin_even.json')
    mfx = mfx.apply_cutoffs(efo=90e3, cfr=0.6)
    loc2png(mfx['loc'], pixelsize=2e-9, method='nonlin', fn='spectrin_even_w_first_filter_V2', filtersize = 2.11932e-9)
    mfx = mfxdta.from_json('spectrin_odd.json')
    mfx = mfx.apply_cutoffs(efo=90e3, cfr=0.6)
    loc2png(mfx['loc'], pixelsize=2e-9, method='nonlin', fn='spectrin_odd_w_first_filter_V2', filtersize = 2.11932e-9)