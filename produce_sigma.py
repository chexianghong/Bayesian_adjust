'''
File: produce_sigma.py
Author: Xianhgong Che
Version: 1.0
Create: 2020-03-14
Description: using different sigma of MODIS PSF to calculate the R2 between spatial degraded Landsat and MODIS NBAR
'''


def cal_psf(_rad):
    import math
    import numpy as np
    _num = _rad * 15
    ps = {}
    
    for _sig in range(100, 710, 10):
        m_psf = np.zeros((_num, _num), dtype = np.float32)
        _sig1 = _sig * _sig
        for i in range(_num):
            for j in range(_num):
                m_psf[i,j]  = 1.0 / (math.sqrt(2 * math.pi) * _sig) * math.exp(- (900*((_num//2-i)**2+(_num//2-j)**2))/(2.0*_sig1))
        ps[_sig] = m_psf / np.sum(m_psf)
    return ps
   
        
def cal_landsat(_dat, _row, _col, val_m, _sigs, pt, _rad):
    # degrading Landsat data
    import numpy as np
    ls = []
    _rad = _rad // 2
   
    _row_n = _row + _rad
    _col_n = _col + _rad
    _dat1 = _dat[(_row_n - _rad) * 15 :(_row_n + _rad + 1) *15, (_col_n - _rad) * 15 : (_col_n + _rad + 1) *15]
    
    if np.where((_dat1 == -9999) | (_dat1 == 20000))[0].size != 0:
        return
        
    for _sig, _psf in _sigs:
        pt[_sig].append((val_m, int(np.sum(_dat1 * _psf))))

def _task1(_rad, _b_id, _rows_m, _cols_m):
   
    from gio import geo_raster as ge
    from gio import geo_raster_ex as gx
    from gio import geo_base as gb
    import linear_regress
    import collections
    import os
    import numpy as np
    
    
    # the location of the test area 
    _row_s = 1500
    _col_s = 350
    
    f_out = '/mnt/data1/xhche/data/NSF/south_dokata/results/p027r030_20180813_sig_b%s_rad%s.csv' % (_b_id, _rad) 
    # if os.path.exists(f_out):
        # return
        
    _bnd_idx = {1:3,2:4,3:1,4:2,5:6,7:7}
    
    # _bnd_idx = {2:3,3:4,4:1,5:2,6:6,7:7}
     
     # Landsat SR directory
    _d_landsat = '/mnt/data1/xhche/data/NSF/south_dokata/landsat_unzip'
    f_modis = '/mnt/data1/xhche/data/NSF/south_dokata/MCD43A4.A2001274.h11v04.006.2016118231306.hdf'
    
    _bnd_m = ge.geo_raster.open(f_modis).get_subdataset('Reflectance_Band%s' % _bnd_idx[_b_id]).get_band()
    geo = _bnd_m.geo_transform
    _geo_m = [geo[0] + _col_s * geo[1], geo[1], 0, geo[3] + _row_s * geo[5], 0 , geo[5]]
    _geo_l = [geo[0] + _col_s * geo[1], geo[1]/15, 0, geo[3] + _row_s * geo[5], 0 , geo[5]/15]
    
    ps = cal_psf(_rad)
    _sigs = list(ps.keys())
    _sigs.sort()
    
    _tmp_l = ge.geo_band_info(_geo_l, _cols_m * 15, _rows_m * 15, gx.modis_projection())
    _tmp_m = ge.geo_band_info(_geo_m, _cols_m, _rows_m, gx.modis_projection())

    _bnd_qc = ge.geo_raster.open(os.path.join(_d_landsat, 'LT05_L1TP_027030_20011001_20160917_01_T1_pixel_qa.tif')).get_band().read_block(_tmp_l)
  
    ls = ['sig,slope,intercept,r2,num']
    pt = collections.defaultdict(lambda:[])
    
    _bnd_l = ge.geo_raster.open(os.path.join(_d_landsat, 'LT05_L1TP_027030_20011001_20160917_01_T1_sr_band%s.tif' % _b_id)).get_band().read_block(_tmp_l)
    _bnd_l_dat = _bnd_l.data
    _bnd_l_dat[_bnd_qc.data & 41 > 0] = -9999
    
    _radus = (_rad // 2) * 15
    _bnd_l_dat = np.pad(_bnd_l_dat, ((_radus, _radus),(_radus, _radus)),'constant', constant_values = (-9999,-9999))

    _bnd_m = _bnd_m.read_block(_tmp_m)

    for _row in range(_rows_m):
        for _col in range(_cols_m):
            val_m = _bnd_m.data[_row, _col]
            if val_m == 32767:
                continue
            cal_landsat(_bnd_l_dat, _row, _col, val_m, list(ps.items()), pt, _rad)
    
    # output the R2 with the different sigma
    for _sig in _sigs:
        _vvs = pt[_sig]
        ls_m, ls_l= zip(*_vvs)
        _b0, _b1, _r2 = linear_regress.linear_regress_rma(ls_m, ls_l)
        ls.append('%s,%s,%s,%s,%s' % (_sig, _b1, _b0, _r2, len(ls_m)))

    with open(f_out, 'w') as _fo:
        _fo.write('\n'.join(ls))   
        
    print ('processing band%s with radius %s done' % (_b_id, _rad))    
    
def main(opts):
    import os 
    from gio import multi_task
    import re

    _rows_m = 160
    _cols_m = 160
    _ts = []

    for _b in [1,2,3,4,5,7]:
    # for _b in [2,3,4,5,6,7]:
        _ts.append([3, _b])

    # _task(3, 1, _rows_m, _cols_m)
    multi_task.run(_task1, [(_t[0], _t[1],_rows_m, _cols_m) for _t in multi_task.load(_ts, opts)], opts)   
    
def usage():
    _p = environ_mag.usage(True)
    # _p.add_argument('-i', '--input', dest='input', required=True)
    # _p.add_argument('-o', '--output', dest='output', required=True)

    return _p

if __name__ == '__main__':
    from gio import environ_mag
    environ_mag.init_path()
    environ_mag.run(main, [environ_mag.config(usage())])

