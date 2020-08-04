'''
File: bayesain_adjust.py
Author: Xianhgong Che
Version: 1.0
Create: 2020-04-02
Description: adjusting the Landsat data base on the optimal sigma using beyesian theory
'''


def conjugate_gradient(_bnd_m, _bnd_l, _rad, _sig, _b_id, _d_out):
	import numpy as np
	import math
	from scipy.ndimage import filters
	from scipy import ndimage
	import time
	import os

	# for _scale in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2]:
	# for _scale in [0.2,0.4,0.6,0.8,1.0,1.2,1.4]:
	# for _scale in [0.6,0.8,1.0]:
	for _scale in [0.4]:
		_Z0 = _bnd_l.data
		_z = _bnd_m.data
		
		_z = np.ma.masked_equal(_z, 32767)
		#_Z0 = np.ma.masked_equal(_Z0, -9999)
		_locs = np.where(_Z0 == 20000)
		_Z0 = np.ma.masked_where((_Z0 == -9999)|(_Z0 == 20000), _Z0)
	
		#_Cn = (np.ma.amin(_z) + (np.ma.amax(_z) - np.ma.amin(_z)) * _scale) ** 2
		_C = np.ma.var(_Z0) 
		_Cn = _C * _scale
		
		_tru = (((_rad * 15 - 1)/2)-0.5)/_sig
		_Z0_f = _Z0.filled(0)

		_HZ0 = ndimage.zoom(filters.gaussian_filter(_Z0_f, _sig, mode = 'nearest', truncate = _tru), 1.0/15, order=0)     # degrading 
		
		_dat_md = ((_z - _HZ0) / _Cn).filled(0)
			
		_dat_md = ndimage.zoom(_dat_md, 15, order=0)
		_dat_r0 = filters.gaussian_filter(_dat_md, _sig, mode = 'nearest', truncate = _tru)
		_dat_p0 = _dat_r0
		
		# lt = []

		for i in range(20):
			_ap_h = ndimage.zoom(filters.gaussian_filter(_dat_p0, _sig, mode = 'nearest', truncate = _tru),1.0/15, order=0)
			_ap = filters.gaussian_filter(ndimage.zoom(_ap_h, 15, order=0), _sig, mode = 'nearest', truncate = _tru) / _Cn + _dat_p0 / _C     # convolving Landsat with MODIS PSF 
		
			_ak = np.sum(_dat_r0 * _dat_r0) / np.sum(_dat_p0 * _ap)
			_Z1 = _Z0 + _ak * _dat_p0
			
			_dat_r1 = _dat_r0 - _ak * _ap
			# if np.ma.average(_dat_r1) < _diff:
				# break
			# lt.append(str(np.average(np.abs(_dat_r1))))
			
			lt.append(str(np.ma.average(np.abs(_Z1-_Z0))))
			
			_rk = np.sum(_dat_r1 * _dat_r1) / np.sum(_dat_r0 * _dat_r0)
			_dat_p1 = _dat_r1 + _rk * _dat_p0 
			
			_dat_r0 = _dat_r1
			_dat_p0 = _dat_p1
			_Z0 = _Z1
			

		_f_out = os.path.join(_d_out, 'L5_p124r044_20001123_b%s_%s_pred.tif' % (_b_id, _scale))
		_Z1.filled(-9999)[_locs] = 20000
		_bnd_l.from_grid(_Z1.astype(np.int16)).save(_f_out, opts=['compress=lzw'])
		
		# f_out = '/mnt/data1/xhche/data/NSF/china_data/results/L5_p124r044_20001123_b%s_diff.csv' % _b_id
		# with open(f_out, 'w') as _fo:
			# _fo.write('\n'.join(lt))  	

def _task(_b_id, _rad, _sig,_rows_m, _cols_m):
	from gio import geo_raster as ge
	from gio import geo_raster_ex as gx
	from gio import geo_base as gb
	import os
	import math
	
	_row_s = 1580
	_col_s = 270
	
	_bnd_idx = {1:3,2:4,3:1,4:2,5:6,7:7}
	
	# _bnd_idx = {2:3,3:4,4:1,5:2,6:6,7:7}
	_d_landsat = '/mnt/data1/xhche/data/NSF/china_data/landsat_unzip'
	_f_modis = '/mnt/data1/xhche/data/NSF/china_data/MCD43A4.A2000304.h28v06.006.2016110011945.hdf'
	_d_out = '/mnt/data1/xhche/data/NSF/china_data/results/'

	_bnd_m = ge.geo_raster.open(_f_modis).get_subdataset('Reflectance_Band%s' % _bnd_idx[_b_id]).get_band()
	geo = _bnd_m.geo_transform
	_geo_m = [geo[0] + _col_s * geo[1], geo[1], 0, geo[3] + _row_s * geo[5], 0 , geo[5]]
	_geo_l = [geo[0] + _col_s * geo[1], geo[1]/15, 0, geo[3] + _row_s * geo[5], 0 , geo[5]/15]
	
	_tmp_l = ge.geo_band_info(_geo_l, _cols_m * 15, _rows_m * 15, gx.modis_projection())
	_tmp_m = ge.geo_band_info(_geo_m, _cols_m, _rows_m, gx.modis_projection())
	
	_bnd_l = ge.geo_raster.open(os.path.join(_d_landsat, 'LT05_L1TP_124044_20001123_20161213_01_T1_sr_band%s.tif' % _b_id)).get_band().read_block(_tmp_l)
	
	# _bnd_l.save('/mnt/data1/xhche/data/NSF/france_data/results/TOA/L7_p124r044_20001030_landsat_b%s.tif' % _b_id, opts=['compress=lzw'])
	
	_bnd_m = _bnd_m.read_block(_tmp_m)
	# _bnd_m.save('/mnt/data1/xhche/data/NSF/france_data/results/TOA/L7_p124r044_20001030_modis_b%s.tif' % _b_id, opts=['compress=lzw'])
	# return

	conjugate_gradient(_bnd_m, _bnd_l, _rad, _sig, _b_id, _d_out)
	
def main(opts):
	from gio import geo_raster_ex as gx
	from gio import geo_raster as ge
	import os
	import re
	import numpy as np
	from gio import multi_task
	import time
	
	_rows_m = 160
	_cols_m = 160
	
	ps = {1:(3,290), 2:(3,270), 3:(3,280), 4:(3, 270), 5:(3,280),7:(3,320)}
	# ps = {2:(3,280), 3:(3,240), 4:(3,260), 5:(3, 250), 6:(3,260),7:(3,290)}
	_ts = []

	for i in [1,2,3,4,5,7]:
	# for i in [2,3,4,5,6,7]:
		_ts.append([i, ps[i]])
	
	#_task(4, 1, 150, _rows_l, _cols_l,_rows_m, _cols_m)	
	multi_task.run(_task, [(_t[0], _t[1][0], _t[1][1], _rows_m, _cols_m) for _t in multi_task.load(_ts, opts)], opts)	
	
def usage():
	_p = environ_mag.usage(True)

	return _p

if __name__ == '__main__':
	from gio import environ_mag
	environ_mag.init_path()
	environ_mag.run(main, [environ_mag.config(usage())])	
		

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	