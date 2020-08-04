'''
File: rma_adjust.py
Author: Xianhgong Che
Version: 1.0
Create: 2020-05-10
Description: adjusting the Landsat data base on the optimal sigma using the linear regression method
'''

def cal_block_difference(_bnd, _dis, _row, _col, _rows, _cols):
	ls = []		
	for _rr in range(max(0, _row -_dis), min(_row + _dis + 1, _rows)):
		for _cc in range(max(0, _col -_dis), min(_col + _dis +1, _cols)):
			_v = _bnd.data[_rr, _cc]
			if _v == _bnd.nodata:
				continue
			else:
				ls.append(_v)
				
	if len(ls):
		return max(ls) - min(ls)
		
	return None
	
def cal_landsat_difference(_bnd, _bnd_qc_dis, _row, _col, _rows, _cols):
	ls = []		
	for _rr in range(max(0, _row -_dis), min(_row + _dis + 1, _rows)):
		for _cc in range(max(0, _col -_dis), min(_col + _dis +1, _cols)):
			_v = _bnd.data[_rr, _cc]
			if _v == _bnd.nodata:
				continue
			else:
				ls.append(_v)
				
	if len(ls):
		return max(ls) - min(ls)
		
	return None
	
	
def cal_mean(_dat, _dat_qc):
	ls = []
	for i in range(15):
		for j in range(15):
			_val = _dat[i,j]
			if _val in [-9999, 20000]:
				return None
			if _dat_qc[i,j] & 41 >0:
				return None
			ls.append(_val)
	if len(ls):
		return [sum(ls)/len(ls), max(ls) - min(ls)]
	else:
		return None


def _task(_bnd_id):
	# building linear relation with PSF aggregation
	from gio import geo_raster as ge
	from gio import geo_raster_ex as gx
	from gio import geo_base as gb
	import os
	import math
	import numpy as np
	import linear_regress
	import cal_metrics
	from scipy.ndimage import filters
	from scipy import ndimage
	
	_rows_m = 160
	_cols_m = 160
	

	d_out = '/mnt/data1/xhche/data/NSF/china_data/results/'
	
	ls_x = []
	ls_y = []
	
	# the parameters of MODIS PSF
	ps = {2:(3,280), 3:(3,240), 4:(3,260), 5:(3, 250), 6:(3,260),7:(3,290)}
	_info = 'L8_p124r044_20190925'
	_ff_l = os.path.join(d_out, '%s_landsat_b%s.tif' % (_info, _bnd_id))
	_ff_m = os.path.join(d_out, '%s_modis_b%s.tif'%(_info, _bnd_id))
	_f_out = os.path.join(d_out, '%s_b%s_rma.tif'%(_info, _bnd_id))
	
	_bnd_l = ge.geo_raster.open(_ff_l).get_band().cache()
	_bnd_m = ge.geo_raster.open(_ff_m).get_band().cache()
	
	_bnd_qc = ge.geo_raster.open('/mnt/data1/xhche/data/NSF/china_data/landsat_unzip/LC08_L1TP_124044_20190925_20191017_01_T1_pixel_qa.tif').get_band().read_block(_bnd_l)
	
	_rad = ps[_bnd_id][0]
	_sig = ps[_bnd_id][1]
	_tru = (((_rad * 15 - 1)/2)-0.5)/_sig
	
	# degrading Landsat data
	_dat_l = ndimage.zoom(filters.gaussian_filter(_bnd_l.data, _sig, mode = 'nearest', truncate = _tru), 1.0/15, order=0)

	for _rr in range(1,159):
		for _cc in range(1,159):

			val_m = _bnd_m.data[_rr, _cc]
			if val_m == _bnd_m.nodata:
				continue
				
			_dat_b = _bnd_l.data[(_rr-1)*15:(_rr+2)*15, (_cc-1)*15: (_cc+2)*15]
			if np.where((_dat_b == -9999)|(_dat_b == 20000))[0].size > 0:
				continue

			_dat_b_qc = _bnd_qc.data[(_rr-1)*15:(_rr+2)*15, (_cc-1)*15: (_cc+2)*15]
			if np.where(_dat_b_qc & 41 > 0)[0].size > 0:
				continue

			ls_x.append(_dat_l[_rr, _cc])
			ls_y.append(val_m)
		
	bias, gain, r2 = linear_regress.linear_regress_rma(ls_x, ls_y)
	RMSD,MSEs,MSEu = cal_metrics.cal_RMSDs(ls_x, ls_y, gain, bias)
	print (_bnd_id, len(ls_x), bias, gain, r2)
	
	_locs = np.where(_bnd_l.data == 20000)
	_dat = np.ma.masked_where((_bnd_l.data == -9999)|(_bnd_l.data == 20000), _bnd_l.data)
	_dat = _dat * gain + bias
	# _dat[np.ma.where(_dat > 10000)] = 10000
	# _dat[np.ma.where(_dat < 0)] = 0
	_dat.filled(-9999)[_locs] = 20000
	_bnd = _bnd_l.from_ma_grid(_dat.astype(np.int16), nodata=-9999)
	_bnd.save(_f_out, opts=['compress=lzw'])
	

def main(opts):

	import os
	import re
	import datetime
	from gio import geo_raster as ge
	_ts = []
	
	for i in [2,3,4,5,6,7]:
	# for i in [1,2,3,4,5,7]:
		_ts.append(i)
			
	from gio import multi_task
	multi_task.run(_task, [_r for _r in multi_task.load(_ts, opts)], opts)

def usage():
	_p = environ_mag.usage(True)
	return _p

if __name__ == '__main__':
	from gio import environ_mag
	environ_mag.init_path()
	environ_mag.run(main, [environ_mag.config(usage())])
	
		