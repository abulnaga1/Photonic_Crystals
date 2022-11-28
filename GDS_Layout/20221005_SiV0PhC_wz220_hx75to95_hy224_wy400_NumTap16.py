'''
Draw gds pattern from python code.

DH:
Version 1: 2018-01-31

Version 2: 2018-09-12

Version 3: 2019-05-30

Options to change pinch point location
Text options

Version 4: 2020-07-30
Options for negative tone resist

Version 5: 2022-19-05
Introduce a column of mirror-hole only devices to aid with diagnostics

'''

import numpy
import copy
import math
import gdspy

# layers
# 1: beam
# 2: holes
# 3: outerbox

def write_beam_single(cell,pltdata,layer=1):

	cell.add(gdspy.Polygon([
			(pltdata.xpos - pltdata.dx/2.0, pltdata.ypos - pltdata.dy/2.0),
			(pltdata.xpos - pltdata.dx/2.0, pltdata.ypos + pltdata.dy/2.0),
			(pltdata.xpos + pltdata.dx/2.0, pltdata.ypos + pltdata.dy/2.0),
			(pltdata.xpos + pltdata.dx/2.0, pltdata.ypos - pltdata.dy/2.0),
			(pltdata.xpos - pltdata.dx/2.0, pltdata.ypos - pltdata.dy/2.0)
			], layer=layer))

def write_beam_array(cell,first_beam,width_list=None):
	global num_beam, beam_spacing

	if width_list is None:
		_dy_list = first_beam.dy*numpy.ones(num_beam, dtype=int)
	else:
		_dy_list = width_list

	_beam_write = copy.copy(first_beam)
	_initial_ypos = first_beam.ypos - _dy_list[0]/2.0
	for i in range(num_beam):
		_beam_write.dy = _dy_list[i]
		_beam_write.ypos = _initial_ypos + i*beam_spacing + numpy.sum(_dy_list[:i]) + _beam_write.dy/2
		write_beam_single(cell,_beam_write)

	return _dy_list

def write_hole_single(cell,pltdata,layer=2):
	global num_circ_points, shot_pitch

	_philist = numpy.linspace(0,2*numpy.pi,num_circ_points)

	_circ_pts = [((round(pltdata.xpos/shot_pitch) + round(pltdata.dx/2 * numpy.cos(phi)/shot_pitch)) * shot_pitch,
		(round(pltdata.ypos/shot_pitch) + round(pltdata.dy/2 * numpy.sin(phi)/shot_pitch)) * shot_pitch) for phi in _philist]

	return gdspy.Polygon(_circ_pts,layer=layer)

def write_hole_1D(cell,beamdata,holedata,num_taper_hole,num_mir_hole_L,num_mir_hole_R,acav,amir,
	end_taper_L=False,end_taper_R=False,num_end_taper=0,reverse_tone=False,layer=2):

	_one_side_holes = (num_taper_hole-1)/2.0
	_idx_list = numpy.linspace(-_one_side_holes,_one_side_holes,num_taper_hole)
	_aper_a = (amir - acav)/(_one_side_holes*_one_side_holes - 0.25)
	_aper_c = (amir - 4*acav*_one_side_holes*_one_side_holes)/(1-4*_one_side_holes*_one_side_holes)
	_acav_list = [_aper_a*x*x + _aper_c for x in _idx_list]
	_amir_list_L = [amir for x in range(int(num_mir_hole_L+num_end_taper*end_taper_L))]
	_amir_list_R = [amir for x in range(int(num_mir_hole_R+num_end_taper*end_taper_R))]
	_aper_list = copy.copy(_amir_list_L)
	_aper_list.extend(_acav_list)
	_aper_list.extend(_amir_list_R)

	_hole_write = copy.copy(holedata)
	_hole_write.xpos = _hole_write.xpos - numpy.sum(numpy.array(_acav_list))/2.0 - numpy.sum(numpy.array(_amir_list_L))
	_taper_scale_list = []
	if num_end_taper > 0:
		_taper_scale_list = numpy.linspace(0,1.0,num_end_taper+2)
		_taper_scale_list_L = _taper_scale_list[1:-1]
		_taper_scale_list_R = numpy.flipud(_taper_scale_list[1:-1])

	for i in range(len(_aper_list)):
		_hole_write.xpos = _hole_write.xpos + _aper_list[i]/2.0
		if i < num_end_taper*end_taper_L:
			_hole_write.dx = holedata.dx * _taper_scale_list_L[i]
			_hole_write.dy = holedata.dy * _taper_scale_list_L[i]
		else:
			_hole_write.dx = holedata.dx
			_hole_write.dy = holedata.dy
		if i >= (len(_aper_list) - num_end_taper*end_taper_R):
			_hole_write.dx = holedata.dx * _taper_scale_list_R[i-len(_aper_list)+num_end_taper]
			_hole_write.dy = holedata.dy * _taper_scale_list_R[i-len(_aper_list)+num_end_taper]
		_single_hole_polygon = write_hole_single(cell,_hole_write)

		if i == 0:
			_hole_polygon_combined = _single_hole_polygon
		else:
			_hole_polygon_combined = gdspy.fast_boolean(_hole_polygon_combined,_single_hole_polygon, 'or', max_points=0, layer = layer)

		_hole_write.xpos = _hole_write.xpos + _aper_list[i]/2.0

	if reverse_tone is True:
		_tmp_beam_polygon = gdspy.Polygon([
				(beamdata.xpos - beamdata.dx/2.0, beamdata.ypos - beamdata.dy/2.0),
				(beamdata.xpos - beamdata.dx/2.0, beamdata.ypos + beamdata.dy/2.0),
				(beamdata.xpos + beamdata.dx/2.0, beamdata.ypos + beamdata.dy/2.0),
				(beamdata.xpos + beamdata.dx/2.0, beamdata.ypos - beamdata.dy/2.0),
				(beamdata.xpos - beamdata.dx/2.0, beamdata.ypos - beamdata.dy/2.0)
				], layer=layer)
		_tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon,_hole_polygon_combined, 'not', max_points=0, layer = layer)
		cell.add(_tmp_beam_polygon)
	else:
		cell.add(_hole_polygon_combined)

def write_hole_1D_mir_sweep(cell,beamdata,holedata,num_taper_hole,num_mir_hole_L,num_mir_hole_R,acav,amir,
	end_taper_L=False,end_taper_R=False,num_end_taper=0,reverse_tone=False,layer=2):

	_one_side_holes = 0
	_amir_list_L = [amir for x in range(int(num_mir_hole_L+num_end_taper*end_taper_L))]
	_amir_list_R = [amir for x in range(int(num_mir_hole_R+num_end_taper*end_taper_R))]
	_aper_list = copy.copy(_amir_list_L)
	_aper_list.extend(_amir_list_R)
	if len(_aper_list)>0:
		_hole_write = copy.copy(holedata)
		_hole_write.xpos = _hole_write.xpos - numpy.sum(numpy.array(_amir_list_L))
		_taper_scale_list = []
		if num_end_taper > 0:
			_taper_scale_list = numpy.linspace(0,1.0,num_end_taper+2)
			_taper_scale_list_L = _taper_scale_list[1:-1]
			_taper_scale_list_R = numpy.flipud(_taper_scale_list[1:-1])

		for i in range(len(_aper_list)):
			_hole_write.xpos = _hole_write.xpos + _aper_list[i]/2.0
			if i < num_end_taper*end_taper_L:
				_hole_write.dx = holedata.dx * _taper_scale_list_L[i]
				_hole_write.dy = holedata.dy * _taper_scale_list_L[i]
			else:
				_hole_write.dx = holedata.dx
				_hole_write.dy = holedata.dy
			if i >= (len(_aper_list) - num_end_taper*end_taper_R):
				_hole_write.dx = holedata.dx * _taper_scale_list_R[i-len(_aper_list)+num_end_taper]
				_hole_write.dy = holedata.dy * _taper_scale_list_R[i-len(_aper_list)+num_end_taper]
			_single_hole_polygon = write_hole_single(cell,_hole_write)

			if i == 0:
				_hole_polygon_combined = _single_hole_polygon
			else:
				_hole_polygon_combined = gdspy.fast_boolean(_hole_polygon_combined,_single_hole_polygon, 'or', max_points=0, layer = layer)

			_hole_write.xpos = _hole_write.xpos + _aper_list[i]/2.0

		if reverse_tone is True:
			_tmp_beam_polygon = gdspy.Polygon([
					(beamdata.xpos - beamdata.dx/2.0, beamdata.ypos - beamdata.dy/2.0),
					(beamdata.xpos - beamdata.dx/2.0, beamdata.ypos + beamdata.dy/2.0),
					(beamdata.xpos + beamdata.dx/2.0, beamdata.ypos + beamdata.dy/2.0),
					(beamdata.xpos + beamdata.dx/2.0, beamdata.ypos - beamdata.dy/2.0),
					(beamdata.xpos - beamdata.dx/2.0, beamdata.ypos - beamdata.dy/2.0)], layer=layer)
			_tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon,_hole_polygon_combined, 'not', max_points=0, layer = layer)
			cell.add(_tmp_beam_polygon)
		else:
			cell.add(_hole_polygon_combined)
	else:
		_tmp_beam_polygon = gdspy.Polygon([
				(beamdata.xpos - beamdata.dx/2.0, beamdata.ypos - beamdata.dy/2.0),
				(beamdata.xpos - beamdata.dx/2.0, beamdata.ypos + beamdata.dy/2.0),
				(beamdata.xpos + beamdata.dx/2.0, beamdata.ypos + beamdata.dy/2.0),
				(beamdata.xpos + beamdata.dx/2.0, beamdata.ypos - beamdata.dy/2.0),
				(beamdata.xpos - beamdata.dx/2.0, beamdata.ypos - beamdata.dy/2.0)], layer=layer)
		cell.add(_tmp_beam_polygon)
        
def write_hole_2D(cell,beamdata,holedata,beam_dy_list,num_taper_hole,num_mir_hole_L,num_mir_hole_R,acav,amir,
	end_taper_L=False,end_taper_R=False,num_end_taper=0,reverse_tone=False):
	global num_beam, beam_spacing

	_initial_ypos = holedata.ypos - beam_dy_list[0]/2.0
	_tmp_beamdata = copy.copy(beamdata)
	_tmp_beamdata.ypos = _tmp_beamdata.ypos - beam_dy_list[0]/2.0

	for i in range(num_beam):
		holedata.ypos = _initial_ypos + i*beam_spacing + numpy.sum(beam_dy_list[:i]) + beam_dy_list[i]/2.0
		_tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i]/2.0
		_tmp_beamdata.dy = beam_dy_list[i]
		write_hole_1D(cell,_tmp_beamdata,holedata,num_taper_hole,num_mir_hole_L,num_mir_hole_R,acav,amir,
			end_taper_L,end_taper_R,num_end_taper,reverse_tone=reverse_tone)
		_tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i]/2.0 + beam_spacing

def write_hole_2D_mir_sweep(cell,beamdata,holedata,beam_dy_list,num_taper_hole,num_mir_hole_L,num_mir_hole_R,acav,amir,
	end_taper_L=False,end_taper_R=False,num_end_taper=0,reverse_tone=False):
	global num_beam, beam_spacing

	_initial_ypos = holedata.ypos - beam_dy_list[0]/2.0
	_tmp_beamdata = copy.copy(beamdata)
	_tmp_beamdata.ypos = _tmp_beamdata.ypos - beam_dy_list[0]/2.0

	for i in range(num_beam):
		num_mir = i
		holedata.ypos = _initial_ypos + i*beam_spacing + numpy.sum(beam_dy_list[:i]) + beam_dy_list[i]/2.0
		_tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i]/2.0
		_tmp_beamdata.dy = beam_dy_list[i]
		write_hole_1D_mir_sweep(cell,_tmp_beamdata,holedata,num_taper_hole,num_mir,num_mir,acav,amir,
			end_taper_L,end_taper_R,num_end_taper,reverse_tone=reverse_tone)
		_tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i]/2.0 + beam_spacing

def write_alignment_mark(cell,alignment_xpos,alignment_ypos,layer=3):
	global mark_thickness, mark_length

	cell.add(gdspy.Polygon([
		(alignment_xpos - mark_length/2.0, alignment_ypos - mark_thickness/2.0),
		(alignment_xpos - mark_length/2.0, alignment_ypos + mark_thickness/2.0),
		(alignment_xpos + mark_length/2.0, alignment_ypos + mark_thickness/2.0),
		(alignment_xpos + mark_length/2.0, alignment_ypos - mark_thickness/2.0),
		(alignment_xpos - mark_length/2.0, alignment_ypos - mark_thickness/2.0)
		], layer=layer))
	cell.add(gdspy.Polygon([
		(alignment_xpos - mark_thickness/2.0, alignment_ypos - mark_length/2.0),
		(alignment_xpos - mark_thickness/2.0, alignment_ypos + mark_length/2.0),
		(alignment_xpos + mark_thickness/2.0, alignment_ypos + mark_length/2.0),
		(alignment_xpos + mark_thickness/2.0, alignment_ypos - mark_length/2.0),
		(alignment_xpos - mark_thickness/2.0, alignment_ypos - mark_length/2.0)
		], layer=layer))

def pinch_pt_polygon(beam_dy,pinch_pt_size,pinch_pt_xloc,pinch_pt_yloc,layer=3):
	_upper_triangle = gdspy.Polygon([
		(pinch_pt_xloc, pinch_pt_yloc + pinch_pt_size/2.0),
		(pinch_pt_xloc - (beam_dy-pinch_pt_size)*math.sqrt(3)/2.0, pinch_pt_yloc + beam_dy/2.0),
		(pinch_pt_xloc + (beam_dy-pinch_pt_size)*math.sqrt(3)/2.0, pinch_pt_yloc + beam_dy/2.0),
		(pinch_pt_xloc, pinch_pt_yloc + pinch_pt_size/2.0)
		], layer=layer)
	_lower_triangle = gdspy.Polygon([
		(pinch_pt_xloc, pinch_pt_yloc - pinch_pt_size/2.0),
		(pinch_pt_xloc - (beam_dy-pinch_pt_size)*math.sqrt(3)/2.0, pinch_pt_yloc - beam_dy/2.0),
		(pinch_pt_xloc + (beam_dy-pinch_pt_size)*math.sqrt(3)/2.0, pinch_pt_yloc - beam_dy/2.0),
		(pinch_pt_xloc, pinch_pt_yloc - pinch_pt_size/2.0)
		], layer=layer)

	return _upper_triangle, _lower_triangle

def pinch_pt_polygon_vertical(beam_dy,pinch_pt_size,pinch_pt_xloc,pinch_pt_yloc,layer=3):
	_left_triangle = gdspy.Polygon([
		(pinch_pt_xloc - pinch_pt_size/2.0, pinch_pt_yloc),
		(pinch_pt_xloc - beam_dy/2.0, pinch_pt_yloc - (beam_dy-pinch_pt_size)*math.sqrt(3)/2.0),
		(pinch_pt_xloc - beam_dy/2.0, pinch_pt_yloc + (beam_dy-pinch_pt_size)*math.sqrt(3)/2.0),
		(pinch_pt_xloc - pinch_pt_size/2.0, pinch_pt_yloc)
		], layer=layer)
	_right_triangle = gdspy.Polygon([
		(pinch_pt_xloc + pinch_pt_size/2.0, pinch_pt_yloc),
		(pinch_pt_xloc + beam_dy/2.0, pinch_pt_yloc - (beam_dy-pinch_pt_size)*math.sqrt(3)/2.0),
		(pinch_pt_xloc + beam_dy/2.0, pinch_pt_yloc + (beam_dy-pinch_pt_size)*math.sqrt(3)/2.0),
		(pinch_pt_xloc + pinch_pt_size/2.0, pinch_pt_yloc)
		], layer=layer)

	return _left_triangle, _right_triangle

def linker_polygon(pltdata,layer=3):
	_linker = gdspy.Polygon([
		(pltdata.xpos - pltdata.dx/2.0, pltdata.ypos - pltdata.dy/2.0),
		(pltdata.xpos - pltdata.dx/2.0, pltdata.ypos + pltdata.dy/2.0),
		(pltdata.xpos + pltdata.dx/2.0, pltdata.ypos + pltdata.dy/2.0),
		(pltdata.xpos + pltdata.dx/2.0, pltdata.ypos - pltdata.dy/2.0),
		(pltdata.xpos - pltdata.dx/2.0, pltdata.ypos - pltdata.dy/2.0)
		], layer=layer)

	return _linker

def write_linker_region(beamdata,xmin,xmax,ymin,ymax,round_corner=False,layer=3):
	global corner_bend_rad, corner_bend_pts, linker_edgeoffset, linker_width, linker_connector_width
	global linker_notch_size, linker_beam_size, linker_xnotches, linker_ynotches

	_tmp_beamdata = copy.copy(beamdata)
	_linkerdata_inner = PLTdata()
	_linkerdata_inner.xpos = _tmp_beamdata.xpos
	_linkerdata_inner.ypos = (ymin+ymax)/2.0
	_linkerdata_inner.dx = _tmp_beamdata.dx
	_linkerdata_inner.dy = ymax - ymin - 2.0*linker_edgeoffset - 2.0*linker_connector_width
	_tmp_linker_inner = linker_polygon(_linkerdata_inner,layer)

	_linkerdata_outer = PLTdata()
	_linkerdata_outer.xpos = _tmp_beamdata.xpos
	_linkerdata_outer.ypos = (ymin+ymax)/2.0
	_linkerdata_outer.dx = _tmp_beamdata.dx + 2.0*linker_width
	_linkerdata_outer.dy = ymax - ymin - 2.0*linker_edgeoffset
	_tmp_linker_outer = linker_polygon(_linkerdata_outer,layer)

	if round_corner is True:
		_tmp_linker_inner = _tmp_linker_inner.fillet(corner_bend_rad, points_per_2pi=corner_bend_pts)
		_tmp_linker_outer = _tmp_linker_outer.fillet(corner_bend_rad, points_per_2pi=corner_bend_pts)

	_tmp_notch_yarray = (_linkerdata_outer.dy - corner_bend_rad*2 - linker_beam_size) * numpy.linspace(0,1.0,num=linker_ynotches)
	for i in range(linker_ynotches):
		_tmp_notch_ypos = ymin + linker_edgeoffset + corner_bend_rad + linker_beam_size/2.0 + _tmp_notch_yarray[i]
		_tmp_beam_polygon_L = gdspy.Polygon([
			(xmin, _tmp_notch_ypos - linker_beam_size/2.0),
			(xmin, _tmp_notch_ypos + linker_beam_size/2.0),
			(xmin + linker_edgeoffset, _tmp_notch_ypos + linker_beam_size/2.0),
			(xmin + linker_edgeoffset, _tmp_notch_ypos - linker_beam_size/2.0),
			(xmin, _tmp_notch_ypos - linker_beam_size/2.0)
			], layer=layer)
		_tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(linker_beam_size, linker_notch_size,
			xmin + linker_edgeoffset/2.0, _tmp_notch_ypos,layer)
		_tmp_beam_polygon_L = gdspy.fast_boolean(_tmp_beam_polygon_L,_tmp_upper_triangle, 'not', max_points=0, layer = layer)
		_tmp_beam_polygon_L = gdspy.fast_boolean(_tmp_beam_polygon_L,_tmp_lower_triangle, 'not', max_points=0, layer = layer)

		_tmp_beam_polygon_R = gdspy.Polygon([
			(xmax - linker_edgeoffset, _tmp_notch_ypos - linker_beam_size/2.0),
			(xmax - linker_edgeoffset, _tmp_notch_ypos + linker_beam_size/2.0),
			(xmax, _tmp_notch_ypos + linker_beam_size/2.0),
			(xmax, _tmp_notch_ypos - linker_beam_size/2.0),
			(xmax - linker_edgeoffset, _tmp_notch_ypos - linker_beam_size/2.0)
			], layer=layer)
		_tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(linker_beam_size, linker_notch_size,
			xmax - linker_edgeoffset/2.0, _tmp_notch_ypos,layer)
		_tmp_beam_polygon_R = gdspy.fast_boolean(_tmp_beam_polygon_R,_tmp_upper_triangle, 'not', max_points=0, layer = layer)
		_tmp_beam_polygon_R = gdspy.fast_boolean(_tmp_beam_polygon_R,_tmp_lower_triangle, 'not', max_points=0, layer = layer)

		_tmp_linker_outer = gdspy.fast_boolean(_tmp_linker_outer,_tmp_beam_polygon_L, 'or', max_points=0, layer = layer)
		_tmp_linker_outer = gdspy.fast_boolean(_tmp_linker_outer,_tmp_beam_polygon_R, 'or', max_points=0, layer = layer)


	_tmp_notch_xarray = (_linkerdata_outer.dx - corner_bend_rad*2 - linker_beam_size) * numpy.linspace(0,1.0,num=linker_xnotches)
	for i in range(linker_xnotches):
		_tmp_notch_xpos = xmin + linker_edgeoffset + corner_bend_rad + linker_beam_size/2.0 + _tmp_notch_xarray[i]
		_tmp_beam_polygon_T = gdspy.Polygon([
			(_tmp_notch_xpos - linker_beam_size/2.0, ymax),
			(_tmp_notch_xpos - linker_beam_size/2.0, ymax - linker_edgeoffset),
			(_tmp_notch_xpos + linker_beam_size/2.0, ymax - linker_edgeoffset),
			(_tmp_notch_xpos + linker_beam_size/2.0, ymax),
			(_tmp_notch_xpos - linker_beam_size/2.0, ymax)
			], layer=layer)
		_tmp_left_triangle, _tmp_right_triangle = pinch_pt_polygon_vertical(linker_beam_size, linker_notch_size,
			_tmp_notch_xpos, ymax - linker_edgeoffset/2.0,layer)
		_tmp_beam_polygon_T = gdspy.fast_boolean(_tmp_beam_polygon_T,_tmp_left_triangle, 'not', max_points=0, layer = layer)
		_tmp_beam_polygon_T = gdspy.fast_boolean(_tmp_beam_polygon_T,_tmp_right_triangle, 'not', max_points=0, layer = layer)

		_tmp_beam_polygon_B = gdspy.Polygon([
			(_tmp_notch_xpos - linker_beam_size/2.0, ymin),
			(_tmp_notch_xpos - linker_beam_size/2.0, ymin + linker_edgeoffset),
			(_tmp_notch_xpos + linker_beam_size/2.0, ymin + linker_edgeoffset),
			(_tmp_notch_xpos + linker_beam_size/2.0, ymin),
			(_tmp_notch_xpos - linker_beam_size/2.0, ymin)
			], layer=layer)
		_tmp_left_triangle, _tmp_right_triangle = pinch_pt_polygon_vertical(linker_beam_size, linker_notch_size,
			_tmp_notch_xpos, ymin + linker_edgeoffset/2.0,layer)
		_tmp_beam_polygon_B = gdspy.fast_boolean(_tmp_beam_polygon_B,_tmp_left_triangle, 'not', max_points=0, layer = layer)
		_tmp_beam_polygon_B = gdspy.fast_boolean(_tmp_beam_polygon_B,_tmp_right_triangle, 'not', max_points=0, layer = layer)


		_tmp_linker_outer = gdspy.fast_boolean(_tmp_linker_outer,_tmp_beam_polygon_T, 'or', max_points=0, layer = layer)
		_tmp_linker_outer = gdspy.fast_boolean(_tmp_linker_outer,_tmp_beam_polygon_B, 'or', max_points=0, layer = layer)

	_linker_combined = gdspy.fast_boolean(_tmp_linker_outer,_tmp_linker_inner, 'not', max_points=0, layer = layer)

	return _linker_combined

def write_left_circ_grating(beamdata, layer=3):
	global shot_pitch, num_circ_grating_points, grating_spacing, grating_linewidth, num_grating
	global circ_grating_base, grating_angle, odd_support_angle, even_support_angle, support_angle_width

	_philist_L = numpy.linspace(numpy.pi/2.0,numpy.pi*3.0/2.0,num_circ_grating_points-1)
	_philist_L = numpy.append(_philist_L,numpy.pi/2.0)
	_tmp_beamdata = copy.copy(beamdata)
	_ini_pt = [(_tmp_beamdata.xpos - _tmp_beamdata.dx/2, _tmp_beamdata.ypos)]

	_radius_inner = circ_grating_base

	for i in range(num_grating+1):
		_left_grating_inner = gdspy.Polygon([
			((round((_tmp_beamdata.xpos-_tmp_beamdata.dx/2)/shot_pitch) + round(_radius_inner * numpy.cos(phi)/shot_pitch)) * shot_pitch,
			(round(_tmp_beamdata.ypos/shot_pitch) + round(_radius_inner * numpy.sin(phi)/shot_pitch)) * shot_pitch)\
			for phi in _philist_L],layer=layer)

		_radius_outer = _radius_inner + grating_spacing
		_left_grating_outer = gdspy.Polygon([
			((round((_tmp_beamdata.xpos-_tmp_beamdata.dx/2)/shot_pitch) + round(_radius_outer * numpy.cos(phi)/shot_pitch)) * shot_pitch,
			(round(_tmp_beamdata.ypos/shot_pitch) + round(_radius_outer * numpy.sin(phi)/shot_pitch)) * shot_pitch)\
			for phi in _philist_L],layer=layer)

		_left_grating_tmp = gdspy.fast_boolean(_left_grating_outer,_left_grating_inner, 'not', max_points=0, layer = layer)

		if (i % 2 == 0):
			_radius_outer = _radius_outer + 10
			_philist_support = numpy.linspace(numpy.pi/2.0 + odd_support_angle - support_angle_width/2.0,
				numpy.pi/2.0 + odd_support_angle + support_angle_width/2.0,
				num_circ_grating_points-1)
			_support_pt = [
				((round((_tmp_beamdata.xpos-_tmp_beamdata.dx/2)/shot_pitch) + round(_radius_outer * numpy.cos(phi)/shot_pitch)) * shot_pitch,
				(round(_tmp_beamdata.ypos/shot_pitch) + round(_radius_outer * numpy.sin(phi)/shot_pitch)) * shot_pitch)\
				for phi in _philist_support]
			_support_pt_combined = copy.copy(_ini_pt)
			_support_pt_combined.extend(_support_pt)
			_support_pt_combined.extend(_ini_pt)
			_support_frame = gdspy.Polygon(_support_pt_combined,layer=layer)
			_left_grating_tmp = gdspy.fast_boolean(_left_grating_tmp,_support_frame, 'not', max_points=0, layer = layer)

			_philist_support = numpy.linspace(numpy.pi*3.0/2.0 - odd_support_angle - support_angle_width/2.0,
				numpy.pi*3.0/2.0 - odd_support_angle + support_angle_width/2.0,
				num_circ_grating_points-1)
			_support_pt = [
				((round((_tmp_beamdata.xpos-_tmp_beamdata.dx/2)/shot_pitch) + round(_radius_outer * numpy.cos(phi)/shot_pitch)) * shot_pitch,
				(round(_tmp_beamdata.ypos/shot_pitch) + round(_radius_outer * numpy.sin(phi)/shot_pitch)) * shot_pitch)\
				for phi in _philist_support]
			_support_pt_combined = copy.copy(_ini_pt)
			_support_pt_combined.extend(_support_pt)
			_support_pt_combined.extend(_ini_pt)
			_support_frame = gdspy.Polygon(_support_pt_combined,layer=layer)
			_left_grating_tmp = gdspy.fast_boolean(_left_grating_tmp,_support_frame, 'not', max_points=0, layer = layer)
		else:
			_radius_outer = _radius_outer + 10
			_philist_support = numpy.linspace(numpy.pi/2.0 + even_support_angle - support_angle_width/2.0,
				numpy.pi/2.0 + even_support_angle + support_angle_width/2.0,
				num_circ_grating_points-1)
			_support_pt = [
				((round((_tmp_beamdata.xpos-_tmp_beamdata.dx/2)/shot_pitch) + round(_radius_outer * numpy.cos(phi)/shot_pitch)) * shot_pitch,
				(round(_tmp_beamdata.ypos/shot_pitch) + round(_radius_outer * numpy.sin(phi)/shot_pitch)) * shot_pitch)\
				for phi in _philist_support]
			_support_pt_combined = copy.copy(_ini_pt)
			_support_pt_combined.extend(_support_pt)
			_support_pt_combined.extend(_ini_pt)
			_support_frame = gdspy.Polygon(_support_pt_combined,layer=layer)
			_left_grating_tmp = gdspy.fast_boolean(_left_grating_tmp,_support_frame, 'not', max_points=0, layer = layer)

		if i == 0:
			_left_grating = _left_grating_tmp
		else:
			_left_grating = gdspy.fast_boolean(_left_grating,_left_grating_tmp, 'or', max_points=0, layer = layer)

		_radius_inner = _radius_outer + grating_linewidth

	_philist_frame = numpy.linspace(numpy.pi/2.0+grating_angle,numpy.pi*3.0/2.0-grating_angle,num_circ_grating_points-1)
	_grating_frame_pt = [
		((round((_tmp_beamdata.xpos-_tmp_beamdata.dx/2)/shot_pitch) + round(_radius_outer * numpy.cos(phi)/shot_pitch)) * shot_pitch,
		(round(_tmp_beamdata.ypos/shot_pitch) + round(_radius_outer * numpy.sin(phi)/shot_pitch)) * shot_pitch)\
		for phi in _philist_frame]

	_grating_frame_pt_combined = copy.copy(_ini_pt)
	_grating_frame_pt_combined.extend(_grating_frame_pt)
	_grating_frame_pt_combined.extend(_ini_pt)

	_grating_frame = gdspy.Polygon(_grating_frame_pt_combined,layer=layer)

	_left_grating = gdspy.fast_boolean(_left_grating,_grating_frame, 'and', max_points=0, layer = layer)

	return _left_grating

def write_right_circ_grating(beamdata, layer=3):
	global shot_pitch, num_circ_grating_points, grating_spacing, grating_linewidth, num_grating
	global circ_grating_base, grating_angle, odd_support_angle, even_support_angle, support_angle_width

	_philist_R = numpy.linspace(-1*numpy.pi/2.0,numpy.pi/2.0,num_circ_grating_points-1)
	_philist_R = numpy.append(_philist_R,-1*numpy.pi/2.0)
	_tmp_beamdata = copy.copy(beamdata)
	_ini_pt = [(_tmp_beamdata.xpos + _tmp_beamdata.dx/2, _tmp_beamdata.ypos)]

	_radius_inner = circ_grating_base

	for i in range(num_grating+1):
		_right_grating_inner = gdspy.Polygon([
			((round((_tmp_beamdata.xpos+_tmp_beamdata.dx/2)/shot_pitch) + round(_radius_inner * numpy.cos(phi)/shot_pitch)) * shot_pitch,
			(round(_tmp_beamdata.ypos/shot_pitch) + round(_radius_inner * numpy.sin(phi)/shot_pitch)) * shot_pitch)\
			for phi in _philist_R],layer=layer)

		_radius_outer = _radius_inner + grating_spacing
		_right_grating_outer = gdspy.Polygon([
			((round((_tmp_beamdata.xpos+_tmp_beamdata.dx/2)/shot_pitch) + round(_radius_outer * numpy.cos(phi)/shot_pitch)) * shot_pitch,
			(round(_tmp_beamdata.ypos/shot_pitch) + round(_radius_outer * numpy.sin(phi)/shot_pitch)) * shot_pitch)\
			for phi in _philist_R],layer=layer)

		_right_grating_tmp = gdspy.fast_boolean(_right_grating_outer,_right_grating_inner, 'not', max_points=0, layer = layer)

		if (i % 2 == 0):
			_radius_outer = _radius_outer + 10
			_philist_support = numpy.linspace(-1*numpy.pi/2.0 + odd_support_angle - support_angle_width/2.0,
				-1*numpy.pi/2.0 + odd_support_angle + support_angle_width/2.0,
				num_circ_grating_points-1)
			_support_pt = [
				((round((_tmp_beamdata.xpos+_tmp_beamdata.dx/2)/shot_pitch) + round(_radius_outer * numpy.cos(phi)/shot_pitch)) * shot_pitch,
				(round(_tmp_beamdata.ypos/shot_pitch) + round(_radius_outer * numpy.sin(phi)/shot_pitch)) * shot_pitch)\
				for phi in _philist_support]
			_support_pt_combined = copy.copy(_ini_pt)
			_support_pt_combined.extend(_support_pt)
			_support_pt_combined.extend(_ini_pt)
			_support_frame = gdspy.Polygon(_support_pt_combined,layer=layer)
			_right_grating_tmp = gdspy.fast_boolean(_right_grating_tmp,_support_frame, 'not', max_points=0, layer = layer)

			_philist_support = numpy.linspace(numpy.pi/2.0 - odd_support_angle - support_angle_width/2.0,
				numpy.pi/2.0 - odd_support_angle + support_angle_width/2.0,
				num_circ_grating_points-1)
			_support_pt = [
				((round((_tmp_beamdata.xpos+_tmp_beamdata.dx/2)/shot_pitch) + round(_radius_outer * numpy.cos(phi)/shot_pitch)) * shot_pitch,
				(round(_tmp_beamdata.ypos/shot_pitch) + round(_radius_outer * numpy.sin(phi)/shot_pitch)) * shot_pitch)\
				for phi in _philist_support]
			_support_pt_combined = copy.copy(_ini_pt)
			_support_pt_combined.extend(_support_pt)
			_support_pt_combined.extend(_ini_pt)
			_support_frame = gdspy.Polygon(_support_pt_combined,layer=layer)
			_right_grating_tmp = gdspy.fast_boolean(_right_grating_tmp,_support_frame, 'not', max_points=0, layer = layer)
		else:
			_radius_outer = _radius_outer + 10
			_philist_support = numpy.linspace(-1*numpy.pi/2.0 + even_support_angle - support_angle_width/2.0,
				-1*numpy.pi/2.0 + even_support_angle + support_angle_width/2.0,
				num_circ_grating_points-1)
			_support_pt = [
				((round((_tmp_beamdata.xpos+_tmp_beamdata.dx/2)/shot_pitch) + round(_radius_outer * numpy.cos(phi)/shot_pitch)) * shot_pitch,
				(round(_tmp_beamdata.ypos/shot_pitch) + round(_radius_outer * numpy.sin(phi)/shot_pitch)) * shot_pitch)\
				for phi in _philist_support]
			_support_pt_combined = copy.copy(_ini_pt)
			_support_pt_combined.extend(_support_pt)
			_support_pt_combined.extend(_ini_pt)
			_support_frame = gdspy.Polygon(_support_pt_combined,layer=layer)
			_right_grating_tmp = gdspy.fast_boolean(_right_grating_tmp,_support_frame, 'not', max_points=0, layer = layer)

		if i == 0:
			_right_grating = _right_grating_tmp
		else:
			_right_grating = gdspy.fast_boolean(_right_grating,_right_grating_tmp, 'or', max_points=0, layer = layer)

		_radius_inner = _radius_outer + grating_linewidth

	_philist_frame = numpy.linspace(-1*numpy.pi/2.0+grating_angle,numpy.pi/2.0-grating_angle,num_circ_grating_points-1)
	_grating_frame_pt = [
		((round((_tmp_beamdata.xpos+_tmp_beamdata.dx/2)/shot_pitch) + round(_radius_outer * numpy.cos(phi)/shot_pitch)) * shot_pitch,
		(round(_tmp_beamdata.ypos/shot_pitch) + round(_radius_outer * numpy.sin(phi)/shot_pitch)) * shot_pitch)\
		for phi in _philist_frame]

	_grating_frame_pt_combined = copy.copy(_ini_pt)
	_grating_frame_pt_combined.extend(_grating_frame_pt)
	_grating_frame_pt_combined.extend(_ini_pt)

	_grating_frame = gdspy.Polygon(_grating_frame_pt_combined,layer=layer)

	_right_grating = gdspy.fast_boolean(_right_grating,_grating_frame, 'and', max_points=0, layer = layer)

	return _right_grating

def write_circ_grating(beamdata,beam_dy_list,circ_grating_support=None,layer=3):
	global num_beam, beam_spacing

	_tmp_beamdata = copy.copy(beamdata)
	_tmp_beamdata.ypos = _tmp_beamdata.ypos - beam_dy_list[0]/2.0
	for i in range(num_beam):
		_tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i]/2.0
		_tmp_beamdata.dy = beam_dy_list[i]

		_circ_grating_L = write_left_circ_grating(_tmp_beamdata, layer = layer)

		_circ_grating_R = write_right_circ_grating(_tmp_beamdata, layer = layer)

		if i == 0:
			_circ_grating_combined = _circ_grating_L
			_circ_grating_combined = gdspy.fast_boolean(_circ_grating_combined,_circ_grating_R, 'or', max_points=0, layer = layer)
		else:
			_circ_grating_combined = gdspy.fast_boolean(_circ_grating_combined,_circ_grating_L, 'or', max_points=0, layer = layer)
			_circ_grating_combined = gdspy.fast_boolean(_circ_grating_combined,_circ_grating_R, 'or', max_points=0, layer = layer)

		_tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i]/2.0 + beam_spacing

	return _circ_grating_combined

def write_pattern_number(cell,xloc,yloc,pattern_number,layer=3):
	global textheight, textwidth, textsep

	text_to_write = str(pattern_number)
	if pattern_number == 0:
		num_digit=1
	else:
		num_digit = numpy.ceil(math.log10(pattern_number))+1

	_pattern_number_to_write = gdspy.Text(
		text_to_write, textheight,
		(xloc - textwidth/2.0 - textsep*(num_digit-1)/4.0, yloc - textheight/2.0),
		layer=layer)

	return _pattern_number_to_write

def write_outer_box(cell,beamdata,beam_dy_list,round_corner=False,direct_write_area=False,alignment_mark=False,
	write_linker=False,pinch_pt_L=False,pinch_pt_R=False,pinch_pt_L_offset=0,pinch_pt_R_offset=0,
	pinch_pt_size=0,circ_grating=False,grating_support_size=None,pattern_number=None,reverse_tone=False,layer=3):

	global num_beam, beam_spacing, edge_offset, corner_bend_rad, corner_bend_pts, linker_edgeoffset
	global linker_width, text_dist_to_top

	_xmin = beamdata.xpos - beamdata.dx/2.0 - int(write_linker)*(linker_width+linker_edgeoffset)
	_xmax = beamdata.xpos + beamdata.dx/2.0 + int(write_linker)*(linker_width+linker_edgeoffset)
	_ymin = beamdata.ypos - beam_dy_list[0]/2.0 - edge_offset
	_ymax = _ymin + (num_beam-1)*beam_spacing + numpy.sum(beam_dy_list) + edge_offset*2

	_outer_box = gdspy.Polygon([
		(_xmin, _ymin), (_xmin, _ymax), (_xmax, _ymax), (_xmax, _ymin), (_xmin, _ymin)
		], layer=layer)

	if round_corner is True:
		_outer_box = _outer_box.fillet(corner_bend_rad, points_per_2pi=corner_bend_pts)

	if direct_write_area is True:
		_tmp_beamdata = copy.copy(beamdata)
		_tmp_beamdata.ypos = _tmp_beamdata.ypos - beam_dy_list[0]/2.0
		for i in range(num_beam):
			_tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i]/2.0
			_tmp_beamdata.dy = beam_dy_list[i]
			_tmp_beam_polygon = gdspy.Polygon([
				(_tmp_beamdata.xpos - _tmp_beamdata.dx/2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy/2.0),
				(_tmp_beamdata.xpos - _tmp_beamdata.dx/2.0, _tmp_beamdata.ypos + _tmp_beamdata.dy/2.0),
				(_tmp_beamdata.xpos + _tmp_beamdata.dx/2.0, _tmp_beamdata.ypos + _tmp_beamdata.dy/2.0),
				(_tmp_beamdata.xpos + _tmp_beamdata.dx/2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy/2.0),
				(_tmp_beamdata.xpos - _tmp_beamdata.dx/2.0, _tmp_beamdata.ypos - _tmp_beamdata.dy/2.0)
				], layer=layer)
			if pinch_pt_L is True:
				_tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(_tmp_beamdata.dy, pinch_pt_size,
					_tmp_beamdata.xpos - _tmp_beamdata.dx/2.0 + pinch_pt_L_offset, _tmp_beamdata.ypos,layer)
				_tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon,_tmp_upper_triangle, 'not', max_points=0, layer = layer)
				_tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon,_tmp_lower_triangle, 'not', max_points=0, layer = layer)
			if pinch_pt_R is True:
				_tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(_tmp_beamdata.dy, pinch_pt_size,
					_tmp_beamdata.xpos + _tmp_beamdata.dx/2.0 - pinch_pt_R_offset, _tmp_beamdata.ypos,layer)
				_tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon,_tmp_upper_triangle, 'not', max_points=0, layer = layer)
				_tmp_beam_polygon = gdspy.fast_boolean(_tmp_beam_polygon,_tmp_lower_triangle, 'not', max_points=0, layer = layer)

			_tmp_beamdata.ypos = _tmp_beamdata.ypos + beam_dy_list[i]/2.0 + beam_spacing
			if i == 0:
				_tmp_beam_polygon_combined = _tmp_beam_polygon
			else:
				_tmp_beam_polygon_combined = gdspy.fast_boolean(_tmp_beam_polygon_combined,_tmp_beam_polygon, 'or', max_points=0, layer = layer)

		_box_write_area = gdspy.fast_boolean(_outer_box,_tmp_beam_polygon_combined, 'not', max_points=0, layer = layer)

		if write_linker is True:
			_linker_combined = write_linker_region(beamdata,_xmin,_xmax,_ymin,_ymax,round_corner,layer = layer)
			_box_write_area = gdspy.fast_boolean(_box_write_area,_linker_combined, 'not', max_points=0, layer = layer)

		if circ_grating is True:
			_circ_grating_combined = write_circ_grating(beamdata,beam_dy_list,grating_support_size,layer=3)
			_box_write_area = gdspy.fast_boolean(_box_write_area,_circ_grating_combined, 'or', max_points=0, layer = layer)

		if pattern_number is not None:
			_left_pattern_number = write_pattern_number(cell,_xmin + linker_edgeoffset + linker_width/2.0,
				_ymax - linker_edgeoffset - text_dist_to_top, pattern_number)
			_right_pattern_number = write_pattern_number(cell,_xmax - linker_edgeoffset - linker_width/2.0,
				_ymax - linker_edgeoffset - text_dist_to_top, pattern_number)
			_pattern_number_combined = gdspy.fast_boolean(_left_pattern_number,_right_pattern_number, 'or', max_points=0, layer = layer)
			_box_write_area = gdspy.fast_boolean(_box_write_area,_pattern_number_combined, 'xor', max_points=0, layer = layer)

	if reverse_tone is True:
		_box_write_area = gdspy.fast_boolean(_box_write_area,_tmp_beam_polygon_combined, 'or', max_points=0, layer = layer)
		_box_write_area_reverse = gdspy.fast_boolean(_outer_box,_box_write_area, 'not', max_points=0, layer = layer)
		cell.add(_box_write_area_reverse)

	else:
		cell.add(_box_write_area)

	if alignment_mark:
		_yoffset = 10000
		_alignment_x = beamdata.xpos
		_alignment_y = _ymax + _yoffset
		write_alignment_mark(cell,_alignment_x,_alignment_y)

def write_support_region(innerframe,layer=4):
	global corner_bend_rad, support_connector_width, support_notch_size, support_beam_size
	global num_xsupport, num_ysupport, write_field_x_size, write_field_y_size

	_support_box_outer = copy.copy(innerframe)
	_ymin = _support_box_outer.ypos - _support_box_outer.dy/2.0
	_ymax = _support_box_outer.ypos + _support_box_outer.dy/2.0
	_xmin = _support_box_outer.xpos - _support_box_outer.dx/2.0
	_xmax = _support_box_outer.xpos + _support_box_outer.dx/2.0
	_tmp_support_outer = linker_polygon(_support_box_outer,layer)

	_support_box_inner = copy.copy(_support_box_outer)
	_support_box_inner.dx = _support_box_inner.dx - 2*overlap_width
	_support_box_inner.dy = _support_box_inner.dy - 2*overlap_width
	_tmp_support_inner = linker_polygon(_support_box_inner,layer)

	_write_field = copy.copy(_support_box_outer)
	_write_field.dx = write_field_x_size
	_write_field.dy = write_field_y_size
	_tmp_write_field = linker_polygon(_write_field,layer)

	_support_outline = copy.copy(_support_box_outer)
	_support_outline.dx = _support_outline.dx + 2*support_connector_width
	_support_outline.dy = _support_outline.dy + 2*support_connector_width
	_tmp_support_outline = linker_polygon(_support_outline,layer)

	_box_write_area = gdspy.fast_boolean(_tmp_support_outer,_tmp_support_inner, 'not', max_points=0, layer = layer)

	_tmp_notch_yarray = (_support_box_outer.dy - corner_bend_rad*2 - support_beam_size) * numpy.linspace(0,1.0,num=num_ysupport)
	for i in range(num_ysupport):
		_tmp_notch_ypos = _ymin + corner_bend_rad + support_beam_size/2.0 + _tmp_notch_yarray[i]
		_tmp_beam_polygon_L = gdspy.Polygon([
			(_xmin - support_connector_width, _tmp_notch_ypos - support_beam_size/2.0),
			(_xmin - support_connector_width, _tmp_notch_ypos + support_beam_size/2.0),
			(_xmin, _tmp_notch_ypos + support_beam_size/2.0),
			(_xmin, _tmp_notch_ypos - support_beam_size/2.0),
			(_xmin - support_connector_width, _tmp_notch_ypos - support_beam_size/2.0)
			], layer=layer)
		_tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(support_beam_size, support_notch_size,
			_xmin - support_connector_width/2.0, _tmp_notch_ypos,layer)
		_tmp_beam_polygon_L = gdspy.fast_boolean(_tmp_beam_polygon_L,_tmp_upper_triangle, 'not', max_points=0, layer = layer)
		_tmp_beam_polygon_L = gdspy.fast_boolean(_tmp_beam_polygon_L,_tmp_lower_triangle, 'not', max_points=0, layer = layer)

		_tmp_beam_polygon_R = gdspy.Polygon([
			(_xmax, _tmp_notch_ypos - support_beam_size/2.0),
			(_xmax, _tmp_notch_ypos + support_beam_size/2.0),
			(_xmax + support_connector_width, _tmp_notch_ypos + support_beam_size/2.0),
			(_xmax + support_connector_width, _tmp_notch_ypos - support_beam_size/2.0),
			(_xmax, _tmp_notch_ypos - support_beam_size/2.0)
			], layer=layer)
		_tmp_upper_triangle, _tmp_lower_triangle = pinch_pt_polygon(support_beam_size, support_notch_size,
			_xmax + support_connector_width/2.0, _tmp_notch_ypos,layer)
		_tmp_beam_polygon_R = gdspy.fast_boolean(_tmp_beam_polygon_R,_tmp_upper_triangle, 'not', max_points=0, layer = layer)
		_tmp_beam_polygon_R = gdspy.fast_boolean(_tmp_beam_polygon_R,_tmp_lower_triangle, 'not', max_points=0, layer = layer)

		_box_write_area = gdspy.fast_boolean(_box_write_area,_tmp_beam_polygon_L, 'or', max_points=0, layer = layer)
		_box_write_area = gdspy.fast_boolean(_box_write_area,_tmp_beam_polygon_R, 'or', max_points=0, layer = layer)

	_tmp_notch_xarray = (_support_box_outer.dx - corner_bend_rad*2 - support_beam_size) * numpy.linspace(0,1.0,num=num_xsupport)
	for i in range(num_xsupport):
		_tmp_notch_xpos = _xmin + corner_bend_rad + support_beam_size/2.0 + _tmp_notch_xarray[i]
		_tmp_beam_polygon_T = gdspy.Polygon([
			(_tmp_notch_xpos - support_beam_size/2.0, _ymax),
			(_tmp_notch_xpos - support_beam_size/2.0, _ymax + support_connector_width),
			(_tmp_notch_xpos + support_beam_size/2.0, _ymax + support_connector_width),
			(_tmp_notch_xpos + support_beam_size/2.0, _ymax),
			(_tmp_notch_xpos - support_beam_size/2.0, _ymax)
			], layer=layer)
		_tmp_left_triangle, _tmp_right_triangle = pinch_pt_polygon_vertical(support_beam_size, support_notch_size,
			_tmp_notch_xpos, _ymax + support_connector_width/2.0,layer)
		_tmp_beam_polygon_T = gdspy.fast_boolean(_tmp_beam_polygon_T,_tmp_left_triangle, 'not', max_points=0, layer = layer)
		_tmp_beam_polygon_T = gdspy.fast_boolean(_tmp_beam_polygon_T,_tmp_right_triangle, 'not', max_points=0, layer = layer)

		_tmp_beam_polygon_B = gdspy.Polygon([
			(_tmp_notch_xpos - support_beam_size/2.0, _ymin - support_connector_width),
			(_tmp_notch_xpos - support_beam_size/2.0, _ymin),
			(_tmp_notch_xpos + support_beam_size/2.0, _ymin),
			(_tmp_notch_xpos + support_beam_size/2.0, _ymin - support_connector_width),
			(_tmp_notch_xpos - support_beam_size/2.0, _ymin - support_connector_width)
			], layer=layer)
		_tmp_left_triangle, _tmp_right_triangle = pinch_pt_polygon_vertical(support_beam_size, support_notch_size,
			_tmp_notch_xpos, _ymin - support_connector_width/2.0,layer)
		_tmp_beam_polygon_B = gdspy.fast_boolean(_tmp_beam_polygon_B,_tmp_left_triangle, 'not', max_points=0, layer = layer)
		_tmp_beam_polygon_B = gdspy.fast_boolean(_tmp_beam_polygon_B,_tmp_right_triangle, 'not', max_points=0, layer = layer)


		_box_write_area = gdspy.fast_boolean(_box_write_area,_tmp_beam_polygon_T, 'or', max_points=0, layer = layer)
		_box_write_area = gdspy.fast_boolean(_box_write_area,_tmp_beam_polygon_B, 'or', max_points=0, layer = layer)


	_support_combined = gdspy.fast_boolean(_tmp_write_field,_tmp_support_outline, 'not', max_points=0, layer = layer)
	_support_combined = gdspy.fast_boolean(_support_combined,_box_write_area, 'or', max_points=0, layer = layer)

	return _support_combined


def write_outer_frame(cell,beamdata,beam_dy_list,pattern_number=None,reverse_tone=False,layer=4):

	global num_beam, beam_spacing, edge_offset, linker_edgeoffset, linker_width, text_dist_to_top

	_ymin = beamdata.ypos - beam_dy_list[0]/2.0 - edge_offset
	_ymax = _ymin + (num_beam-1)*beam_spacing + numpy.sum(beam_dy_list) + edge_offset*2

	_tmp_beamdata = copy.copy(beamdata)
	_support_inner = PLTdata()
	_support_inner.xpos = _tmp_beamdata.xpos
	_support_inner.ypos = (_ymin+_ymax)/2.0
	_support_inner.dx = beamdata.dx + 2*linker_width + 2*linker_edgeoffset
	_support_inner.dy = _ymax - _ymin

	_frame_write_area = write_support_region(_support_inner,layer)

	if pattern_number is not None:
		_pattern_number_to_write = write_pattern_number(cell,_support_inner.xpos,
			_support_inner.ypos + _support_inner.dy/2.0 + 5e3, pattern_number,layer)

		_frame_write_area = gdspy.fast_boolean(_frame_write_area,_pattern_number_to_write, 'not', max_points=0, layer = layer)

	cell.add(_frame_write_area)

class PLTdata:
	def __init__(self):
		self.xpos = None
		self.ypos = None
		self.dx = None
		self.dy = None

beams = gdspy.Cell('beams')

# Unit for all parameters: nm

# Parameters for exposure
write_field_size = 100e3
write_field_x_size = 50e3
write_field_y_size = 90e3
shot_pitch = 0.1
pattern_sep = 100e3

# Parameters for nanobeam
num_beam = 11
pinch_pt_size = [100]
beam_spacing = 5e3
num_circ_points = 181
hole_pos_offset = 0
edge_offset = 5e3
corner_bend_rad = 3e3
corner_bend_pts = 41
acav = 165
amir = 177
wy_list = [400]
wy_end_taper_list = []
hx_list = [75,80,85,90,95]
hy_list = [224]
mirror_list = [2, 4, 8, 10, 12, 16]
taper_list = [16]

# Parameters for alignment marks
mark_thickness = 1e3
mark_length = 5e3

# Parameters for the linker
linker_edgeoffset = 0
linker_width = 6e3
linker_connector_width = 1e3
linker_notch_size = 500
linker_beam_size = 500
linker_xnotches = 0
linker_ynotches = 0

# Parameter for the support
support_connector_width = 1e3
support_notch_size = 400
support_beam_size = 1e3
num_xsupport = 2
num_ysupport = 10
overlap_width = 0

# Parameters for circular grating structure
circ_grating_support = [180]
num_circ_grating_points = 181
grating_spacing = 280
grating_linewidth = 120 # design width = 120nm
num_grating = 4
circ_grating_base = 640
grating_angle = 15*numpy.pi/180 # 15 degree in radian
odd_support_angle = 55*numpy.pi/180 # 55 degree in radian
even_support_angle = numpy.pi/2.0 # 90 degree in radian
support_angle_width = 10*numpy.pi/180 # 10 degree in radian

# Parameters for the text
textheight = 5e3
textwidth = textheight*5.0/9.0
textsep = textheight*8.0/9.0
text_dist_to_top = 3e3


matrix_y_size = len(hy_list)
matrix_x_size = len(mirror_list)
block_y_size = len(hx_list)
block_x_size = len(taper_list)
block_x_sep = pattern_sep * matrix_x_size
block_y_sep = (150e3) * matrix_y_size
k = 0

origin_x = 0
origin_y = 0

for m in range(block_y_size):
    origin_y = origin_y + m * block_y_sep
    for n in range(block_x_size):
        origin_x = origin_x + n * block_x_sep
        for i in range(matrix_y_size):
            for j in range(matrix_x_size):
                num_tap = taper_list[n]
                num_mirr = mirror_list[j]
                rowicolj = PLTdata()
                rowicolj.dx = 10e3
                rowicolj.dy = wy_list[0]
                rowicolj.xpos = write_field_size/2 + pattern_sep*j + origin_x
                rowicolj.ypos = write_field_size/2 - (beam_spacing+rowicolj.dy)*(num_beam-1)/2 + (150e3)*i +origin_y
                dy_list = write_beam_array(beams,rowicolj)

                circ_rowicolj = PLTdata()
                circ_rowicolj.xpos = rowicolj.xpos - hole_pos_offset
                circ_rowicolj.ypos = rowicolj.ypos
                circ_rowicolj.dx = hx_list[m]
                circ_rowicolj.dy = hy_list[i]
                write_hole_2D(beams,rowicolj,circ_rowicolj,dy_list,num_taper_hole=num_tap,num_mir_hole_L=num_mirr/2,
                    num_mir_hole_R=num_mirr/2,acav=acav,amir=amir,num_end_taper=3,end_taper_L=False,end_taper_R=False,reverse_tone=True)

                write_outer_box(beams,rowicolj,dy_list,round_corner=False,direct_write_area=True,write_linker=True,
                        circ_grating=True,reverse_tone=True)
                write_outer_frame(beams,rowicolj,dy_list,k)

                k = k+1

        origin_x = 0
    origin_y = 0
    
    
#Write mirror sweep patterns
cavity_period = acav
mirror_period = amir
num_rows = matrix_y_size * block_y_size

hole_spacings = numpy.linspace(cavity_period, mirror_period, num=num_rows)

n = 0
for a in hole_spacings:
    
    
    rowicolj = PLTdata()
    rowicolj.dx = 10e3
    rowicolj.dy = wy_list[0]
    rowicolj.xpos = write_field_size/2 + pattern_sep*matrix_x_size
    rowicolj.ypos = write_field_size/2 - (beam_spacing+rowicolj.dy)*(num_beam-1)/2 + (150e3)*n + origin_y
    dy_list = write_beam_array(beams,rowicolj)

    circ_rowicolj = PLTdata()
    circ_rowicolj.xpos = rowicolj.xpos - hole_pos_offset
    circ_rowicolj.ypos = rowicolj.ypos
    circ_rowicolj.dx = hx_list[2]
    circ_rowicolj.dy = hy_list[0]
                
    write_hole_2D_mir_sweep(beams,rowicolj,circ_rowicolj,dy_list,num_taper_hole=0,num_mir_hole_L=num_mirr/2,
                        num_mir_hole_R=num_mirr/2,acav=cavity_period,amir=a,num_end_taper=0,end_taper_L=False,end_taper_R=False,reverse_tone=True)
    write_outer_box(beams,rowicolj,dy_list,round_corner=False,direct_write_area=True,write_linker=True,
                circ_grating=True,reverse_tone=True)
    write_outer_frame(beams,rowicolj,dy_list,k)
    
    k = k+1
    n = n+1

gdspy.write_gds('20221005_SiV0PhC_wz220_hx75to95_hy224_wy400_NumTap16.gds', unit=1.0e-9, precision=1.0e-11)
