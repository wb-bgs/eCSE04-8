/*--------------------------------------------------------------------
 *
 *	Copyright (c) 1991-2021 by the GMT Team (https://www.generic-mapping-tools.org/team.html)
 *	See LICENSE.TXT file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU Lesser General Public License as published by
 *	the Free Software Foundation; version 3 or any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU Lesser General Public License for more details.
 *
 *	Contact info: www.generic-mapping-tools.org
 *--------------------------------------------------------------------*/
/*
 *
 *			G M T _ P L O T . C
 *
 *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * GMT_plot.c contains code related to plotting maps.  These functions requires
 * we pass both the GMT and PSL control structure pointers.
 *
 * Author:	Paul Wessel
 * Date:	1-JAN-2010
 * Version:	5
 *
 * PUBLIC Functions include (39):
 *
 *	gmt_draw_map_scale 	 : Plot map scale
 *	gmt_draw_map_rose 	 : Plot map rose
 *	gmt_draw_map_panel 	 : Plot map panel
 *	gmt_draw_front 		 : Draw a front line
 *	gmt_geo_line 		 : Plots line in lon/lat on maps, takes care of periodicity jumps
 *	gmt_plot_geo_ellipse 	 : Plots ellipse in lon/lat on maps, takes care of periodicity jumps
 *	gmt_geo_polygons 	 : Plots polygon in lon/lat on maps, takes care of periodicity jumps
 *	gmt_geo_rectangle 	 : Plots rectangle in lat/lon on maps, takes care of periodicity jumps
 *	gmt_map_basemap 	 : Generic basemap function
 *	gmt_map_clip_off 	 : Deactivate map region clip path
 *	gmt_map_clip_on 	 : Activate map region clip path
 *	gmt_BB_clip_on		 : Activate Bounding Box clip path
 *	gmt_plane_perspective	 : Adds PS matrix to simulate perspective plotting
 *	gmt_plot_line 		 : Plots path (in projected coordinates), takes care of boundary crossings
 *	gmt_vertical_axis 	 : Draw 3-D vertical axes
 *	gmt_xy_axis 		 : Draw x or y axis
 *	gmt_linearx_grid 	 : Draw linear x grid lines
 *	gmt_setfill              :
 *	gmt_setfont              :
 *	gmt_draw_map_inset       :
 *	gmt_setpen               :
 *	gmt_draw_custom_symbol   :
 *	gmt_add_label_record     :
 *	gmt_contlabel_save_begin :
 *	gmt_contlabel_save_end   :
 *	gmt_textpath_init        :
 *	gmt_contlabel_plot       :
 *	gmt_export2proj4         :
 *	gmt_plotinit             :
 *	gmt_plotcanvas           :
 *	gmt_plotend              :
 *	gmt_geo_polarcap_segment :
 *	gmt_geo_vector           :
 *	gmtlib_create_ps            :
 *	gmtlib_free_ps_ptr          :
 *	gmtlib_free_ps              :
 *	gmtlib_read_ps              :
 *	gmtlib_write_ps             :
 *	gmtlib_duplicate_ps         :
 *
 */

#include "gmt_dev.h"
#include "gmt_internals.h"
#ifdef _WIN32
#include <windows.h>
#include <tlhelp32.h>
#endif

#define gmt_M_axis_is_geo_strict(C,axis) (((axis) == GMT_X && gmt_M_type (C, GMT_IN, axis) & GMT_IS_LON) || ((axis) == GMT_Y && gmt_M_type (C, GMT_IN, axis) & GMT_IS_LAT))
#define PSL_IZ(PSL,z) ((int)lrint ((z) * PSL->internal.dpu))

/* Local variables to this file */

static unsigned int GMT_n_annotations_skip[4] = {0, 0, 0, 0};
static size_t GMT_n_annotations[4] = {0, 0, 0, 0};
static size_t GMT_alloc_annotations[4] = {0, 0, 0, 0};
static double *GMT_x_annotation[4] = {NULL, NULL, NULL, NULL}, *GMT_y_annotation[4] = {NULL, NULL, NULL, NULL};

/* THese functions are public but used in a static function so declared here to avoid resorting */
void gmt_linearx_grid (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, double dval);

/* Get bitmapped 600 dpi GMT glyph for timestamp.  The glyph is a 90 x 220 pixel 1-bit image
	and it is here represented as ceil (220 / 8) * 90 = 2520 bytes */

static unsigned char GMT_glyph[2520] = {
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x7f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x07, 0xff, 0xff, 0x80, 0x00,
0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x1f, 0xff, 0xff, 0xf8, 0x00, 0x00, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x3f, 0xff, 0xff, 0xff, 0xc0, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x7f,
0xff, 0xff, 0xff, 0xfc, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
0xe0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff,
0xc0, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3f, 0xff, 0xfc, 0x00, 0x00, 0x01,
0xff, 0xc0, 0x00, 0x00, 0x01, 0xff, 0xe0, 0xff, 0xff, 0xff, 0xff, 0xff, 0xf8, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0xff, 0xff, 0xff, 0x00, 0x00, 0x01, 0xff, 0xe0, 0x00, 0x00,
0x01, 0xff, 0xf0, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00,
0x00, 0x0f, 0xff, 0xff, 0xff, 0x80, 0x00, 0x03, 0xff, 0xe0, 0x00, 0x00, 0x00, 0xff, 0xf8, 0x7f,
0xff, 0xff, 0xff, 0xff, 0xff, 0x80, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3f, 0xff, 0xff,
0xff, 0xc0, 0x00, 0x03, 0xff, 0xf0, 0x00, 0x00, 0x00, 0xff, 0xfc, 0x3f, 0xff, 0xff, 0xff, 0xff,
0xff, 0xe0, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x01, 0xff, 0xff, 0xff, 0xff, 0xc0, 0x00, 0x03,
0xff, 0xf0, 0x00, 0x00, 0x00, 0xff, 0xfc, 0x1f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xf0, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x00, 0x07, 0xff, 0xff, 0xff, 0xff, 0xe0, 0x00, 0x07, 0xff, 0xf8, 0x00, 0x00,
0x00, 0xff, 0xfe, 0x0f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfc, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00,
0x1f, 0xff, 0xfc, 0x3f, 0xff, 0xe0, 0x00, 0x07, 0xff, 0xf8, 0x00, 0x00, 0x00, 0xff, 0xff, 0x07,
0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x7f, 0xff, 0x80, 0x07,
0xff, 0xf0, 0x00, 0x07, 0xff, 0xfc, 0x00, 0x00, 0x00, 0xff, 0xff, 0x07, 0x80, 0x00, 0xff, 0xfc,
0x00, 0x0f, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x01, 0xff, 0xfc, 0x00, 0x01, 0xff, 0xf8, 0x00, 0x07,
0xff, 0xfc, 0x00, 0x00, 0x00, 0xff, 0xff, 0x80, 0x00, 0x00, 0x1f, 0xfc, 0x00, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x07, 0xff, 0xe0, 0x00, 0x00, 0x7f, 0xf8, 0x00, 0x07, 0xff, 0xfc, 0x00, 0x00,
0x00, 0xff, 0xff, 0xc0, 0x00, 0x00, 0x0f, 0xfe, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x1f,
0xff, 0x00, 0x00, 0x00, 0x3f, 0xf8, 0x00, 0x0f, 0xff, 0xfe, 0x00, 0x00, 0x00, 0xff, 0xff, 0xc0,
0x00, 0x00, 0x07, 0xff, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x7f, 0xfc, 0x00, 0x00, 0x00,
0x3f, 0xf8, 0x00, 0x0f, 0xff, 0xfe, 0x00, 0x00, 0x00, 0xff, 0xff, 0xe0, 0x00, 0x00, 0x03, 0xff,
0x80, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0xff, 0xf8, 0x00, 0x00, 0x00, 0x1f, 0xf8, 0x00, 0x0f,
0xff, 0xff, 0x00, 0x00, 0x00, 0xff, 0xff, 0xe0, 0x00, 0x00, 0x01, 0xff, 0xc0, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x03, 0xff, 0xe0, 0x00, 0x00, 0x00, 0x1f, 0xf8, 0x00, 0x0f, 0xff, 0xff, 0x00, 0x00,
0x00, 0xff, 0xbf, 0xf0, 0x00, 0x00, 0x00, 0xff, 0xc0, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x07, 0xff,
0x80, 0x00, 0x00, 0x00, 0x1f, 0xf8, 0x00, 0x0f, 0xfd, 0xff, 0x80, 0x00, 0x01, 0xff, 0xbf, 0xf0,
0x00, 0x00, 0x00, 0x7f, 0xe0, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x0f, 0xff, 0x00, 0x00, 0x00, 0x00,
0x1f, 0xf8, 0x00, 0x0f, 0xfd, 0xff, 0x80, 0x00, 0x01, 0xff, 0x9f, 0xf8, 0x00, 0x00, 0x00, 0x7f,
0xf0, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x1f, 0xfc, 0x00, 0x00, 0x00, 0x00, 0x1f, 0xf8, 0x00, 0x0f,
0xfd, 0xff, 0x80, 0x00, 0x01, 0xff, 0x9f, 0xf8, 0x00, 0x00, 0x00, 0x3f, 0xf8, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x3f, 0xf8, 0x00, 0x00, 0x00, 0x00, 0x1f, 0xf8, 0x00, 0x1f, 0xfc, 0xff, 0xc0, 0x00,
0x01, 0xff, 0x9f, 0xf8, 0x00, 0x00, 0x00, 0x1f, 0xf8, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x7f, 0xf0,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x1f, 0xf8, 0xff, 0xc0, 0x00, 0x01, 0xff, 0x8f, 0xfc,
0x00, 0x00, 0x00, 0x1f, 0xfc, 0x00, 0x00, 0x0f, 0x00, 0x00, 0xff, 0xe0, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x1f, 0xf8, 0xff, 0xe0, 0x00, 0x01, 0xff, 0x8f, 0xfc, 0x00, 0x00, 0x00, 0x0f,
0xfc, 0x00, 0x00, 0x0f, 0x00, 0x01, 0xff, 0xc0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x1f,
0xf8, 0x7f, 0xe0, 0x00, 0x03, 0xff, 0x8f, 0xfc, 0x00, 0x00, 0x00, 0x0f, 0xfe, 0x00, 0x00, 0x0f,
0x00, 0x01, 0xff, 0x80, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x1f, 0xf8, 0x7f, 0xe0, 0x00,
0x03, 0xff, 0x8f, 0xfc, 0x00, 0x00, 0x00, 0x07, 0xfe, 0x00, 0x00, 0x0f, 0x00, 0x03, 0xff, 0x80,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x1f, 0xf8, 0x7f, 0xf0, 0x00, 0x03, 0xff, 0x0f, 0xfe,
0x00, 0x00, 0x00, 0x07, 0xff, 0x00, 0x00, 0x0f, 0x00, 0x07, 0xff, 0x80, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x1f, 0xf8, 0x3f, 0xf0, 0x00, 0x03, 0xff, 0x07, 0xfe, 0x00, 0x00, 0x00, 0x07,
0xff, 0x00, 0x00, 0x0f, 0x00, 0x07, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x1f,
0xf8, 0x3f, 0xf8, 0x00, 0x07, 0xff, 0x07, 0xfe, 0x00, 0x00, 0x00, 0x07, 0xff, 0x00, 0x00, 0x0f,
0x00, 0x07, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x1f, 0xf8, 0x3f, 0xf8, 0x00,
0x07, 0xff, 0x07, 0xfe, 0x00, 0x00, 0x00, 0x03, 0xff, 0x00, 0x00, 0x0f, 0x00, 0x07, 0xff, 0x00,
0x00, 0x03, 0xff, 0xff, 0xff, 0xf8, 0x00, 0x1f, 0xf8, 0x1f, 0xf8, 0x00, 0x07, 0xfe, 0x07, 0xfe,
0x00, 0x00, 0x00, 0x03, 0xff, 0x00, 0x00, 0x0f, 0x00, 0x07, 0xff, 0x00, 0x00, 0x03, 0xff, 0xff,
0xff, 0xf8, 0x00, 0x1f, 0xf8, 0x1f, 0xfc, 0x00, 0x0f, 0xfe, 0x07, 0xfe, 0x00, 0x00, 0x00, 0x03,
0xff, 0x00, 0x00, 0x0f, 0x00, 0x07, 0xff, 0x00, 0x00, 0x03, 0xff, 0xff, 0xff, 0xf8, 0x00, 0x1f,
0xf8, 0x1f, 0xfc, 0x00, 0x0f, 0xfe, 0x07, 0xfe, 0x00, 0x00, 0x00, 0x03, 0xff, 0x00, 0x00, 0x0f,
0x00, 0x07, 0xff, 0x00, 0x00, 0x03, 0xff, 0xff, 0xff, 0xf8, 0x00, 0x1f, 0xf8, 0x0f, 0xfc, 0x00,
0x1f, 0xfc, 0x07, 0xfe, 0x00, 0x00, 0x00, 0x03, 0xff, 0x00, 0x00, 0x0f, 0x00, 0x07, 0xff, 0x00,
0x00, 0x03, 0xff, 0xff, 0xff, 0xfc, 0x00, 0x1f, 0xf8, 0x0f, 0xfe, 0x00, 0x1f, 0xfc, 0x07, 0xfe,
0x00, 0x00, 0x00, 0x03, 0xff, 0x00, 0x00, 0x0f, 0x00, 0x07, 0xff, 0x00, 0x00, 0x03, 0xff, 0xff,
0xff, 0xfc, 0x00, 0x1f, 0xf8, 0x0f, 0xfe, 0x00, 0x1f, 0xf8, 0x07, 0xfe, 0x00, 0x00, 0x00, 0x07,
0xff, 0x00, 0x00, 0x0f, 0x00, 0x03, 0xff, 0x80, 0x00, 0x03, 0xff, 0xff, 0xff, 0xfc, 0x00, 0x1f,
0xf8, 0x07, 0xfe, 0x00, 0x3f, 0xf8, 0x0f, 0xfe, 0x00, 0x00, 0x00, 0x07, 0xff, 0x00, 0x00, 0x0f,
0x00, 0x03, 0xff, 0x80, 0x00, 0x00, 0x03, 0xf8, 0x0f, 0xfc, 0x00, 0x1f, 0xf8, 0x07, 0xff, 0x00,
0x3f, 0xf0, 0x0f, 0xfc, 0x00, 0x00, 0x00, 0x07, 0xfe, 0x00, 0x00, 0x0f, 0x00, 0x03, 0xff, 0x80,
0x00, 0x00, 0x00, 0x00, 0x0f, 0xfc, 0x00, 0x1f, 0xf8, 0x07, 0xff, 0x00, 0x7f, 0xf0, 0x0f, 0xfc,
0x00, 0x00, 0x00, 0x0f, 0xfe, 0x00, 0x00, 0x0f, 0x00, 0x01, 0xff, 0xc0, 0x00, 0x00, 0x00, 0x00,
0x0f, 0xfc, 0x00, 0x1f, 0xf8, 0x03, 0xff, 0x00, 0x7f, 0xe0, 0x0f, 0xfc, 0x00, 0x00, 0x00, 0x0f,
0xfe, 0x00, 0x00, 0x0f, 0x00, 0x00, 0xff, 0xe0, 0x00, 0x00, 0x00, 0x00, 0x0f, 0xfc, 0x00, 0x1f,
0xf8, 0x03, 0xff, 0x80, 0xff, 0xe0, 0x0f, 0xfc, 0x00, 0x00, 0x00, 0x1f, 0xfc, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x7f, 0xf0, 0x00, 0x00, 0x00, 0x00, 0x0f, 0xfe, 0x00, 0x1f, 0xfc, 0x03, 0xff, 0x80,
0xff, 0xc0, 0x1f, 0xf8, 0x00, 0x00, 0x00, 0x1f, 0xfc, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x3f, 0xf0,
0x00, 0x00, 0x00, 0x00, 0x0f, 0xfe, 0x00, 0x1f, 0xfc, 0x01, 0xff, 0x81, 0xff, 0xc0, 0x1f, 0xf8,
0x00, 0x00, 0x00, 0x3f, 0xf8, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x1f, 0xf8, 0x00, 0x00, 0x00, 0x00,
0x07, 0xfe, 0x00, 0x0f, 0xfc, 0x01, 0xff, 0x81, 0xff, 0x80, 0x1f, 0xf8, 0x00, 0x00, 0x00, 0x3f,
0xf0, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x0f, 0xfe, 0x00, 0x00, 0x00, 0x00, 0x07, 0xfe, 0x00, 0x0f,
0xfc, 0x00, 0xff, 0xc3, 0xff, 0x80, 0x3f, 0xf0, 0x00, 0x00, 0x00, 0x7f, 0xf0, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x07, 0xff, 0x00, 0x00, 0x00, 0x00, 0x07, 0xff, 0x00, 0x0f, 0xfc, 0x00, 0xff, 0xc3,
0xff, 0x00, 0x3f, 0xf0, 0x00, 0x00, 0x00, 0xff, 0xe0, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x03, 0xff,
0xc0, 0x00, 0x00, 0x00, 0x0f, 0xff, 0x00, 0x0f, 0xfc, 0x00, 0xff, 0xc7, 0xfe, 0x00, 0x3f, 0xe0,
0x00, 0x00, 0x01, 0xff, 0xc0, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x01, 0xff, 0xe0, 0x00, 0x00, 0x00,
0x0f, 0xff, 0x00, 0x0f, 0xfc, 0x00, 0x7f, 0xcf, 0xfe, 0x00, 0x7f, 0xe0, 0x00, 0x00, 0x03, 0xff,
0x80, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x7f, 0xf8, 0x00, 0x00, 0x00, 0x0f, 0xff, 0x80, 0x0f,
0xfc, 0x00, 0x7f, 0xef, 0xfc, 0x00, 0x7f, 0xc0, 0x00, 0x00, 0x07, 0xff, 0x00, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x1f, 0xfe, 0x00, 0x00, 0x00, 0x1f, 0xff, 0x80, 0x07, 0xfc, 0x00, 0x7f, 0xff,
0xf8, 0x00, 0xff, 0xc0, 0x00, 0x00, 0x0f, 0xfe, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x0f,
0xff, 0x80, 0x00, 0x00, 0x1f, 0xff, 0x80, 0x07, 0xfc, 0x00, 0x7f, 0xff, 0xf0, 0x00, 0xff, 0x80,
0x00, 0x00, 0x0f, 0xfc, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x03, 0xff, 0xf0, 0x00, 0x00,
0x3f, 0xff, 0xc0, 0x07, 0xfc, 0x00, 0x3f, 0xff, 0xf0, 0x01, 0xff, 0x80, 0x00, 0x00, 0x3f, 0xf8,
0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0xff, 0xfe, 0x00, 0x00, 0xff, 0xff, 0xc0, 0x07,
0xfc, 0x00, 0x3f, 0xff, 0xe0, 0x01, 0xff, 0x00, 0x00, 0x00, 0x7f, 0xf0, 0x00, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x00, 0x3f, 0xff, 0xe0, 0x03, 0xff, 0xff, 0xe0, 0x07, 0xfc, 0x00, 0x3f, 0xff,
0xc0, 0x03, 0xfe, 0x00, 0x00, 0x00, 0xff, 0xc0, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00,
0x0f, 0xff, 0xff, 0xff, 0xff, 0xdf, 0xe0, 0x03, 0xfe, 0x00, 0x1f, 0xff, 0x80, 0x07, 0xfe, 0x00,
0x00, 0x01, 0xff, 0x80, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x01, 0xff, 0xff, 0xff,
0xff, 0xcf, 0xf0, 0x03, 0xfe, 0x00, 0x1f, 0xff, 0x00, 0x07, 0xfc, 0x00, 0x00, 0x03, 0xff, 0x00,
0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x7f, 0xff, 0xff, 0xff, 0xc7, 0xf0, 0x03,
0xfe, 0x00, 0x1f, 0xff, 0x00, 0x0f, 0xf8, 0x00, 0x00, 0x0f, 0xfc, 0x00, 0x00, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0xff, 0xff, 0xff, 0x83, 0xf0, 0x03, 0xfe, 0x00, 0x0f, 0xfe,
0x00, 0x1f, 0xf0, 0x00, 0x00, 0x1f, 0xf8, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00,
0x00, 0x01, 0xff, 0xff, 0xff, 0x81, 0xf8, 0x01, 0xfe, 0x00, 0x0f, 0xfc, 0x00, 0x1f, 0xe0, 0x00,
0x00, 0x3f, 0xe0, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3f, 0xff,
0xff, 0x00, 0xf8, 0x01, 0xfe, 0x00, 0x0f, 0xf8, 0x00, 0x3f, 0xe0, 0x00, 0x00, 0x7f, 0xc0, 0x00,
0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0xff, 0xfc, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x07, 0xc0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0f,
};

struct GMT_CIRCLE {	/* Helper variables needed to draw great or small circle heads */
	double lon[2], lat[2];	/* Coordinates of arc end points */
	double A[3], B[3];	/* Cartesian vector of arc end points */
	double P[3];		/* Cartesian vector of the pole */
	bool longway;		/* True if the arc > 180 degrees */
	double r0;		/* Arc length in degrees */
	double r;		/* Will be 180 less if longway is true, otherwise r == r0 */
	double colat;		/* Colatitude of circle relative to pole */
	double rot;		/* Full opening angle of vector arc */
};

/* Local functions */

/* Converting LaTeX strings to EPS files.  This is based on the discussion we had at
 * https://github.com/GenericMappingTools/gmt/issues/4563#issuecomment-743374160.
 * As long as the user has latex, dvips and required fonts installed this should work
 * for everybody. P. Wessel, Dec 11, 2020.
 */

bool gmt_text_is_latex (struct GMT_CTRL *GMT, const char *string) {
	/* Detect if string contains LaTeX commands, i.e., "....@[LaTeX...@[ ..." or  "....<math>LaTeX...</math> ..." */
	char *p;
	gmt_M_unused (GMT);
	if (string == NULL || string[0] == '\0') return false;
	if ((p = strstr (string, "@[")) && strstr (&p[1], "@[")) return true;
	if ((p = strstr (string, "<math>")) && strstr (&p[1], "</math>")) return true;
	return false;
}

GMT_LOCAL bool gmtplot_has_title_breaks (struct GMT_CTRL *GMT, const char *string) {
	/* Returns true if string has line-break escape sequences in it */
	gmt_M_unused (GMT);
	if (string == NULL || string[0] == '\0') return false;
	return ((strstr (string, "<break>") || strstr (string, "@^")) ? true : false);
}

GMT_LOCAL unsigned char * gmtplot_latex_eps (struct GMT_CTRL *GMT, struct GMT_FONT *F, const char *string, struct imageinfo *h) {
	/* Convert a string containing LaTeX syntax to an EPS image */
	unsigned int i, o;
	int error = 0;
	char *text = NULL, *font, *code;
	char tmpdir[PATH_MAX] = {""}, here[PATH_MAX] = {""}, file[PATH_MAX] = {""}, cmd[PATH_MAX] = {""};
	unsigned char *picture = NULL;
	FILE *fp = NULL;
	struct GMTAPI_CTRL *API = GMT->parent;

	if (gmtplot_has_title_breaks (GMT, string)) {
		GMT_Report (API, GMT_MSG_ERROR, "LaTeX expressions are only allowed in single-line strings\n");
		return NULL;
	}

	if (gmt_check_executable (GMT, "latex", "--version", NULL, NULL)) {
		GMT_Report (API, GMT_MSG_DEBUG, "latex found.\n");
	}
	else {
		GMT_Report (API, GMT_MSG_ERROR, "latex is not installed or not in your executable path - cannot process LaTeX to DVI.\n");
		return NULL;
	}
	if (gmt_check_executable (GMT, "dvips", "--version", NULL, NULL)) {
		GMT_Report (API, GMT_MSG_DEBUG, "dvips found.\n");
	}
	else {
		GMT_Report (API, GMT_MSG_ERROR, "dvips is not installed or not in your executable path - cannot convert DVI to EPS.\n");
		return NULL;
	}

	/* Create unique directory for outputs, stored in tmpdir */

	if (gmt_get_tempname (API, "gmt_latex", NULL, tmpdir))
		return NULL;
	if (gmt_mkdir (tmpdir)) {
		GMT_Report (API, GMT_MSG_ERROR, "Unable to create directory %s - exiting.\n", tmpdir);
		return NULL;
	}
	/* Remember where we are */
	if (getcwd (here, PATH_MAX) == NULL) {
		GMT_Report (API, GMT_MSG_ERROR, "Unable to determine current working directory - exiting.\n");
		return NULL;
	}
	gmt_replace_backslash_in_path (here);
	/* Use tmpdir as the current directory */
	if (chdir (tmpdir)) {
		GMT_Report (API, GMT_MSG_ERROR, "Unable to change directory to %s - exiting.\n", tmpdir);
		return NULL;
	}
	/* Create LaTeX file */
	sprintf (file, "gmt_eq.tex");
	if ((fp = fopen (file, "w")) == NULL) {
		GMT_Report (API, GMT_MSG_ERROR, "gmtplot_latex_eps: Could not create LaTeX file %s.\n", file);
		return NULL;
	}

	/* Replace any @[ or <math> ... </math> with $ */
	text = strdup (string);
	for (i = o = 0; i < strlen (string); i++) {
		if (string[i] == '@' && string[i+1] == '[')
			text[o++] = '$', i++;
		else if (!strncmp (&string[i], "<math>", 6U))
			text[o++] = '$', i += 5;
		else if (!strncmp (&string[i], "</math>", 7U))
			text[o++] = '$', i += 6;
		else
			text[o++] = string[i];
	}
	text[o] = '\0'; /* Terminate the now shortened string */

	/* Check title font selection and pick corresponding fontpackagename and fontcode, if possible */
	switch (F->id) {
		case 0:  case 1:  case 2:   case 3: font = "helvet";    code = "phv";	break;
		case 4:  case 5:  case 6:   case 7: font = "mathptmx";  code = "ptm";	break;
		case 8:  case 9:  case 10: case 11: font = "courier";   code = "pcr";	break;
		case 12: font = "symbol";   code = "psy";   break;
		case 13: case 14: case 15: case 16: font = "avantgar";  code = "pag";	break;
		case 17: case 18: case 19: case 20: font = "bookman";   code = "pbk";	break;
		case 21: case 22: case 23: case 24: font = "helvet";    code = "phv";	break;
		case 25: case 26: case 27: case 28: font = "newcent";   code = "pnc";	break;
		case 29: case 30: case 31: case 32: font = "mathpazo";  code = "ppl";	break;
		case 33: font = "zapfchan"; code = "pzc";	break;
		case 34: font = "zapfding"; code = "pzd";	break;
		default: font = code = NULL;    /* Go with default */
	}
	/* Write LaTeX file content */
	fprintf (fp, "\\documentclass{article}\n"); /* Default to 10p font size */
	if (font) { /* Impose a selected font family, otherwise take default Computer Modern */
		GMT_Report (API, GMT_MSG_DEBUG, "gmtplot_latex_eps: Selecting font %s [%s].\n", font, code);
		fprintf (fp, "\\usepackage[T1]{fontenc}\\usepackage[utf8]{inputenc}\\usepackage{%s}\n", font);
	}
	fprintf (fp, "\\begin{document}\n\\thispagestyle{empty}\n");    /* No page number */
	if (code) /* Select font */
		fprintf (fp, "\\fontfamily{%s}\\selectfont\n", code);
	fprintf (fp, "%s\n\\end{document}\n", text);
	fclose (fp);
	gmt_M_str_free (text);

	/* Make script file for running latex and dvips */
#ifdef _WIN32
	sprintf (file, "gmt_eq.bat");
	sprintf (cmd, "start /B gmt_eq.bat");
#else
	sprintf (file, "gmt_eq.sh");
	sprintf (cmd, "sh gmt_eq.sh");
#endif
	if ((fp = fopen (file, "w")) == NULL) {
		GMT_Report (API, GMT_MSG_ERROR, "gmtplot_latex_eps: Could not create script file %s.\n", file);
		return NULL;
	}
#ifdef _WIN32
	fprintf (fp, "latex -interaction=nonstopmode gmt_eq.tex > NUL\ndvips -q -E gmt_eq.dvi -o equation.eps\n");
#else
	fprintf (fp, "latex -interaction=nonstopmode gmt_eq.tex > /dev/null\ndvips -q -E gmt_eq.dvi -o equation.eps\n");
#endif
	fclose (fp);

	/* Run the script via a system call */
	if ((error = system (cmd))) {
		GMT_Report (API, GMT_MSG_ERROR, "gmtplot_latex_eps: Running \"%s\" returned error %d.\n", cmd, error);
		GMT_Report (API, GMT_MSG_ERROR, "gmtplot_latex_eps: Please run it manually to learn what LaTeX packages you are missing.\n", cmd, error);
		GMT_Report (API, GMT_MSG_ERROR, "gmtplot_latex_eps: The script and logs can be found here: %s\n", tmpdir);
		return NULL;
	}
	else {  /* Success, now remove the temp files but not worry about the return code here */
#ifdef _WIN32
		system ("del gmt_eq.*");
#else
		system ("rm -f gmt_eq.*");
#endif
	}
	/* Retrieve the EPS code */
	gmt_M_memset (h, 1U, struct imageinfo); /* Initialize information struct */
	if (PSL_loadeps (GMT->PSL, "equation.eps", h, &picture)) {
		GMT_Report (API, GMT_MSG_ERROR, "gmtplot_latex_eps: Unable to load EPS file equation.eps\n");
		return NULL;
	}
	/* Clean up */
	if (gmt_remove_dir (API, tmpdir, false)) {
		GMT_Report (API, GMT_MSG_ERROR, "gmtplot_latex_eps: Unable to remove temporary directory %s!\n", tmpdir);
		PSL_free (picture);
		return NULL;
	}
	/* Change back to original directory */
	if (chdir (here)) {
		GMT_Report (API, GMT_MSG_ERROR, "gmtplot_latex_eps: Unable to change directory back to %s - exiting.\n", here);
		PSL_free (picture);
		return NULL;
	}

	/* Return the EPS data and size information via header */
	return picture;
}

/*	GMT_LINEAR PROJECTION MAP BOUNDARY	*/

GMT_LOCAL void gmtplot_linear_map_boundary (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	unsigned int cap = PSL->internal.line_cap;
	double x_length, y_length;

	x_length = GMT->current.proj.rect[XHI] - GMT->current.proj.rect[XLO];
	y_length = GMT->current.proj.rect[YHI] - GMT->current.proj.rect[YLO];

	PSL_command (PSL, "/PSL_slant_y 0 def /PSL_slant_x 0 def\n");	/* Unless x-annotations are slanted there is no adjustment. PSL_slant_y may be revised in gmt_xy_axis */

	if (GMT->current.map.frame.draw) {
		/* Temporarily change to square cap so rectangular frames have neat corners */
		PSL_setlinecap (PSL, PSL_SQUARE_CAP);

		if (GMT->current.map.frame.side[W_SIDE]) gmt_xy_axis (GMT, GMT->current.proj.rect[XLO], GMT->current.proj.rect[YLO], y_length, s, n,
			&GMT->current.map.frame.axis[GMT_Y], true, GMT->current.map.frame.side[W_SIDE]);	/* West or left y-axis */
		if (GMT->current.map.frame.side[E_SIDE]) gmt_xy_axis (GMT, GMT->current.proj.rect[XHI], GMT->current.proj.rect[YLO], y_length, s, n,
			&GMT->current.map.frame.axis[GMT_Y], false, GMT->current.map.frame.side[E_SIDE]);	/* East or right y-axis */
		if (GMT->current.map.frame.side[S_SIDE]) gmt_xy_axis (GMT, GMT->current.proj.rect[XLO], GMT->current.proj.rect[YLO], x_length, w, e,
			&GMT->current.map.frame.axis[GMT_X], true, GMT->current.map.frame.side[S_SIDE]);	/* South or lower x-axis */
		if (GMT->current.map.frame.side[N_SIDE]) gmt_xy_axis (GMT, GMT->current.proj.rect[XLO], GMT->current.proj.rect[YHI], x_length, w, e,
			&GMT->current.map.frame.axis[GMT_X], false, GMT->current.map.frame.side[N_SIDE]);	/* North or upper x-axis */

		PSL_setlinecap (PSL, cap);	/* Reset back to default */
	}
	if (!GMT->current.map.frame.header[0] || GMT->current.map.frame.plotted_header) return;	/* No title (and optional subtitle) today */

	PSL_comment (PSL, "Placing plot title\n");

	if (!GMT->current.map.frame.draw || GMT->current.map.frame.side[N_SIDE] <= GMT_AXIS_DRAW || GMT->current.setting.map_frame_type == GMT_IS_INSIDE)
		PSL_defunits (PSL, "PSL_H_y", GMT->current.setting.map_title_offset);	/* No ticks or annotations, offset by map_title_offset only */
	else
		PSL_command (PSL, "/PSL_H_y PSL_L_y PSL_LH add %d add def\n", PSL_IZ (PSL, GMT->current.setting.map_title_offset));	/* For title adjustment */

	PSL_command (PSL, "%d %d PSL_H_y add PSL_slant_y add M\n", PSL_IZ (PSL, 0.5 * x_length), PSL_IZ (PSL, y_length));

	gmt_map_title (GMT, 0.0, 0.0);
}

GMT_LOCAL unsigned int gmtplot_get_primary_annot (struct GMT_PLOT_AXIS *A) {
	/* Return the primary annotation item number [== GMT_ANNOT_UPPER if there are no unit set]*/

	unsigned int i, no[2] = {GMT_ANNOT_UPPER, GMT_ANNOT_LOWER};
	double val[2], s;

	for (i = 0; i < 2; i++) {
		val[i] = 0.0;
		if (!A->item[no[i]].active || A->item[no[i]].type == 'i' || A->item[no[i]].type == 'I') continue;
		switch (A->item[no[i]].unit) {
			case 'Y': case 'y':
				s = GMT_DAY2SEC_F * 365.25;
				break;
			case 'O': case 'o':
				s = GMT_DAY2SEC_F * 30.5;
				break;
			case 'U': case 'u':
				s = GMT_DAY2SEC_F * 7.0;
				break;
			case 'K': case 'k':
			case 'D': case 'd':
				s = GMT_DAY2SEC_F;
				break;
			case 'H': case 'h':
				s = GMT_HR2SEC_F;
				break;
			case 'M': case 'm':
				s = GMT_MIN2SEC_F;
				break;
			case 'C': case 'c':
				s = 1.0;
				break;
			default:	/* No unit specified - probably not a time axis */
				s = 1.0;
				break;
		}
		val[i] = A->item[no[i]].interval * s;
	}
	if (A->item[GMT_ANNOT_UPPER].special) return GMT_ANNOT_UPPER;
	return ((val[0] > val[1]) ? GMT_ANNOT_UPPER : GMT_ANNOT_LOWER);
}

GMT_LOCAL bool gmtplot_skip_second_annot (unsigned int item, double x, double x2[], unsigned int n, unsigned int primary) {
	unsigned int i;
	bool found;
	double small;

	if (n < 2) return (false);	/* Need at least two points so no need to skip */
	if (item == primary) return (false);		/* Not working on secondary annotation */
	if (!x2) return (false);			/* None given */

	small = (x2[1] - x2[0]) * GMT_CONV4_LIMIT;
	for (i = 0, found = false; !found && i < n; i++)
		found = (fabs (x2[i] - x) < small);
	return (found);
}

GMT_LOCAL void gmtplot_map_latline (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double lat, double west, double east)		/* Draws a line of constant latitude */ {
	uint64_t nn;
	double *llon = NULL, *llat = NULL;
#ifdef DEBUG
	uint64_t k;
	FILE *fp = NULL;
	char name[GMT_LEN64] = {""};
	if (GMT->hidden.gridline_debug && (GMT->hidden.gridline_kind == 'x' || GMT->hidden.gridline_val != lat)) return;
#endif
	nn = gmtlib_latpath (GMT, lat, west, east, &llon, &llat);
#ifdef DEBUG
	if (GMT->hidden.gridline_debug) {
		snprintf (name, GMT_LEN64, "gridline_y_%g_ll.txt", lat);
		fp = fopen (name, "w");
		for (k = 0; k < nn; k++) fprintf (fp, "%g\t%g\n", llon[k], llat[k]);
		fclose (fp);
	}
#endif
	GMT->current.plot.n = gmt_geo_to_xy_line (GMT, llon, llat, nn);
#ifdef DEBUG
	if (GMT->hidden.gridline_debug) {
		snprintf (name, GMT_LEN64, "gridline_y_%g_xy.txt", lat);
		fp = fopen (name, "w");
		for (k = 0; k < GMT->current.plot.n; k++) fprintf (fp, "%g\t%g\t%d\n", GMT->current.plot.x[k], GMT->current.plot.y[k], GMT->current.plot.pen[k]);
		fclose (fp);
	}
#endif

	if (GMT->current.plot.n > 1) {	/* Need at least 2 points for a line */
		PSL_comment (PSL, "Lat = %g\n", lat);
		if (GMT->current.map.parallel_straight) {	/* Simplify to a 2-point straight line */
			GMT->current.plot.x[1] = GMT->current.plot.x[GMT->current.plot.n-1];
			GMT->current.plot.y[1] = GMT->current.plot.y[GMT->current.plot.n-1];
			GMT->current.plot.pen[1] = GMT->current.plot.pen[GMT->current.plot.n-1];
			GMT->current.plot.n = 2;
		}
		gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);
	}
	gmt_M_free (GMT, llon);
	gmt_M_free (GMT, llat);
}

GMT_LOCAL void gmtplot_map_lonline (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double lon, double south, double north)	/* Draws a line of constant longitude */ {
	uint64_t nn;
	double *llon = NULL, *llat = NULL;
#ifdef DEBUG
	uint64_t k;
	FILE *fp = NULL;
	char name[GMT_LEN64] = {""};
	if (GMT->hidden.gridline_debug && (GMT->hidden.gridline_kind == 'y' || GMT->hidden.gridline_val != lon)) return;
#endif
	nn = gmtlib_lonpath (GMT, lon, south, north, &llon, &llat);
#ifdef DEBUG
	if (GMT->hidden.gridline_debug) {
		snprintf (name, GMT_LEN64, "gridline_x_%g_ll.txt", lon);
		fp = fopen (name, "w");
		for (k = 0; k < nn; k++) fprintf (fp, "%g\t%g\n", llon[k], llat[k]);
		fclose (fp);
	}
#endif
	GMT->current.plot.n = gmt_geo_to_xy_line (GMT, llon, llat, nn);
#ifdef DEBUG
	if (GMT->hidden.gridline_debug) {
		snprintf (name, GMT_LEN64, "gridline_x_%g_xy.txt", lon);
		fp = fopen (name, "w");
		for (k = 0; k < GMT->current.plot.n; k++) fprintf (fp, "%g\t%g\t%d\n", GMT->current.plot.x[k], GMT->current.plot.y[k], GMT->current.plot.pen[k]);
		fclose (fp);
	}
#endif

	if (GMT->current.plot.n > 1) {	/* Need at least 2 points for a line */
		PSL_comment (PSL, "Lon = %g\n", lon);
		if (GMT->current.map.meridian_straight) {	/* Simplify to a 2-point straight line */
			GMT->current.plot.x[1] = GMT->current.plot.x[GMT->current.plot.n-1];
			GMT->current.plot.y[1] = GMT->current.plot.y[GMT->current.plot.n-1];
			GMT->current.plot.pen[1] = GMT->current.plot.pen[GMT->current.plot.n-1];
			GMT->current.plot.n = 2;
		}
		gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);
	}
	gmt_M_free (GMT, llon);
	gmt_M_free (GMT, llat);
}

GMT_LOCAL void gmtplot_linearx_oblgrid (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, double dval) {
	/* x gridlines in oblique coordinates for all but the Oblique Mercator projection [which already is oblique] */
	double *x = NULL, *lon = NULL, *lat = NULL, *lat_obl = NULL, tval, p_cap, s_cap;
	unsigned int idup = 0, i, j, k, nx, np, nc1 = 0, nc2, npc, np1;
	bool cap = false;
	gmt_M_unused(PSL); gmt_M_unused(w); gmt_M_unused(e);

	/* Ideally we should determine the w/e/s/n of the oblique coordinates but here we will simply
	 * create oblique coordinates for the full 0/360/-90/90, convert to regular coordinates and
	 * then truncate points outside the actual w/e/s/n */

	/* Do we have duplicate e and w boundaries ? */
	p_cap = fabs (GMT->current.setting.map_polar_cap[0]);
	s_cap = -p_cap;
	idup = (gmt_M_is_azimuthal(GMT)) ? 1 : 0;
	cap = !doubleAlmostEqual (p_cap, 90.0);	/* true if we have a polar cap specified */
	tval = (n - s) * GMT->current.setting.map_line_step / GMT->current.map.height;

	nx = gmtlib_linear_array (GMT, 0.0, TWO_PI, D2R * dval, D2R * GMT->current.map.frame.axis[GMT_X].phase, &x);
	np = gmtlib_linear_array (GMT, -90.0, 90.0, tval, 0.0, &lat_obl);
	np1 = nc2 = np - 1;	/* Nominal number of points in path */
	lon = gmt_M_memory (GMT, NULL, np+2, double);	/* Allow 2 more slots for possibly inserted cap-latitudes */
	lat = gmt_M_memory (GMT, NULL, np+2, double);
	for (i = 0; i < nx - idup; i++) {	/* For each oblique meridian to draw */
		/* Create lon,lat arrays of oblique coordinates for this meridian */
		for (k = j = 0; k < np; k++, j++) {
			gmtlib_iobl (GMT, &lon[j], &lat[j], x[i], D2R * lat_obl[k]);	/* Get regular coordinates of this point */
			lon[j] *= R2D;	lat[j] *= R2D;	/* Convert back to degrees */
			if (lat_obl[k] < s_cap && k < np1 && lat_obl[k+1] > s_cap)	{	/* Must insert S pole cap latitude point */
				j++;	gmtlib_iobl (GMT, &lon[j], &lat[j], x[i], D2R * s_cap);
				lon[j] *= R2D;	lat[j] *= R2D;	/* Back to degrees */
				nc1 = j;
			}
			else if (lat_obl[k] < p_cap && k < np1 && lat_obl[k+1] > p_cap) {	/* Must insert N pole cap latitude point */
				j++; gmtlib_iobl (GMT, &lon[j], &lat[j], x[i], D2R * p_cap);
				lon[j] *= R2D;	lat[j] *= R2D;	/* Back to degrees */
				nc2 = j;
			}
		}
		if (cap) {	/* Only plot the line between the two polar caps */
			npc = nc2 - nc1 + 1;	/* Number of points along meridian bounded by caps */
			if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, &lon[nc1], &lat[nc1], npc)) == 0) continue;
		}
		else {		/* No polar cap in effect, plot entire meridian */
			if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, lon, lat, j)) == 0) continue;
		}
		gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);
	}
	if (nx) gmt_M_free (GMT, x);
	if (cap) {	/* Draw the polar cap(s) meridians with a separate lon spacing */
		nx = gmtlib_linear_array (GMT, 0.0, TWO_PI, D2R * GMT->current.setting.map_polar_cap[1], D2R * GMT->current.map.frame.axis[GMT_X].phase, &x);
		for (i = 0; i < nx - idup; i++) {
			for (k = j = 0; k < np; k++, j++) {
				gmtlib_iobl (GMT, &lon[j], &lat[j], x[i], D2R * lat_obl[k]);	/* Get regular coordinates of this point */
				lon[j] *= R2D;	lat[j] *= R2D;	/* Back to degrees */
				if (lat_obl[k] < s_cap && k < np1 && lat_obl[k+1] > s_cap)	{	/* Must insert S pole cap latitude point */
					j++;	gmtlib_iobl (GMT, &lon[j], &lat[j], x[i], D2R * s_cap);
					lon[j] *= R2D;	lat[j] *= R2D;	/* Back to degrees */
					nc1 = j;
				}
				else if (lat_obl[k] < p_cap && k < np1 && lat_obl[k+1] > p_cap) {	/* Must insert N pole cap latitude point */
					j++; gmtlib_iobl (GMT, &lon[j], &lat[j], x[i], D2R * p_cap);
					lon[j] *= R2D;	lat[j] *= R2D;	/* Back to degrees */
					nc2 = j;
				}
			}
			if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, lon, lat, nc1+1)) > 0)
				gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);
			if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, &lon[nc2], &lat[nc2], j-nc2)) > 0)
				gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);
		}
		if (nx) gmt_M_free (GMT, x);
	}
	gmt_M_free (GMT, lat_obl);
	gmt_M_free (GMT, lon);
	gmt_M_free (GMT, lat);
}

GMT_LOCAL void gmtplot_lineary_grid (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, double dval) {
	double *y = NULL;
	char *type = (gmt_M_y_is_lat (GMT, GMT_IN)) ? "parallel" : "y";
	unsigned int i, ny = 0;

	if (GMT->current.proj.z_down) {
		if (GMT->current.proj.z_down == GMT_ZDOWN_Z) /* z = n - r */
			ny = gmtlib_linear_array (GMT, 0.0, GMT->current.proj.z_radius-s, dval, GMT->current.map.frame.axis[GMT_Y].phase, &y);
		else if (GMT->current.proj.z_down == GMT_ZDOWN_ZP) /* z = n - r */
			ny = gmtlib_linear_array (GMT, GMT->current.proj.z_radius-n, GMT->current.proj.z_radius-s, dval, GMT->current.map.frame.axis[GMT_Y].phase, &y);
		for (i = 0; i < ny; i++) {
			if (GMT->current.proj.z_down == GMT_ZDOWN_ZP)
				y[i] = GMT->current.proj.z_radius - y[i];	/* These are the z values needed for positioning */
			else
				y[i] = GMT->common.R.wesn[YHI] - y[i];	/* These are the radial values needed for positioning */
		}
	}
	else
		ny = gmtlib_linear_array (GMT, s, n, dval, GMT->current.map.frame.axis[GMT_Y].phase, &y);
	for (i = 0; i < ny; i++) {
		GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Draw %s = %g from %g to %g\n", type, y[i], w, e);
		gmtplot_map_latline (GMT, PSL, y[i], w, e);
	}
	if (ny) gmt_M_free (GMT, y);

}

GMT_LOCAL void gmtplot_x_grid (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double s, double n, double *x, unsigned int nx) {
	unsigned int i;
	double x1, y1, x2, y2;

	for (i = 0; i < nx; i++) {
		if (gmt_M_x_is_lon (GMT, GMT_IN))
			gmtplot_map_lonline (GMT, PSL, x[i], s, n);
		else {
			gmt_geo_to_xy (GMT, x[i], s, &x1, &y1);
			gmt_geo_to_xy (GMT, x[i], n, &x2, &y2);
			PSL_plotsegment (PSL, x1, y1, x2, y2);
		}
	}
}

GMT_LOCAL void gmtplot_lineary_oblgrid (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, double dval) {
	/* y gridlines in oblique coordinates for all but the Oblique Mercator projection [which already is oblique] */

	double *y = NULL, *lon = NULL, *lon_obl = NULL, *lat = NULL, tval, p_cap;
	bool cap;
	unsigned int i, k, ny, np;
	gmt_M_unused(PSL); gmt_M_unused(s); gmt_M_unused(n);

	/* Ideally we should determine the w/e/s/n of the oblique coordinates but here we will simply
	 * create oblique coordinates for the full 0/360/-90/90, convert to regular coordinates and
	 * then truncate points outside the actual w/e/s/n */

	ny = gmtlib_linear_array (GMT, -M_PI_2, M_PI_2, D2R * dval, D2R * GMT->current.map.frame.axis[GMT_Y].phase, &y);
	tval = (e - w) * GMT->current.setting.map_line_step / GMT->current.map.width;
	np = gmtlib_linear_array (GMT, 0.0, TWO_PI, D2R * tval, 0.0, &lon_obl);
	lon = gmt_M_memory (GMT, NULL, np+2, double);	/* Allow 2 more slots for possibly inserted cap-latitudes */
	lat = gmt_M_memory (GMT, NULL, np+2, double);
	for (i = 0; i < ny; i++) {
		for (k = 0; k < np; k++) {
			gmtlib_iobl (GMT, &lon[k], &lat[k], lon_obl[k], y[i]);	/* Get regular coordinates of this point */
			lon[k] *= R2D;	lat[k] *= R2D;	/* Convert to degrees */
		}
		if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, lon, lat, np)) == 0) continue;
		gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);
	}
	if (ny) gmt_M_free (GMT, y);
	p_cap = fabs (GMT->current.setting.map_polar_cap[0]);
	cap = !doubleAlmostEqual (p_cap, 90.0);	/* true if we have a polar cap specified */
	if (cap) {	/* Draw the polar cap(s) with a separate spacing */
		p_cap = D2R * GMT->current.setting.map_polar_cap[0];
		for (k = 0; k < np; k++) {	/* S polar cap */
			gmtlib_iobl (GMT, &lon[k], &lat[k], lon_obl[k], -p_cap);	/* Get regular coordinates of this point */
			lon[k] *= R2D;	lat[k] *= R2D;	/* Convert to degrees */
		}
		if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, lon, lat, np)) > 0)
			gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);
		for (k = 0; k < np; k++) {	/* N polar cap */
			gmtlib_iobl (GMT, &lon[k], &lat[k], lon_obl[k], p_cap);	/* Get regular coordinates of this point */
			lon[k] *= R2D;	lat[k] *= R2D;	/* Convert to degrees */
		}
		if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, lon, lat, np)) > 0)
			gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);
	}
	gmt_M_free (GMT, lon_obl);
	gmt_M_free (GMT, lon);
	gmt_M_free (GMT, lat);
}

GMT_LOCAL void gmtplot_y_grid (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double *y, unsigned int ny) {
	unsigned int i;
	double x1, y1, x2, y2;

	for (i = 0; i < ny; i++) {
		if (gmt_M_y_is_lat (GMT, GMT_IN))
			gmtplot_map_latline (GMT, PSL, y[i], w, e);
		else {
			gmt_geo_to_xy (GMT, w, y[i], &x1, &y1);
			gmt_geo_to_xy (GMT, e, y[i], &x2, &y2);
			PSL_plotsegment (PSL, x1, y1, x2, y2);
		}
	}
}

void plot_timex_grid (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, unsigned int item) {
	unsigned int nx;
	double *x = NULL;

	nx = gmtlib_time_array (GMT, w, e, &GMT->current.map.frame.axis[GMT_X].item[item], &x);
	gmtplot_x_grid (GMT, PSL, s, n, x, nx);
	if (x) gmt_M_free (GMT, x);
}

GMT_LOCAL void gmtplot_timey_grid (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, unsigned int item) {
	unsigned int ny;
	double *y = NULL;

	ny = gmtlib_time_array (GMT, s, n, &GMT->current.map.frame.axis[GMT_Y].item[item], &y);
	gmtplot_y_grid (GMT, PSL, w, e, y, ny);
	if (y) gmt_M_free (GMT, y);
}

GMT_LOCAL void gmtplot_logx_grid (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, double dval) {
	unsigned int nx;
	double *x = NULL;

	nx = gmtlib_log_array (GMT, w, e, dval, &x);
	gmtplot_x_grid (GMT, PSL, s, n, x, nx);
	if (x) gmt_M_free (GMT, x);
}

GMT_LOCAL void gmtplot_logy_grid (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, double dval) {
	unsigned int ny;
	double *y = NULL;

	ny = gmtlib_log_array (GMT, s, n, dval, &y);
	gmtplot_y_grid (GMT, PSL, w, e, y, ny);
	if (y) gmt_M_free (GMT, y);
}

GMT_LOCAL void gmtplot_powx_grid (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, double dval) {
	unsigned int nx;
	double *x = NULL;

	nx = gmtlib_pow_array (GMT, w, e, dval, 0, &x);
	gmtplot_x_grid (GMT, PSL, s, n, x, nx);
	if (x) gmt_M_free (GMT, x);
}

GMT_LOCAL void gmtplot_powy_grid (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, double dval) {
	unsigned int ny;
	double *y = NULL;

	ny = gmtlib_pow_array (GMT, s, n, dval, 1, &y);
	gmtplot_y_grid (GMT, PSL, w, e, y, ny);
	if (y) gmt_M_free (GMT, y);
}

GMT_LOCAL double gmtplot_shift_gridline (struct GMT_CTRL *GMT, double val, unsigned int type) {
	/* Only for oblique projections: If any of the corners are exactly multiples of annotation
	 * or tick intervals then the gridline intersection may fail (tangent or slightly outside
	 * due to round-off).  We determine which gridlines go through the corners and shift them
	 * a tiny bit to the inside to ensure crossings */
	double shift = 0.0;
	if (!GMT->common.R.oblique) return shift;	/* Return zero if not an oblique projection */

	if (type == GMT_X) {
		if (gmt_M_360_range (val, GMT->common.R.wesn_orig[XLO])) val = GMT->common.R.wesn_orig[XLO];
		else if (gmt_M_360_range (val, GMT->common.R.wesn_orig[XHI])) val = GMT->common.R.wesn_orig[XHI];
		if (doubleAlmostEqualZero (val, GMT->common.R.wesn_orig[XLO])) shift = +GMT_CONV4_LIMIT * fabs (GMT->common.R.wesn_orig[XHI] - GMT->common.R.wesn_orig[XLO]);	/* Add this to lon to get a slightly larger longitude to ensure crossing */
		else if (doubleAlmostEqualZero (val, GMT->common.R.wesn_orig[XHI])) shift = -GMT_CONV4_LIMIT * fabs (GMT->common.R.wesn_orig[XHI] - GMT->common.R.wesn_orig[XLO]);	/* Add this to lon to get a slightly smaller longitude to ensure crossing */
	}
	else {
		if (doubleAlmostEqualZero (val, GMT->common.R.wesn_orig[YLO])) shift = +GMT_CONV4_LIMIT * fabs (GMT->common.R.wesn_orig[YHI] - GMT->common.R.wesn_orig[YLO]);	/* Add this to lon to get a slightly larger longitude to ensure crossing */
		else if (doubleAlmostEqualZero (val, GMT->common.R.wesn_orig[YHI])) shift = -GMT_CONV4_LIMIT * fabs (GMT->common.R.wesn_orig[YHI] - GMT->common.R.wesn_orig[YLO]);	/* Add this to lon to get a slightly smaller longitude to ensure crossing */

	}
	if (shift != 0.0) GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "Adjusted argument %g by %g\n", val, shift);
	return shift;
}

/*	FANCY RECTANGULAR PROJECTION MAP BOUNDARY	*/

GMT_LOCAL void gmtplot_fancy_frame_offset (struct GMT_CTRL *GMT, double angle, double shift[2]) {
	/* Given the angle of the axis, return the coordinate adjustments needed to
	 * shift in order to plot the outer 1-2 parallel frame lines (shift[0|1] */

	double s, c;
	sincos (angle, &s, &c);
	shift[0] =  GMT->current.setting.map_frame_width * s;
	shift[1] = -GMT->current.setting.map_frame_width * c;
}

GMT_LOCAL void gmtplot_fancy_frame_straightlon_checkers (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, bool secondary_too) {
	/* Plot checkers along straight longitude boundaries */
	int i, k, nx;
	unsigned int shade, item[2] = {GMT_TICK_UPPER, GMT_TICK_LOWER};
	double dx, w1, val, v1, v2, x1, x2, y1, y2, shift_s[2], shift_n[2], scale[2];
	struct GMT_PLOT_AXIS_ITEM *T = NULL;

	scale[0] = (secondary_too) ? 0.5 : 1.0;
	scale[1] = 1.5;

	gmt_geo_to_xy (GMT, e, s, &x1, &y1);
	gmt_geo_to_xy (GMT, w, s, &x2, &y2);
	gmtplot_fancy_frame_offset (GMT, d_atan2 (y2 - y1, x2 - x1), shift_s);

	gmt_geo_to_xy (GMT, w, n, &x1, &y1);
	gmt_geo_to_xy (GMT, e, n, &x2, &y2);
	gmtplot_fancy_frame_offset (GMT, d_atan2 (y2 - y1, x2 - x1), shift_n);

	for (k = 0; k < 1 + secondary_too; k++) {
		T = &GMT->current.map.frame.axis[GMT_X].item[item[k]];
		if (T->active) {
			dx = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_X].type, T);
			shade = (urint (floor ((w - GMT->current.map.frame.axis[GMT_X].phase)/ dx))) % 2;
			w1 = floor ((w - GMT->current.map.frame.axis[GMT_X].phase)/ dx) * dx + GMT->current.map.frame.axis[GMT_X].phase;
			nx = (w1 > e) ? -1 : irint (((e - w1) / dx + GMT_CONV4_LIMIT));
			for (i = 0; i <= nx; i++) {
				shade = !shade;
				val = w1 + i * dx;
				v1 = MAX (val, w);
				v2 = MIN (val + dx, e);
				if (v2 - v1 < GMT_CONV8_LIMIT) continue;
				PSL_setcolor (PSL, shade ? GMT->current.setting.map_frame_pen.rgb : GMT->PSL->init.page_rgb, PSL_IS_STROKE);
				if (GMT->current.map.frame.side[S_SIDE] & GMT_AXIS_TICK) {
					gmt_geo_to_xy (GMT, v1, s, &x1, &y1);
					gmt_geo_to_xy (GMT, v2, s, &x2, &y2);
					PSL_plotsegment (PSL, x1-0.5*scale[k]*shift_s[0], y1-0.5*scale[k]*shift_s[1], x2-0.5*scale[k]*shift_s[0], y2-0.5*scale[k]*shift_s[1]);
				}
				if (GMT->current.map.frame.side[N_SIDE] & GMT_AXIS_TICK) {
					gmt_geo_to_xy (GMT, v1, n, &x1, &y1);
					gmt_geo_to_xy (GMT, v2, n, &x2, &y2);
					PSL_plotsegment (PSL, x1-0.5*scale[k]*shift_n[0], y1-0.5*scale[k]*shift_n[1], x2-0.5*scale[k]*shift_n[0], y2-0.5*scale[k]*shift_n[1]);
				}
			}
		}
	}
}

GMT_LOCAL void gmtplot_fancy_frame_curvedlon_checkers (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, bool secondary_too) {
	/* Plot checkers along curved longitude boundaries */
	int i, k, nx;
	unsigned int shade, item[2] = {GMT_TICK_UPPER, GMT_TICK_LOWER};
	double dx, w1, v1, v2, val, x1, x2, y1, y2, az1, az2, dr, scale[2], radius_s, radius_n;
	struct GMT_PLOT_AXIS_ITEM *T;

	scale[0] = (secondary_too) ? 0.5 : 1.0;
	scale[1] = 1.5;
	dr = 0.5 * GMT->current.setting.map_frame_width;
	gmt_geo_to_xy (GMT, w, s, &x1, &y1);
	gmt_geo_to_xy (GMT, w, n, &x2, &y2);
	radius_s = hypot (x1 - GMT->current.proj.c_x0, y1 - GMT->current.proj.c_y0);
	radius_n = hypot (x2 - GMT->current.proj.c_x0, y2 - GMT->current.proj.c_y0);

	for (k = 0; k < 1 + secondary_too; k++) {
		T = &GMT->current.map.frame.axis[GMT_X].item[item[k]];
		if (T->active) {
			dx = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_X].type, T);
			shade = urint (floor ((w - GMT->current.map.frame.axis[GMT_X].phase) / dx)) % 2;
			w1 = floor ((w - GMT->current.map.frame.axis[GMT_X].phase)/dx) * dx + GMT->current.map.frame.axis[GMT_X].phase;
			nx = (w1 > e) ? -1 : irint ((e-w1) / dx + GMT_CONV4_LIMIT);
			for (i = 0; i <= nx; i++) {
				shade = !shade;
				val = w1 + i * dx;
				v1 = MAX (val, w);
				v2 = MIN (val + dx, e);
				if (v2 - v1 < GMT_CONV8_LIMIT) continue;
				PSL_setcolor (PSL, shade ? GMT->current.setting.map_frame_pen.rgb : GMT->PSL->init.page_rgb, PSL_IS_STROKE);
				if (GMT->current.map.frame.side[S_SIDE] & GMT_AXIS_TICK) {
					gmt_geo_to_xy (GMT, v2, s, &x1, &y1);
					gmt_geo_to_xy (GMT, v1, s, &x2, &y2);
					az1 = d_atan2d (y1 - GMT->current.proj.c_y0, x1 - GMT->current.proj.c_x0);
					az2 = d_atan2d (y2 - GMT->current.proj.c_y0, x2 - GMT->current.proj.c_x0);
					if (GMT->current.proj.north_pole) {
						if (az1 < az2) az1 += 360.0;
						PSL_plotarc (PSL, GMT->current.proj.c_x0, GMT->current.proj.c_y0, radius_s+scale[k]*dr, az2, az1, PSL_MOVE|PSL_STROKE);
					}
					else {
						if (az2 < az1) az2 += 360.0;
						PSL_plotarc (PSL, GMT->current.proj.c_x0, GMT->current.proj.c_y0, radius_s-scale[k]*dr, az1, az2, PSL_MOVE|PSL_STROKE);
					}
				}
				if (GMT->current.map.frame.side[N_SIDE] & GMT_AXIS_TICK) {
					gmt_geo_to_xy (GMT, v2, n, &x1, &y1);
					gmt_geo_to_xy (GMT, v1, n, &x2, &y2);
					az1 = d_atan2d (y1 - GMT->current.proj.c_y0, x1 - GMT->current.proj.c_x0);
					az2 = d_atan2d (y2 - GMT->current.proj.c_y0, x2 - GMT->current.proj.c_x0);
					if (GMT->current.proj.north_pole) {
						if (az1 < az2) az1 += 360.0;
						PSL_plotarc (PSL, GMT->current.proj.c_x0, GMT->current.proj.c_y0, radius_n-scale[k]*dr, az2, az1, PSL_MOVE|PSL_STROKE);
					}
					else {
						if (az2 < az1) az2 += 360.0;
						PSL_plotarc (PSL, GMT->current.proj.c_x0, GMT->current.proj.c_y0, radius_n+scale[k]*dr, az1, az2, PSL_MOVE|PSL_STROKE);
					}
				}
			}
		}
	}
}

GMT_LOCAL void gmtplot_fancy_frame_straightlat_checkers (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, bool secondary_too) {
	/* Plot checkers along straight latitude boundaries */
	int i, k, ny;
	unsigned int shade, item[2] = {GMT_TICK_UPPER, GMT_TICK_LOWER};
	double dy, s1, val, v1, v2, x1, x2, y1, y2, shift_w[2], shift_e[2], scale[2];
	struct GMT_PLOT_AXIS_ITEM *T = NULL;

	scale[0] = (secondary_too) ? 0.5 : 1.0;
	scale[1] = 1.5;

	gmt_geo_to_xy (GMT, w, s, &x1, &y1);
	gmt_geo_to_xy (GMT, w, n, &x2, &y2);
	gmtplot_fancy_frame_offset (GMT, d_atan2 (y2 - y1, x2 - x1), shift_w);

	gmt_geo_to_xy (GMT, e, s, &x1, &y1);
	gmt_geo_to_xy (GMT, e, n, &x2, &y2);
	gmtplot_fancy_frame_offset (GMT, d_atan2 (y2 - y1, x2 - x1), shift_e);

	/* Tick S-N axes */

	for (k = 0; k < 1 + secondary_too; k++) {
		T = &GMT->current.map.frame.axis[GMT_Y].item[item[k]];
		if (T->active) {
			dy = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_Y].type, T);
			shade = urint (floor ((s - GMT->current.map.frame.axis[GMT_Y].phase) / dy)) % 2;
			s1 = floor((s - GMT->current.map.frame.axis[GMT_Y].phase)/dy) * dy + GMT->current.map.frame.axis[GMT_Y].phase;
			ny = (s1 > n) ? -1 : irint ((n-s1) / dy + GMT_CONV4_LIMIT);
			for (i = 0; i <= ny; i++) {
				shade = !shade;
				val = s1 + i * dy;
				v1 = MAX (val, s);
				v2 = MIN (val + dy, n);
				if (v2 - v1 < GMT_CONV8_LIMIT) continue;
				PSL_setcolor (PSL, shade ? GMT->current.setting.map_frame_pen.rgb : GMT->PSL->init.page_rgb, PSL_IS_STROKE);
				if (GMT->current.map.frame.side[W_SIDE] & GMT_AXIS_TICK) {
					gmt_geo_to_xy (GMT, w, v1, &x1, &y1);
					gmt_geo_to_xy (GMT, w, v2, &x2, &y2);
					PSL_plotsegment (PSL, x1-0.5*scale[k]*shift_w[0], y1-0.5*scale[k]*shift_w[1], x2-0.5*scale[k]*shift_w[0], y2-0.5*scale[k]*shift_w[1]);
				}
				if (GMT->current.map.frame.side[E_SIDE] & GMT_AXIS_TICK) {
					gmt_geo_to_xy (GMT, e, v1, &x1, &y1);
					gmt_geo_to_xy (GMT, e, v2, &x2, &y2);
					PSL_plotsegment (PSL, x1+0.5*scale[k]*shift_e[0], y1+0.5*scale[k]*shift_e[1], x2+0.5*scale[k]*shift_e[0], y2+0.5*scale[k]*shift_e[1]);
				}
			}
		}
	}
}

GMT_LOCAL void gmtplot_fancy_frame_straight_outline (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double lonA, double latA, double lonB, double latB, unsigned int side, bool secondary_too) {
	unsigned int k, kn = 1, axis;
	double scale = 1.0, x[2], y[2], angle, s, c, dx, dy, Ldx, Ldy, f0 = 1.0, f1 = 1.0;
	struct GMT_PLOT_AXIS_ITEM *T = NULL;

	if (!GMT->current.map.frame.side[side]) return;	/* Do not draw this frame side */
	switch (side) {	/* Determine if there are missing partner sides that invalidates extensions */
		case S_SIDE:
			if (!GMT->current.map.frame.side[W_SIDE]) f0 = 0;
			if (!GMT->current.map.frame.side[E_SIDE]) f1 = 0;
			break;
		case E_SIDE:
			if (!GMT->current.map.frame.side[S_SIDE]) f0 = 0;
			if (!GMT->current.map.frame.side[N_SIDE]) f1 = 0;
			break;
		case N_SIDE:
			if (!GMT->current.map.frame.side[E_SIDE]) f0 = 0;
			if (!GMT->current.map.frame.side[W_SIDE]) f1 = 0;
			break;
		case W_SIDE:
			if (!GMT->current.map.frame.side[N_SIDE]) f0 = 0;
			if (!GMT->current.map.frame.side[S_SIDE]) f1 = 0;
			break;
	}

	if (secondary_too) {
		scale = 0.5;
		++kn;
	}
	axis = side % 2;	/* Gives 0 for GMT_X and 1 for GMT_Y */
	T = &GMT->current.map.frame.axis[axis].item[GMT_TICK_UPPER];
	if (!T->active) return;

	gmt_geo_to_xy (GMT, lonA, latA, &x[0], &y[0]);
	gmt_geo_to_xy (GMT, lonB, latB, &x[1], &y[1]);
	angle = d_atan2 (y[1] - y[0], x[1] - x[0]);
	sincos (angle, &s, &c);
	/* Ldx/dy is the components of the extension of the fancy frame into neighboring side frames */
	Ldx = (GMT->current.setting.map_frame_type == GMT_IS_ROUNDED) ? 0.0 : GMT->current.setting.map_frame_width * c;
	Ldy = (GMT->current.setting.map_frame_type == GMT_IS_ROUNDED) ? 0.0 : GMT->current.setting.map_frame_width * s;
	/* dx,dy is the outward shift components to draw the outside (and possibly half-way) parallel frame outline */
	dx =  GMT->current.setting.map_frame_width * s;
	dy = -GMT->current.setting.map_frame_width * c;
	PSL_plotsegment (PSL, x[0]-f0*Ldx, y[0]-f0*Ldy, x[1]+f1*Ldx, y[1]+f1*Ldy);
	for (k = 0; k < kn; k++) {
		x[0] += scale*dx;
		y[0] += scale*dy;
		x[1] += scale*dx;
		y[1] += scale*dy;
		PSL_plotsegment (PSL, x[0]-f0*Ldx, y[0]-f0*Ldy, x[1]+f1*Ldx, y[1]+f1*Ldy);
	}
}

GMT_LOCAL void gmtplot_fancy_frame_curved_outline (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double lonA, double latA, double lonB, double latB, unsigned int side, bool secondary_too) {
	double scale[2] = {1.0, 1.0}, escl, x1, x2, y1, y2, radius, r_inc, az1, az2, da0, da, width, s, fw, fe;

	if (!GMT->current.map.frame.side[side]) return;	/* This side is inactive */
	fw = (GMT->current.map.frame.side[W_SIDE]) ? 1.0 : 0.0;	/* Only extend if W side is plotted */
	fe = (GMT->current.map.frame.side[E_SIDE]) ? 1.0 : 0.0;	/* Only extend if E side is plotted */
	if (secondary_too) scale[0] = scale[1] = 0.5;
	width = GMT->current.setting.map_frame_width;
	escl = (GMT->current.setting.map_frame_type == GMT_IS_ROUNDED) ? 0.0 : 1.0;	/* Want rounded corners */
	gmt_geo_to_xy (GMT, lonA, latA, &x1, &y1);
	gmt_geo_to_xy (GMT, lonB, latB, &x2, &y2);
	radius = hypot (x1 - GMT->current.proj.c_x0, y1 - GMT->current.proj.c_y0);
	s = ((GMT->current.proj.north_pole && side == 2) || (!GMT->current.proj.north_pole && side == 0)) ? -1.0 : +1.0;	/* North: needs shorter radius.  South: Needs longer radius (opposite in S hemi) */
	r_inc = s*scale[0] * width;
	if (gmt_M_is_azimuthal(GMT) && gmt_M_360_range (lonA, lonB)) {	/* Full 360-degree circle */
		PSL_plotarc (PSL, GMT->current.proj.c_x0, GMT->current.proj.c_y0, radius, 0.0, 360.0, PSL_MOVE|PSL_STROKE);
		PSL_plotarc (PSL, GMT->current.proj.c_x0, GMT->current.proj.c_y0, radius + r_inc, 0.0, 360.0, PSL_MOVE|PSL_STROKE);
		if (secondary_too) PSL_plotarc (PSL, GMT->current.proj.c_x0, GMT->current.proj.c_y0, radius + 2.0 * r_inc, 0.0, 360.0, PSL_MOVE|PSL_STROKE);
	}
	else {
		az1 = d_atan2d (y1 - GMT->current.proj.c_y0, x1 - GMT->current.proj.c_x0);
		az2 = d_atan2d (y2 - GMT->current.proj.c_y0, x2 - GMT->current.proj.c_x0);
		if (!GMT->current.proj.north_pole) gmt_M_double_swap (az1, az2);	/* In S hemisphere, must draw in opposite direction */
		while (az1 < 0.0) az1 += 360.0;	/* Wind az1 to be in the 0-360 range */
		while (az2 < az1) az2 += 360.0;	/* Likewise ensure az1 > az1 and is now in the 0-720 range */
		da0 = R2D * escl * width /radius;
		da  = R2D * escl * width / (radius + r_inc);
		PSL_plotarc (PSL, GMT->current.proj.c_x0, GMT->current.proj.c_y0, radius, az1-fw*da0, az2+fe*da0, PSL_MOVE|PSL_STROKE);
		PSL_plotarc (PSL, GMT->current.proj.c_x0, GMT->current.proj.c_y0, radius + r_inc, az1-fw*da, az2+fe*da, PSL_MOVE|PSL_STROKE);
		if (secondary_too) {
			r_inc *= 2.0;
			da = R2D * escl * width / (radius + r_inc);
			PSL_plotarc (PSL, GMT->current.proj.c_x0, GMT->current.proj.c_y0, radius + r_inc, az1-fw*da, az2+fe*da, PSL_MOVE|PSL_STROKE);
		}
	}
}

GMT_LOCAL void gmtplot_rounded_framecorners (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, bool secondary_too) {
	unsigned int k, kn, item[2] = {GMT_TICK_UPPER, GMT_TICK_LOWER};
	double x1, y1, x2, y2, anglew, anglee, x, y, width;
	struct GMT_PLOT_AXIS_ITEM *Tx = NULL, *Ty = NULL;

	if (GMT->current.setting.map_frame_type != GMT_IS_ROUNDED) return;	/* Only do this if rounded corners are requested */

	gmt_geo_to_xy (GMT, w, n, &x1, &y1);
	gmt_geo_to_xy (GMT, w, s, &x2, &y2);
	anglew = d_atan2d (y2 - y1, x2 - x1);

	gmt_geo_to_xy (GMT, e, s, &x1, &y1);
	gmt_geo_to_xy (GMT, e, n, &x2, &y2);
	anglee = d_atan2d (y2 - y1, x2 - x1);

	width = ((secondary_too) ? 0.5 : 1.0) * fabs (GMT->current.setting.map_frame_width);
	kn = (secondary_too) ? 2 : 1;
	for (k = 0; k < kn; k++) {
		Tx = &GMT->current.map.frame.axis[GMT_X].item[item[k]];
		Ty = &GMT->current.map.frame.axis[GMT_Y].item[item[k]];
		if (!(Tx->active && Ty->active)) continue;
		if (GMT->current.map.frame.side[S_SIDE] && GMT->current.map.frame.side[E_SIDE]) {
			gmt_geo_to_xy (GMT, e, s, &x, &y);
			PSL_plotarc (PSL, x, y, (k+1)*width, 180.0+anglee, 270.0+anglee, PSL_MOVE|PSL_STROKE);
		}
		if (GMT->current.map.frame.side[E_SIDE] && GMT->current.map.frame.side[N_SIDE]) {
			gmt_geo_to_xy (GMT, e, n, &x, &y);
			PSL_plotarc (PSL, x, y, (k+1)*width, 270.0+anglee, 360.0+anglee, PSL_MOVE|PSL_STROKE);
		}
		if (GMT->current.map.frame.side[N_SIDE] && GMT->current.map.frame.side[W_SIDE]) {
			gmt_geo_to_xy (GMT, w, n, &x, &y);
			PSL_plotarc (PSL, x, y, (k+1)*width, 180.0+anglew, 270.0+anglew, PSL_MOVE|PSL_STROKE);
		}
		if (GMT->current.map.frame.side[W_SIDE] && GMT->current.map.frame.side[S_SIDE]) {
			gmt_geo_to_xy (GMT, w, s, &x, &y);
			PSL_plotarc (PSL, x, y, (k+1)*width, 270.0+anglew, 360.0+anglew, PSL_MOVE|PSL_STROKE);
		}
		if ((gmt_M_is_azimuthal(GMT) || gmt_M_is_conical(GMT)) && GMT->current.map.frame.side[W_SIDE] && GMT->current.map.frame.side[E_SIDE]) {	/* Round off the pointy head? */
			if (doubleAlmostEqual (GMT->common.R.wesn[YHI], 90.0)) {
				gmt_geo_to_xy (GMT, w, n, &x, &y);
				PSL_plotarc (PSL, x, y, (k+1)*width, anglee, 180.0+anglew, PSL_MOVE|PSL_STROKE);
			}
			else if (doubleAlmostEqual (GMT->common.R.wesn[YLO], -90.0)) {
				gmt_geo_to_xy (GMT, w, s, &x, &y);
				PSL_plotarc (PSL, x, y, (k+1)*width, anglew-90.0, anglee-90.0, PSL_MOVE|PSL_STROKE);
			}
		}
	}
}

#if 0
/* Nov-11-2014 PW: For reference until we know there are no side effects with the new one below */
GMT_LOCAL void gmtplot_wesn_map_boundary_old (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	uint64_t i, np = 0;
	double *xx = NULL, *yy = NULL;

	gmt_setpen (GMT, &GMT->current.setting.map_frame_pen);

	if (GMT->current.map.frame.side[W_SIDE] && GMT->current.map.frame.side[E_SIDE] && GMT->current.map.frame.side[S_SIDE] && GMT->current.map.frame.side[N_SIDE]) {
		/* Want the entire boundary so use a single path to void notches at corners */
		np = gmt_graticule_path (GMT, &xx, &yy, 1, false, w, e, s, n);
		for (i = 0; i < np; i++)
			gmt_geo_to_xy (GMT, xx[i], yy[i], &xx[i], &yy[i]);
		PSL_plotline (PSL, xx, yy, (int)np, PSL_MOVE + PSL_STROKE + PSL_CLOSE);
		gmt_M_free (GMT, xx);
		gmt_M_free (GMT, yy);
		return;
	}

	/* Just do the sides that were requested */

	if (GMT->current.map.frame.side[W_SIDE]) {	/* West */
		np = gmtlib_map_path (GMT, w, s, w, n, &xx, &yy);
		for (i = 0; i < np; i++)
			gmt_geo_to_xy (GMT, xx[i], yy[i], &xx[i], &yy[i]);
		PSL_plotline (PSL, xx, yy, (int)np, PSL_MOVE + PSL_STROKE);
		gmt_M_free (GMT, xx);
		gmt_M_free (GMT, yy);
	}
	if (GMT->current.map.frame.side[E_SIDE]) {	/* East */
		np = gmtlib_map_path (GMT, e, s, e, n, &xx, &yy);
		for (i = 0; i < np; i++)
			gmt_geo_to_xy (GMT, xx[i], yy[i], &xx[i], &yy[i]);
		PSL_plotline (PSL, xx, yy, (int)np, PSL_MOVE + PSL_STROKE);
		gmt_M_free (GMT, xx);
		gmt_M_free (GMT, yy);
	}
	if (GMT->current.map.frame.side[S_SIDE]) {	/* South */
		np = gmtlib_map_path (GMT, w, s, e, s, &xx, &yy);
		for (i = 0; i < np; i++)
			gmt_geo_to_xy (GMT, xx[i], yy[i], &xx[i], &yy[i]);
		PSL_plotline (PSL, xx, yy, (int)np, PSL_MOVE + PSL_STROKE);
		gmt_M_free (GMT, xx);
		gmt_M_free (GMT, yy);
	}
	if (GMT->current.map.frame.side[N_SIDE]) {	/* North */
		np = gmtlib_map_path (GMT, w, n, e, n, &xx, &yy);
		for (i = 0; i < np; i++)
			gmt_geo_to_xy (GMT, xx[i], yy[i], &xx[i], &yy[i]);
		PSL_plotline (PSL, xx, yy, (int)np, PSL_MOVE + PSL_STROKE);
		gmt_M_free (GMT, xx);
		gmt_M_free (GMT, yy);
	}
}
#endif

GMT_LOCAL void gmtplot_wesn_map_boundary (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	/* Draw 0-4 boundary sides.  If more than 1 then ensure we draw a continuous line to
	 * avoid notches at sharp corners. */
	uint64_t i, n_sides = 0, np = 0, n_set = 0;
	int this, next = GMT_NOTSET, flag;
	double lonstart[4] = {w, e, e, w}, lonstop[4] = {e, e, w, w};
	double latstart[4] = {s, s, n, n}, latstop[4] = {s, n, n, s};
	double *xx = NULL, *yy = NULL;

	gmt_setpen (GMT, &GMT->current.setting.map_frame_pen);

	/* Determine how many sides are requested and find first side to be skipped (if any) */

	for (this = S_SIDE; this <= W_SIDE; this++) {
		if (GMT->current.map.frame.side[this]) n_sides++;
		else if (next == GMT_NOTSET) next = this;
	}
	if (n_sides == 0) return;	/* Nuthin' to do */

	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "gmtplot_wesn_map_boundary n_sides = %" PRIu64 "\n", n_sides);
	this = next;		/* First side to be skipped (== GMT_NOTSET if none to be skipped) */
	flag = PSL_MOVE;	/* Need to move to the start of the line the first time */
	while (n_set < n_sides) {	/* While more sides need to be plotted we loop counter-clockwise */
		this++;		/* Go to next side (this will be S_SIDE if all 4 sides should be plotted) */
		if (this > W_SIDE) this = S_SIDE;	/* Wrap around */
		if (!GMT->current.map.frame.side[this]) {	/* Side to be skipped */
			flag = PSL_MOVE;	/* Need to move to the start of the new line after the gap */
			continue;
		}
		/* Get coordinates for this border */
		np = gmtlib_map_path (GMT, lonstart[this], latstart[this], lonstop[this], latstop[this], &xx, &yy);
		for (i = 0; i < np; i++) /* Project to inches */
			gmt_geo_to_xy (GMT, xx[i], yy[i], &xx[i], &yy[i]);
		next = this + 1;	/* Find ID of next side */
		if (next > W_SIDE) next = S_SIDE;	/* Wrap around */
		if (!GMT->current.map.frame.side[next]) flag |= PSL_STROKE;	/* Must stroke line since a gap follows */
		else if (n_sides == 4 && this == W_SIDE) flag |= (PSL_STROKE+PSL_CLOSE);	/* Must close and stroke line since no gaps */
		GMT_Report (GMT->parent, GMT_MSG_DEBUG, "gmtplot_wesn_map_boundary doing side %d with flag = %d\n", this, flag);

		PSL_plotline (PSL, xx, yy, (int)np, flag);
		gmt_M_free (GMT, xx);
		gmt_M_free (GMT, yy);
		n_set++;
		flag = 0;
	}
}

GMT_LOCAL double gmtplot_fancy_fat_pen_width (struct GMT_CTRL *GMT) {
	/* Return the pen width of the black/white checker fancy frame in points */
	return (fabs (0.01 * GMT->current.setting.map_frame_percent * GMT->current.setting.map_frame_width) * GMT->session.u2u[GMT_INCH][GMT_PT]);
}

GMT_LOCAL double gmtplot_fancy_thin_pen_width (struct GMT_CTRL *GMT, double width) {
	/* Return the pen width of the enclosing checker fancy frame in points as 10% of fat pen width */
	return (0.1 * width / ( GMT->current.setting.map_frame_percent * 0.01 ));
}

GMT_LOCAL void gmtplot_fancy_map_boundary (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	double fat_pen, thin_pen;
	bool dual = false;
	unsigned int cap = PSL->internal.line_cap;

	if (!(GMT->current.setting.map_frame_type & GMT_IS_FANCY)) {	/* Draw plain boundary and return */
		gmtplot_wesn_map_boundary (GMT, PSL, w, e, s, n);
		return;
	}

	PSL_setcolor (PSL, GMT->current.setting.map_frame_pen.rgb, PSL_IS_STROKE);

	fat_pen = gmtplot_fancy_fat_pen_width (GMT);
	if (GMT->current.map.frame.axis[GMT_Y].item[GMT_TICK_LOWER].active) {	/* Need two-layer frame */
		fat_pen *= 0.5;
		dual = true;
	}
	thin_pen = gmtplot_fancy_thin_pen_width (GMT, fat_pen);

	/* Draw frame checkers */
	/* This needs to be done with BUTT cap since checker segments are drawn as heavy lines */

	PSL_setlinewidth (PSL, fat_pen);
	PSL_setlinecap (PSL, PSL_BUTT_CAP);

	gmtplot_fancy_frame_straightlat_checkers (GMT, PSL, w, e, s, n, dual);
	gmtplot_fancy_frame_straightlon_checkers (GMT, PSL, w, e, s, n, dual);

	/* Draw the outline on top of the checkers */
	/* Reset line cap, etc. */

	PSL_setlinecap (PSL, cap);
	PSL_setcolor (PSL, GMT->current.setting.map_frame_pen.rgb, PSL_IS_STROKE);
	PSL_setlinewidth (PSL, thin_pen);

	gmtplot_fancy_frame_straight_outline (GMT, PSL, w, s, e, s, S_SIDE, dual);
	gmtplot_fancy_frame_straight_outline (GMT, PSL, e, s, e, n, E_SIDE, dual);
	gmtplot_fancy_frame_straight_outline (GMT, PSL, e, n, w, n, N_SIDE, dual);
	gmtplot_fancy_frame_straight_outline (GMT, PSL, w, n, w, s, W_SIDE, dual);

	gmtplot_rounded_framecorners (GMT, PSL, w, e, s, n, dual);
}

GMT_LOCAL void gmtplot_rect_map_boundary (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double x0, double y0, double x1, double y1) {
	unsigned int cap = PSL->internal.line_cap;

	gmt_setpen (GMT, &GMT->current.setting.map_frame_pen);
	/* Temporarily change to square cap so rectangular frames have neat corners */
	PSL_setlinecap (PSL, PSL_SQUARE_CAP);

	if (GMT->current.map.frame.side[W_SIDE] & GMT_AXIS_DRAW) PSL_plotsegment (PSL, x0, y0, x0, y1);	/* West */
	if (GMT->current.map.frame.side[E_SIDE] & GMT_AXIS_DRAW) PSL_plotsegment (PSL, x1, y0, x1, y1);	/* East */
	if (GMT->current.map.frame.side[S_SIDE] & GMT_AXIS_DRAW) PSL_plotsegment (PSL, x0, y0, x1, y0);	/* South */
	if (GMT->current.map.frame.side[N_SIDE] & GMT_AXIS_DRAW) PSL_plotsegment (PSL, x0, y1, x1, y1);	/* North */
	PSL_setlinecap (PSL, cap);	/* Reset back to default */
}

/*	GMT_POLAR (S or N) PROJECTION MAP BOUNDARY	*/

GMT_LOCAL void gmtplot_polar_map_boundary (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	bool dual = false;
	unsigned int cap = PSL->internal.line_cap;
	double thin_pen, fat_pen, lon_range = e - w;
	double x0, x1, x2, y0, y1, y2, start, stop;

	if (GMT->common.R.oblique) { /* Draw rectangular boundary and return */
		gmtplot_rect_map_boundary (GMT, PSL, 0.0, 0.0, GMT->current.proj.rect[XHI], GMT->current.proj.rect[YHI]);
		return;
	}

	if (!GMT->current.proj.north_pole && gmt_M_is_Spole (s)) /* Cannot have southern boundary */
		GMT->current.map.frame.side[S_SIDE] = GMT_AXIS_NONE;
	if (GMT->current.proj.north_pole && gmt_M_is_Npole (n)) /* Cannot have northern boundary */
		GMT->current.map.frame.side[N_SIDE] = GMT_AXIS_NONE;
	if (gmt_M_360_range (w, e) || doubleAlmostEqualZero (e, w))
		GMT->current.map.frame.side[E_SIDE] = GMT->current.map.frame.side[W_SIDE] = GMT_AXIS_NONE;

	if (!(GMT->current.setting.map_frame_type & GMT_IS_FANCY)) {	/* Draw plain boundary and return */
		gmtplot_wesn_map_boundary (GMT, PSL, w, e, s, n);
		return;
	}

	/* Here draw fancy map boundary */

	fat_pen = gmtplot_fancy_fat_pen_width (GMT);
	if (GMT->current.map.frame.axis[GMT_Y].item[GMT_TICK_LOWER].active) {	/* Need two-layer frame */
		fat_pen *= 0.5;
		dual = true;
	}
	thin_pen = gmtplot_fancy_thin_pen_width (GMT, fat_pen);

	/* Draw frame checkers */
	/* This needs to be done with BUTT cap since checker segments are drawn as heavy lines */

	PSL_setlinewidth (PSL, fat_pen);
	PSL_setlinecap (PSL, PSL_BUTT_CAP);

	gmtplot_fancy_frame_straightlat_checkers (GMT, PSL, w, e, s, n, dual);
	gmtplot_fancy_frame_curvedlon_checkers (GMT, PSL, w, e, s, n, dual);

	/* Draw the outline on top of the checkers */
	/* Reset line cap, etc. */

	PSL_setlinecap (PSL, cap);
	PSL_setcolor (PSL, GMT->current.setting.map_frame_pen.rgb, PSL_IS_STROKE);
	PSL_setlinewidth (PSL, thin_pen);

	gmtplot_fancy_frame_curved_outline (GMT, PSL, w, s, e, s, S_SIDE, dual);
	gmtplot_fancy_frame_straight_outline (GMT, PSL, e, s, e, n, E_SIDE, dual);
	gmtplot_fancy_frame_curved_outline (GMT, PSL, w, n, e, n, N_SIDE, dual);
	gmtplot_fancy_frame_straight_outline (GMT, PSL, w, n, w, s, W_SIDE, dual);

	gmtplot_rounded_framecorners (GMT, PSL, w, e, s, n, dual);

	if (GMT->current.proj.north_pole && gmt_M_is_Npole (n) && !(doubleAlmostEqual (lon_range, 180.0) || doubleAlmostEqual (lon_range, 360.0))) {
		/* Connect the outer straight borders with an arc using map frame pen */
		double radius = fat_pen * GMT->session.u2u[GMT_PT][GMT_INCH];	/* Get radius in inches */
		PSL_setlinewidth (PSL, thin_pen);
		gmt_geo_to_xy (GMT, w, n, &x0, &y0);
		gmt_geo_to_xy (GMT, w, s, &x1, &y1);
		stop = d_atan2d (y1 - y0, x1 - x0) - 90;
		gmt_geo_to_xy (GMT, e, s, &x2, &y2);
		start = d_atan2d (y2 - y0, x2 - x0) + 90;
		if (stop < start) stop += 360.0;
		if (stop > start) PSL_plotarc (PSL, x0, y0, radius, start, stop, PSL_MOVE|PSL_STROKE);
	}
	else if (!GMT->current.proj.north_pole && gmt_M_is_Spole (s) && !(doubleAlmostEqual (lon_range, 180.0) || doubleAlmostEqual (lon_range, 360.0))) {
		/* Connect the outer straight borders with an arc using map frame pen */
		double radius = fat_pen * GMT->session.u2u[GMT_PT][GMT_INCH];	/* Get radius in inches */
		PSL_setlinewidth (PSL, thin_pen);
		gmt_geo_to_xy (GMT, e, s, &x0, &y0);
		gmt_geo_to_xy (GMT, e, n, &x1, &y1);
		stop = d_atan2d (y1 - y0, x1 - x0) - 90;
		gmt_geo_to_xy (GMT, w, n, &x2, &y2);
		start = d_atan2d (y2 - y0, x2 - x0) + 90;
		if (stop < start) stop += 360.0;
		if (stop > start) PSL_plotarc (PSL, x0, y0, radius, start, stop, PSL_MOVE|PSL_STROKE);
	}
}

/*	CONIC PROJECTION MAP BOUNDARY	*/

GMT_LOCAL void gmtplot_conic_map_boundary (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	bool dual = false;
	unsigned int cap = PSL->internal.line_cap;
	double thin_pen, fat_pen;

	if (GMT->common.R.oblique) { /* Draw rectangular boundary and return */
		gmtplot_rect_map_boundary (GMT, PSL, 0.0, 0.0, GMT->current.proj.rect[XHI], GMT->current.proj.rect[YHI]);
		return;
	}

	if (!(GMT->current.setting.map_frame_type & GMT_IS_FANCY)) {	/* Draw plain boundary and return */
		gmtplot_wesn_map_boundary (GMT, PSL, w, e, s, n);
		return;
	}

	/* Here draw fancy map boundary */

	if (!GMT->current.proj.north_pole && gmt_M_is_Spole (s)) /* Cannot have southern boundary */
		GMT->current.map.frame.side[S_SIDE] = GMT_AXIS_NONE;
	if (GMT->current.proj.north_pole && gmt_M_is_Npole (n)) /* Cannot have northern boundary */
		GMT->current.map.frame.side[N_SIDE] = GMT_AXIS_NONE;

	fat_pen = gmtplot_fancy_fat_pen_width (GMT);
	if (GMT->current.map.frame.axis[GMT_Y].item[GMT_TICK_LOWER].active) {	/* Need two-layer frame */
		fat_pen *= 0.5;
		dual = true;
	}
	thin_pen = gmtplot_fancy_thin_pen_width (GMT, fat_pen);

	/* Draw frame checkers */
	/* This needs to be done with BUTT cap since checker segments are drawn as heavy lines */

	PSL_setlinewidth (PSL, fat_pen);
	PSL_setlinecap (PSL, PSL_BUTT_CAP);

	gmtplot_fancy_frame_straightlat_checkers (GMT, PSL, w, e, s, n, dual);
	gmtplot_fancy_frame_curvedlon_checkers (GMT, PSL, w, e, s, n, dual);

	/* Draw the outline on top of the checkers */
	/* Reset line cap, etc. */

	PSL_setlinecap (PSL, cap);
	PSL_setcolor (PSL, GMT->current.setting.map_frame_pen.rgb, PSL_IS_STROKE);
	PSL_setlinewidth (PSL, thin_pen);

	gmtplot_fancy_frame_curved_outline (GMT, PSL, w, s, e, s, S_SIDE, dual);
	gmtplot_fancy_frame_straight_outline (GMT, PSL, e, s, e, n, E_SIDE, dual);
	gmtplot_fancy_frame_curved_outline (GMT, PSL, w, n, e, n, N_SIDE, dual);
	gmtplot_fancy_frame_straight_outline (GMT, PSL, w, n, w, s, W_SIDE, dual);

	gmtplot_rounded_framecorners (GMT, PSL, w, e, s, n, dual);
}

/*	OBLIQUE MERCATOR PROJECTION MAP FUNCTIONS	*/

GMT_LOCAL void gmtplot_oblmrc_map_boundary (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	gmt_M_unused(w); gmt_M_unused(e); gmt_M_unused(s); gmt_M_unused(n);
	gmtplot_rect_map_boundary (GMT, PSL, 0.0, 0.0, GMT->current.proj.rect[XHI], GMT->current.proj.rect[YHI]);
}

/*	MOLLWEIDE and HAMMER-AITOFF EQUAL AREA PROJECTION MAP FUNCTIONS	*/

GMT_LOCAL void gmtplot_ellipse_map_boundary (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	if (GMT->common.R.oblique) { /* Draw rectangular boundary and return */
		gmtplot_rect_map_boundary (GMT, PSL, 0.0, 0.0, GMT->current.proj.rect[XHI], GMT->current.proj.rect[YHI]);
		return;
	}
	gmtplot_wesn_map_boundary (GMT, PSL, w, e, s, n);	/* Draw outline first, then turn off non-existent sides */
	if (gmt_M_is_Spole (GMT->common.R.wesn[YLO])) /* Cannot have southern boundary */
		GMT->current.map.frame.side[S_SIDE] = GMT_AXIS_NONE;
	if (gmt_M_is_Npole (GMT->common.R.wesn[YHI])) /* Cannot have northern boundary */
		GMT->current.map.frame.side[N_SIDE] = GMT_AXIS_NONE;
}

GMT_LOCAL void gmtplot_basic_map_boundary (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	if (GMT->common.R.oblique) { /* Draw rectangular boundary and return */
		gmtplot_rect_map_boundary (GMT, PSL, 0.0, 0.0, GMT->current.proj.rect[XHI], GMT->current.proj.rect[YHI]);
		return;
	}
	gmtplot_wesn_map_boundary (GMT, PSL, w, e, s, n);
}

/*
 *	GENERIC MAP PLOTTING FUNCTIONS
 */

GMT_LOCAL int gmtplot_genper_map_boundary (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	uint64_t nr;
	gmt_M_unused(w); gmt_M_unused(e); gmt_M_unused(s); gmt_M_unused(n);

	if (GMT->common.R.oblique) {	/* Draw rectangular boundary and return */
		gmtplot_rect_map_boundary (GMT, PSL, 0.0, 0.0, GMT->current.proj.rect[XHI], GMT->current.proj.rect[YHI]);
		return 0;
	}

	gmt_setpen (GMT, &GMT->current.setting.map_frame_pen);

	nr = GMT->current.map.n_lon_nodes + GMT->current.map.n_lat_nodes;
	if (nr >= GMT->current.plot.n_alloc) gmt_get_plot_array (GMT);

	if (GMT->current.proj.g_debug > 1) GMT_Report (GMT->parent, GMT_MSG_DEBUG, "genper_map_boundary nr = %" PRIu64 "\n", nr);

	gmtlib_genper_map_clip_path (GMT, nr, GMT->current.plot.x, GMT->current.plot.y);

	PSL_plotline (PSL, GMT->current.plot.x, GMT->current.plot.y, (int)nr, PSL_MOVE|PSL_STROKE);

	return 0;
}

GMT_LOCAL void gmtplot_circle_map_boundary (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	gmt_M_unused(w); gmt_M_unused(e); gmt_M_unused(s); gmt_M_unused(n);
	if (GMT->common.R.oblique) { /* Draw rectangular boundary and return */
		gmtplot_rect_map_boundary (GMT, PSL, 0.0, 0.0, GMT->current.proj.rect[XHI], GMT->current.proj.rect[YHI]);
		return;
	}

	gmt_setpen (GMT, &GMT->current.setting.map_frame_pen);

	PSL_plotarc (PSL, GMT->current.proj.r, GMT->current.proj.r, GMT->current.proj.r, 0.0, 360.0, PSL_MOVE|PSL_STROKE);
}

GMT_LOCAL void gmtplot_theta_r_map_boundary (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	bool circles = false;
	unsigned int flag = PSL_CLOSE;
	uint64_t i, k = 0, nr, n_max_path;
	double a, da;

	/* This boundary may be a circle, a donut, a wedge, or a donut with a missing sector.  In all cases we wish
	 * to draw the outline in one go to ensure the corners are joined properly and not leaving small notches. */

	gmt_setpen (GMT, &GMT->current.setting.map_frame_pen);

	if (GMT->current.proj.flip) {	/* Flipped, so N is the inside circle; check if radius matches n */
		if (doubleAlmostEqual (n, GMT->current.proj.flip_radius) && gmt_M_is_zero (GMT->current.proj.radial_offset))
			GMT->current.map.frame.side[N_SIDE] = GMT_AXIS_NONE;	/* No donuts, please */
	}
	else {	/* No flip, so s is inside circle. Check if s is zero and there is not an offset */
		if (gmt_M_is_zero (s) && gmt_M_is_zero (GMT->current.proj.radial_offset))
			GMT->current.map.frame.side[S_SIDE] = GMT_AXIS_NONE;		/* No donuts, please */
	}
	if (gmt_M_360_range (w, e) || doubleAlmostEqualZero (e, w)) {	/* Draw a full 360 circle so no E/W sides will be drawn */
		GMT->current.map.frame.side[E_SIDE] = GMT->current.map.frame.side[W_SIDE] = GMT_AXIS_NONE;
		circles = true;
	}

	nr = GMT->current.map.n_lon_nodes;	/* Points needed to draw a full circle */
	n_max_path = 2 * (nr + 2);	/* Max length of boundary */
	while (n_max_path > GMT->current.plot.n_alloc) gmt_get_plot_array (GMT);	/* Ensure we have enough plot memory */
	da = fabs (GMT->common.R.wesn[XHI] - GMT->common.R.wesn[XLO]) / (nr - 1);	/* Steps in azimuth along the curved boundary */
	if (GMT->current.map.frame.side[N_SIDE] & GMT_AXIS_DRAW) {	/* Must draw the N circular boundary from W (XLO) to E (XHI), all at YHI */
		for (i = 0; i < nr; i++, k++) {
			a = GMT->common.R.wesn[XLO] + i * da;
			gmt_geo_to_xy (GMT, a, GMT->common.R.wesn[YHI], &GMT->current.plot.x[k], &GMT->current.plot.y[k]);
		}
		if (circles) {	/* Nothing to connect to, so plot this circle we have as is */
			PSL_plotline (PSL, GMT->current.plot.x, GMT->current.plot.y, (int)nr, PSL_MOVE|PSL_STROKE|PSL_CLOSE);
			k = 0;	/* Start all over in case another circle is needed */
		}
	}
	/* Now at E (XHI, YHI).  If we need to add a radial E boundary then add it now to the array, ending at XHI, YLO */
	if (GMT->current.map.frame.side[E_SIDE] & GMT_AXIS_DRAW) {
		gmt_geo_to_xy (GMT, GMT->common.R.wesn[XHI], GMT->common.R.wesn[YHI], &GMT->current.plot.x[k], &GMT->current.plot.y[k]);	k++;
		gmt_geo_to_xy (GMT, GMT->common.R.wesn[XHI], GMT->common.R.wesn[YLO], &GMT->current.plot.x[k], &GMT->current.plot.y[k]);	k++;
	}
	/* Now at E (XHI, YLO).  Do we hook in the other partial circle? */
	if (GMT->current.map.frame.side[S_SIDE] & GMT_AXIS_DRAW) {	/* Must draw the S circular boundary from E (XHI) backwards to W (XLO), all at YLO */
		for (i = 0; i < nr; i++, k++) {
			a = GMT->common.R.wesn[XHI] - i * da;
			gmt_geo_to_xy (GMT, a, GMT->common.R.wesn[YLO], &GMT->current.plot.x[k], &GMT->current.plot.y[k]);
		}
		if (circles) {	/* Nothing to connect to, so plot this circle as is */
			PSL_plotline (PSL, GMT->current.plot.x, GMT->current.plot.y, (int)nr, PSL_MOVE|PSL_STROKE);
			k = 0;	/* Start all over */
		}
	}
	else if (k) {   /* Must plot now */
		PSL_plotline (PSL, GMT->current.plot.x, GMT->current.plot.y, (int)k, PSL_MOVE|PSL_STROKE);
		k = 0;
	}
	/* Now at W (XLO, YLO).  If need to add a radial W boundary add it now to the array, ending at XLO, YHI (where we started) */
	if (GMT->current.map.frame.side[W_SIDE] & GMT_AXIS_DRAW) {
		gmt_geo_to_xy (GMT, GMT->common.R.wesn[XLO], GMT->common.R.wesn[YLO], &GMT->current.plot.x[k], &GMT->current.plot.y[k]);	k++;
		gmt_geo_to_xy (GMT, GMT->common.R.wesn[XLO], GMT->common.R.wesn[YHI], &GMT->current.plot.x[k], &GMT->current.plot.y[k]);	k++;
	}
	for (i = S_SIDE; i <= W_SIDE; i++) if ((GMT->current.map.frame.side[i] & GMT_AXIS_DRAW) == 0) flag = 0;
	if (k)	/* Finally draw the strange donut/pacman boundary as one piece */
		PSL_plotline (PSL, GMT->current.plot.x, GMT->current.plot.y, (int)k, PSL_MOVE|PSL_STROKE|flag);
}

GMT_LOCAL void gmtplot_map_tick (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double *xx, double *yy, unsigned int *sides, double *angles, unsigned int nx, unsigned int type, double len) {
	double angle, xl, yl, c, s, tick_length;
	unsigned int i;
	bool set_angle;

	/* The set_angle bit trieds to handle the fact that the given angles from crossing may need an adjustment depending on
	 * which side of a rectangular box it occurs.  There are exception for round maps, etc.  It is a bit nebulous and could
	 * need a better explanation.  For instance, I commented out the Gnomonic case which is needed for annotations but not here, apparently */
	set_angle = ((!GMT->common.R.oblique && !(gmt_M_is_azimuthal(GMT) || gmt_M_is_conical(GMT))) || GMT->common.R.oblique);
	if (!GMT->common.R.oblique && (GMT->current.proj.projection_GMT == GMT_GENPER || GMT->current.proj.projection_GMT == GMT_POLYCONIC)) set_angle = true;

	for (i = 0; i < nx; i++) {
		if (!GMT->current.proj.edge[sides[i]]) continue;
		if ((GMT->current.map.frame.side[sides[i]] & GMT_AXIS_TICK) == 0) continue;
		if (!(GMT->current.setting.map_annot_oblique & GMT_OBL_ANNOT_ANYWHERE) && ((type == 0 && (sides[i] % 2)) || (type == 1 && !(sides[i] % 2)))) continue;
		angle = ((GMT->current.setting.map_annot_oblique & GMT_OBL_ANNOT_NORMAL_TICKS) ? (sides[i] - 1) * 90.0 : angles[i]);
		if (set_angle) {	/* Adjust angle to fit the range of angles relative to each side */
			if (sides[i] == 0 && angle < 180.0) angle -= 180.0;
			if (sides[i] == 1 && (angle > 90.0 && angle < 270.0)) angle -= 180.0;
			if (sides[i] == 2 && angle > 180.0) angle -= 180.0;
			if (sides[i] == 3 && (angle < 90.0 || angle > 270.0)) angle -= 180.0;
		}
		sincosd (angle, &s, &c);
		tick_length = len;
		if (GMT->current.setting.map_annot_oblique & GMT_OBL_ANNOT_EXTEND_TICKS) {
			if (sides[i] % 2) {
				/* if (fabs (c) > cosd (GMT->current.setting.map_annot_min_angle)) continue; */
				if (fabs (c) < sind (GMT->current.setting.map_annot_min_angle)) continue;
				tick_length /= fabs(c);
			}
			else {
				if (fabs (s) < sind (GMT->current.setting.map_annot_min_angle)) continue;
				tick_length /= fabs(s);
			}
		}
		xl = tick_length * c;
		yl = tick_length * s;
		PSL_plotsegment (PSL, xx[i], yy[i], xx[i]+xl, yy[i]+yl);
	}
}

GMT_LOCAL void gmtplot_map_lontick (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double lon, double south, double north, double len) {
	unsigned int i, nc;
	struct GMT_XINGS *xings = NULL;

	nc = gmtlib_map_loncross (GMT, lon, south, north, &xings);
	for (i = 0; i < nc; i++)
		gmtplot_map_tick (GMT, PSL, xings[i].xx, xings[i].yy, xings[i].sides, xings[i].angle, xings[i].nx, 0, len);
	if (nc) gmt_M_free (GMT, xings);
}

GMT_LOCAL void gmtplot_map_lattick (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double lat, double west, double east, double len) {
	unsigned int i, nc;

	struct GMT_XINGS *xings = NULL;

	nc = gmtlib_map_latcross (GMT, lat, west, east, &xings);
	for (i = 0; i < nc; i++)
		gmtplot_map_tick (GMT, PSL, xings[i].xx, xings[i].yy, xings[i].sides, xings[i].angle, xings[i].nx, 1, len);
	if (nc) gmt_M_free (GMT, xings);
}

GMT_LOCAL bool gmtplot_annot_too_crowded (struct GMT_CTRL *GMT, double x, double y, unsigned int side) {
	/* Checks if the proposed annotation is too close to a previously plotted annotation */
	unsigned int i;
	double d_min;

	if (GMT->current.setting.map_annot_min_spacing <= 0.0) return (false);

	for (i = 0, d_min = DBL_MAX; i < GMT_n_annotations[side]; i++)
		d_min = MIN (d_min, hypot (GMT_x_annotation[side][i] - x, GMT_y_annotation[side][i] - y));
	if (d_min < GMT->current.setting.map_annot_min_spacing) {
		GMT_n_annotations_skip[side]++;
		return (true);
	}

	/* OK to plot and add to list */

	if (GMT_n_annotations[side] == GMT_alloc_annotations[side]) gmt_M_malloc2 (GMT, GMT_x_annotation[side], GMT_y_annotation[side], GMT_n_annotations[side], &(GMT_alloc_annotations[side]), double);
	GMT_x_annotation[side][GMT_n_annotations[side]] = x, GMT_y_annotation[side][GMT_n_annotations[side]] = y, GMT_n_annotations[side]++;

	return (false);
}

/*! . */
GMT_LOCAL double gmtplot_get_annot_offset (struct GMT_CTRL *GMT, bool *flip, unsigned int level) {
	/* Return offset in inches for text annotation.  If annotation
	 * is to be placed 'inside' the map, set flip to true */

	double a = GMT->current.setting.map_annot_offset[level];
	if (GMT->current.setting.map_frame_type & GMT_IS_INSIDE) {	/* Inside annotation */
		a = -fabs (a);
		a -= fabs (GMT->current.setting.map_tick_length[GMT_PRIMARY]);
		*flip = true;
	}
	else if (a >= 0.0) {	/* Outside annotation */
		double dist = GMT->current.setting.map_tick_length[GMT_PRIMARY];	/* Length of tickmark (could be negative) */
		/* For fancy frame we must consider that the frame width might exceed the ticklength */
		if (GMT->current.setting.map_frame_type & GMT_IS_FANCY && GMT->current.setting.map_frame_width > dist) dist = GMT->current.setting.map_frame_width;
		if (dist > 0.0) a += dist;
		*flip = false;
	}
	else {		/* Inside annotation via negative tick length */
		if (GMT->current.setting.map_tick_length[GMT_PRIMARY] < 0.0) a += GMT->current.setting.map_tick_length[GMT_PRIMARY];
		*flip = true;
	}

	return (a);
}

GMT_LOCAL double gmtplot_curved_boundary_offset (struct GMT_CTRL *GMT, double lon, double lat, int type, unsigned int level) {
	/* For global curved maps like Mollweide etc, the parallel labeling occurs against an increasingly
	 * sloping map border as we approach the poles.  Here we try to compensate for that by adding a bit
	 * of extra annotation offset boost that depends on the latitude and size of annotation font */
	double lat1, lat2, dlat, x1, x2, y1, y2, h, boost;
	if (type == 0) return 0.0;	/* No boost for any longitude annotations at top/bottom of these global projections */
	if (!gmt_M_is_misc (GMT)) return 0.0;	/* Not one of the global misc projections */
	if (GMT->common.R.oblique) return 0.0;	/* Not a w/e/s/n curved frame */

	/* OK, here we have Hammer, Robinson, etc. */

	dlat = (GMT->common.R.wesn[YHI] - GMT->common.R.wesn[YLO]) / 180.0;	/* So 1 degree for a global map */
	if (doubleAlmostEqual (lat, -90)) lat1 = lat, lat2 = lat + dlat;	/* South pole point */
	else if (doubleAlmostEqual (lat, 90)) lat1 = lat - dlat, lat2 = lat;	/* North pole point */
	else lat1 = lat - dlat,	lat2 = lat + dlat;
	gmt_geo_to_xy (GMT, lon, lat1, &x1, &y1);
	gmt_geo_to_xy (GMT, lon, lat2, &x2, &y2);
	h = GMT->current.setting.font_annot[level].size * GMT->session.u2u[GMT_PT][GMT_INCH] * GMT_LET_HEIGHT;	/* Approximate height of annotations */
	boost = 0.375 * h * fabs (x2 - x1) / (y2 - y1);	/* Started with 0.5, 0.25 was too small so split the difference */
	return boost;
}

GMT_LOCAL bool gmtplot_skip_end_annotation (struct GMT_CTRL *GMT, unsigned int axis, double x, double length) {
	/* Returns true if this is an annotation at the end of the axis and we are either in GMT_IS_INSIDE mode or user added modifier +e[l|u] to the axis annotation settings */
	if (GMT->current.map.frame.axis[axis].skip[0] && fabs (x) < GMT_CONV4_LIMIT) return true;   /* Skip annotation if it falls ~exactly at the start of the axis */
	if (GMT->current.map.frame.axis[axis].skip[1] && fabs (x - length) < GMT_CONV4_LIMIT) return true;   /* Skip annotation if t falls ~exactly at the end of the axis */
	if (GMT->current.setting.map_frame_type & GMT_IS_INSIDE && (fabs (x) < GMT_CONV4_LIMIT || fabs (x - length) < GMT_CONV4_LIMIT)) return true;   /* Skip annotation on edges when MAP_FRAME_TYPE = inside */
	return false;	/* No, let's plot it */
}

GMT_LOCAL void gmtplot_map_symbol (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, struct GMT_XINGS *xings, char *label, unsigned int type, bool annot, unsigned int level, unsigned int form) {
	/* type = 0 for lon and 1 for lat */

	double line_angle, text_angle, div, tick_length, o_len, len, ca, sa, boost;
	double *xx = xings->xx, *yy = xings->yy, *line_angles = xings->angle;
	unsigned int i, annot_type, justify, *sides = xings->sides, nx = xings->nx;
	bool flip;
	len = gmtplot_get_annot_offset (GMT, &flip, level);	/* Get annotation offset, and flip justification if "inside" */
	annot_type = 2 << type;		/* 2 = NS, 4 = EW */

	for (i = 0; i < nx; i++) {
		if (!(GMT->current.setting.map_annot_oblique & GMT_OBL_ANNOT_ANYWHERE) && ((type == 0 && (sides[i] % 2)) || (type == 1 && !(sides[i] % 2)))) continue;

		if (gmtlib_prepare_label (GMT, line_angles[i], sides[i], xx[i], yy[i], type, &line_angle, &text_angle, &justify)) continue;

		boost = gmtplot_curved_boundary_offset (GMT, xings->lon[i], xings->lat[i], type, level);
		sincosd (line_angle, &sa, &ca);
		tick_length = GMT->current.setting.map_tick_length[GMT_PRIMARY];
		o_len = len;
		if (!flip && GMT->current.setting.map_annot_oblique & annot_type) o_len = tick_length;
		if (GMT->current.setting.map_annot_oblique & GMT_OBL_ANNOT_EXTEND_TICKS) {
			div = ((sides[i] % 2) ? fabs (ca) : fabs (sa));
			o_len /= div;
		}
		else
			o_len += copysign (boost, o_len);
		xx[i] += o_len * ca;
		yy[i] += o_len * sa;
		if (!flip && (GMT->current.setting.map_annot_oblique & annot_type) && GMT->current.setting.map_annot_offset[level] > 0.0) {
			if (sides[i] % 2)
				xx[i] += (sides[i] == 1) ? GMT->current.setting.map_annot_offset[level] : -GMT->current.setting.map_annot_offset[level];
			else
				yy[i] += (sides[i] == 2) ? GMT->current.setting.map_annot_offset[level] : -GMT->current.setting.map_annot_offset[level];
		}

		if (annot) {
			double length = (type == GMT_X) ? GMT->current.map.width : GMT->current.map.height;
			double pos = (type == GMT_X) ? xx[i] : yy[i];
			if (gmtplot_annot_too_crowded (GMT, xx[i], yy[i], sides[i])) continue;
			if (gmtplot_skip_end_annotation (GMT, type, pos, length)) continue;	/* Don't want annotations exactly at one or both ends of the axis */
			if (GMT->current.proj.three_D && GMT->current.proj.z_project.cos_az > 0) {	/* Rotate annotation when seen "from North" */
				if (!flip) justify = gmt_flip_justify (GMT, justify);
				text_angle += 180.0;
			}
			else
				if (flip) justify = gmt_flip_justify (GMT, justify);
			PSL_plottext (PSL, xx[i], yy[i], GMT->current.setting.font_annot[level].size, label, text_angle, justify, form);
		}
	}
}

GMT_LOCAL void gmtplot_map_symbol_ew (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double lat, char *label, double west, double east, bool annot, unsigned int level, unsigned int form) {
	unsigned int i, nc;
	struct GMT_XINGS *xings = NULL;

	nc = gmtlib_map_latcross (GMT, lat, west, east, &xings);
	for (i = 0; i < nc; i++)
		//gmtplot_map_symbol (GMT, PSL, xings[i].xx, xings[i].yy, xings[i].sides, xings[i].angle, label, xings[i].nx, 1, annot, level, form);
		gmtplot_map_symbol (GMT, PSL, &(xings[i]), label, 1, annot, level, form);
	if (nc) gmt_M_free (GMT, xings);
}

GMT_LOCAL void gmtplot_map_symbol_ns (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double lon, char *label, double south, double north, bool annot, unsigned int level, unsigned int form) {
	unsigned int i, k, nc;
	struct GMT_XINGS *xings = NULL;
	bool flip = (gmt_M_type (GMT, GMT_IN, GMT_X) == GMT_IS_LON && gmt_M_type (GMT, GMT_IN, GMT_Y) != GMT_IS_LAT && GMT->current.proj.scale[GMT_Y] < 0.0);
	/* flip deals with the problem when x is lon and geographic annotation machinery is used but y is Cartesian and upside down */
	nc = gmtlib_map_loncross (GMT, lon, south, north, &xings);
	for (i = 0; i < nc; i++) {
		if (flip) for (k = 0; k < xings[i].nx; k++) {	/* Must turn sides 0 and 2 into sides 2 and 0 */
			if ((xings[i].sides[k] % 2) == 0) xings[i].sides[k] = 2 - xings[i].sides[k];	/* Flip up and down sides */
		}
		//gmtplot_map_symbol (GMT, PSL, xings[i].xx, xings[i].yy, xings[i].sides, xings[i].angle, label, xings[i].nx, 0, annot, level, form);
		gmtplot_map_symbol (GMT, PSL, &(xings[i]), label, 0, annot, level, form);
	}
	if (nc) gmt_M_free (GMT, xings);
}

GMT_LOCAL void gmtplot_z_gridlines (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double zmin, double zmax, int plane, unsigned int mode3d, int quadrant) {
	unsigned int k, i, nz, n, item[2] = {GMT_GRID_UPPER, GMT_GRID_LOWER};
	int qplane = 1 - plane;
	double dz, zz, dd, min, max, z0, z1, xx, dummy, value, *z = NULL, *d = NULL;
	char *plane_name[2] = {"y-z", "x-z"};

	if (qplane == GMT_X) {
		value = (quadrant == 4) ? GMT->common.R.wesn[YHI] : GMT->common.R.wesn[YLO];
		min = GMT->current.proj.rect[XLO];	max = GMT->current.proj.rect[XHI];
	}
	else {
		value = (quadrant == 3) ? GMT->common.R.wesn[XHI] : GMT->common.R.wesn[XLO];
		min = GMT->current.proj.rect[YLO];	max = GMT->current.proj.rect[YHI];
	}

	for (k = 0; k < 2; k++) {
		if (fabs (GMT->current.setting.map_grid_cross_size[k]) > 0.0) continue;

		dz = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_Z].type, &GMT->current.map.frame.axis[GMT_Z].item[item[k]]);	/* Gridline spacing in z */

		if (!GMT->current.map.frame.axis[GMT_Z].item[item[k]].active || fabs(dz) == 0.0) continue;

		PSL_comment (PSL, "%s gridlines %s\n", plane_name[plane], k ? "(secondary)" : "(primary)");

		gmt_setpen (GMT, &GMT->current.setting.map_grid_pen[k]);

		nz = gmtlib_coordinate_array (GMT, zmin, zmax, &GMT->current.map.frame.axis[GMT_Z].item[item[k]], &z, NULL);
		for (i = 0; i < nz; i++) {	/* Here z acts as y and x|y acts as x */
			/* Draw one horizontal line */
			zz = gmt_z_to_zz (GMT, z[i]);
			PSL_plotsegment (PSL, min, zz, max, zz);
		}
		z0 = gmt_z_to_zz (GMT, zmin);	/* These are the projected min and max z values, i.e., the ends of vertical gridlines on the back walls */
		z1 = gmt_z_to_zz (GMT, zmax);

		dd = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[qplane].type, &GMT->current.map.frame.axis[qplane].item[item[k]]);
		if (!GMT->current.map.frame.axis[qplane].item[item[k]].active || fabs(dd) == 0.0) continue;
		n = gmtlib_coordinate_array (GMT, GMT->common.R.wesn[2*qplane], GMT->common.R.wesn[2*qplane+1], &GMT->current.map.frame.axis[qplane].item[item[k]], &d, NULL);
		for (i = 0; i < n; i++) {	/* Here z acts as y and x|y acts as x */
			/* Draw one horizontal line */
			if (qplane == GMT_X)	/* The "x" are projection x or longitude */
				gmt_geo_to_xy (GMT, d[i], value, &xx, &dummy);
			else	/* The "x" are projection y or latitude */
				gmt_geo_to_xy (GMT, value, d[i], &dummy, &xx);

			PSL_plotsegment (PSL, xx, z0, xx, z1);
		}
		gmt_M_free (GMT, d);

		PSL_setdash (PSL, NULL, 0);
	}
	if (!GMT->current.map.frame.draw_wall && z && (mode3d & GMT_3D_BOX) == 0) {
		PSL_plotsegment (PSL, min, z0, min, z1);
		PSL_plotsegment (PSL, max, z0, max, z1);
	}
	gmt_M_free (GMT, z);
}

GMT_LOCAL int gmtplot_save_current_gridlines (struct GMT_CTRL *GMT) {
	/* If only primary gridlines are drawn, we save information to gmt.history */

	if (!(GMT->current.map.frame.axis[GMT_X].item[GMT_GRID_UPPER].active || GMT->current.map.frame.axis[GMT_Y].item[GMT_GRID_UPPER].active)) return (GMT_NOERROR);	/* Primary gridlines not selected, so bail */
	if (GMT->current.map.frame.axis[GMT_X].item[GMT_GRID_LOWER].active || GMT->current.map.frame.axis[GMT_Y].item[GMT_GRID_LOWER].active) return (GMT_NOERROR);		/* Secondary gridlines selected, so bail */

	GMT->current.plot.gridline_spacing[GMT_X] = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_X].type, &GMT->current.map.frame.axis[GMT_X].item[GMT_GRID_UPPER]);
	GMT->current.plot.gridline_spacing[GMT_Y] = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_Y].type, &GMT->current.map.frame.axis[GMT_Y].item[GMT_GRID_UPPER]);
	GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "Save current gridline information to gmt.history\n");
	return (GMT_NOERROR);
}

GMT_LOCAL int gmtplot_get_current_gridlines (struct GMT_CTRL *GMT, double *dx, double *dy) {
	/* Obtain the previous gridline intervals, if nonzero */

	if (gmt_M_is_zero (GMT->current.plot.gridline_spacing[GMT_X]) && gmt_M_is_zero (GMT->current.plot.gridline_spacing[GMT_Y])) return (GMT_NOERROR);	/* No gridlines drawn yet */
	*dx = GMT->current.plot.gridline_spacing[GMT_X];
	*dy = GMT->current.plot.gridline_spacing[GMT_Y];
	return (GMT_NOERROR);
}

GMT_LOCAL void gmtplot_map_gridlines (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	unsigned int k, i, np, item[2] = {GMT_GRID_UPPER, GMT_GRID_LOWER};
	double dx, dy, *v = NULL;
	bool reset = false;
	struct GMT_PLOT_AXIS *A[2] = {NULL, NULL};

	for (k = i = 0; k < 2; k++) {	/* First check if any gridlines are requested */
		if (fabs (GMT->current.setting.map_grid_cross_size[k]) > 0.0) continue;
		if (GMT->current.setting.map_grid_cross_type[k] > GMT_CROSS_NORMAL) continue;
		if (!(GMT->current.map.frame.axis[GMT_X].item[item[k]].active || GMT->current.map.frame.axis[GMT_Y].item[item[k]].active)) continue;
		i++;
	}
	if (i == 0) return;	/* No gridlines requested */

	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Entering gmtplot_map_gridlines\n");
	reset = gmtlib_genper_reset (GMT, reset);

	if (gmtplot_save_current_gridlines (GMT)) return;	/* Save a file with primary only gridline interval used under modern mode only; return if error */

	for (k = 0; k < 2; k++) {
		if (fabs (GMT->current.setting.map_grid_cross_size[k]) > 0.0) continue;	/* Must have a size of zero to be a grid line, else it is a cross or tick */

		A[GMT_X] = &GMT->current.map.frame.axis[GMT_X];	/* Short-hand for x-axis */
		A[GMT_Y] = &GMT->current.map.frame.axis[GMT_Y];	/* Short-hand for y-axis */
		dx = gmtlib_get_map_interval (GMT, A[GMT_X]->type, &A[GMT_X]->item[item[k]]);	/* x grid spacing; will be 0 if custom intervals */
		dy = gmtlib_get_map_interval (GMT, A[GMT_Y]->type, &A[GMT_Y]->item[item[k]]);	/* y grid spacing; will be 0 if custom intervals */

		if (!(A[GMT_X]->item[item[k]].active || A[GMT_Y]->item[item[k]].active)) continue;	/* Neither is active */

		PSL_comment (PSL, "%s\n", k ? "Map gridlines (secondary)" : "Map gridlines (primary)");

		gmt_setpen (GMT, &GMT->current.setting.map_grid_pen[k]);

		if (A[GMT_X]->item[k].special && (np = gmtlib_load_custom_annot (GMT, &GMT->current.map.frame.axis[GMT_X], 'g', &v, NULL))) {
			gmtplot_x_grid (GMT, PSL, s, n, v, np);
			gmt_M_free (GMT, v);
		}
		else if (!A[GMT_X]->item[item[k]].active || fabs(dx) == 0.0) { /* Nothing */ }
		else if (GMT->current.proj.xyz_projection[GMT_X] == GMT_TIME)
			plot_timex_grid (GMT, PSL, w, e, s, n, item[k]);
		else if (GMT->current.proj.xyz_projection[GMT_X] == GMT_LOG10)
			gmtplot_logx_grid (GMT, PSL, w, e, s, n, dx);
		else if (GMT->current.proj.xyz_projection[GMT_X] == GMT_POW)
			gmtplot_powx_grid (GMT, PSL, w, e, s, n, dx);
		else if (GMT->current.map.frame.obl_grid)	/* Draw oblique grid lines that go S to N */
			gmtplot_linearx_oblgrid (GMT, PSL, w, e, s, n, dx);
		else	/* Draw grid lines that go S to N */
			gmt_linearx_grid (GMT, PSL, w, e, s, n, dx);

		if (A[GMT_Y]->item[k].special && (np = gmtlib_load_custom_annot (GMT, &GMT->current.map.frame.axis[GMT_Y], 'g', &v, NULL))) {
			gmtplot_y_grid (GMT, PSL, w, e, v, np);
			gmt_M_free (GMT, v);
		}
		else if (!A[GMT_Y]->item[item[k]].active || fabs(dy) == 0.0) { /* Nothing */ }
		else if (GMT->current.proj.xyz_projection[GMT_Y] == GMT_TIME)
			gmtplot_timey_grid (GMT, PSL, w, e, s, n, item[k]);
		else if (GMT->current.proj.xyz_projection[GMT_Y] == GMT_LOG10)
			gmtplot_logy_grid (GMT, PSL, w, e, s, n, dy);
		else if (GMT->current.proj.xyz_projection[GMT_Y] == GMT_POW)
			gmtplot_powy_grid (GMT, PSL, w, e, s, n, dy);
		else if (GMT->current.map.frame.obl_grid)	/* Draw oblique grid lines that go S to N */
			gmtplot_lineary_oblgrid (GMT, PSL, w, e, s, n, dx);
		else	/* Draw grid lines that go E to W */
			gmtplot_lineary_grid (GMT, PSL, w, e, s, n, dy);

		PSL_setdash (PSL, NULL, 0);
	}
	gmtlib_genper_reset (GMT, reset);
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Exiting gmtplot_map_gridlines\n");
}

GMT_LOCAL double gmtplot_cross_angle (struct GMT_CTRL *GMT, double lon, double lat, double dx, double dy, unsigned int type) {
	/* (lon,alt) is the point and we want a tangent angle there in radians */
	if (gmt_M_is_geographic (GMT, GMT_IN)) {	/* General case with arbitrary orientations */
		double x0, y0, x1, y1;
		if (type == GMT_X) {	/* Compute angle from d/dx */
			gmt_geo_to_xy (GMT, lon - dx, lat, &x0, &y0);
			gmt_geo_to_xy (GMT, lon + dx, lat, &x1, &y1);
		}
		else {	/* Get d/dy */
			gmt_geo_to_xy (GMT, lon, ((lat - dy) < -90.0) ? -90.0 : lat - dy, &x0, &y0);
			gmt_geo_to_xy (GMT, lon, ((lat + dy) >  90.0) ?  90.0 : lat + dy, &x1, &y1);
		}
		return (d_atan2 (y1-y0, x1-x0));
	}
	/* Cartesian is always x horizontal, y 90 degrees up */
	return ((type == GMT_X) ? 0.0 : M_PI_2);
}

GMT_LOCAL void gmtplot_map_gridcross (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	unsigned int i, j, k, nx, ny, item[2] = {GMT_GRID_UPPER, GMT_GRID_LOWER};
	double x0, y0, xa, xb, ya, yb, xi, yj, *x = NULL, *y = NULL;
	double angle, Ca, Sa, L, dx, dy;
	struct GMT_PLOT_AXIS *A[2] = {NULL, NULL};

	for (k = i = 0; k < 2; k++)
		if (GMT->current.setting.map_grid_cross_size[k] > 0.0 && GMT->current.setting.map_grid_cross_type[k] == GMT_CROSS_NORMAL) i++;

	if (i == 0) return;	/* No grid crosses requested */

	gmt_map_clip_on (GMT, GMT->session.no_rgb, 3);

	for (k = 0; k < 2; k++) {
		if (gmt_M_is_zero (GMT->current.setting.map_grid_cross_size[k])) continue;

		PSL_comment (PSL, "%s\n", k ? "Map gridcrosses (secondary)" : "Map gridcrosses (primary)");

		gmt_setpen (GMT, &GMT->current.setting.map_grid_pen[k]);

		A[GMT_X] = &GMT->current.map.frame.axis[GMT_X];	/* Short-hand for x-axis */
		A[GMT_Y] = &GMT->current.map.frame.axis[GMT_Y];	/* Short-hand for y-axis */
		nx = gmtlib_coordinate_array (GMT, w, e, &A[GMT_X]->item[item[k]], &x, NULL);
		ny = gmtlib_coordinate_array (GMT, s, n, &A[GMT_Y]->item[item[k]], &y, NULL);
		dy = (ny > 1) ? y[1] - y[0] : GMT->current.map.dlat;
		dx = (nx > 1) ? x[1] - x[0] : GMT->current.map.dlon;

		L = 0.5 * GMT->current.setting.map_grid_cross_size[k];

		for (j = 0; j < ny; j++) {
			yj = y[j];
			for (i = 0; i < nx; i++) {
				if (gmt_M_pole_is_point(GMT) && doubleAlmostEqual (fabs (yj), 90.0)) {	/* Only place one grid cross at the poles for maps where the poles are points */
					xi = GMT->current.proj.central_meridian;
					i = nx;	/* This ends the loop for this particular latitude */
				}
				else
					xi = x[i];

				if (gmt_map_outside (GMT, xi, yj)) continue;	/* Outside map */

				gmt_geo_to_xy (GMT, xi, yj, &x0, &y0);	/* Grid crossing center */
				angle = gmtplot_cross_angle (GMT, xi, yj, dx, dy, GMT_X);
				sincos (angle, &Sa, &Ca);
				xa = x0 - L * Ca;
				xb = x0 + L * Ca;
				ya = y0 - L * Sa;
				yb = y0 + L * Sa;
				PSL_plotsegment (PSL, xa, ya, xb, yb);

				angle += M_PI_2;	/* Since crosses are orthogonal */
				sincos (angle, &Sa, &Ca);
				xa = x0 - L * Ca;
				xb = x0 + L * Ca;
				ya = y0 - L * Sa;
				yb = y0 + L * Sa;
				PSL_plotsegment (PSL, xa, ya, xb, yb);
			}
		}
		if (nx) gmt_M_free (GMT, x);
		if (ny) gmt_M_free (GMT, y);

		PSL_setdash (PSL, NULL, 0);

	}
	gmt_map_clip_off (GMT);
}

GMT_LOCAL void gmtplot_set_gridcross_limbs (struct GMT_CTRL *GMT, unsigned int axis, unsigned int kind, double value, unsigned int *B, unsigned int *E) {
	*B = *E = 0;	/* Default is symmetric tick */
	if (GMT->current.setting.map_grid_cross_type[kind] == GMT_CROSS_SYMM) return;	/* Symmetric tick */
	if (gmt_M_is_zero (value)) return;	/* Symmetrical as well for zero */
	if (gmt_M_type (GMT, GMT_IN, axis) == GMT_IS_LON) {	/* Worry about longitudes */
		if (doubleAlmostEqualZero (fabs (value), 180.0)) return;	/* Symmetrical at both 0 and 180 longitude */
		if (value > 0.0 && value < 180.0) *B = 1; else *E = 1;	/* One-sided */
	}
	else {
		if (value > 0.0) *B = 1; else *E = 1;	/* One-sided */
	}
}

GMT_LOCAL void gmtplot_map_gridticks (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	/* Draw symmetric or asymmetric grid ticks along the courser gridlines for parallels and meridians */
	bool single, point_point = false;
	unsigned int i, j, k, nx, ny, G_nx, G_ny, G_i, G_j, B = 0, E = 0, item[2] = {GMT_GRID_UPPER, GMT_GRID_LOWER};
	double x0, y0, xa, xb, ya, yb, xi, yj, *x = NULL, *y = NULL, *G_x = NULL, *G_y = NULL;
	double angle, Ca, Sa, L, sgn[2] = {1.0, 0.0}, G_dx, G_dy, dx, dy, ys, yn;

	for (k = i = 0; k < 2; k++)
		if (GMT->current.setting.map_grid_cross_type[k] > GMT_CROSS_NORMAL) i++;

	if (i == 0) return;	/* No gridline ticks requested */

	if (gmtplot_get_current_gridlines (GMT, &G_dx, &G_dy)) {
		GMT_Report (GMT->parent, GMT_MSG_ERROR, "Gridline tick embellishments specified but no gridlines have been laid down first\n");
		return;
	}

	if (gmt_M_pole_is_point (GMT)) {	/* Might have two separate domains of gridlines */
		point_point = true;
		if (GMT->current.proj.projection_GMT == GMT_POLAR) {	/* Different for polar graphs since "lat" = 0 is at the center */
			ys =  90.0 - GMT->current.setting.map_polar_cap[0];
			yn = n;
		}
		else {
			ys = MAX (s, -GMT->current.setting.map_polar_cap[0]);
			yn = MIN (n,  GMT->current.setting.map_polar_cap[0]);
		}
	}
	else {	/* No polar cap to worry about */
		ys = s;
		yn = n;
	}
	/* Get the course gridline spacings */
	G_nx = gmtlib_linear_array (GMT, w, e, G_dx, 0.0, &G_x);
	G_ny = gmtlib_linear_array (GMT, ys, yn, G_dy, 0.0, &G_y);

	gmt_map_clip_on (GMT, GMT->session.no_rgb, 3);

	for (k = 0; k < 2; k++) {
		if (gmt_M_is_zero (GMT->current.setting.map_grid_cross_size[k])) continue;

		PSL_comment (PSL, "%s\n", k ? "Map gridticks (secondary)" : "Map gridticks (primary)");

		gmt_setpen (GMT, &GMT->current.setting.map_grid_pen[k]);

		/* Get the detailed gridline tick spacings */
		nx = gmtlib_coordinate_array (GMT, w, e, &GMT->current.map.frame.axis[GMT_X].item[item[k]], &x, NULL);
		ny = gmtlib_coordinate_array (GMT, ys, yn, &GMT->current.map.frame.axis[GMT_Y].item[item[k]], &y, NULL);
		/* Get a small increment in x and y for computing local angle */
		dy = (ny > 1) ? y[1] - y[0] : GMT->current.map.dlat;
		dx = (nx > 1) ? x[1] - x[0] : GMT->current.map.dlon;

		L = 0.5 * fabs (GMT->current.setting.map_grid_cross_size[k]);

		for (G_j = 0; G_j < G_ny; G_j++) {	/* For each gridline parallel */
			yj = G_y[G_j];	/* Current parallel */
			single = (gmt_M_pole_is_point(GMT) && doubleAlmostEqualZero (fabs (yj), 90.0));	/* Only place one grid tick at the poles for maps where the poles are points */
			if (gmt_M_pole_is_point(GMT) && doubleAlmostEqualZero (fabs (yj), 90.0)) continue;	/* No grid tick at single pole points */
			gmtplot_set_gridcross_limbs (GMT, GMT_Y, k, yj, &B, &E);
			for (i = 0; i < nx; i++) {	/* Going along this parallel to place ticks */
				if (single) {
					xi = GMT->current.proj.central_meridian;
					i = nx;	/* This ends the loop for this particular parallel */
				}
				else
					xi = x[i];	/* Current longitude */
				if (gmt_M_is_zero (fmod (xi, G_dx))) continue;	/* Not draw on top of gridlines */
				if (gmt_map_outside (GMT, xi, yj)) continue;	/* Outside map */
				gmt_geo_to_xy (GMT, xi, yj, &x0, &y0);
				angle = gmtplot_cross_angle (GMT, xi, yj, dx, dy, GMT_Y);
				sincos (angle, &Sa, &Ca);
				xa = x0 - L * Ca * sgn[B];
				xb = x0 + L * Ca * sgn[E];
				ya = y0 - L * Sa * sgn[B];
				yb = y0 + L * Sa * sgn[E];
				PSL_plotsegment (PSL, xa, ya, xb, yb);
			}
		}

		for (G_i = 0; G_i < G_nx; G_i++) {	/* For each gridline meridian */
			xi = G_x[G_i];	/* Current meridian */
			gmtplot_set_gridcross_limbs (GMT, GMT_X, k, xi, &B, &E);
			for (j = 0; j < ny; j++) {
				yj = y[j];
				if (point_point && doubleAlmostEqualZero (fabs (yj), GMT->current.setting.map_polar_cap[0])) continue;	/* No latitude ticks long polar circles */
				if (gmt_map_outside (GMT, xi, yj)) continue;	/* Outside map */
				if (gmt_M_pole_is_point(GMT) && doubleAlmostEqualZero (fabs (yj), 90.0)) continue; 	/* No grid tick at single pole points */
				if (gmt_M_is_zero (fmod (yj, G_dy))) continue;	/* Not draw on top of gridlines */
				gmt_geo_to_xy (GMT, xi, yj, &x0, &y0);
				angle = gmtplot_cross_angle (GMT, xi, yj, dx, dy, GMT_X);
				sincos (angle, &Sa, &Ca);
				xa = x0 - L * Ca * sgn[B];
				xb = x0 + L * Ca * sgn[E];
				ya = y0 - L * Sa * sgn[B];
				yb = y0 + L * Sa * sgn[E];
				PSL_plotsegment (PSL, xa, ya, xb, yb);
			}
		}

		if (point_point && k == GMT_SECONDARY && MAX (fabs(s), fabs(n)) > GMT->current.setting.map_polar_cap[0]) {
			/* At least one of the southern or northern polar cap parallels are plotted */
			unsigned int nx_p, ny_p, cap, m;
			double limit[2] = {s, n}, sign[2] = {-1.0, +1.0}, *xp = NULL, *yp = NULL;

			/* Do coarser meridional lines inside polar cap */
			nx_p = gmtlib_linear_array (GMT, w, e, GMT->current.setting.map_polar_cap[1], GMT->current.map.frame.axis[GMT_X].phase, &xp);
			for (cap = 0; cap < 2; cap++) {	/* For south and north cap */
				if (GMT->current.setting.map_polar_cap[0] < fabs (limit[cap])) {	/* Tick this polar cap */
					yj = sign[cap] * GMT->current.setting.map_polar_cap[0];
					gmtplot_set_gridcross_limbs (GMT, GMT_Y, k, yj, &B, &E);
					for (i = 0; i < nx; i++) {	/* Going along this parallel to place ticks */
						xi = x[i];	/* Current longitude */
						if (gmt_M_is_zero (fmod (xi, GMT->current.setting.map_polar_cap[1]))) continue;	/* Not draw on top of coarse gridlines */
						if (gmt_map_outside (GMT, xi, yj)) continue;	/* Outside map */
						gmt_geo_to_xy (GMT, xi, yj, &x0, &y0);
						angle = gmtplot_cross_angle (GMT, xi, yj, dx, dy, GMT_Y);
						sincos (angle, &Sa, &Ca);
						xa = x0 - L * Ca * sgn[B];
						xb = x0 + L * Ca * sgn[E];
						ya = y0 - L * Sa * sgn[B];
						yb = y0 + L * Sa * sgn[E];
						PSL_plotsegment (PSL, xa, ya, xb, yb);
					}
					if (cap == 0) /* South cap */
						ny_p = gmtlib_coordinate_array (GMT, -90.0, ys, &GMT->current.map.frame.axis[GMT_Y].item[item[k]], &yp, NULL);
					else
						ny_p = gmtlib_coordinate_array (GMT, yn, 90.0, &GMT->current.map.frame.axis[GMT_Y].item[item[k]], &yp, NULL);
					for (m = 0; m < nx_p; m++) {	/* For the coarse meridional polar cap lines */
						xi = xp[m];	/* Current meridian */
						gmtplot_set_gridcross_limbs (GMT, GMT_X, k, xi, &B, &E);
						for (j = 0; j < ny_p; j++) {
							yj = yp[j];
							if (point_point && doubleAlmostEqualZero (fabs (yj), GMT->current.setting.map_polar_cap[0])) continue;	/* No latitude ticks long polar circles */
							if (gmt_map_outside (GMT, xi, yj)) continue;	/* Outside map */
							if (gmt_M_pole_is_point(GMT) && doubleAlmostEqualZero (fabs (yj), 90.0)) continue; 	/* No grid tick at single pole points */
							gmt_geo_to_xy (GMT, xi, yj, &x0, &y0);
							angle = gmtplot_cross_angle (GMT, xi, yj, dx, dy, GMT_X);
							sincos (angle, &Sa, &Ca);
							xa = x0 - L * Ca * sgn[B];
							xb = x0 + L * Ca * sgn[E];
							ya = y0 - L * Sa * sgn[B];
							yb = y0 + L * Sa * sgn[E];
							PSL_plotsegment (PSL, xa, ya, xb, yb);
						}
					}
					if (ny_p) gmt_M_free (GMT, yp);
				}
			}
			if (nx_p) gmt_M_free (GMT, xp);
		}

		if (nx) gmt_M_free (GMT, x);
		if (ny) gmt_M_free (GMT, y);

		PSL_setdash (PSL, NULL, 0);
	}
	if (G_nx) gmt_M_free (GMT, G_x);
	if (G_ny) gmt_M_free (GMT, G_y);

	gmt_map_clip_off (GMT);
}

GMT_LOCAL bool gmtplot_skip_polar_apex_annotation (struct GMT_CTRL *GMT, unsigned int i, double *val, unsigned int ny) {
	/* Determine if the W and E annotations on a polar basemap can be placed when the radius from the center to
	 * the annotation location is zero.  This depends on w-e range and a few special cases */
	double annot_height, annot_offset, alpha, beta;
	if (GMT->current.proj.projection_GMT != GMT_POLAR) return false;	/* No action unless the polar -JP|p projection */
	if (!(GMT->current.map.frame.side[W_SIDE] == GMT_AXIS_ALL && GMT->current.map.frame.side[E_SIDE] == GMT_AXIS_ALL)) return false;	/* Requires -BWE to have overprinting issues */
	if (GMT->current.proj.flip) {	/* Here the north value is at the potential apex */
		//if (i < (ny-1) || !doubleAlmostEqualZero (GMT->common.R.wesn[YHI], val[i])) return false;	/* Not apex label */
		if (i < (ny-1) || !doubleAlmostEqualZero (GMT->current.proj.flip_radius, val[i])) return false;	/* Not apex label */
	}
	else {	/* Here radius 0 is at the potential apex */
		if (i > 0 || !gmt_M_is_zero (val[i])) return false;	/*  Not apex label  */
	}
	/* OK, here the label for W or E is right at the apex.  Only plot if overprinting is unlikely */
	if (doubleAlmostEqualZero (GMT->common.R.wesn[XHI] - GMT->common.R.wesn[XLO], 180.0)) return false;	/* No apex when angular range is exactly 180 */
	/* Determine maximum apex angle that still avoids overlap - this depends on fontsize and ticklength */
	annot_height  = GMT_LET_HEIGHT * GMT->current.setting.font_annot[GMT_PRIMARY].size / PSL_POINTS_PER_INCH;	/* Approximate height of annotation in inches */
	annot_offset = MAX (0.0, GMT->current.setting.map_annot_offset[GMT_PRIMARY]) + MAX (0.0, GMT->current.setting.map_tick_length[GMT_ANNOT_UPPER]);	/* Add up distance from axis to annotation */
	beta = 0.5 * (180 + GMT->common.R.wesn[XLO] - GMT->common.R.wesn[XHI]);	/* half-angle between the two gridline directions at the apex at the W and E sides */
	alpha = (gmt_M_is_zero (annot_offset)) ? 90.0 : R2D * atan (0.5 * annot_height / annot_offset);	/* Angle between gridline direction and line from apex to nearest textbox corner */
	if (alpha > beta) return true;	/* Apex angle not large enough to avoid overprinting */
	//if (range > 135.0) return true;	/* Apex angle not large enough to avoid overprinting */
	return false;	/* OK to place the label */
}

GMT_LOCAL void gmtplot_map_tickitem (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, unsigned int item) {
	unsigned int i, nx = 0, ny = 0;
	bool do_x, do_y;
	double dx, dy, *val = NULL, len, shift = 0.0;

	if (! (GMT->current.map.frame.axis[GMT_X].item[item].active || GMT->current.map.frame.axis[GMT_Y].item[item].active)) return;

	dx = (GMT->current.map.frame.axis[GMT_X].file_custom) ? 1.0 : gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_X].type, &GMT->current.map.frame.axis[GMT_X].item[item]);
	dy = (GMT->current.map.frame.axis[GMT_Y].file_custom) ? 1.0 : gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_Y].type, &GMT->current.map.frame.axis[GMT_Y].item[item]);

	if (dx <= 0.0 && dy <= 0.0) return;

	do_x = dx > 0.0 && GMT->current.map.frame.axis[GMT_X].item[item].active && (item == GMT_ANNOT_UPPER ||
		(item == GMT_TICK_UPPER && dx != gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_X].type, &GMT->current.map.frame.axis[GMT_X].item[GMT_ANNOT_UPPER])) ||
		(item == GMT_TICK_LOWER && dx != gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_X].type, &GMT->current.map.frame.axis[GMT_X].item[GMT_ANNOT_LOWER])));
	do_y = dy > 0.0 && GMT->current.map.frame.axis[GMT_Y].item[item].active && (item == GMT_ANNOT_UPPER ||
		(item == GMT_TICK_UPPER && dy != gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_Y].type, &GMT->current.map.frame.axis[GMT_Y].item[GMT_ANNOT_UPPER])) ||
		(item == GMT_TICK_LOWER && dy != gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_Y].type, &GMT->current.map.frame.axis[GMT_Y].item[GMT_ANNOT_LOWER])));
	len = GMT->current.setting.map_tick_length[item];
	if (GMT->current.setting.map_frame_type & GMT_IS_INSIDE) len = -fabs (len);	/* Negative to become inside */

	GMT->current.map.on_border_is_outside = true;	/* Temporarily, points on the border are outside */

	if (do_x) {	/* Draw grid lines that go E to W */
		if (GMT->current.map.frame.axis[GMT_X].file_custom)
			nx = gmtlib_coordinate_array (GMT, w, e, &GMT->current.map.frame.axis[GMT_X].item[item], &val, NULL);
		else
			nx = gmtlib_linear_array (GMT, w, e, dx, GMT->current.map.frame.axis[GMT_X].phase, &val);
		for (i = 0; i < nx; i++)  {
			shift = gmtplot_shift_gridline (GMT, val[i], GMT_X);
			gmtplot_map_lontick (GMT, PSL, val[i] + shift, s, n, len);
		}
		if (nx) gmt_M_free (GMT, val);
	}

	if (do_y) {	/* Draw grid lines that go S to N */
		if (GMT->current.proj.z_down) {
			if (GMT->current.map.frame.axis[GMT_Y].file_custom)
				ny = gmtlib_coordinate_array (GMT, 0.0, GMT->current.proj.z_radius-s, &GMT->current.map.frame.axis[GMT_Y].item[item], &val, NULL);
			else if (GMT->current.proj.z_down == GMT_ZDOWN_Z) /* z = n - r */
				ny = gmtlib_linear_array (GMT, 0.0, GMT->current.proj.z_radius-s, dy, GMT->current.map.frame.axis[GMT_Y].phase, &val);
			else if (GMT->current.proj.z_down == GMT_ZDOWN_ZP) /* r = rp - z */
				ny = gmtlib_linear_array (GMT, GMT->current.proj.z_radius-n, GMT->current.proj.z_radius-s, dy, GMT->current.map.frame.axis[GMT_Y].phase, &val);
			for (i = 0; i < ny; i++) {
				if (GMT->current.proj.z_down == GMT_ZDOWN_ZP)
					val[i] = GMT->current.proj.z_radius - val[i];	/* These are the z values needed for positioning */
				else
					val[i] = GMT->common.R.wesn[YHI] - val[i];	/* These are the radial values needed for positioning */
			}
		}
		else {
			if (GMT->current.map.frame.axis[GMT_Y].file_custom)
				ny = gmtlib_coordinate_array (GMT, s, n, &GMT->current.map.frame.axis[GMT_Y].item[item], &val, NULL);
			else
				ny = gmtlib_linear_array (GMT, s, n, dy, GMT->current.map.frame.axis[GMT_Y].phase, &val);
		}
		for (i = 0; i < ny; i++) {
			if (gmtplot_skip_polar_apex_annotation (GMT, i, val, ny)) continue;
			shift = gmtplot_shift_gridline (GMT, val[i], GMT_Y);
			gmtplot_map_lattick (GMT, PSL, val[i] + shift, w, e, len);
		}
		if (ny) gmt_M_free (GMT, val);
	}

	GMT->current.map.on_border_is_outside = false;	/* Reset back to default */
}

GMT_LOCAL void gmtplot_map_tickmarks (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	/* Tickmarks at annotation interval has already been done except when annotations were not desired */

	if (!(gmt_M_is_geographic (GMT, GMT_IN) || GMT->current.proj.projection_GMT == GMT_POLAR)) return;	/* Tickmarks already done by linear axis */

	PSL_comment (PSL, "Map tickmarks\n");

	gmt_setpen (GMT, &GMT->current.setting.map_tick_pen[GMT_PRIMARY]);
	gmtplot_map_tickitem (GMT, PSL, w, e, s, n, GMT_ANNOT_UPPER);
	if (!(GMT->current.setting.map_frame_type & GMT_IS_FANCY)) {	/* Draw plain boundary and return */
		gmtplot_map_tickitem (GMT, PSL, w, e, s, n, GMT_TICK_UPPER);
		gmt_setpen (GMT, &GMT->current.setting.map_tick_pen[GMT_SECONDARY]);
		gmtplot_map_tickitem (GMT, PSL, w, e, s, n, GMT_TICK_LOWER);
	}

	PSL_setdash (PSL, NULL, 0);
}

bool gmtlib_set_do_seconds (struct GMT_CTRL *GMT, double inc) {
	/* Determines if seconds are to be labelled based on size of increment */
	if (GMT->current.plot.calclock.geo.order[2] == -1) return (false);			/* Seconds not requested by format */
	if (GMT->current.plot.calclock.geo.n_sec_decimals > 0) return (true);			/* If asked for ss.xxx annotations */
	if (fabs (60.0 * fmod (fmod (inc, 1.0) * 60.0, 1.0)) >= 1.0) return (true);	/* Multiples of >= 1 sec intervals */
	return (false);
}

GMT_LOCAL void gmtplot_label_trim (char *label, int stage) {
	/* Used to shorten secondary annotations by eliminating the leading digits. E.g. if the
	 * primary annotation is 30 degrees and the secondary is 30:05, 30:10, etc then we remove
	 * the leading 1 (degrees) or 2 (degrees and minutes) part of the annotation.
	 * This can only work if we are annotating with the ddd:mm:sss[.xx] format as for ddd.xxxxx
	 * we cannot do so.  So if primary is 30:20 and secondary is 30:20:02, 30:20:04 etc we
	 * end upremoving the leading 30:20 (stage = 2).
	 */
	size_t i;
	if (!label) return;	/* No label given */
	if (!stage) return;	/* Not asked to do anything */
	/* Must remove leading stuff (e.g., ddd<degree_sign>) for 2ndary annotations */
	for (i = 0; stage && label[i]; i++)
		if (!strchr ("0123456789-+", label[i])) stage--;
	while (label[i])
		label[stage++] = label[i++];	/* Copy over the later part of the label to the beginning */
	label[stage] = '\0';
}

GMT_LOCAL void gmtplot_radial_annot_setup (struct GMT_CTRL *GMT, double angle, unsigned int *justify, double *text_angle, double *line_angle, double *dc, double *ds) {
	if (angle < 0.0) angle += 360.0;
	if (GMT->current.proj.got_azimuths)
		angle = 90.0 - angle;
	if (GMT->current.proj.projection_GMT == GMT_POLAR)	/* Shift by 90 so it works like N polar */
		gmtlib_polar_prepare_label (GMT, angle-90-GMT->current.proj.p_base_angle, E_SIDE, line_angle, text_angle, justify);
	else
		gmtlib_polar_prepare_label (GMT, angle, E_SIDE, line_angle, text_angle, justify);

	if (GMT->current.proj.projection_GMT != GMT_POLAR && !GMT->current.proj.north_pole) {
		/* Adjust for the fact that for S polar maps, things are a bit backwards */
		*text_angle = -*text_angle;
		*justify = (*justify == PSL_ML) ? PSL_MR : PSL_ML;
		*line_angle = 190.0 - *line_angle;
	}
	if (GMT->current.proj.projection_GMT != GMT_POLAR) {
		*line_angle += 180.0;
		*justify = (*justify == PSL_ML) ? PSL_BR : PSL_TL;
	}
	else
		*justify = (*justify == PSL_ML) ? PSL_TL : PSL_BR;
	*line_angle -= 45;
	/* Shift annotation by the amount implied by GMT_ANNOT_OFFSET */
	sincosd (*line_angle, ds, dc);
	*ds *= MAX (GMT->current.setting.map_annot_offset[GMT_PRIMARY], 0.0);
	*dc *= MAX (GMT->current.setting.map_annot_offset[GMT_PRIMARY], 0.0);

	//fprintf (stderr, "Just = %d Angle = %g Line = %g dx = %g dy = %g\n", *justify, *text_angle, *line_angle, *dc, *ds);
}

GMT_LOCAL void gmtplot_consider_internal_annotations (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	/* Special annotation routine for annotating latitude (or radius) along a specified meridian.  This allows us to annotate
	 * maps that do not have the axes we wish to annotate, such as annotating latitudes for a 0/360 azimuthal map. */
	unsigned int i, nx = 0, ny = 0, first = 0, form, lonlat, justify;
	bool do_minutes, do_seconds, do_grid = false;
	char label[GMT_LEN256] = {""}, format[GMT_LEN64] = {""}, **label_c = NULL;
	double *val = NULL, *tval = NULL, dx, dy, x0, y0, x1, y1, L, ds, dc, line_angle, text_angle, angle, dyg, dxg, xa, xb, ya, yb, Sa, Ca;

	if (GMT->current.map.frame.internal_annot == 0) return;	/* Not requested */

	form = gmt_setfont (GMT, &GMT->current.setting.font_annot[GMT_PRIMARY]);
	PSL_comment (PSL, "Map annotations (Internal)\n");
	PSL_settextmode (PSL, PSL_TXTMODE_MINUS);	/* Replace hyphens with minus signs */

	if (GMT->current.map.frame.internal_annot == 1) {	/* Placement of latitude or radial annotations along selected meridian */
		dy = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_Y].type, &GMT->current.map.frame.axis[GMT_Y].item[GMT_ANNOT_UPPER]);
		if (gmt_M_y_is_lat (GMT, GMT_IN)) {
			do_minutes = (fabs (fmod (dy, 1.0)) > GMT_CONV4_LIMIT);
			do_seconds = gmtlib_set_do_seconds (GMT, dy);
			lonlat = 1;
		}
		else {	/* Also, we know that GMT->current.setting.format_geo_out = -1 in this case */
			do_minutes = do_seconds = 0;
			lonlat = 2;
		}
		if (GMT->current.plot.r_theta_annot) {	/* Make format for radial term */
			char tmp[GMT_LEN64] = {""};
			strcpy (format, GMT->current.setting.format_float_map);
			gmt_get_format (GMT, dy, NULL, NULL, tmp);
			strncpy (GMT->current.setting.format_float_map, tmp, GMT_LEN64-1);
		}

		if (GMT->current.proj.z_down) {	/* Want to annotate depth rather than radius */
			if (GMT->current.map.frame.axis[GMT_Y].file_custom)
				ny = gmtlib_coordinate_array (GMT, 0.0, GMT->current.proj.z_radius-s, &GMT->current.map.frame.axis[GMT_Y].item[GMT_ANNOT_UPPER], &tval, &label_c);
			else if (GMT->current.proj.z_down == GMT_ZDOWN_Z) /* z = n - r */
				ny = gmtlib_linear_array (GMT, 0.0, GMT->current.proj.z_radius-s, dy, GMT->current.map.frame.axis[GMT_Y].phase, &tval);
			else if (GMT->current.proj.z_down == GMT_ZDOWN_ZP) /* z = n - r */
				ny = gmtlib_linear_array (GMT, GMT->current.proj.z_radius-n, GMT->current.proj.z_radius-s, dy, GMT->current.map.frame.axis[GMT_Y].phase, &tval);
			val = gmt_M_memory (GMT, NULL, ny, double);
			for (i = 0; i < ny; i++) {
				if (GMT->current.proj.z_down == GMT_ZDOWN_ZP)
					val[i] = GMT->current.proj.z_radius - tval[i];	/* These are the z values needed for positioning */
				else
					val[i] = GMT->common.R.wesn[YHI] - tval[i];	/* These are the radial values needed for positioning */
			}
		}
		else {				/* Annotate radius */
			if (GMT->current.map.frame.axis[GMT_Y].file_custom)
				ny = gmtlib_coordinate_array (GMT, s, n, &GMT->current.map.frame.axis[GMT_Y].item[GMT_ANNOT_UPPER], &val, &label_c);
			else
				ny = gmtlib_linear_array (GMT, s, n, dy, GMT->current.map.frame.axis[GMT_Y].phase, &val);
			tval = val;	/* Here they are the same thing */
		}

		/* Lot of exceptions to check for that affects if first or last annotatino should be skipped */
		if (GMT->current.proj.z_down == GMT_ZDOWN_Z && ny && tval[ny-1] < GMT->current.proj.z_radius)
			ny--;
		else if (GMT->current.proj.z_down == GMT_ZDOWN_ZP && ny) {
			if ((tval[0] - GMT->current.proj.flip_radius + n) < 0.5*dy)
			first = 1;
		}
		else if (GMT->current.proj.projection_GMT == GMT_POLAR && ny && GMT->current.proj.flip && GMT->current.proj.radial_offset > 0.0)
			ny--;
		else if (GMT->current.proj.projection_GMT == GMT_POLAR && ny && (GMT->current.proj.radial_offset > 0.0 || (!GMT->current.proj.flip && tval[0] > 0.0)))
			first = 1;
		else if (GMT->current.proj.got_elevations && n < 90.0)
			ny--;
		else if (GMT->current.proj.flip && ny && tval[ny-1] < GMT->current.proj.flip_radius)
			ny--;

		gmtplot_radial_annot_setup (GMT, GMT->current.map.frame.internal_arg, &justify, &text_angle, &line_angle, &dc, &ds);
		/* Shall we place grid crosses or have user selected gridlines? */
		if ((dyg = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_Y].type, &GMT->current.map.frame.axis[GMT_Y].item[GMT_GRID_UPPER])) > 0.0) {
			if ((dx = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_X].type, &GMT->current.map.frame.axis[GMT_X].item[GMT_GRID_UPPER])) > 0.0) {
				if (fabs (fmod (GMT->current.map.frame.internal_arg, dx)) > GMT_CONV4_LIMIT)
					do_grid = true;
				else if (dyg > dy)
					do_grid = true;
			}
			else
				do_grid = true;
		}
		else
			do_grid = true;

		if (do_grid) {
			/* Either pick current grid-cross size or use half the annotation size as backup */
			L = 0.5 * fabs ((GMT->current.setting.map_grid_cross_size[GMT_PRIMARY] > 0.0) ? GMT->current.setting.map_grid_cross_size[GMT_PRIMARY] : 0.5 * GMT->current.setting.font_annot[GMT_PRIMARY].size / PSL_POINTS_PER_INCH);
			gmt_setpen (GMT, &GMT->current.setting.map_grid_pen[GMT_PRIMARY]);
		}

		for (i = first; i < ny; i++) {
			if (label_c && label_c[i] && label_c[i][0])
				strncpy (label, label_c[i], GMT_LEN256-1);
			else
				gmtlib_get_annot_label (GMT, tval[i], label, do_minutes, do_seconds, 1, lonlat, GMT->current.map.is_world);
			gmt_geo_to_xy (GMT, GMT->current.map.frame.internal_arg, val[i], &x0, &y0);
			PSL_plottext (PSL, x0+dc, y0+ds, GMT->current.setting.font_annot[GMT_PRIMARY].size, label, text_angle, justify, form);
			if (do_grid) {
				gmt_geo_to_xy (GMT, GMT->current.map.frame.internal_arg, tval[i] - copysign (GMT->current.map.dlat, tval[i]), &x1, &y1);
				angle = d_atan2 (y1-y0, x1-x0);
				sincos (angle, &Sa, &Ca);
				xa = x0 - L * Ca;	xb = x0 + L * Ca;
				ya = y0 - L * Sa;	yb = y0 + L * Sa;
				PSL_plotsegment (PSL, xa, ya, xb, yb);
				angle += M_PI_2;
				sincos (angle, &Sa, &Ca);
				xa = x0 - L * Ca;	xb = x0 + L * Ca;
				ya = y0 - L * Sa;	yb = y0 + L * Sa;
				PSL_plotsegment (PSL, xa, ya, xb, yb);
			}
		}
		if (ny) gmt_M_free (GMT, val);
		if (label_c) {
			for (i = 0; i < ny; i++) gmt_M_str_free (label_c[i]);
			gmt_M_free (GMT, label_c);
		}
		if (GMT->current.proj.z_down) gmt_M_free (GMT, tval);
		if (GMT->current.plot.r_theta_annot)	/* Restore the format */
			strcpy (GMT->current.setting.format_float_map, format);
	}

	if (GMT->current.map.frame.internal_annot == 2) {	/* Placement of longitude annotations along selected parallel */
		dx = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_X].type, &GMT->current.map.frame.axis[GMT_X].item[GMT_ANNOT_UPPER]);
		do_minutes = (fabs (fmod (dx, 1.0)) > GMT_CONV4_LIMIT);
		do_seconds = gmtlib_set_do_seconds (GMT, dx);

		if (GMT->current.map.frame.axis[GMT_X].file_custom)
			nx = gmtlib_coordinate_array (GMT, w, e, &GMT->current.map.frame.axis[GMT_X].item[GMT_ANNOT_UPPER], &val, &label_c);
		else
			nx = gmtlib_linear_array (GMT, w, e, dx, GMT->current.map.frame.axis[GMT_X].phase, &val);
		//if (nx && doubleAlmostEqualZero (val[0], w)) first = 1;
		if (nx && doubleAlmostEqualZero (val[nx-1], e)) nx--;

		if ((dyg = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_X].type, &GMT->current.map.frame.axis[GMT_X].item[GMT_GRID_UPPER])) > 0.0) {
			if ((dxg = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_X].type, &GMT->current.map.frame.axis[GMT_X].item[GMT_GRID_UPPER])) > 0.0) {
				if (fabs (fmod (GMT->current.map.frame.internal_arg, dxg)) > GMT_CONV4_LIMIT)
					do_grid = true;
				else if (dxg > dx)
					do_grid = true;
			}
			else
				do_grid = true;
		}
		else
			do_grid = true;

		if (do_grid) {
			/* Either pick current grid-cross size or use half the annotation size as backup */
			L = 0.5 * fabs ((GMT->current.setting.map_grid_cross_size[GMT_PRIMARY] > 0.0) ? GMT->current.setting.map_grid_cross_size[GMT_PRIMARY] : 0.5 * GMT->current.setting.font_annot[GMT_PRIMARY].size / PSL_POINTS_PER_INCH);
			gmt_setpen (GMT, &GMT->current.setting.map_grid_pen[GMT_PRIMARY]);
		}

		for (i = first; i < nx; i++) {
			if (label_c && label_c[i] && label_c[i][0])
				strncpy (label, label_c[i], GMT_LEN256-1);
			else
				gmtlib_get_annot_label (GMT, val[i], label, do_minutes, do_seconds, 1, 0, GMT->current.map.is_world);
			gmt_geo_to_xy (GMT, val[i], GMT->current.map.frame.internal_arg, &x0, &y0);
			gmt_geo_to_xy (GMT, val[i] + GMT->current.map.dlon, GMT->current.map.frame.internal_arg, &x1, &y1);
			angle = d_atan2 (y1-y0, x1-x0);	/* Angle at longitude annotation */
			sincos (angle+M_PI_4, &dc, &ds);	/* Add 45 to get distance to BL corner of boundingbox for annotation */
			ds *= MAX (GMT->current.setting.map_annot_offset[GMT_PRIMARY], 0.0);
			dc *= MAX (GMT->current.setting.map_annot_offset[GMT_PRIMARY], 0.0);
			PSL_plottext (PSL, x0+dc, y0+ds, GMT->current.setting.font_annot[GMT_PRIMARY].size, label, R2D*angle, PSL_BL, form);
			if (do_grid) {
				sincos (angle, &Sa, &Ca);
				xa = x0 - L * Ca;	xb = x0 + L * Ca;
				ya = y0 - L * Sa;	yb = y0 + L * Sa;
				PSL_plotsegment (PSL, xa, ya, xb, yb);
				angle += M_PI_2;
				sincos (angle, &Sa, &Ca);
				xa = x0 - L * Ca;	xb = x0 + L * Ca;
				ya = y0 - L * Sa;	yb = y0 + L * Sa;
				PSL_plotsegment (PSL, xa, ya, xb, yb);
			}
		}
		if (nx) gmt_M_free (GMT, val);
		if (label_c) {
			for (i = 0; i < nx; i++) gmt_M_str_free (label_c[i]);
			gmt_M_free (GMT, label_c);
		}
	}
	if (val) gmt_M_free (GMT, val);
	PSL_settextmode (PSL, PSL_TXTMODE_HYPHEN);	/* Back to leave as is */
}

GMT_LOCAL void gmtplot_map_label (struct GMT_CTRL *GMT, double x, double y, char *label, double angle, int just, unsigned int axis, bool below) {
	/* Function to use to set axis labels for Cartesian basemaps and colorbars */
	struct PSL_CTRL *PSL= GMT->PSL;
	bool pos_set = (gmt_M_is_zero (x) && gmt_M_is_zero (y));	/* If the current point has already been placed */

	if (gmt_text_is_latex (GMT, label)) {	/* Detected LaTeX commands, i.e., "....@[LaTeX...@[ ..." or  "....<math>LaTeX...</math> ..." */
		bool set_L_off = (axis == GMT_X && !below);	/* May need to ensure extra offset for title */
		double w, h;
		unsigned char *eps = NULL;
		struct imageinfo header;

		if ((eps = gmtplot_latex_eps (GMT, &GMT->current.setting.font_label, label, &header)) == NULL) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "gmtplot_map_label: Conversion of LaTeX label to EPS failed\n");
			return;	/* Done */
		}
		/* Scale up EPS dimensions by the ratio of label font size to LaTeX default size of 10p */
		w = (header.width / 72.0)  * (GMT->current.setting.font_label.size / 10.0);
		h = (header.height / 72.0) * (GMT->current.setting.font_label.size / 10.0);
		/* Place EPS file as label, then free eps */
		GMT_Report (GMT->parent, GMT_MSG_DEBUG, "gmtplot_map_label: Conversion of LaTeX label gave dimensions %g x %g\n", w, h);
		PSL_command (PSL, "V\n");	/* Keep the relative changes inside a save/restore block */
		/* If we plot label below the axis then we must adjust for the fact that the base y-coordinate is based on label font height of "M",
		 * but now we have an EPS image of given (and presumably larger) height.  So we adjust by the difference in those two values */
		if (pos_set) {	/* Translate origin to currentpoint since already set by calling function, and possibly rotate by +/- 90 degrees */
			double sgn[2] = {-1.0, 1.0};
			if (fabs (angle) > 0.0) PSL_command (PSL, "currentpoint T %g R\n", angle); else PSL_command (PSL, "currentpoint T\n");
			if (below) PSL_command (PSL, "0 %d PSL_LH sub neg M currentpoint T\n", (int)lrint (h * PSL->internal.y2iy));
			/* I am somehow missing the label offset for y-axes so I need to account for it here in order to get correct result */
			if (axis == GMT_Y) PSL_command (PSL, "0 %d M currentpoint T\n", (int)lrint (sgn[below]*GMT->current.setting.map_label_offset[GMT_Y] * PSL->internal.y2iy));
		}
		// PSL_command (PSL, "V currentpoint -2000 0 G 4000 0 D S U\n");	/* Debug line for base of label; keep for future debugging */
		PSL_plotlatexeps (PSL, x, y, w, h, PSL_BC, eps, GMT->current.setting.font_label.fill.rgb, &header);	/* Place the EPS plot */
		PSL_command (PSL, "U\n");	/* Close up the block */
		if (set_L_off)	/* Must reset the value of PSL_LH to the EPS image height so the title baseline can be adjusted */
			PSL_command (PSL, "/PSL_LH %d def\n", (int)lrint (h * PSL->internal.y2iy));
		PSL_free (eps);
		return;	/* Done on this end */
	}
	else {	/* Regular text label */
		int form = gmt_setfont (GMT, &GMT->current.setting.font_label);
		double s = (pos_set) ? -1.0 : +1.0;
		PSL_plottext (PSL, x, y, s * GMT->current.setting.font_label.size, label, angle, just, form);
	}
}

bool gmtplot_skip_pole_lat_annotation (struct GMT_CTRL *GMT, double lat) {
	/* Here, latitude is -90 or +90 and the question is do we annotation that latitude? */
	if (GMT->current.proj.projection_GMT == GMT_ECKERT4   && fabs (lat) > 85.0) return true;	/* Slope too gentle near poles to adjust via boost */
	if (GMT->current.proj.projection_GMT == GMT_HAMMER    && fabs (lat) > 85.0) return true;	/* Slope too gentle near poles to adjust via boost */
	if (GMT->current.proj.projection_GMT == GMT_MOLLWEIDE && fabs (lat) > 88.0) return true;	/* Slope too gentle near poles to adjust via boost */
	if (fabs (lat) < 90.0) return false;	/* Not at the poles so OK to annotate */
	if (GMT->current.proj.polar) return true;	/* Cannot since inside map */
	return false;
}

static char *gmtplot_revise_angles_just (struct GMT_CTRL *GMT, double line_angle, unsigned side, unsigned int def_just[], double *t_angle, unsigned int *just) {
	/* Given the side (W or E) and angle of the boundary, determine the justification of the text and possibly flip it 180 */
	static char *flip[2] = {"", "neg "};
	unsigned int k = GMT->current.proj.got_azimuths ? 1 : 0, pside;
	if (side == W_SIDE) side = 0; else side = 1;	/* Turn W_SIDE to 0 and E_SIDE to 1 for use with flip and def_just arrays */
	pside = GMT->current.proj.got_azimuths ? 1-side : side;	/* Since azimuths go the other way */
	if ((line_angle < -90.00 || line_angle > 90.0)) {	/* Must flip the text angle 180 degrees and also flip the justification */
		*t_angle = 180.0;
		*just = def_just[pside];
	}
	else {	/* Normally the annotation is just plotted along the radial which is the x-axis after the rotation */
		*t_angle = 0.0;
		*just = def_just[1-pside];
	}
	if (side == 1) k = (k + 1) % 2;	/* We want positive direction to be outwards */
	// fprintf (stderr, "Side %c Input: line-angle = %g, just = %d  Output: t-angle = %g, just = %d neg = %s\n", name[side], line_angle, def_just[0], *t_angle, *just, flip[k]);
	return flip[k];
}

GMT_LOCAL void gmtplot_map_annotate (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n) {
	unsigned int i, k, nx = 0, ny = 0, last, form, remove[2] = {0,0}, trim, add, x_kind = 0;
	bool do_minutes, do_seconds, done_Greenwich, done_Dateline, check_edges;
	bool full_lat_range, proj_A, proj_B, annot_0_and_360 = false, dual[2], is_dual, annot, is_world_save, lon_wrap_save;
	char label[GMT_LEN256] = {""};
	char **label_c = NULL;
	double *val = NULL, dx[2], dy[2], w2, s2, del, shift = 0.0;

	if (!(gmt_M_x_is_lon (GMT, GMT_IN) || gmt_M_y_is_lat (GMT, GMT_IN) || GMT->current.proj.projection_GMT == GMT_POLAR)) return;	/* Annotations and header already done by gmtplot_linear_map_boundary */

	is_world_save = GMT->current.map.is_world;
	lon_wrap_save = GMT->current.map.lon_wrap;

	if (GMT->current.proj.projection_GMT == GMT_POLAR) {	/* r-theta requires more info to know what theta is */
		if (GMT->current.proj.angle_kind == GMT_IS_LON) x_kind = 0;
		else if (GMT->current.proj.angle_kind == GMT_IS_LAT) x_kind = 3;
		else x_kind = 4;
	}

	if (GMT->current.map.frame.header[0] && !GMT->current.map.frame.plotted_header) {	/* Make plot header for geographic maps */
		if (gmt_M_is_geographic (GMT, GMT_IN) || GMT->current.map.frame.side[N_SIDE] & GMT_AXIS_ANNOT) {
			PSL_setfont (PSL, GMT->current.setting.font_annot[GMT_PRIMARY].id);
			PSL_command (PSL, "/PSL_H_y %d ", PSL_IZ (PSL, GMT->current.setting.map_tick_length[GMT_PRIMARY] + GMT->current.setting.map_annot_offset[GMT_PRIMARY] + GMT->current.setting.map_title_offset));
			if (GMT->current.setting.map_degree_symbol == gmt_none)
				PSL_deftextdim (PSL, "-h", GMT->current.setting.font_annot[GMT_PRIMARY].size, "100");
			else {
				snprintf (label, GMT_LEN16, "100%c",  (int)GMT->current.setting.ps_encoding.code[GMT->current.setting.map_degree_symbol]);
				PSL_deftextdim (PSL, "-h", GMT->current.setting.font_annot[GMT_PRIMARY].size, label);
			}
			PSL_command (PSL, "add def\n");
		}
		else
			PSL_defunits (PSL, "PSL_H_y", GMT->current.setting.map_title_offset + GMT->current.setting.map_tick_length[GMT_PRIMARY]);

		PSL_command (PSL, "%d %d PSL_H_y add M\n", PSL_IZ (PSL, GMT->current.proj.rect[XHI] * 0.5), PSL_IZ (PSL, GMT->current.proj.rect[YHI]));
		gmt_map_title (GMT, 0.0, 0.0);
	}

	gmtplot_consider_internal_annotations (GMT, PSL, w, e, s, n);	/* Handle any special case of internal annotations */

	if (GMT->current.proj.edge[S_SIDE] || GMT->current.proj.edge[N_SIDE]) {
		dx[0] = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_X].type, &GMT->current.map.frame.axis[GMT_X].item[GMT_ANNOT_UPPER]);
		dx[1] = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_X].type, &GMT->current.map.frame.axis[GMT_X].item[GMT_ANNOT_LOWER]);
		/* Determine if we should annotate both 0 and 360 degrees */

		full_lat_range = (fabs (180.0 - fabs (GMT->common.R.wesn[YHI] - GMT->common.R.wesn[YLO])) < GMT_CONV4_LIMIT);
		proj_A = (GMT->current.proj.projection_GMT == GMT_MERCATOR || GMT->current.proj.projection_GMT == GMT_OBLIQUE_MERC ||
			GMT->current.proj.projection_GMT == GMT_WINKEL || GMT->current.proj.projection_GMT == GMT_ECKERT4 || GMT->current.proj.projection_GMT == GMT_ECKERT6 ||
			GMT->current.proj.projection_GMT == GMT_ROBINSON || GMT->current.proj.projection_GMT == GMT_CYL_EQ || GMT->current.proj.projection_GMT == GMT_CYL_STEREO ||
			GMT->current.proj.projection_GMT == GMT_CYL_EQDIST || GMT->current.proj.projection_GMT == GMT_MILLER || GMT->current.proj.projection_GMT == GMT_LINEAR);
		proj_B = (GMT->current.proj.projection_GMT == GMT_HAMMER || GMT->current.proj.projection_GMT == GMT_MOLLWEIDE ||
			GMT->current.proj.projection_GMT == GMT_SINUSOIDAL);
		if (gmt_M_is_conical (GMT) && gmt_M_360_range (w, e))	/* Special case since 360 longitudes do not form a circle but a pacman shape */
			annot_0_and_360 = true;
		else
			annot_0_and_360 = (is_world_save && (proj_A || (!full_lat_range && proj_B)));
	}
	else
		dx[0] = dx[1] = 0.0;
	if (GMT->current.proj.edge[E_SIDE] || GMT->current.proj.edge[W_SIDE]) {
		dy[0] = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_Y].type, &GMT->current.map.frame.axis[GMT_Y].item[GMT_ANNOT_UPPER]);
		dy[1] = gmtlib_get_map_interval (GMT, GMT->current.map.frame.axis[GMT_Y].type, &GMT->current.map.frame.axis[GMT_Y].item[GMT_ANNOT_LOWER]);
	}
	else
		dy[0] = dy[1] = 0.0;

	if (GMT->current.map.frame.axis[GMT_X].file_custom) dx[0] = 1.0;	/* To pass checks below */
	if (GMT->current.map.frame.axis[GMT_Y].file_custom) dy[0] = 1.0;	/* To pass checks below */

	if (dx[0] <= 0.0 && dy[0] <= 0.0) return;

	dual[GMT_X] = (dx[1] > 0.0);
	dual[GMT_Y] = (dy[1] > 0.0);
	is_dual = (dual[GMT_X] | dual[GMT_Y]);
	check_edges = (!GMT->common.R.oblique && (GMT->current.setting.map_frame_type & GMT_IS_INSIDE));

	PSL_comment (PSL, "Map annotations (external)\n");
	PSL_settextmode (PSL, PSL_TXTMODE_MINUS);	/* Replace hyphens with minus signs */

	form = gmt_setfont (GMT, &GMT->current.setting.font_annot[GMT_PRIMARY]);

	GMT->current.map.on_border_is_outside = true;	/* Temporarily, points on the border are outside */
	if (!GMT->common.R.oblique) {
		GMT->current.map.is_world = false;
		if (!(GMT->current.proj.projection_GMT == GMT_GENPER || GMT->current.proj.projection_GMT == GMT_GNOMONIC)) GMT->current.map.lon_wrap = false;
	}

	w2 = (dx[1] > 0.0) ? floor (w / dx[1]) * dx[1] : 0.0;
	s2 = (dy[1] > 0.0) ? floor (s / dy[1]) * dy[1] : 0.0;

	if (dual[GMT_X]) remove[GMT_X] = (dx[0] < (1.0/60.0)) ? 2 : 1;
	if (dual[GMT_Y]) remove[GMT_Y] = (dy[0] < (1.0/60.0)) ? 2 : 1;
	add = (is_dual) ? 1 : 0;
	for (k = 0; k < 1 + add; k++) {
		if (dx[k] > 0.0 && (gmt_M_x_is_lon (GMT, GMT_IN) || GMT->current.proj.projection_GMT == GMT_POLAR)) {	/* Annotate the S and N boundaries */
			done_Greenwich = done_Dateline = false;
			do_minutes = (fabs (fmod (dx[k], 1.0)) > GMT_CONV4_LIMIT);
			do_seconds = gmtlib_set_do_seconds (GMT, dx[k]);

			if (GMT->current.map.frame.axis[GMT_X].file_custom)
				nx = gmtlib_coordinate_array (GMT, w, e, &GMT->current.map.frame.axis[GMT_X].item[GMT_ANNOT_UPPER], &val, &label_c);
			else
				nx = gmtlib_linear_array (GMT, w, e, dx[k], GMT->current.map.frame.axis[GMT_X].phase, &val);
			last = nx - 1;
			for (i = 0; i < nx; i++) {	/* Worry that we do not try to plot 0 and 360 OR -180 and +180 on top of each other */
				if (check_edges && ((i == 0 && val[i] == w) || (i == last && val[i] == e)))
					continue;	/* To avoid/limit clipping of annotations */
				if (gmt_M_is_zero (val[i]))
					done_Greenwich = true;		/* OK, want to plot 0 */
				if (doubleAlmostEqual (val[i], -180.0))
					done_Dateline = true;	/* OK, want to plot -180 */
				/* Only annotate val[i] if
				 *	(1) projection is such that 0/360 or -180/180 are in different x/y locations, OR
				 *	(2) Plot 360 if 0 hasn't been plotted, OR
				 *	(3) plot +180 if -180 hasn't been plotted
				 */
				annot = annot_0_and_360 || !((done_Greenwich && doubleAlmostEqual (val[i], 360.0)) || (done_Dateline && doubleAlmostEqual (val[i], 180.0)));
				trim = 0;
				if (dual[GMT_X] && k == 0) {
					del = fmod (val[i] - w2, dx[1]);
					if (gmt_M_is_zero (del) || doubleAlmostEqual (del, dx[1]))
						annot = false;
					else
						trim = remove[GMT_X];
				}
				if (label_c && label_c[i] && label_c[i][0])
					strncpy (label, label_c[i], GMT_LEN256-1);
				else
					gmtlib_get_annot_label (GMT, val[i], label, do_minutes, do_seconds, !trim, x_kind, is_world_save);
				gmtplot_label_trim (label, trim);
				shift = gmtplot_shift_gridline (GMT, val[i], GMT_X);
				gmtplot_map_symbol_ns (GMT, PSL, val[i]+shift, label, s, n, annot, k, form);
			}
			if (nx) gmt_M_free (GMT, val);
			if (label_c) {
				for (i = 0; i < nx; i++) gmt_M_str_free (label_c[i]);
				gmt_M_free (GMT, label_c);
			}
		}

		if (dy[k] > 0.0 && (gmt_M_y_is_lat (GMT, GMT_IN) || GMT->current.proj.projection_GMT == GMT_POLAR)) {	/* Annotate W and E boundaries */
			unsigned int lonlat;
			double *tval = NULL;
			char format[GMT_LEN64] = {""};

			if (gmt_M_y_is_lat (GMT, GMT_IN)) {
				do_minutes = (fabs (fmod (dy[k], 1.0)) > GMT_CONV4_LIMIT);
				do_seconds = gmtlib_set_do_seconds (GMT, dy[k]);
				lonlat = 1;
			}
			else {	/* Also, we know that GMT->current.setting.format_geo_out = -1 in this case */
				do_minutes = do_seconds = 0;
				lonlat = 2;
			}
			if (GMT->current.plot.r_theta_annot) {	/* Make format for radial term */
				char tmp[GMT_LEN64] = {""};
				strcpy (format, GMT->current.setting.format_float_map);
				gmt_get_format (GMT, dy[k], NULL, NULL, tmp);
				strncpy (GMT->current.setting.format_float_map, tmp, GMT_LEN64-1);
			}

			if (GMT->current.proj.z_down) {	/* Want to annotate depth rather than radius */
				if (GMT->current.map.frame.axis[GMT_Y].file_custom)
					ny = gmtlib_coordinate_array (GMT, 0.0, GMT->current.proj.z_radius-s, &GMT->current.map.frame.axis[GMT_Y].item[GMT_ANNOT_UPPER], &tval, &label_c);
				else if (GMT->current.proj.z_down == GMT_ZDOWN_Z) /* z = n - r */
					ny = gmtlib_linear_array (GMT, 0.0, GMT->current.proj.z_radius-s, dy[k], GMT->current.map.frame.axis[GMT_Y].phase, &tval);
				else if (GMT->current.proj.z_down == GMT_ZDOWN_ZP) /* z = n - r */
					ny = gmtlib_linear_array (GMT, GMT->current.proj.z_radius-n, GMT->current.proj.z_radius-s, dy[k], GMT->current.map.frame.axis[GMT_Y].phase, &tval);
				val = gmt_M_memory (GMT, NULL, ny, double);
				for (i = 0; i < ny; i++) {
					if (GMT->current.proj.z_down == GMT_ZDOWN_ZP)
						val[i] = GMT->current.proj.z_radius - tval[i];	/* These are the z values needed for positioning */
					else
						val[i] = GMT->common.R.wesn[YHI] - tval[i];	/* These are the radial values needed for positioning */
				}
			}
			else {				/* Annotate radius */
				if (GMT->current.map.frame.axis[GMT_Y].file_custom)
					ny = gmtlib_coordinate_array (GMT, s, n, &GMT->current.map.frame.axis[GMT_Y].item[GMT_ANNOT_UPPER], &val, &label_c);
				else
					ny = gmtlib_linear_array (GMT, s, n, dy[k], GMT->current.map.frame.axis[GMT_Y].phase, &val);
				tval = val;	/* Here they are the same thing */
			}
			PSL_command (PSL, "/PSL_AH1 0\n");	/* Initiate annotation width  */
			last = ny - 1;
			for (i = 0; i < ny; i++) {
				if (gmtplot_skip_pole_lat_annotation (GMT, val[i]))
					continue;	/* Cannot place S or N 90 degree annotation for this projection */
				if (gmtplot_skip_polar_apex_annotation (GMT, i, val, ny)) continue;
				annot = true, trim = 0;
				if (check_edges && ((i == 0 && val[i] == s) || (i == last && val[i] == n)))
					continue;	/* To avoid/limit clipping of annotations */

				if (dual[GMT_Y] && k == 0) {
					del = fmod (val[i] - s2, dy[1]);
					if (gmt_M_is_zero (del) || doubleAlmostEqual (del, dy[1]))
						annot = false;
					else
						trim = remove[GMT_Y];
				}
				if (label_c && label_c[i] && label_c[i][0])
					strncpy (label, label_c[i], GMT_LEN256-1);
				else
					gmtlib_get_annot_label (GMT, tval[i], label, do_minutes, do_seconds, !trim, lonlat, is_world_save);
				gmtplot_label_trim (label, trim);
				shift = gmtplot_shift_gridline (GMT, val[i], GMT_Y);
				gmtplot_map_symbol_ew (GMT, PSL, val[i]+shift, label, w, e, annot, k, form);
				PSL_deftextdim (PSL, "-w", GMT->current.setting.font_annot[GMT_PRIMARY].size, label);
				PSL_command (PSL, "mx\n");		/* Update the longest annotation stored in PSL_AH1 */
			}
			PSL_command (PSL, "def\n");	/* Finalize the definition of longest annotation found */
			if (ny) gmt_M_free (GMT, val);
			if (label_c) {
				for (i = 0; i < ny; i++) gmt_M_str_free (label_c[i]);
				gmt_M_free (GMT, label_c);
			}
			if (GMT->current.proj.z_down) gmt_M_free (GMT, tval);
			if (GMT->current.plot.r_theta_annot)	/* Restore the format */
				strcpy (GMT->current.setting.format_float_map, format);
		}
	}

	if (GMT->current.proj.projection_GMT == GMT_POLAR && GMT->current.map.frame.axis[GMT_Y].label[0] && !gmt_M_360_range (GMT->common.R.wesn[XLO], GMT->common.R.wesn[XHI])) {
		/* May need to label radial axis */
		struct GMT_PLOT_AXIS *A = &GMT->current.map.frame.axis[GMT_Y];
		double x0, x1, y0, y1, xm, ym, line_angle, t_angle;
		double off = GMT->current.setting.map_label_offset[GMT_Y] + MAX (0.0, GMT->current.setting.map_annot_offset[GMT_PRIMARY]) + MAX (0.0, GMT->current.setting.map_tick_length[GMT_PRIMARY]);
		unsigned int just, we_side[2] = {W_SIDE, E_SIDE}, def_just[2] = {PSL_TC, PSL_BC};
		char *neg;
		if (GMT->current.map.frame.side[we_side[0]] & GMT_AXIS_ANNOT) {
			/* Compute x, y, angles etc */
			gmt_geo_to_xy (GMT, GMT->common.R.wesn[XHI], GMT->common.R.wesn[YLO], &x0, &y0);
			gmt_geo_to_xy (GMT, GMT->common.R.wesn[XHI], GMT->common.R.wesn[YHI], &x1, &y1);
			xm = 0.5 * (x0 + x1);	ym = 0.5 * (y0 + y1);	line_angle = d_atan2d (y1 - y0, x1 - x0);
			neg = gmtplot_revise_angles_just (GMT, line_angle, we_side[0], def_just, &t_angle, &just);
			/* Center origin on mid-axis and rotate plot x-axis to coincide with line_angle direction */
			PSL_setorigin (PSL, xm, ym, line_angle, PSL_FWD);
			PSL_setcurrentpoint (PSL, 0.0, 0.0);	/* Make the center point of axis the current point */
			PSL_command (PSL, "0 PSL_AH1 %d add %sG\n", PSL_IZ (PSL, off), neg);	/* Move outwards in the x-direction by the annotation width and offset */
			gmtplot_map_label (GMT, 0.0, 0.0, A->label, t_angle, just, GMT_Y, true);	/* Place the label */
			PSL_setorigin (PSL, -xm, -ym, -line_angle, PSL_INV);	/* Undo the coordinate system change */
		}
		if (GMT->current.map.frame.side[we_side[1]] & GMT_AXIS_ANNOT) {
			/* Compute x, y, angles etc */
			gmt_geo_to_xy (GMT, GMT->common.R.wesn[XLO], GMT->common.R.wesn[YLO], &x0, &y0);
			gmt_geo_to_xy (GMT, GMT->common.R.wesn[XLO], GMT->common.R.wesn[YHI], &x1, &y1);
			xm = 0.5 * (x0 + x1);	ym = 0.5 * (y0 + y1);	line_angle = d_atan2d (y1 - y0, x1 - x0);
			neg = gmtplot_revise_angles_just (GMT, line_angle, we_side[1], def_just,  &t_angle, &just);
			/* Center origin on mid-axis and rotate plot x-axis to coincide with line_angle direction */
			PSL_setorigin (PSL, xm, ym, line_angle, PSL_FWD);
			PSL_setcurrentpoint (PSL, 0.0, 0.0);	/* Make the center point of axis the current point */
			PSL_command (PSL, "0 PSL_AH1 %d add %sG\n", PSL_IZ (PSL, off), neg);	/* Move outwards in the x-direction by the annotation width and offset */
			gmtplot_map_label (GMT, 0.0, 0.0, A->label, t_angle, just, GMT_Y, false);
			PSL_setorigin (PSL, -xm, -ym, -line_angle, PSL_INV);	/* Undo the coordinate system change */
		}
	}

	GMT->current.map.on_border_is_outside = false;	/* Reset back to default */
	GMT->current.map.is_world = is_world_save;
	GMT->current.map.lon_wrap = lon_wrap_save;
	PSL_settextmode (PSL, PSL_TXTMODE_HYPHEN);	/* Back to leave as is */
}

GMT_LOCAL void gmtplot_map_boundary (struct GMT_CTRL *GMT) {
	double w, e, s, n;
	struct PSL_CTRL *PSL= GMT->PSL;

	if (!GMT->current.map.frame.draw && GMT->current.proj.projection_GMT != GMT_LINEAR) return;	/* We have a separate check in linear_map_boundary */
	if (GMT->current.map.frame.no_frame) return;	/* Specifically did not want frame */

	if (GMT->current.map.frame.order == GMT_BASEMAP_BEFORE && GMT->current.map.frame.basemap_flag & GMT_BASEMAP_FRAME_AFTER) return;	/* Wrong order */
	if (GMT->current.map.frame.order == GMT_BASEMAP_AFTER  && !(GMT->current.map.frame.basemap_flag & GMT_BASEMAP_FRAME_AFTER)) return;	/* Wrong order */

	w = GMT->common.R.wesn[XLO], e = GMT->common.R.wesn[XHI], s = GMT->common.R.wesn[YLO], n = GMT->common.R.wesn[YHI];

	PSL_comment (PSL, "Start of map frame\n");

	gmt_setpen (GMT, &GMT->current.setting.map_frame_pen);
	PSL_setcolor (PSL, GMT->current.setting.map_frame_pen.rgb, PSL_IS_STROKE);

	switch (GMT->current.proj.projection_GMT) {
		case GMT_LINEAR:
			if (gmt_M_is_geographic (GMT, GMT_IN))	/* xy is lonlat */
				gmtplot_fancy_map_boundary (GMT, PSL, w, e, s, n);
			else
				gmtplot_linear_map_boundary (GMT, PSL, w, e, s, n);
			break;
		case GMT_POLAR:
			gmtplot_theta_r_map_boundary (GMT, PSL, w, e, s, n);
			break;
		case GMT_MERCATOR:
		case GMT_CYL_EQ:
		case GMT_CYL_EQDIST:
		case GMT_CYL_STEREO:
		case GMT_MILLER:
			gmtplot_fancy_map_boundary (GMT, PSL, w, e, s, n);
			break;
		case GMT_ALBERS:
		case GMT_ECONIC:
		case GMT_LAMBERT:
		case GMT_POLYCONIC:
			gmtplot_conic_map_boundary (GMT, PSL, w, e, s, n);
			break;
		case GMT_OBLIQUE_MERC:
			gmtplot_oblmrc_map_boundary (GMT, PSL, w, e, s, n);
			break;
		case GMT_GENPER:
			gmtplot_genper_map_boundary (GMT, PSL, w, e, s, n);
			break;
		case GMT_STEREO:
		case GMT_ORTHO:
		case GMT_LAMB_AZ_EQ:
		case GMT_AZ_EQDIST:
		case GMT_GNOMONIC:
			if (GMT->current.proj.polar)
				gmtplot_polar_map_boundary (GMT, PSL, w, e, s, n);
			else
				gmtplot_circle_map_boundary (GMT, PSL, w, e, s, n);
			break;
		case GMT_HAMMER:
		case GMT_MOLLWEIDE:
		case GMT_SINUSOIDAL:
			gmtplot_ellipse_map_boundary (GMT, PSL, w, e, s, n);
			break;
		case GMT_TM:
		case GMT_UTM:
		case GMT_CASSINI:
		case GMT_WINKEL:
		case GMT_ECKERT4:
		case GMT_ECKERT6:
		case GMT_ROBINSON:
		case GMT_VANGRINTEN:
			gmtplot_basic_map_boundary (GMT, PSL, w, e, s, n);
			break;
	}
	PSL_comment (PSL, "End of map frame\n");
}

/* gmt_map_basemap will create a basemap for the given area.
 * Scaling and wesn are assumed to be passed through the GMT->current.proj-structure (see GMT_project.h)
 * Tickmark info are passed through the GMT->current.map.frame-structure
 */

GMT_LOCAL bool gmtplot_is_fancy_boundary (struct GMT_CTRL *GMT) {
	switch (GMT->current.proj.projection_GMT) {
		case GMT_LINEAR:
			return (gmt_M_is_geographic (GMT, GMT_IN));
			break;
		case GMT_MERCATOR:
		case GMT_CYL_EQ:
		case GMT_CYL_EQDIST:
		case GMT_CYL_STEREO:
		case GMT_MILLER:
			return (true);
			break;
		case GMT_ALBERS:
		case GMT_ECONIC:
		case GMT_LAMBERT:
			return (!GMT->common.R.oblique);
			break;
		case GMT_STEREO:
		case GMT_ORTHO:
		case GMT_GENPER:
		case GMT_LAMB_AZ_EQ:
		case GMT_AZ_EQDIST:
		case GMT_GNOMONIC:
		case GMT_VANGRINTEN:
			return (GMT->current.proj.polar && !GMT->common.R.oblique);
			break;
		case GMT_POLAR:
		case GMT_OBLIQUE_MERC:
		case GMT_HAMMER:
		case GMT_MOLLWEIDE:
		case GMT_SINUSOIDAL:
		case GMT_TM:
		case GMT_UTM:
		case GMT_CASSINI:
		case GMT_WINKEL:
		case GMT_ECKERT4:
		case GMT_ECKERT6:
		case GMT_ROBINSON:
			return (false);
			break;
		default:
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Failure in gmtplot_is_fancy_boundary - notify developers\n");
			return (false);
	}
}

GMT_LOCAL void gmtplot_vertical_wall (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, int quadrant, double *nesw, bool back, unsigned int mode3d) {
	/* Draws the two vertical walls in the back of the 3-D plot */
	int plane = (quadrant + 1) % 2;	/* Black magic that gives us GMT_X (for yz plane) or GMT_Y (for xz plane) */
	double xx[4], yy[4];
	gmt_plane_perspective (GMT, plane, nesw[quadrant % 4]);
	/* Assign the coordinates of the four corners of the wall (hence assumed to be rectangular) */
	xx[0] = xx[1] = nesw[(quadrant+1)%4];	xx[2] = xx[3] = nesw[(quadrant+3)%4];
	yy[0] = yy[3] = GMT->current.proj.zmin;	yy[1] = yy[2] = GMT->current.proj.zmax;
	if (mode3d & GMT_3D_WALL && back) {	/* Called for one of the back walls and we wish to draw it */
		int wplane = 1 - plane;	/* Get integer for the correct plain side */
		if (GMT->current.map.frame.paint[wplane]) {	/* First paint the back wall, with no outline drawn */
			PSL_setfill (PSL, GMT->current.map.frame.fill[wplane].rgb, 0);
			PSL_plotbox (PSL, nesw[(quadrant+1)%4], GMT->current.proj.zmin, nesw[(quadrant+3)%4], GMT->current.proj.zmax);
		}
		if (back)	/* Gridlines was requested, so draw those now */
			gmtplot_z_gridlines (GMT, PSL, GMT->common.R.wesn[ZLO], GMT->common.R.wesn[ZHI], plane, mode3d, quadrant);

		if (GMT->current.map.frame.draw_wall) {	/* We wanted the outline of the back-wall drawn, do now to overwrite any perimeter gridlines */
			gmt_setpen (GMT, &GMT->current.map.frame.pen);
			PSL_plotline (PSL, xx, yy, 4, PSL_MOVE|PSL_STROKE);
		}
	}
	else if (back)	/* Get here if no back-walls were filled our outlined */
		gmtplot_z_gridlines (GMT, PSL, GMT->common.R.wesn[ZLO], GMT->common.R.wesn[ZHI], plane, mode3d, quadrant);
}

GMT_LOCAL void gmtplot_cube_box (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, int quadrant, double *nesw) {
	/* Draws the 3-D cube lines */
	int plane = (quadrant + 1) % 2;	/* Black magic that gives us GMT_X (for yz plane) or GMT_Y (for xz plane) */
	double xx[4], yy[4];
	gmt_plane_perspective (GMT, plane, nesw[quadrant % 4]);
	/* Assign the coordinates of the four corners of the wall (hence assumed to be rectangular) */
	xx[0] = xx[1] = nesw[(quadrant+1)%4];	xx[2] = xx[3] = nesw[(quadrant+3)%4];
	yy[0] = yy[3] = GMT->current.proj.zmin;	yy[1] = yy[2] = GMT->current.proj.zmax;
	gmt_setpen (GMT, &GMT->current.map.frame.pen);
	PSL_plotline (PSL, xx, yy, 4, PSL_MOVE|PSL_STROKE);
}

GMT_LOCAL void gmtplot_timestamp (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double x, double y, unsigned int justify, char *U_label) {
	/* x, y = location of the time stamp box
	 * justify indicates the corner of the box that (x,y) refers to, see below
	 * U_label = label to be plotted to the right of the box
	 *
	 *   9        10       11
	 *   |------------------|
	 *   5         6        7
	 *   |------------------|
	 *   1         2        3
	 */

	time_t right_now;
	char label[GMT_LEN512] = {""}, text[GMT_LEN256] = {""};
	double dim[3] = {0.365, 0.15, 0.032};	/* Predefined dimensions in inches */
	double unset_rgb[4] = {-1.0, -1.0, -1.0, 0.0};

	/* Plot time string in format defined by format_time_stamp */

	right_now = time ((time_t *)0);
	strftime (text, sizeof(text), GMT->current.setting.format_time_stamp, localtime (&right_now));
	snprintf (label, GMT_LEN256, "  %s  ", text);

	PSL_command (PSL, "%% Begin GMT time-stamp\nV\n");
	PSL_setorigin (PSL, x, y, 0.0, PSL_FWD);
	PSL_setlinewidth (PSL, 0.25);
	PSL_setfont (PSL, GMT->current.setting.font_logo.id);
	PSL_defunits (PSL, "PSL_g_w", dim[0]);	/* Size of the black [GMT] box */
	PSL_defunits (PSL, "PSL_g_h", dim[1]);
	PSL_deftextdim (PSL, "PSL_b", 8.0, label);	/* Size of the white [timestamp] box (use only length) */

	/* When justification is not BL (justify == 1), add some PostScript code to move to the
		location where the lower left corner of the time stamp box is to be drawn */

	switch ((justify + 3) % 4) {
		case 1:	/* Center */
			PSL_command (PSL, "PSL_g_w PSL_b_w add 2 div neg 0 T\n"); break;
		case 2:	/* Right justify */
			PSL_command (PSL, "PSL_g_w PSL_b_w add neg 0 T\n"); break;
	}
	switch (justify / 4) {
		case 1: /* Middle */
			PSL_command (PSL, "0 PSL_g_h 2 div neg T\n"); break;
		case 2: /* Top justify */
			PSL_command (PSL, "0 PSL_g_h neg T\n"); break;
	}

	/* Now draw black box with GMT logo, and white box with time stamp */

	PSL_setfill (PSL, GMT->current.setting.map_default_pen.rgb, 1);
	PSL_plotsymbol (PSL, 0.5*dim[0], 0.5*dim[1], dim, PSL_RECT);
	PSL_plotcolorimage (PSL, 0.0, 0.0, dim[0], dim[1], PSL_BL, GMT_glyph, 220, 90, 1);
	PSL_setfill (PSL, GMT->PSL->init.page_rgb, 1);
	PSL_command (PSL, "PSL_g_h PSL_b_w PSL_g_w 0 Sb\n");
	PSL_plottext (PSL, dim[0], dim[2], 8.0, label, 0.0, PSL_BL, 0);

	/* Optionally, add additional label to the right of the box */

	if (U_label && U_label[0]) {
		snprintf (label, GMT_LEN512, "   %s", U_label);
		PSL_plottext (PSL, 0.0, 0.0, -7.0, label, 0.0, PSL_BL, 0);
	}

	PSL_command (PSL, "U\n%% End GMT time-stamp\n");

	/* Reset fill style to empty and no outline and reset linewidth */
	PSL_setfill (PSL, unset_rgb, 0);
	PSL->current.linewidth = GMT_NOTSET;
}

GMT_LOCAL void gmtplot_echo_command (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, struct GMT_OPTION *options) {
	/* This routine will echo the command and its arguments to the
	 * PostScript output file so that the user can see what scales
	 * etc was used to produce this plot.  Any options with arguments
	 * containing spaces will be enclosed in single quotes.
	 */
	size_t length = 0;
	char outstring[GMT_LEN1024] = {""}, not_used[GMT_LEN32] = {""};
	struct GMT_OPTION *opt = NULL;

	if (GMT->current.setting.run_mode == GMT_MODERN)
		PSL_command (PSL, "\n%% PostScript produced by:\n%%@GMT: gmt %s", gmt_current_name (GMT->init.module_name, not_used));
	else
		PSL_command (PSL, "\n%% PostScript produced by:\n%%@GMT: gmt %s", GMT->init.module_name);
	for (opt = options; opt; opt = opt->next) {
		if (length >= GMT_LEN512) {
			PSL_command (PSL, "%s \\\n%%@GMT:+", outstring);
			length = 0;
			gmt_M_memset (outstring, GMT_LEN1024, char);
		}
		strcat (outstring, " ");	length++;
		if (!(opt->option == GMT_OPT_INFILE || opt->option == GMT_OPT_OUTFILE)) {
			if (strchr (opt->arg, ' ')) outstring[length++] = '\'';
			outstring[length++] = '-';
			outstring[length++] = opt->option;
		}
		if ((strlen (opt->arg) + length) < GMT_LEN1024) strcat (outstring, opt->arg);
		length += strlen (opt->arg);
		if (strchr (opt->arg, ' ')) outstring[length++] = '\'';
	}
	PSL_command (PSL, "%s\n", outstring);
}

GMT_LOCAL void gmtplot_NaN_pen_up (double x[], double y[], unsigned int pen[], uint64_t n) {
	/* Ensure that if there are NaNs we set pen = PSL_MOVE */

	uint64_t i, n1;

	for (i = 0, n1 = n - 1; i < n; i++) {
		if (gmt_M_is_dnan (x[i]) || gmt_M_is_dnan (y[i])) {
			pen[i] = PSL_MOVE;
			if (i < n1) pen[i+1] = PSL_MOVE;	/* Since the next point must become the new anchor */
		}
	}
}

GMT_LOCAL void gmtplot_northstar (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double x0, double y0, double r) {
	/* Draw a fancy 5-pointed North star */
	unsigned int a;
	double r2, x[4], y[4], dir, dir2, s, c;

	r2 = r * 0.3;
	for (a = 0; a <= 360; a += 72) {	/* Azimuth of the 5 points on the star */
		/* Solid half */
		x[0] = x[3] = x0, y[0] = y[3] = y0;
		dir = 90.0 - (double)a;
		sincosd (dir, &s, &c);
		x[1] = x0 + r * c;
		y[1] = y0 + r * s;
		dir2 = dir - 36.0;
		sincosd (dir2, &s, &c);
		x[2] = x0 + r2 * c;
		y[2] = y0 + r2 * s;
		PSL_setfill (PSL, GMT->current.setting.map_default_pen.rgb, 1);
		PSL_plotpolygon (PSL, x, y, 4);
		/* Hollow half */
		x[0] = x[3] = x0, y[0] = y[3] = y0;
		sincosd (dir, &s, &c);
		x[1] = x0 + r * c;
		y[1] = y0 + r * s;
		dir2 = dir + 36.0;
		sincosd (dir2, &s, &c);
		x[2] = x0 + r2 * c;
		y[2] = y0 + r2 * s;
		PSL_setfill (PSL, GMT->PSL->init.page_rgb, 1);
		PSL_plotpolygon (PSL, x, y, 4);
	}
}

#define M_VW	0.005
#define M_HL	0.075
#define M_HW	0.025

#define DIST_TO_2ND_POINT 1.0

#define F_VW	0.01
#define F_HL	0.15
#define F_HW	0.05

GMT_LOCAL void gmtplot_draw_mag_rose (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, struct GMT_MAP_ROSE *mr) {
	/* Magnetic compass rose */
	unsigned int i, k, level, just, ljust[4] = {PSL_TC, PSL_ML, PSL_BC, PSL_MR}, n_tick = 0, form;
	double ew_angle, angle, R[2], tlen[3], L, s, c, lon, lat, x[5], y[5], xp[5], yp[5];
	double offset, t_angle, scale[2], base, v_angle, *val = NULL, dim[PSL_MAX_DIMS];
	char label[GMT_LEN16], *type[2] = {"inner", "outer"};
	struct GMT_FILL f;

	gmt_xy_to_geo (GMT, &lon, &lat, mr->refpoint->x, mr->refpoint->y);

	/* Initialize fill structure */
	gmt_init_fill (GMT, &f, GMT->current.setting.color_patch[GMT_BGD][0], GMT->current.setting.color_patch[GMT_BGD][1], GMT->current.setting.color_patch[GMT_BGD][2]);
	ew_angle = gmt_azim_to_angle (GMT, lon, lat, DIST_TO_2ND_POINT, 90.0);	/* Get angle of E-W direction at this location */

	R[GMT_ROSE_PRIMARY] = 0.75 * 0.5 * mr->size;
	R[GMT_ROSE_SECONDARY] = 0.5 * mr->size;
	tlen[0] = GMT->current.setting.map_tick_length[GMT_TICK_UPPER];
	tlen[1] = GMT->current.setting.map_tick_length[GMT_ANNOT_UPPER];
	tlen[2] = 1.5 * GMT->current.setting.map_tick_length[GMT_ANNOT_UPPER];
	scale[GMT_ROSE_PRIMARY] = 0.85;
	scale[GMT_ROSE_SECONDARY] = 1.0;
	GMT->current.plot.r_theta_annot = false;	/* Just in case it was turned on in gmt_map.c */

	PSL_settextmode (PSL, PSL_TXTMODE_MINUS);	/* Replace hyphens with minus signs */

	for (level = 0; level < 2; level++) {	/* Inner (0) and outer (1) angles */
		if (level == GMT_ROSE_PRIMARY && mr->kind != 2) continue;	/* Sorry, not magnetic directions */
		if (mr->draw_circle[level]) {
			gmt_setfill (GMT, NULL, 1);
			PSL_comment (PSL, "Draw magnetic rose %s circle\n", type[level]);
			gmt_setpen (GMT, &mr->pen[level]);
			s = 2.0 * R[level];
			PSL_plotsymbol (PSL, mr->refpoint->x, mr->refpoint->y, &s, PSL_CIRCLE);
		}
		offset = (level == GMT_ROSE_PRIMARY) ? mr->declination : 0.0;
		gmt_setpen (GMT, &GMT->current.setting.map_tick_pen[level]);
		n_tick = gmtlib_linear_array (GMT, 0.0, 360.0, mr->g_int[level], 0.0, &val);
		PSL_comment (PSL, "Draw %d tickmarks for magnetic rose %s circle\n", n_tick, type[level]);
		for (i = 0; i < n_tick - 1; i++) {	/* Increments of fine tickmarks (-1 to avoid repeating 360) */
			angle = 90.0 - (offset + val[i]);	/* Since val is azimuth */
			k = (gmt_M_is_zero (fmod (val[i], mr->a_int[level]))) ? 2 : ((gmt_M_is_zero (fmod (val[i], mr->f_int[level]))) ? 1 : 0);
			sincosd (ew_angle + angle, &s, &c);
			x[0] = mr->refpoint->x + R[level] * c, y[0] = mr->refpoint->y + R[level] * s;
			x[1] = mr->refpoint->x + (R[level] - scale[level]*tlen[k]) * c, y[1] = mr->refpoint->y + (R[level] - scale[level]*tlen[k]) * s;
			PSL_plotsegment (PSL, x[0], y[0], x[1], y[1]);
		}
		gmt_M_free (GMT, val);

		form = gmt_setfont (GMT, &GMT->current.setting.font_annot[level]);
		n_tick = gmtlib_linear_array (GMT, 0.0, 360.0, mr->a_int[level], 0.0, &val);
		PSL_comment (PSL, "Draw %d tickmarks and annotations for magnetic rose %s circle\n", n_tick, type[level]);
		for (i = 0; i < n_tick - 1; i++) {	/* Increments of annotations (-1 to avoid repeating 360) */
			angle = 90.0 - (offset + val[i]);	/* Since val is azimuth */
			sincosd (ew_angle + angle, &s, &c);
			x[0] = mr->refpoint->x + (R[level] + GMT->current.setting.map_annot_offset[level]) * c, y[0] =
					 mr->refpoint->y + (R[level] + GMT->current.setting.map_annot_offset[level]) * s;
			if (GMT->current.setting.map_degree_symbol == gmt_none)
				snprintf (label, GMT_LEN16, "%ld", lrint (val[i]));
			else
				snprintf (label, GMT_LEN16, "%ld%c", lrint (val[i]), (int)GMT->current.setting.ps_encoding.code[GMT->current.setting.map_degree_symbol]);
			t_angle = fmod ((double)(ew_angle - val[i] - offset) + 360.0, 360.0);	/* Now in 0-360 range */
			if (t_angle > 180.0) t_angle -= 180.0;	/* Now in -180/180 range */
			if (t_angle > 90.0 || t_angle < -90.0) t_angle -= copysign (180.0, t_angle);
			just = (y[0] <= mr->refpoint->y) ? PSL_TC : PSL_BC;
			v_angle = val[i] - ew_angle;
			if (level == GMT_ROSE_SECONDARY && doubleAlmostEqual (v_angle, 90.0))
				t_angle = -90.0, just = PSL_BC;
			if (level == GMT_ROSE_SECONDARY && doubleAlmostEqual (v_angle, 270.0))
				t_angle = 90.0, just = PSL_BC;
			PSL_plottext (PSL, x[0], y[0], GMT->current.setting.font_annot[level].size, label, t_angle, just, form);
		}
		gmt_M_free (GMT, val);
	}

	/* Draw extra tick for the 4 main compass directions */
	gmt_setpen (GMT, &GMT->current.setting.map_tick_pen[GMT_SECONDARY]);
	base = R[GMT_ROSE_SECONDARY] + GMT->current.setting.map_annot_offset[GMT_SECONDARY] + GMT->current.setting.font_annot[GMT_ROSE_SECONDARY].size / PSL_POINTS_PER_INCH;
	PSL_comment (PSL, "Draw 4 tickmarks and annotations for cardinal directions of magnetic rose\n", n_tick);
	for (i = 0, k = 1; i < 360; i += 90, k++) {	/* 90-degree increments of tickmarks */
		angle = (double)i;
		sincosd (ew_angle + angle, &s, &c);
		x[0] = mr->refpoint->x + R[GMT_ROSE_SECONDARY] * c, y[0] = mr->refpoint->y + R[GMT_ROSE_SECONDARY] * s;
		x[1] = mr->refpoint->x + (R[GMT_ROSE_SECONDARY] + tlen[0]) * c, y[1] = mr->refpoint->y + (R[GMT_ROSE_SECONDARY] + tlen[0]) * s;
		PSL_plotsegment (PSL, x[0], y[0], x[1], y[1]);
		if (k == 4) k = 0;
		if (!mr->label[k][0]) continue;	/* No label desired */
		x[0] = mr->refpoint->x + base * c, y[0] = mr->refpoint->y + base * s;
		x[1] = mr->refpoint->x + (base + 2.0 * tlen[2]) * c, y[1] = mr->refpoint->y + (base + 2.0 * tlen[2]) * s;
		PSL_plotsegment (PSL, x[0], y[0], x[1], y[1]);
		if (k == 2 && mr->label[2][0] == '*') {	/* Wanted '*' instead of N */
			x[0] = mr->refpoint->x + (base + 2.0*tlen[2] + GMT->current.setting.map_title_offset + 0.025*mr->size) * c, y[0] = mr->refpoint->y + (base + 2.0*tlen[2] + GMT->current.setting.map_title_offset + 0.025*mr->size) * s;
			gmtplot_northstar (GMT, PSL, x[0], y[0], 0.1*mr->size);
		}
		else {
			x[0] = mr->refpoint->x + (base + 2.0*tlen[2] + GMT->current.setting.map_title_offset) * c, y[0] = mr->refpoint->y + (base + 2.0*tlen[2] + GMT->current.setting.map_title_offset) * s;
			form = gmt_setfont (GMT, &GMT->current.setting.font_title);
			PSL_plottext (PSL, x[0], y[0], GMT->current.setting.font_title.size, mr->label[k], ew_angle, ljust[k], form);
			gmt_setpen (GMT, &GMT->current.setting.map_tick_pen[GMT_SECONDARY]);
		}
	}

	gmt_M_memset (dim, PSL_MAX_DIMS, double);
	dim[11] = 0.5 * GMT->current.setting.map_default_pen.width;
	PSL_defpen (PSL, "PSL_vecheadpen", GMT->current.setting.map_default_pen.width, GMT->current.setting.map_default_pen.style, GMT->current.setting.map_default_pen.offset, GMT->current.setting.map_default_pen.rgb);
	if (mr->kind == 2) {	/* Compass needle and label */
		char tmpstring[GMT_LEN64] = {""};
		PSL_comment (PSL, "Draw magnetic rose declination arrow and optional label\n");
		sincosd (ew_angle + (90.0 - mr->declination), &s, &c);
		L = R[GMT_ROSE_PRIMARY] - 2.0 * tlen[2];
		x[0] = mr->refpoint->x - L * c, y[0] = mr->refpoint->y - L * s;
		x[1] = mr->refpoint->x + L * c, y[1] = mr->refpoint->y + L * s;
		dim[PSL_VEC_XTIP]        = x[1];
		dim[PSL_VEC_YTIP]        = y[1];
		dim[PSL_VEC_TAIL_WIDTH]  = M_VW * mr->size;
		dim[PSL_VEC_HEAD_LENGTH] = M_HL * mr->size;
		dim[PSL_VEC_HEAD_WIDTH]  = M_HW * mr->size,
		dim[PSL_VEC_HEAD_SHAPE]  = GMT->current.setting.map_vector_shape;
		dim[PSL_VEC_STATUS]      = PSL_VEC_END | PSL_VEC_FILL;
		gmt_setpen (GMT, &GMT->current.setting.map_default_pen);
		gmt_setfill (GMT, &f, 1);
		PSL_plotsymbol (PSL, x[0], y[0], dim, PSL_VECTOR);
		t_angle = fmod (ew_angle + 90.0 - mr->declination + 360.0, 360.0);	/* Now in 0-360 range */
		if (fabs (t_angle) > 90.0) t_angle -= copysign (180.0, t_angle);
		sincosd (t_angle, &s, &c);
		x[0] = mr->refpoint->x - 2.0 * M_VW * mr->size * s, y[0] = mr->refpoint->y + 2.0 * M_VW * mr->size * c;
		if (strcmp (mr->dlabel, "-")) {	/* Want declination labeling unless when giving "-" */
			if (mr->dlabel[0] == 0) {	/* Want default label */
				gmtlib_get_annot_label (GMT, mr->declination, tmpstring, true, false, true, 0, GMT->current.map.is_world);
				snprintf (mr->dlabel, GMT_LEN256, "@~d@~ = %s", tmpstring);
			}
			form = gmt_setfont (GMT, &GMT->current.setting.font_label);
			PSL_plottext (PSL, x[0], y[0], GMT->current.setting.font_label.size, mr->dlabel, t_angle, PSL_BC, form);
		}
	}
	else {			/* Just geographic directions and a centered arrow */
		L = mr->size - 4.0*tlen[2];
		x[0] = x[1] = x[4] = 0.0,	x[2] = -0.25 * mr->size,		x[3] = -x[2];
		y[0] = -0.5 * L,		y[1] = -y[0], y[2] = y[3] = 0.0,	y[4] = y[1] + GMT->current.setting.map_annot_offset[GMT_PRIMARY];
		gmtlib_rotate2D (GMT, x, y, 5, mr->refpoint->x, mr->refpoint->y, ew_angle, xp, yp);	/* Coordinate transformation and placement of the 4 labels */
		dim[PSL_VEC_XTIP]        = xp[1];
		dim[PSL_VEC_YTIP]        = yp[1];
		dim[PSL_VEC_TAIL_WIDTH]  = F_VW * mr->size;
		dim[PSL_VEC_HEAD_LENGTH] = F_HL * mr->size;dim[PSL_VEC_HEAD_WIDTH] = F_HW * mr->size;
		dim[PSL_VEC_HEAD_SHAPE]  = GMT->current.setting.map_vector_shape;
		dim[PSL_VEC_STATUS]      = PSL_VEC_END | PSL_VEC_FILL;
		gmt_setfill (GMT, &f, 1);
		PSL_defpen (PSL, "PSL_vecheadpen", GMT->current.setting.map_frame_pen.width, "", 0, f.rgb);
		PSL_plotsymbol (PSL, xp[0], yp[0], dim, PSL_VECTOR);
		s = 0.25 * mr->size;
		gmt_init_fill (GMT, &f, -1.0, -1.0, -1.0);
		gmt_setfill (GMT, &f, 1);
		PSL_plotsymbol (PSL, mr->refpoint->x, mr->refpoint->y, &s, PSL_CIRCLE);
		PSL_plotsegment (PSL, xp[2], yp[2], xp[3], yp[3]);
	}

	PSL_settextmode (PSL, PSL_TXTMODE_HYPHEN);	/* Back to leave as is */
}

/* These are used to scale the plain arrow given rose size */

#define ROSE_LENGTH_SCL1	(0.5 * M_SQRT2)
#define ROSE_LENGTH_SCL2	0.5
#define ROSE_WIDTH_SCL1		0.2
#define ROSE_WIDTH_SCL2		0.2
#define ROSE_WIDTH_SCL3		0.2

GMT_LOCAL void gmtplot_draw_dir_rose (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, struct GMT_MAP_ROSE *mr) {
	unsigned int i, kind, form, just[4] = {PSL_TC, PSL_ML, PSL_BC, PSL_MR};
	int k;
	double angle, L[4], R[4], x[PSL_MAX_DIMS], y[8], xp[8], yp[8], tx[3], ty[3];
	double lon, lat, s, c, rot[4] = {0.0, 45.0, 22.5, -22.5};
	struct GMT_FILL f;

	/* Initialize fill structure */
	gmt_init_fill (GMT, &f, GMT->current.setting.color_patch[GMT_BGD][0], GMT->current.setting.color_patch[GMT_BGD][1], GMT->current.setting.color_patch[GMT_BGD][2]);

	gmt_xy_to_geo (GMT, &lon, &lat, mr->refpoint->x, mr->refpoint->y);
	angle = gmt_azim_to_angle (GMT, lon, lat, DIST_TO_2ND_POINT, 90.0);	/* Get angle of E-W direction at this location */

	gmt_setpen (GMT, &GMT->current.setting.map_tick_pen[GMT_PRIMARY]);

	if (mr->type == GMT_ROSE_DIR_FANCY) {	/* Fancy scale */
		PSL_comment (PSL, "Draw fancy directional rose of level %d\n", mr->kind);
		mr->size *= 0.5;	/* Got diameter, use radius for calculations */
		L[0] = mr->size;
		L[1] = ROSE_LENGTH_SCL1 * mr->size;
		L[2] = L[3] = ROSE_LENGTH_SCL2 * mr->size;
		R[0] = ROSE_WIDTH_SCL1 * mr->size;
		R[1] = ROSE_WIDTH_SCL2 * mr->size;
		R[2] = R[3] = ROSE_WIDTH_SCL3 * mr->size;
		mr->kind--;	/* Turn 1-3 into 0-2 */
		if (mr->kind == 2) mr->kind = 3;	/* Trick so that we can draw 8 rather than 4 points */
		for (k = kind = mr->kind; k >= 0; k--, kind--) {
			/* Do 4 blades 90 degrees apart, aligned with main axes & relative to (0,0) */
			x[0] = L[kind], x[1] = x[7] = 0.5 * M_SQRT2 * R[kind], x[2] = x[6] = 0.0;
			y[0] = y[4] = 0.0, y[1] = y[3] = 0.5 * M_SQRT2 * R[kind], y[2] = L[kind];
			x[3] = x[5] = -x[1], x[4] = -x[0];
			y[5] = y[7] = -y[1], y[6] = -y[2];
			gmtlib_rotate2D (GMT, x, y, 8, mr->refpoint->x, mr->refpoint->y, rot[kind] + angle, xp, yp);	/* Coordinate transformation and placement of the 4 labels */
			PSL_setfill (PSL, GMT->PSL->init.page_rgb, 1);
			PSL_plotpolygon (PSL, xp, yp, 8);	/* Outline of 4-pointed star */
			tx[0] = mr->refpoint->x, ty[0] = mr->refpoint->y;
			/* Fill positive halves of the 4-pointed blades */
			PSL_setfill (PSL, GMT->current.setting.map_default_pen.rgb, 1);
			tx[1] = xp[0], ty[1] = yp[0], tx[2] = xp[7], ty[2] = yp[7];
			PSL_plotpolygon (PSL, tx, ty, 3);	/* East */
			tx[1] = xp[1], ty[1] = yp[1], tx[2] = xp[2], ty[2] = yp[2];
			PSL_plotpolygon (PSL, tx, ty, 3);	/* North */
			tx[1] = xp[3], ty[1] = yp[3], tx[2] = xp[4], ty[2] = yp[4];
			PSL_plotpolygon (PSL, tx, ty, 3);	/* West */
			tx[1] = xp[5], ty[1] = yp[5], tx[2] = xp[6], ty[2] = yp[6];
			PSL_plotpolygon (PSL, tx, ty, 3);	/* South */
		}
		sincosd (angle, &s, &c);
		x[0] = x[2] = 0.0, x[1] = L[0] + GMT->current.setting.map_title_offset; x[3] = -x[1];
		y[1] = y[3] = 0.0, y[2] = L[0] + GMT->current.setting.map_title_offset; y[0] = -y[2];
		gmtlib_rotate2D (GMT, x, y, 4, mr->refpoint->x, mr->refpoint->y, angle, xp, yp);	/* Coordinate transformation and placement of the 4 labels */
		form = gmt_setfont (GMT, &GMT->current.setting.font_title);
		for (i = 0; i < 4; i++) PSL_plottext (PSL, xp[i], yp[i], GMT->current.setting.font_title.size, mr->label[i], angle, just[i], form);
	}
	else {			/* Plain North arrow w/circle */
		PSL_comment (PSL, "Draw plain directional rose\n");
		sincosd (angle, &s, &c);
		gmt_M_memset (x, PSL_MAX_DIMS, double);
		x[0] = x[1] = x[4] = 0.0, x[2] = -0.25 * mr->size, x[3] = -x[2];
		y[0] = -0.5 * mr->size, y[1] = -y[0], y[2] = y[3] = 0.0; y[4] = y[1] + GMT->current.setting.map_annot_offset[GMT_PRIMARY];
		gmtlib_rotate2D (GMT, x, y, 5, mr->refpoint->x, mr->refpoint->y, angle, xp, yp);	/* Coordinate transformation and placement of the 4 labels */
		x[0] = xp[1], x[1] = yp[1];
		x[2] = F_VW * mr->size, x[3] = F_HL * mr->size, x[4] = F_HW * mr->size;
		x[5] = GMT->current.setting.map_vector_shape, x[6] = PSL_VEC_END | PSL_VEC_FILL;
		gmt_setfill (GMT, &f, 1);
		PSL_defpen (PSL, "PSL_vecheadpen", GMT->current.setting.map_frame_pen.width, "", 0, f.rgb);
		PSL_plotsymbol (PSL, xp[0], yp[0], x, PSL_VECTOR);
		s = 0.25 * mr->size;
		gmt_init_fill (GMT, &f, -1.0, -1.0, -1.0);
		gmt_setfill (GMT, &f, 1);
		PSL_plotsymbol (PSL, mr->refpoint->x, mr->refpoint->y, &s, PSL_CIRCLE);
		PSL_plotsegment (PSL, xp[2], yp[2], xp[3], yp[3]);
		if (mr->label[2][0]) {	/* Wanted the north label */
			form = gmt_setfont (GMT, &GMT->current.setting.font_title);
			PSL_plottext (PSL, xp[4], yp[4], GMT->current.setting.font_title.size, mr->label[2], angle, PSL_BC, form);
		}
	}
}

GMT_LOCAL void gmtplot_savepen (struct GMT_CTRL *GMT, struct GMT_PEN *pen) {
	/* gmt_getpen retrieves the current pen in PSL. */
	struct PSL_CTRL *PSL = GMT->PSL;
	if (!pen) return;
	pen->width = PSL->current.linewidth;
	pen->offset = PSL->current.offset;
	if (PSL->current.style[0])
		strncpy (pen->style, PSL->current.style, GMT_PEN_LEN-1);
	else
		memset (pen->style, 0, GMT_PEN_LEN);
	gmt_M_rgb_copy (pen->rgb, PSL->current.rgb[PSL_IS_STROKE]);
}

enum plot_operand {
	LEFT_OPERAND = 0,
	RIGHT_OPERAND1 = 1,
	RIGHT_OPERAND2 = 2
};

GMT_LOCAL int gmtplot_custum_failed_bool_test_string (struct GMT_CTRL *GMT, struct GMT_CUSTOM_SYMBOL_ITEM *s, double size[], char *text, bool *retval) {
	unsigned int k;
	int type;
	bool result, delete[2] = {false, false};
	char *arg[2];
	gmt_M_unused (size);	/* Numerical values not used here */
	for (k = 0; k < 2; k++) {	/* Load up the left and right operands */
		type = (s->var[k] >= GMT_VAR_WORD) ? GMT_VAR_WORD : s->var[k];	/* If word then (s->var[k] - GMT_VAR_WORD) is the word number (0, 1, ...) */
		switch (type) {
			case GMT_VAR_STRING:	/* trailing text comparison */
				arg[k] = text;		break;
			case GMT_VAR_WORD:		/* trailing word comparison */
				arg[k] = gmt_get_word (text, " \t", s->var[0] - GMT_VAR_WORD);	delete[k] = true;	break;
			case GMT_CONST_STRING:	/* Constant text comparison */
				arg[k] = s->string;	break;
			default:	/* Should not get here */
				GMT_Report (GMT->parent, GMT_MSG_ERROR, "Unrecognized text variable type (%d) passed to gmtplot_custum_failed_bool_test_string\n", s->var[k]);
				return GMT_NOT_A_VALID_TYPE;		break;
		}
	}
	switch (s->operator) {
		case '<':	/* < */
			result = (strcmp (arg[LEFT_OPERAND], arg[RIGHT_OPERAND1]) < 0);
			break;
		case 'L':	/* <= */
			result = (strcmp (arg[LEFT_OPERAND], arg[RIGHT_OPERAND1]) <= 0);
			break;
		case '=':	/* == */
			result = (strcmp (arg[LEFT_OPERAND], arg[RIGHT_OPERAND1]) == 0);
			break;
		case 'G':	/* >= */
			result = (strcmp (arg[LEFT_OPERAND], arg[RIGHT_OPERAND1]) >= 0);
			break;
		case '>':	/* > */
			result = (strcmp (arg[LEFT_OPERAND], arg[RIGHT_OPERAND1]) > 0);
			break;
		default:
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Unrecognized symbol macro text operator (%d = '%c') passed to gmt_draw_custom_symbol\n", s->operator, (char)s->operator);
			return GMT_PARSE_ERROR;
			break;
	}
	for (k = 0; k < 2; k++)
		if (delete[k]) gmt_M_str_free (arg[k]);	/* Free word that we extracted from text */

	if (s->negate) result = !result;	/* Negate the test since we used a ! operator , e.g., != */
	*retval = !result;		/* Return the opposite of the test result */
	return (GMT_NOERROR);
}

GMT_LOCAL int gmtplot_custum_failed_bool_test (struct GMT_CTRL *GMT, struct GMT_CUSTOM_SYMBOL_ITEM *s, double size[], char *text, bool *retval) {
	unsigned int k;
	int error = GMT_NOERROR;
	bool result;
	double arg[3];

	/* Determine if we have text comparisons to deal with, if so, call the string version of this function */

	for (k = 0; k < 2; k++)
		if (s->var[k] == GMT_VAR_STRING || s->var[k] >= GMT_VAR_WORD) {
			error = gmtplot_custum_failed_bool_test_string (GMT, s, size, text, retval);
			return error;
		}

	/* Here we have numerical comparisons only */

	/* Perform the boolean comparison and return false if test is true */

	for (k = 0; k < 3; k++) {	/* Load up the left and 1-2 right operands */
		switch (s->var[k]) {
			case GMT_VAR_SIZE:	/* Symbol size */
				arg[k] = size[0];	break;
			case GMT_VAR_IS_Y:	/* User y-coordinate */
				arg[k] = GMT->current.io.curr_rec[GMT_Y];	break;
			case GMT_VAR_IS_X:	/* User x-coordinate */
				arg[k] = GMT->current.io.curr_rec[GMT_X];	break;
			case GMT_CONST_VAR:	/* A numeric constant */
				arg[k] = s->const_val[k];	break;
			default:		/* One of the variables 1-n */
				arg[k] = size[s->var[k]];	break;
		}
	}
	switch (s->operator) {
		case '<':	/* < */
			result = (arg[LEFT_OPERAND] < arg[RIGHT_OPERAND1]);
			break;
		case 'L':	/* <= */
			result = (arg[LEFT_OPERAND] <= arg[RIGHT_OPERAND1]);
			break;
		case '=':	/* == */
			result = (arg[LEFT_OPERAND] == arg[RIGHT_OPERAND1]);
			break;
		case 'G':	/* >= */
			result = (arg[LEFT_OPERAND] >= arg[RIGHT_OPERAND1]);
			break;
		case '>':	/* > */
			result = (arg[LEFT_OPERAND] > arg[RIGHT_OPERAND1]);
			break;
		case '%':	/* % */
			result = (!gmt_M_is_zero (fmod (arg[LEFT_OPERAND], arg[RIGHT_OPERAND1])));
			break;
		case 'I':	/* [] inclusive range */
			result = (arg[LEFT_OPERAND] >= arg[RIGHT_OPERAND1] && arg[LEFT_OPERAND] <= arg[RIGHT_OPERAND2]);
			break;
		case 'i':	/* <> exclusive range */
			result = (arg[LEFT_OPERAND] > arg[RIGHT_OPERAND1] && arg[LEFT_OPERAND] < arg[RIGHT_OPERAND2]);
			break;
		case 'l':	/* [> in/ex-clusive range */
			result = (arg[LEFT_OPERAND] >= arg[RIGHT_OPERAND1] && arg[LEFT_OPERAND] < arg[RIGHT_OPERAND2]);
			break;
		case 'r':	/* <] ex/in-clusive range */
			result = (arg[LEFT_OPERAND] > arg[RIGHT_OPERAND1] && arg[LEFT_OPERAND] <= arg[RIGHT_OPERAND2]);
			break;
		case 'E':	/* var == NaN */
			result = gmt_M_is_dnan (arg[LEFT_OPERAND]);
			break;
		default:
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Unrecognized symbol macro numerical operator (%d = '%c') passed to gmt_draw_custom_symbol\n", s->operator, (char)s->operator);
			return GMT_PARSE_ERROR;
			break;
	}
	if (s->negate) result = !result;	/* Negate the test since we used a ! operator , e.g., != */
	*retval = !result;		/* Return the opposite of the test result */
	return (GMT_NOERROR);
}

GMT_LOCAL void gmtplot_flush_symbol_piece (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double *x, double *y, uint64_t *n, struct GMT_PEN *p, struct GMT_FILL *f, int outline, bool *flush) {
	int draw_outline = (outline && p->rgb[0] != -1) ? 1 : 0;

	if (draw_outline) gmt_setpen (GMT, p);
	if (outline == 2) {	/* Stroke path only */
		PSL_plotline (PSL, x, y, (int)*n, PSL_MOVE|PSL_STROKE);
	}
	else {	/* Fill polygon and possibly stroke outline */
		gmt_setfill (GMT, f, draw_outline);
		PSL_plotpolygon (PSL, x, y, (int)*n);
	}
	*flush = false;
	*n = 0;
}

GMT_LOCAL void gmtplot_format_symbol_string (struct GMT_CTRL *GMT, struct GMT_CUSTOM_SYMBOL_ITEM *s, double size[], char *text) {
	/* Returns the [possibly reformatted] string to use for the letter macro.
	 * These are the things that can happen:
	 * 1. Action = GMT_SYMBOL_TEXT means we have a static fixed text string; just copy
	 * 2. s->text is $t.  Then we use the trailing text from the input as the text.
	 *    Optionally, a single word of the trailing text can be selected by appending
	 *    a word integer, $t0 is the first word, $t1 the second, etc.
	 * 3. We have a format statement that contains free-form text with interspersed
	 *    special formatting commands.  These have the syntax
	 *    %X  Add longitude or x using chosen default format.
	 *    %Y  Add latitude or y using chosen default format.
	 *    $<n>[+X|Y|T]  Format the numerical variable $<n>; if
	 *      followed by +X|Y|T we format as lon, lat, or absolute dateTtime,
	 *      else we use FORMAT_FLOAT_OUT.
	 *    Limitation: Currently, $<n> expects <n> to be 0-9 only.
	 */
	unsigned int n;
	if (s->action == GMT_SYMBOL_TEXT)	/* Constant text given, just duplicate */
		strcpy (text, s->string);
	else if (!strcmp (s->string, "$t"))	/* Get entire string from trailing text in the input */
		strcpy (text, GMT->current.io.curr_trailing_text);
	else if (!strncmp (s->string, "$t", 2U) && isdigit (s->string[2])) {	/* Get word number n from trailing text in the input */
		char *word = gmt_get_word (GMT->current.io.curr_trailing_text, " \t", s->var[0]);
		if (word) {	/* Got it */
			strcpy (text, word);
			gmt_M_str_free (word);
		}
		else {
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "No word # %d in the trailing text - return all text instead\n", s->var[0]);
			strcpy (text, GMT->current.io.curr_trailing_text);
		}
	}
	else {	/* Must replace special items within a template string */
		unsigned int n_skip, in, out;
		char tmp[GMT_LEN64] = {""};
		gmt_M_memset (text, GMT_LEN256, char);
		for (in = out = 0; s->string[in]; in++) {
			switch (s->string[in]) {
				case '%':	/* Possibly a special %X, %Y request */
					if (s->string[in+1] == 'X' || s->string[in+1] == 'Y') {	/* Yes it was */
						n = (s->string[in+1] == 'X') ? GMT_X : GMT_Y;
						gmt_ascii_format_col (GMT, tmp, GMT->current.io.curr_rec[n], GMT_IN, n);
						strcat (text, tmp);
						in++;	/* Skip past the X or Y */
						out += (unsigned int)strlen (tmp);
					}
					else /* Just a % sign */
						text[out++] = s->string[in];
					break;
				case '$':	/* Possibly a variable $n */
					if (isdigit (s->string[in+1])) {	/* Yes, it was */
						n = (s->string[in+1] - '0');
						n_skip = 1;
						if (s->string[in+2] == '+' && strchr ("TXY", s->string[in+3])) {	/* Specific formatting requested */
							if (s->string[in+3] == 'X') gmt_ascii_format_col (GMT, tmp, size[n], GMT_IN, GMT_X);
							else if (s->string[in+3] == 'Y') gmt_ascii_format_col (GMT, tmp, size[n], GMT_IN, GMT_Y);
							else if (s->string[in+3] == 'T') gmt_format_abstime_output (GMT, size[n], tmp);
							n_skip += 2;
						}
						else
							snprintf (tmp, GMT_LEN64, GMT->current.setting.format_float_out, size[n]);
						strcat (text, tmp);
						in += n_skip;	/* Skip past the $n[+X|Y|T] */
						out += (unsigned int)strlen (tmp);
					}
					else	/* Just pass regular text along */
						text[out++] = s->string[in];
					break;
				default:	/* Just pass regular text along */
					text[out++] = s->string[in];
					break;
			}
		}
	}
}

GMT_LOCAL void gmtplot_encodefont (struct PSL_CTRL *PSL, int font_no, char *name, unsigned int id) {
	/* Create the custom symbol macro that selects the correct font and size for the symbol item */

	bool encode = (PSL->init.encoding && !PSL->internal.font[font_no].encoded);

	if (PSL->internal.comments) PSL_command (PSL, "%% Set font encoding and size for this custom symbol %s item %d\n", name, id);
	PSL_command (PSL, "/PSL_symbol_%s_setfont_%d {", name, id);
	if (encode) {	/* Re-encode fonts with Standard+ or ISOLatin1[+] encodings */
		PSL_command (PSL, " PSL_font_encode %d get 0 eq {%s_Encoding /%s /%s PSL_reencode PSL_font_encode %d 1 put} if", font_no, PSL->init.encoding, PSL->internal.font[font_no].name, PSL->internal.font[font_no].name, font_no);
		PSL->internal.font[font_no].encoded = true;
	}
	PSL_command (PSL, " F%d } def\n", font_no);
}

GMT_LOCAL void gmtplot_contlabel_debug (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, struct GMT_CONTOUR *G) {
	uint64_t row;
	double size[1] = {0.025};

	/* If called we simply draw the helper lines or points to assist in debug */

	gmt_setpen (GMT, &G->debug_pen);
	if (G->fixed) {	/* Place a small open circle at each fixed point */
		PSL_setfill (PSL, GMT->session.no_rgb, PSL_OUTLINE);
		for (row = 0; row < (uint64_t)G->f_n; row++)
			PSL_plotsymbol (PSL, G->f_xy[0][row], G->f_xy[1][row], size, PSL_CIRCLE);
	}
	else if (G->crossing) {	/* Draw a thin line */
		uint64_t seg;
		unsigned int *pen = NULL;
		struct GMT_DATASEGMENT *S = NULL;
		for (seg = 0; seg < G->X->n_segments; seg++) {
			S = G->X->table[0]->segment[seg];	/* Current segment */
			pen = gmt_M_memory (GMT, NULL, S->n_rows, unsigned int);
			for (row = 1, pen[0] = PSL_MOVE; row < S->n_rows; row++) pen[row] = PSL_DRAW;
			gmt_plot_line (GMT, S->data[GMT_X], S->data[GMT_Y], pen, S->n_rows, PSL_LINEAR);
			gmt_M_free (GMT, pen);
		}
	}
}

GMT_LOCAL void gmtplot_contlabel_drawlines (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, struct GMT_CONTOUR *G, unsigned int mode) {
	uint64_t seg, k;
	unsigned int *pen = NULL;
	struct GMT_CONTOUR_LINE *L = NULL;
	for (seg = 0; seg < G->n_segments; seg++) {
		L = G->segment[seg];	/* Pointer to current segment */
		if (L->annot && mode == 1) continue; /* Annotated lines done with curved text routine */
		gmt_setpen (GMT, &L->pen);
		pen = gmt_M_memory (GMT, NULL, L->n, unsigned int);
		for (k = 1, pen[0] = PSL_MOVE; k < L->n; k++) pen[k] = PSL_DRAW;
		PSL_comment (PSL, "%s: %s\n", G->line_name, L->name);
		gmt_plot_line (GMT, L->x, L->y, pen, L->n, PSL_LINEAR);
		gmt_M_free (GMT, pen);
	}
}

GMT_LOCAL void gmtplot_contlabel_plotlabels (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, struct GMT_CONTOUR *G, unsigned int mode) {
	/* mode controls what takes place:
	 * mode = 1: We place all the PSL variables required to use the text for clipping of painting.
	 * mode = 2: We paint the text that is stored in the PSL variables.
	 * mode = 4: We use the text stored in the PSL variables to set up a clip path.  CLipping is turned ON.
	 * mode = 8: We draw the lines.
	 * mode = 16: We turn clip path OFF.
	 * mode = 32: We want rounded rectangles instead of straight rectangular boxes [straight text only].
	 */
	int justify = 0, form = 0, *just = NULL, *node = NULL, *nlabels_per_segment = NULL, *npoints_per_segment = NULL;
	uint64_t k, m, seg, n_points = 0;
	unsigned int n_segments = 0, n_labels = 0, first_point_in_segment = 0, this_seg;
	double *angle = NULL, *xpath = NULL, *ypath = NULL, *xtxt = NULL, *ytxt = NULL;
	char **txt = NULL, **pen = NULL, **fonts = NULL;
	struct GMT_CONTOUR_LINE *L = NULL;
	void *A1 = NULL, *A2 = NULL;

	form = mode;					/* Which actions to take */
	if (G->box & 4) form |= PSL_TXT_ROUND;		/* Want round box shape */
	if (G->curved_text) form |= PSL_TXT_CURVED;	/* Want text set along curved path */
	if (G->fillbox) form |= PSL_TXT_FILLBOX;	/* Want the box filled */
	if (G->box & 1) form |= PSL_TXT_DRAWBOX;	/* Want box outline */
	if (G->font_label.form & 2)	/* Must supply both font rgb and pen rgb and use charpath and S to fill and outline text */
		form |= (G->font_label.form & 8) ? PSL_TXT_PENFILL : PSL_TXT_FILLPEN;
	if (mode & PSL_TXT_INIT) {	/* Determine and places all PSL attributes */
		char *font = NULL;
		if (G->number_placement && G->n_cont == 1)		/* Special 1-label justification check */
			justify = G->end_just[(G->number_placement+1)/2];	/* Gives index 0 or 1 */
		else
			justify = G->just;

		for (seg = 0; seg < G->n_segments; seg++) {
			L = G->segment[seg];	/* Pointer to current segment */
			n_segments++;			/* Number of segments */
			n_points += (unsigned int)L->n;		/* Total number of points in all segments so far */
			n_labels += L->n_labels;	/* Number of labels so far */
		}

		if (n_labels == 0) return;	/* There are no labels */

		gmt_M_malloc2 (GMT, xpath, ypath, n_points, NULL, double);
		gmt_M_malloc2 (GMT, npoints_per_segment, nlabels_per_segment, n_segments, NULL, int);
		if (G->curved_text) {	/* Must pass node locations of labels */
			node = gmt_M_memory (GMT, NULL, n_labels, int);
			A1 = node;
		}
		else {			/* Must pass x,y locations of labels */
			gmt_M_malloc2 (GMT, xtxt, ytxt, n_points, NULL, double);
			A1 = xtxt;
			A2 = ytxt;
		}
		angle = gmt_M_memory (GMT, NULL, n_labels, double);
		txt   = gmt_M_memory (GMT, NULL, n_labels, char *);
		pen   = gmt_M_memory (GMT, NULL, n_segments, char *);
		fonts = gmt_M_memory (GMT, NULL, n_labels, char *);
		PSL_setfont (PSL, G->font_label.id);
		for (seg = m = this_seg = 0; seg < G->n_segments; seg++) {	/* Process all segments, skip those without labels */
			L = G->segment[seg];	/* Pointer to current segment */
			npoints_per_segment[this_seg] = (int)L->n;	/* Points along this segment path */
			if (this_seg > 0) first_point_in_segment += npoints_per_segment[this_seg-1];		/* First point id in combined path */
			gmt_M_memcpy (&xpath[first_point_in_segment], L->x, L->n, double);	/* Append this segment path to the combined path array */
			gmt_M_memcpy (&ypath[first_point_in_segment], L->y, L->n, double);
			nlabels_per_segment[this_seg] = L->n_labels;		/* Number of labels for this path */
			pen[this_seg] = strdup (PSL_makepen (GMT->PSL, L->pen.width, L->pen.rgb, L->pen.style, L->pen.offset));	/* Get pen PSL setting for this segment */
			for (k = 0; k < L->n_labels; k++, m++) {	/* Process all labels for this segment */
				angle[m] = L->L[k].angle;
				txt[m]   = L->L[k].label;
				if (G->curved_text)	/* Need local node number for text placement */
					node[m]  = (int)L->L[k].node;	/* node is a local index for the relevant segment */
				else {	/* Need coordinate of text placement */
					xtxt[m]    = L->L[k].x;
					ytxt[m]    = L->L[k].y;
				}
				if (G->font_label.form & 2) {	/* Must supply both font rgb and pen rgb and use charpath and S to fill and outline text; see PSL_labels.ps */
					char string[GMT_LEN128] = {""};
					char *F = PSL_makecolor (PSL, L->L[k].rgb);
					char *P = PSL_makepen (PSL, G->font_label.pen.width, G->font_label.pen.rgb, G->font_label.pen.style, G->font_label.pen.offset);
					font = PSL_makefontsize (PSL, G->font_label.size);
					snprintf (string, GMT_LEN128, "{%s} FS %s %s", F, P, font);	/* E.g., "{1 0 0 C} FS 4 W 0 A [] 0 B 667 F5" */
					fonts[m] = strdup (string);
				}
				else {	/* Regular font fill */
					font = PSL_makefont (PSL, G->font_label.size, L->L[k].rgb);	/* e.g., "1 0 0 C 667 F5" */
					fonts[m] = strdup (font);
				}
			}
			this_seg++;
		}

		PSL_comment (PSL, "Store path and label attributes:\n");
		gmt_textpath_init (GMT, &G->pen, G->rgb);
		PSL_comment (PSL, "Store pens used for each line segment:\n");
		psl_set_txt_array (PSL, "path_pen", pen, n_segments);
		/* While all contours & quoted lines have same justification and font (except maybe color), this is not true of pstext items
		 * so we use separate arrays for both font and justify [here justify and possibly font are constants]  */
		just = gmt_M_memory (GMT, NULL, n_labels, int);
		for (k = 0; k < n_labels; k++)
			just[k] = abs (justify);
		PSL_comment (PSL, "Store text justification for each text label:\n");
		psl_set_int_array   (PSL, "label_justify", just, n_labels);
		PSL_comment (PSL, "Store font setting for each text label:\n");
		psl_set_txt_array   (PSL, "label_font", fonts, n_labels);
		gmt_M_free (GMT, just);
		for (k = 0; k < n_labels; k++)
			gmt_M_str_free (fonts[k]);
		gmt_M_free (GMT, fonts);
	}
	PSL_plottextline (PSL, xpath, ypath, npoints_per_segment, n_segments, A1, A2, txt, angle, nlabels_per_segment, G->font_label.size, justify, G->clearance, form);
	if (mode & PSL_TXT_INIT) {	/* Free up the things we allocated above */
		gmt_M_free (GMT, npoints_per_segment);
		gmt_M_free (GMT, nlabels_per_segment);
		for (k = 0; k < n_segments; k++) gmt_M_str_free (pen[k]);
		gmt_M_free (GMT, pen);
		gmt_M_free (GMT, angle);
		gmt_M_free (GMT, txt);
		gmt_M_free (GMT, xpath);
		gmt_M_free (GMT, ypath);
		if (G->curved_text)
			gmt_M_free (GMT, node);
		else {
			gmt_M_free (GMT, xtxt);
			gmt_M_free (GMT, ytxt);
		}
	}
}

#if 0
GMT_LOCAL void gmtplot_ellipsoid_name_convert2 (char *inname, char outname[]) {
	/* Convert the ellipsoid names to the slightly different way that they are called in proj4 */
	if (!strcmp(inname, "WGS84"))
		sprintf(outname, "WGS-84");
	else if (!strcmp(inname, "WGS72"))
		sprintf(outname, "WGS-72");
	else if (!strcmp(inname, "WGS66"))
		sprintf(outname, "WGS-66");
	else if (!strcmp(inname, "WGS60"))
		sprintf(outname, "WGS-60");
	else if (!strcmp(inname, "airy"))
		sprintf(outname, "Airy");
	else if (!strcmp(inname, "mod_airy"))
		sprintf(outname, "Airy-Ireland");
	else if (!strcmp(inname, "andrae"))
		sprintf(outname, "Andrae");
	else if (!strcmp(inname, "APL4.9"))
		sprintf(outname, "APL4.9");
	else if (!strcmp(inname, "aust_SA"))
		sprintf(outname, "Australian");
	else if (!strcmp(inname, "bessel"))
		sprintf(outname, "Bessel");
	else if (!strcmp(inname, "bess_nam"))
		sprintf(outname, "Bessel-Namibia");
	else if (!strcmp(inname, "clark66"))
		sprintf(outname, "Clarke-1866");
	else if (!strcmp(inname, "clark80"))
		sprintf(outname, "Clarke-1880");
	else if (!strcmp(inname, "CPM"))
		sprintf(outname, "CPM");
	else if (!strcmp(inname, "delmbr"))
		sprintf(outname, "Delambre");
	else if (!strcmp(inname, "engelis"))
		sprintf(outname, "Engelis");
	else if (!strcmp(inname, "evrst30"))
		sprintf(outname, "Everest-1830");
	else if (!strcmp(inname, "evrst48"))
		sprintf(outname, "Everest-1830-Kertau");
	else if (!strcmp(inname, "evrst56"))
		sprintf(outname, "Everest-1830-Kalianpur");
	else if (!strcmp(inname, "evrstSS"))
		sprintf(outname, "Everest-1830-Timbalai");
	else if (!strcmp(inname, "fschr60"))
		sprintf(outname, "Fischer-1960");
	else if (!strcmp(inname, "fschr60m"))
		sprintf(outname, "Fischer-1960-SouthAsia");
	else if (!strcmp(inname, "fschr68"))
		sprintf(outname, "Fischer-1968");
	else if (!strcmp(inname, "GRS80"))
		sprintf(outname, "GRS-80");
	else if (!strcmp(inname, "GRS67"))
		sprintf(outname, "GRS-67");
	else if (!strcmp(inname, "helmert"))
		sprintf(outname, "Helmert-1906");
	else if (!strcmp(inname, "hough"))
		sprintf(outname, "Hough");
	else if (!strcmp(inname, "intl"))
		sprintf(outname, "Hayford-1909");
	else if (!strcmp(inname, "new_intl"))
		sprintf(outname, "International-1967");
	else if (!strcmp(inname, "MERIT"))
		sprintf(outname, "MERIT-83");
	else if (!strcmp(inname, "krass"))
		sprintf(outname, "Krassovsky");
	else if (!strcmp(inname, "kaula"))
		sprintf(outname, "Kaula");
	else if (!strcmp(inname, "NWL-9D"))
		sprintf(outname, "NWL9D");
	else if (!strcmp(inname, "IAU76"))
		sprintf(outname, "IAG-75");
	else if (!strcmp(inname, "lerch"))
		sprintf(outname, "Lerch");
	else if (!strcmp(inname, "mprts"))
		sprintf(outname, "Maupertius");
	else if (!strcmp(inname, "SEasia"))
		sprintf(outname, "Modified-Fischer-1960");
	else if (!strcmp(inname, "SGS85"))
		sprintf(outname, "SGS-85");
	else if (!strcmp(inname, "plessis"))
		sprintf(outname, "Plessis");
	else if (!strcmp(inname, "walbeck"))
		sprintf(outname, "Walbeck");
	else if (!strcmp(inname, "sphere"))
		sprintf(outname, "Sphere");
	else if (!strcmp(inname, "sphere"))
		sprintf(outname, "FlatEarth");
	else
		outname[0] = '\0';
}
#endif

void gmtlib_ellipsoid_name_convert (char *inname, char outname[]) {
	/* Convert the ellipsoid names to the slightly different way that they are called in proj4 */
	if (!strcmp(inname, "WGS-84"))
		sprintf(outname, "WGS84");
	else if (!strcmp(inname, "WGS-72"))
		sprintf(outname, "WGS72");
	else if (!strcmp(inname, "WGS-66"))
		sprintf(outname, "WGS66");
	else if (!strcmp(inname, "WGS-60"))
		sprintf(outname, "WGS60");
	else if (!strcmp(inname, "Airy"))
		sprintf(outname, "airy");
	else if (!strcmp(inname, "Airy-Ireland"))
		sprintf(outname, "mod_airy");
	else if (!strcmp(inname, "Andrae"))
		sprintf(outname, "andrae");
	else if (!strcmp(inname, "APL4.9"))
		sprintf(outname, "APL4.9");
	else if (!strcmp(inname, "Australian"))
		sprintf(outname, "aust_SA");
	else if (!strcmp(inname, "Bessel"))
		sprintf(outname, "bessel");
	else if (!strcmp(inname, "Bessel-Namibia"))
		sprintf(outname, "bess_nam");
	else if (!strcmp(inname, "Clarke-1866"))
		sprintf(outname, "clark66");
	else if (!strcmp(inname, "Clarke-1880"))
		sprintf(outname, "clark80");
	else if (!strcmp(inname, "CPM"))
		sprintf(outname, "CPM");
	else if (!strcmp(inname, "Delambre"))
		sprintf(outname, "delmbr");
	else if (!strcmp(inname, "Engelis"))
		sprintf(outname, "engelis");
	else if (!strcmp(inname, "Everest-1830"))
		sprintf(outname, "evrst30");
	else if (!strcmp(inname, "Everest-1830-Kertau"))
		sprintf(outname, "evrst48");
	else if (!strcmp(inname, "Everest-1830-Kalianpur"))
		sprintf(outname, "evrst56");
	else if (!strcmp(inname, "Everest-1830-Timbalai"))
		sprintf(outname, "evrstSS");
	else if (!strcmp(inname, "Fischer-1960"))
		sprintf(outname, "fschr60");
	else if (!strcmp(inname, "Fischer-1960-SouthAsia"))
		sprintf(outname, "fschr60m");
	else if (!strcmp(inname, "Fischer-1968"))
		sprintf(outname, "fschr68");
	else if (!strcmp(inname, "GRS-80"))
		sprintf(outname, "GRS80");
	else if (!strcmp(inname, "GRS-67"))
		sprintf(outname, "GRS67");
	else if (!strcmp(inname, "Helmert-1906"))
		sprintf(outname, "helmert");
	else if (!strcmp(inname, "Hough"))
		sprintf(outname, "hough");
	else if (!strcmp(inname, "Hayford-1909"))
		sprintf(outname, "intl");
	else if (!strcmp(inname, "International-1967"))
		sprintf(outname, "new_intl");
	else if (!strcmp(inname, "MERIT-83"))
		sprintf(outname, "MERIT");
	else if (!strcmp(inname, "Krassovsky"))
		sprintf(outname, "krass");
	else if (!strcmp(inname, "Kaula"))
		sprintf(outname, "kaula");
	else if (!strcmp(inname, "NWL-9D"))
		sprintf(outname, "NWL9D");
	else if (!strcmp(inname, "IAG-75"))
		sprintf(outname, "IAU76");
	else if (!strcmp(inname, "Lerch"))
		sprintf(outname, "lerch");
	else if (!strcmp(inname, "Maupertius"))
		sprintf(outname, "mprts");
	else if (!strcmp(inname, "Modified-Fischer-1960"))
		sprintf(outname, "SEasia");
	else if (!strcmp(inname, "SGS-85"))
		sprintf(outname, "SGS85");
	else if (!strcmp(inname, "Plessis"))
		sprintf(outname, "plessis");
	else if (!strcmp(inname, "Walbeck"))
		sprintf(outname, "walbeck");
	else if (!strcmp(inname, "Sphere"))
		sprintf(outname, "sphere");
	else if (!strcmp(inname, "FlatEarth"))
		sprintf(outname, "sphere");
	else
		sprintf(outname, "unnamed");
}

#if 0
/* Used to dump an array to file for debug */
GMT_LOCAL void gmtplot_dumpfile (struct GMT_CTRL *GMT, double *x, double *y, unsigned int *pen, uint64_t n, char *file) {
	FILE *fp = fopen (file, "w");
	uint64_t k;
	unsigned int ps = 0;
	double w;
	for (k = 0; k < n; k++) {
		w = 2.0 * gmt_half_map_width (GMT, y[k]);
		if (pen) {
			ps += pen[k];
			fprintf (fp, "%.16g\t%.16g\t%d\t%.16g\n", x[k], y[k], ps, w);
		}
		else
			fprintf (fp, "%.16g\t%.16g\t%.16g\n", x[k], y[k], w);
			}
	fclose (fp);
}
#endif

GMT_LOCAL uint64_t gmtplot_geo_polygon (struct GMT_CTRL *GMT, double *lon, double *lat, uint64_t n, bool init, const char *comment) {
	/* When geographic data are plotted, polygons that cross the west map boundary will
	 * sometimes appear on the area bounded by the east map boundary - they "wrap around".
	 * This usually means we have a global map with (east-west) = 360.
	 * This function solves this by determining the polygon outline three times:
	 * First time: Truncate polygon between left and right border
	 * Second time: Find where the polygon jumps and set all the points between jumps to
	 *	       the point on the west boundary at the same latitude.
	 * Third time: Find where the polygon jumps and set all the points between jumps to
	 *	       the point on the east boundary at the same latitude.
	 * In reality it depends on the nature of the first jump in which order we do the
	 * west and east truncation above.
	 * If the polygon is clipped or wraps around at a periodic boundary then we must
	 * be careful how we draw the outline (if selected).  This function only lays down
	 * the paths; filling/outline is controlled by higher powers (gmt_geo_polygons).
	 */

#define JUMP_L 0
#define JUMP_R 1
#define JUMP_B 0
#define JUMP_T 1

	uint64_t total = 0;
	double *xp = NULL, *yp = NULL;
	struct PSL_CTRL *PSL= GMT->PSL;

	if (gmt_M_eq (PSL->current.rgb[PSL_IS_FILL][0], -1.0)) {
		/* Just draw optional outline, no fill, nor pattern */
	}
	else if (gmt_M_is_azimuthal (GMT) || !GMT->current.map.is_world) {	/* Testing without !is_world map to rediscover the original issue */
		/* Because points way outside the map might get close to the antipode we must
		 * clip the polygon first.  The new radial clip handles this by excluding points
		 * beyond the horizon and adding arcs along the boundary between exit points
		 */

		if ((GMT->current.plot.n = gmt_clip_to_map (GMT, lon, lat, n, &xp, &yp)) == 0) return 0;		/* All points are outside region */
		//gmtplot_dumpfile (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, "raw.txt");

		if (init) {
			PSL_comment (PSL, "Temporarily set FO to P for complex polygon building\n");
			PSL_command (PSL, "/FO {P}!\n");		/* Temporarily replace FO so we can build a complex path of closed polygons using {P} */
		}
		PSL_comment (PSL, comment);
		PSL_plotpolygon (PSL, xp, yp, (unsigned int)GMT->current.plot.n);	/* Fill Cartesian polygon and possibly draw outline */
		/* Free the memory we are done with */
		gmt_M_free (GMT, xp);
		gmt_M_free (GMT, yp);
		total = GMT->current.plot.n;
	}
	else if (GMT->current.proj.projection_GMT == GMT_TM && gmt_M_360_range (GMT->common.R.wesn[XLO], GMT->common.R.wesn[XHI])) {	/* Here, any jumps are in the y-direction */
		uint64_t k, first, i;
		int jump_dir = JUMP_B;
		bool jump, plot_main = true;
		double y_on_border[2] = {GMT->current.proj.rect[YHI], GMT->current.proj.rect[YLO]};

		/* Here we come for the periodic Transverse Mercator projection where longitude wrapping happens in the y-direction */

		if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, lon, lat, n)) == 0) return 0;		/* Convert to (x,y,pen) - return if nothing to do */
		if (init) {
			PSL_comment (PSL, "Temporarily set FO to P for complex polygon building\n");
			PSL_command (PSL, "/FO {P}!\n");		/* Temporarily replace FO so we can build a complex path of closed polygons using {P} */
		}
		PSL_comment (PSL, comment);
		//gmtplot_dumpfile (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, "raw.txt");

		/* Check if there are any boundary jumps in the data as evidenced by pen up [PSL_MOVE] */

		jump = (*GMT->current.map.will_it_wrap) (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.n, &first);	/* Polygon does indeed wrap */

		if (!jump) {	/* We happened to avoid the periodic boundary - just paint and return */
			PSL_plotpolygon (PSL, GMT->current.plot.x, GMT->current.plot.y, (unsigned int)GMT->current.plot.n);
			return GMT->current.plot.n;
		}

		GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Polygon wraps in y-direction for TM global projection\n");

		/* Polygon wraps and we will plot it up to three times by truncating the part that would wrap the wrong way.
		 * Here we cannot use the clipped/wrapped polygon to draw outline - that is done at the end, separately.
		 * Note: We may still have problems if a polygon that wraps across the top/bottom boundary also is clipped
		 * by the left or right boundary. */

		/* Temporary array to hold the modified y values */

		yp = gmt_M_memory (GMT, NULL, GMT->current.plot.n, double);

		/* Do the main truncation of bulk of polygon */

		for (i = 0, jump = false; i < GMT->current.plot.n; i++) {
			if (GMT->current.plot.pen[i] & PSL_MOVE && i) {
				jump = !jump;
				jump_dir = (GMT->current.plot.y[i] > GMT->current.map.half_height) ? JUMP_T : JUMP_B;
			}
			yp[i] = (jump) ? y_on_border[jump_dir] : GMT->current.plot.y[i];
		}

		//gmtplot_dumpfile (GMT, GMT->current.plot.x, yp, NULL, GMT->current.plot.n, "main.txt");
		if (plot_main) {
			PSL_plotpolygon (PSL, GMT->current.plot.x, yp, (unsigned int)GMT->current.plot.n);	/* Paint the truncated polygon */
			total = GMT->current.plot.n;
		}

		/* Then do the Bottom truncation since some wrapped pieces might not have been plotted (k > 0 means we found a piece) */

		jump_dir = (GMT->current.plot.y[first] > GMT->current.map.half_height) ? JUMP_B : JUMP_T;	/* Opposite */
		for (i = k = 0, jump = true; i < GMT->current.plot.n; i++) {
			if ((GMT->current.plot.pen[i] & PSL_MOVE) && i) {
				jump = !jump;
				jump_dir = (GMT->current.plot.y[i] > GMT->current.map.half_height) ? JUMP_T : JUMP_B;
			}
			yp[i] = (jump || jump_dir == JUMP_T) ? y_on_border[JUMP_T]: GMT->current.plot.y[i], k++;
		}
		if (k) {
			//gmtplot_dumpfile (GMT, GMT->current.plot.x, yp, NULL, GMT->current.plot.n, "B.txt");
			PSL_plotpolygon (PSL, GMT->current.plot.x, yp, (unsigned int)GMT->current.plot.n);	/* Paint the truncated polygon */
			total += GMT->current.plot.n;
		}

		/* Then do the Top truncation since some wrapped pieces might not have been plotted (k > 0 means we found a piece) */

		jump_dir = (GMT->current.plot.y[first] > GMT->current.map.half_height) ? JUMP_T : JUMP_B;	/* Opposite */
		for (i = k = 0, jump = true; i < GMT->current.plot.n; i++) {
			if ((GMT->current.plot.pen[i] & PSL_MOVE) && i) {
				jump = !jump;
				jump_dir = (GMT->current.plot.y[i] > GMT->current.map.half_height) ? JUMP_T : JUMP_B;
			}
			yp[i] = (jump || jump_dir == JUMP_L) ? y_on_border[JUMP_B] : GMT->current.plot.y[i], k++;
		}
		if (k) {
			//gmtplot_dumpfile (GMT, GMT->current.plot.x, yp, NULL, GMT->current.plot.n, "T.txt");
			PSL_plotpolygon (PSL, GMT->current.plot.x, yp, (unsigned int)GMT->current.plot.n);	/* Paint the truncated polygon */
			total = GMT->current.plot.n;
		}

		/* Free the memory we are done with */
		gmt_M_free (GMT, yp);
	}
	else {
		uint64_t k, first, i;
		int jump_dir = JUMP_L;
		bool jump, plot_main = true;
		double (*x_on_border[2]) (struct GMT_CTRL *, double) = {NULL, NULL};

		/* Here we come for all other non-azimuthal projections */

		if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, lon, lat, n)) == 0) return 0;		/* Convert to (x,y,pen) - return if nothing to do */
		if (init) {
			PSL_comment (PSL, "Temporarily set FO to P for complex polygon building\n");
			PSL_command (PSL, "/FO {P}!\n");		/* Temporarily replace FO so we can build a complex path of closed polygons using {P} */
		}
		PSL_comment (PSL, comment);
		//gmtplot_dumpfile (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, "raw.txt");

		if (gmt_M_is_cartesian (GMT, GMT_IN)) {		/* Not geographic data so there are no periodic boundaries to worry about */
			PSL_plotpolygon (PSL, GMT->current.plot.x, GMT->current.plot.y, (unsigned int)GMT->current.plot.n);
			return GMT->current.plot.n;
		}

		/* Check if there are any boundary jumps in the data as evidenced by pen up [PSL_MOVE] */

		jump = (*GMT->current.map.will_it_wrap) (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.n, &first);	/* Polygon does indeed wrap */

		if (!jump) {	/* We happened to avoid the periodic boundary - just paint and return */
			PSL_plotpolygon (PSL, GMT->current.plot.x, GMT->current.plot.y, (unsigned int)GMT->current.plot.n);
			return GMT->current.plot.n;
		}

		/* Polygon wraps and we will plot it up to three times by truncating the part that would wrap the wrong way.
		 * Here we cannot use the clipped/wrapped polygon to draw outline - that is done at the end, separately */

		/* Temporary array to hold the modified x values */

		xp = gmt_M_memory (GMT, NULL, GMT->current.plot.n, double);

		x_on_border[JUMP_R] = gmtlib_left_boundary;	/* Pointers to functions that supply the x-coordinate of boundary for given y */
		x_on_border[JUMP_L] = gmtlib_right_boundary;

		/* Do the main truncation of bulk of polygon */

		for (i = 0, jump = false; i < GMT->current.plot.n; i++) {
			if (GMT->current.plot.pen[i] & PSL_MOVE && i) {
				jump = !jump;
				jump_dir = (GMT->current.plot.x[i] > GMT->current.map.half_width) ? JUMP_R : JUMP_L;
			}
			xp[i] = (jump) ? (*x_on_border[jump_dir]) (GMT, GMT->current.plot.y[i]) : GMT->current.plot.x[i];
		}
		if (GMT->current.map.is_world && doubleAlmostEqualZero (GMT->current.plot.y[0], GMT->current.plot.y[GMT->current.plot.n-1]) && !(gmt_M_is_zero (GMT->current.plot.y[0]) || doubleAlmostEqualZero (GMT->current.plot.y[0], GMT->current.proj.rect[YHI])) ) {
			/* Watch for E-W line crossing jumps. We have no defense against these erratic near-pole polygons that are split across two periodic boundaries...
			 * The test checks for horizontal lines NOT at top of bottom of plot, for global maps only.
			 * For an example that triggers this, see test psxy/nearpole.sh and comments therein. PW, 5/29/2018 */
			double w = 1.9 * gmt_half_map_width (GMT, GMT->current.plot.y[0]);
			if (fabs (xp[GMT->current.plot.n-1] - xp[0]) > w) {	/* Does the jump exceed 85% of map width? */
				GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "Truncated wrapped polygon still has E-W jump and was skipped. Please report to developers if plot has artifacts.\n");
				plot_main = false;
			}
		}
		//gmtplot_dumpfile (GMT, xp, GMT->current.plot.y, NULL, GMT->current.plot.n, "main.txt");
		if (plot_main) {
			PSL_plotpolygon (PSL, xp, GMT->current.plot.y, (unsigned int)GMT->current.plot.n);	/* Paint the truncated polygon */
			total = GMT->current.plot.n;
		}

		/* Then do the Left truncation since some wrapped pieces might not have been plotted (k > 0 means we found a piece) */

		jump_dir = (GMT->current.plot.x[first] > GMT->current.map.half_width) ? JUMP_L : JUMP_R;	/* Opposite */
		for (i = k = 0, jump = true; i < GMT->current.plot.n; i++) {
			if ((GMT->current.plot.pen[i] & PSL_MOVE) && i) {
				jump = !jump;
				jump_dir = (GMT->current.plot.x[i] > GMT->current.map.half_width) ? JUMP_R : JUMP_L;
			}
			xp[i] = (jump || jump_dir == JUMP_R) ? (*x_on_border[JUMP_R]) (GMT, GMT->current.plot.y[i]) : GMT->current.plot.x[i], k++;
		}
		if (k) {
			//gmtplot_dumpfile (GMT, xp, GMT->current.plot.y, NULL, GMT->current.plot.n, "L.txt");
			PSL_plotpolygon (PSL, xp, GMT->current.plot.y, (unsigned int)GMT->current.plot.n);	/* Paint the truncated polygon */
			total += GMT->current.plot.n;
		}

		/* Then do the R truncation since some wrapped pieces might not have been plotted (k > 0 means we found a piece) */

		jump_dir = (GMT->current.plot.x[first] > GMT->current.map.half_width) ? JUMP_R : JUMP_L;	/* Opposite */
		for (i = k = 0, jump = true; i < GMT->current.plot.n; i++) {
			if ((GMT->current.plot.pen[i] & PSL_MOVE) && i) {
				jump = !jump;
				jump_dir = (GMT->current.plot.x[i] > GMT->current.map.half_width) ? JUMP_R : JUMP_L;
			}
			xp[i] = (jump || jump_dir == JUMP_L) ? (*x_on_border[JUMP_L]) (GMT, GMT->current.plot.y[i]) : GMT->current.plot.x[i], k++;
		}
		if (k) {
			//gmtplot_dumpfile (GMT, xp, GMT->current.plot.y, NULL, GMT->current.plot.n, "R.txt");
			PSL_plotpolygon (PSL, xp, GMT->current.plot.y, (unsigned int)GMT->current.plot.n);	/* Paint the truncated polygon */
			total = GMT->current.plot.n;
		}

		/* Free the memory we are done with */
		gmt_M_free (GMT, xp);
	}
	return (total);
}

GMT_LOCAL void gmtplot_reverse_polygon (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S) {
	uint64_t k, n1 = S->n_rows - 1;
	/* Reverse the direction of this polygon, i.e, swap points n-1-k and k */
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Make polygon to clockwise\n");
	for (k = 0; k < S->n_rows/2; k++) {
		gmt_M_double_swap (S->data[GMT_X][k], S->data[GMT_X][n1-k]);
		gmt_M_double_swap (S->data[GMT_Y][k], S->data[GMT_Y][n1-k]);
	}
}

#if 0
GMT_LOCAL uint64_t gmtplot_geo_polarcap_segment_orig (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S, bool first, const char *comment) {
	/* Special treatment for polar caps since they must add int parts of possibly curved periodic boundaries
	 * from the pole up to the intersection with the cap perimeter.  We handle this case separately here.
	 * This is in response to issue # 852. P. Wessel */
	uint64_t k0, perim_n, n_new, m, n = S->n_rows, k;
	double start_lon, stop_lon, yc = 0.0, dx, pole_lat = 90.0 * S->pole;
	double *x_perim = NULL, *y_perim = NULL, *plon = NULL, *plat = NULL;
	static char *pole = "S N";
	int type;
#if 0
	FILE *fp;
#endif
	/* We want this code to be used for the misc. global projections but also global cylindrical or linear(if degrees) maps */
	if (!(gmt_M_is_misc(GMT) || (GMT->current.map.is_world  && (gmt_M_is_cylindrical(GMT) || (gmt_M_is_linear(GMT) && gmt_M_is_geographic(GMT,GMT_IN)))))) return 0;	/* We are only concerned with the global misc projections here */

	/* Global projection need to handle pole path properly */
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Try to include %c pole in polar cap path\n", pole[S->pole+1]);
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "First longitude = %g.  Last longitude = %g\n", S->data[GMT_X][0], S->data[GMT_X][n-1]);
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "West longitude = %g.  East longitude = %g\n", GMT->common.R.wesn[XLO], GMT->common.R.wesn[XHI]);
	type = gmtlib_determine_pole (GMT, S->data[GMT_X], S->data[GMT_Y], n);
	if (abs(type) == 2) {	/* The algorithm only works for clockwise polygon so anything CCW we simply reverse... */
		gmtplot_reverse_polygon (GMT, S);
		type = (type == -2) ? -1 : +1;	/* Now just going clockwise */
	}
	start_lon = GMT->common.R.wesn[XHI];
	stop_lon  = GMT->common.R.wesn[XLO];

	for (k = 1, k0 = 0; k0 == 0 && k < n; k++) {	/* Determine where the perimeter crossing with the west boundary occurs */
		if (k && (GMT->common.R.wesn[XLO]-S->data[GMT_X][k]) >= 0.0 && (GMT->common.R.wesn[XLO]-S->data[GMT_X][k-1]) <= 0.0) k0 = k;
	}
	/* Determine the latitude of that crossing */
	if (k0) {	/* Occurred somewhere along the perimeter between points k0 and k0-1 */
		gmt_M_set_delta_lon (S->data[GMT_X][k0-1], S->data[GMT_X][k0], dx);	/* Handles the 360 jump cases */
		yc = S->data[GMT_Y][k0-1] - (S->data[GMT_Y][k0] - S->data[GMT_Y][k0-1]) * (S->data[GMT_X][k0-1] - GMT->common.R.wesn[XLO]) / dx;
	}
	else	/* Very first point is at the right longitude */
		yc = S->data[GMT_Y][k0];
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Crossing at %g,%g\n", GMT->common.R.wesn[XLO], yc);
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "k at point closest to lon %g is = %d [n = %d]\n", GMT->common.R.wesn[XLO], (int)k0, (int)n);
	/* Then, add path from pole to start longitude, then copy perimeter path, then add path from stop longitude back to pole */
	/* WIth cylindrical projections we may not go all the way to the pole, so we adjust pole_lat based on -R: */
	if (pole_lat < GMT->common.R.wesn[YLO]) pole_lat = GMT->common.R.wesn[YLO];
	if (pole_lat >GMT->common.R.wesn[YHI])  pole_lat = GMT->common.R.wesn[YHI];
	/* 1. Calculate the path from yc to pole: */
	perim_n = gmtlib_lonpath (GMT, start_lon, pole_lat, yc, &x_perim, &y_perim);
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Created path from %g/%g to %g/%g [%d points]\n", start_lon, pole_lat, start_lon, yc, perim_n);
	/* 2. Allocate enough space for new polar cap polygon */
	n_new = 2 * perim_n + n;
	plon = gmt_M_memory (GMT, NULL, n_new, double);
	plat = gmt_M_memory (GMT, NULL, n_new, double);
	/* Start off with the path from the pole to the crossing */
	gmt_M_memcpy (plon, x_perim, perim_n, double);
	gmt_M_memcpy (plat, y_perim, perim_n, double);
	/* Now walk from k0 to the end of polygon, wrapping around if needed */
	m = perim_n;	/* Index of next output point */
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Add perimeter data from k0->n [%d->%d], then 0->k0 [%d]\n", k0, n, k0);
	for (k = k0; k < n; k++, m++) {
		plon[m] = S->data[GMT_X][k];
		plat[m] = S->data[GMT_Y][k];
	}
	for (k = 0; k < k0; k++, m++) {
		plon[m] = S->data[GMT_X][k];
		plat[m] = S->data[GMT_Y][k];
	}
	/* Then add the opposite path to the pole, switching the longitude to stop_lon */
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Add path from %g/%g to %g/%g [%d points]\n", stop_lon, yc, stop_lon, pole_lat, perim_n);
	for (k = perim_n-1; k > 0; k--, m++) {
		plon[m] = stop_lon;
		plat[m] = y_perim[k];
	}
	/* Finally add the duplicate pole at the end of polygon so it is closed */
	plon[m] = stop_lon;
	plat[m++] = pole_lat;
	gmt_M_free (GMT, x_perim);	gmt_M_free (GMT, y_perim);	/* No longer needed */
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "New path has %d points, we allocated %d points\n", m, n_new);
#if 0
	fp = fopen ("shit.txt", "w");
	for (k = 0; k < m; k++) {
		fprintf (fp, "%lg\t%lg\n", plon[k], plat[k]);
	}
	fclose (fp);
#endif
	k = gmtplot_geo_polygon (GMT, plon, plat, m, first, comment);	/* Plot filled polygon [no outline] */
	gmt_M_free (GMT, plon);	gmt_M_free (GMT, plat);			/* No longer needed */
	return (k);	/* Number of points plotted */
}
#endif

GMT_LOCAL bool gmtplot_at_pole (double *lat, uint64_t n) {
	return (lat[0] == lat[n-1] && fabs (lat[0]) == 90.0);
}

GMT_LOCAL uint64_t gmtplot_geo_polygon_segment (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S, bool add_pole, bool first, const char *comment) {
	/* Handles the laying down of polygons suitable for filling only; outlines are done separately later.
	 * Polar caps need special treatment in that we must add a detour to the pole.
	 * That detour will not be drawn, only used for fill. However, due to the insanity that is called GIS,
	 * some user polygons may already have artificial lines drawn to the pole in order to work in GIS.
	 * Thus, if we detect such a pole as part of the line then we do NOT add yet another detour.  */

	uint64_t n = S->n_rows, k;
	double *plon = S->data[GMT_X], *plat = S->data[GMT_Y], t_lat;	/* Default is to plot incoming array as is via plon,plat pointers */
	bool ap = gmtplot_at_pole (plat, n);	/* Is the first and last point exactly at the pole? */
	bool free_memory = false;
	struct GMT_DATASEGMENT_HIDDEN *SH = gmt_get_DS_hidden (S);
	if (ap) plon[n-1] = plon[0];	/* Just enforce the same longitude at the pole point */
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Polar cap: %d\n", (int)add_pole);
	if (add_pole) {	/* Make sure there is not already a detour in the data as given */
		double p_lat = SH->pole * 90.0;	/* Latitude of the pole in question */
		bool need_detour = true;	/* Until proven otherwise we assume we must add a detour */
		if (GMT->common.R.oblique) {	/* Determine if any of the 4 map corners are inside this polygon */
			double X, Y;
			gmt_set_inside_mode (GMT, NULL, GMT_IOO_SPHERICAL);
			gmt_xy_to_geo (GMT, &X, &Y, GMT->current.proj.rect[XLO], GMT->current.proj.rect[YLO]);
			GMT->current.proj.corner[0] = gmt_inonout (GMT, X, Y, S);
			gmt_xy_to_geo (GMT, &X, &Y, GMT->current.proj.rect[XHI], GMT->current.proj.rect[YLO]);
			GMT->current.proj.corner[1] = gmt_inonout (GMT, X, Y, S);
			gmt_xy_to_geo (GMT, &X, &Y, GMT->current.proj.rect[XHI], GMT->current.proj.rect[YHI]);
			GMT->current.proj.corner[2] = gmt_inonout (GMT, X, Y, S);
			gmt_xy_to_geo (GMT, &X, &Y, GMT->current.proj.rect[XLO], GMT->current.proj.rect[YHI]);
			GMT->current.proj.corner[3] = gmt_inonout (GMT, X, Y, S);
			need_detour = false;	/* Trying a different tack for these cases */
		}
		for (k = 0; need_detour && k < S->n_rows; k++) {	/* Check every point */
			if (doubleAlmostEqual (S->data[GMT_Y][k], p_lat)) {	/* Point is exactly at the pole in question */
				/* We want to distinguish between a path that gently touches the pole and one that has a fake straight detour to the pole.
				 * We assume a fake detour will have the same longitudes for this point and the previous and that they are both either +/-180 or 0. */
				if (k && doubleAlmostEqual (S->data[GMT_X][k], S->data[GMT_X][k-1]) && (doubleAlmostEqual (fabs (S->data[GMT_X][k]), 180.0) || gmt_M_is_zero (S->data[GMT_X][k])))
					need_detour = false;	/* Well, what do you know. Probably arcGIS or some other handicapped program */
			}
		}
		if (!need_detour) {	/* Do not add another detour but process via gmt_geo_polarcap_segment to handle jumps in our polygon */
			GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Path already had a detour to the pole, skip adding another detour\n");
			add_pole = false;
			n = gmt_geo_polarcap_segment (GMT, S, &plon, &plat);
			if (plon != S->data[GMT_X])	/* Pointer changed so we allocated, must delete below */
				free_memory = true;
		}
	}
	if (add_pole) {	/* If we get here then a detour will be needed */
		if ((n = gmt_geo_polarcap_segment (GMT, S, &plon, &plat)) == 0) {	/* Not a global map */
			/* Here we must detour to the N or S pole, then resample the path */
			n = S->n_rows + 2;	/* Add new first and last point to connect to the pole */
			plon = gmt_M_memory (GMT, NULL, n, double);	/* This memory must be freed below */
			plat = gmt_M_memory (GMT, NULL, n, double);
			free_memory = true;
			t_lat = SH->pole * 90.0;	/* This is presumably the correct pole, but could fail if just touching the pole */
			if (S->data[GMT_Y][0] * t_lat < 0.0) t_lat = -t_lat;	/* Well, I'll be damned... */
			plat[0] = plat[n-1] = t_lat;
			plon[0] = S->data[GMT_X][0];
			plon[n-1] = S->data[GMT_X][S->n_rows-1];
			gmt_M_memcpy (&plon[1], S->data[GMT_X], S->n_rows, double);
			gmt_M_memcpy (&plat[1], S->data[GMT_Y], S->n_rows, double);
			if (GMT->current.map.path_mode == GMT_RESAMPLE_PATH && (n = gmt_fix_up_path (GMT, &plon, &plat, n, 0.0, 0)) == 0) {
				gmt_M_free (GMT, plon);
				gmt_M_free (GMT, plat);
				return 0;
			}
		}
		else if (plon != S->data[GMT_X])	/* Pointer changed so we allocated, must delete below */
			free_memory = true;
	}
	k = gmtplot_geo_polygon (GMT, plon, plat, n, first, comment);	/* Plot filled polygon [no outline] */
	if (free_memory) {	/* Delete what we allocated */
		gmt_M_free (GMT, plon);
		gmt_M_free (GMT, plat);
	}
	return (k);	/* Number of points plotted */
}

GMT_LOCAL float gmtplot_inch_to_degree_scale (struct GMT_CTRL *GMT, double lon0, double lat0, double azimuth) {
	/* Used to determine delta radius in degrees for thickness of vector lines to be drawn as small or
	 * great circles.  We must convert pen thickness to some sense of spherical degrees.
	 * Determine the map scale at (lon, lat) in a direction normal to the vector and use that
	 * to scale items in inches to spherical degrees.  We repeat this for all locations except
	 * if the projection is perspective, when we pick the center as the map distortion will be the least here.
	 * This scaling is approximate only but needed to convert geovector head lengths to degrees. */

	double tlon, tlat, x0, y0, dx, x1, y1, length;
	float scale;

	length = 0.001 * (GMT->common.R.wesn[YHI] - GMT->common.R.wesn[YLO]);		/* 0.1 percent of latitude extent is fairly small */
	gmt_geo_to_xy (GMT, lon0, lat0, &x0, &y0);					/* Get map position in inches for close point */
	gmtlib_get_point_from_r_az (GMT, lon0, lat0, length, azimuth-90.0, &tlon, &tlat);	/* ANearby arbitrary 2nd point in direction normal to vector at (lon0,lat0) */
	gmt_geo_to_xy (GMT, tlon, tlat, &x1, &y1);					/* Get map position in inches for close point */
	dx = fabs (x1 - x0);
	if (dx > (0.25 * GMT->current.map.half_width)) dx = GMT->current.map.width - dx;
	scale = (float) (length / hypot (dx, y1 - y0));				/* This scales a length in inches to degrees, approximately */
	return (scale);
}

GMT_LOCAL uint64_t gmtplot_great_circle_arc (struct GMT_CTRL *GMT, double *A, double *B, double step, bool longway, double **xp, double **yp) {
	/* Given vectors A and B, return great circle path sampled every step.  Shorest path is selected unless longway is true */
	/* Determine unit vector pole of great circle or use the one given by small circle pole and its opening rot */
	uint64_t k, n;
	double P[3], X[3], R[3][3], R0[3][3], c, w, *xx = NULL, *yy = NULL;

	gmt_cross3v (GMT, A, B, P);	/* Parallel to rotation pole */
	gmt_normalize3v (GMT, P);		/* Rotation pole unit vector */
	c = d_acosd (gmt_dot3v (GMT, A, B));	/* opening angle in degrees */

	if (longway) {	/* Want to go the long way */
		c = 360.0 - c;
		P[0] = -P[0], P[1] = -P[1], P[2] = -P[2];
	}
	if (gmt_M_is_zero (step)) step = GMT->current.map.path_step;	/* Use default map-step if given as 0 */
	n = lrint (ceil (c / step)) + 1;	/* Number of segments needed for smooth curve from A to B inclusive */
	step = D2R * c / (n - 1);	/* Adjust step for exact fit, convert to radians */
	gmt_M_malloc2 (GMT, xx, yy, n, NULL, double);	/* Allocate space for arrays */
	gmtlib_init_rot_matrix (R0, P);			/* Get partial rotation matrix since no actual angle is applied yet */
	for (k = 0; k < n; k++) {	/* March along the arc */
		w = k * step;					/* Opening angle from A to this point X */
		gmt_M_memcpy (R, R0, 9, double);			/* Get a copy of the "0-angle" rotation matrix */
		gmtlib_load_rot_matrix (w, R, P);			/* Build the actual rotation matrix for this angle */
		gmt_matrix_vect_mult (GMT, 3U, R, A, X);			/* Rotate point A towards B and get X */
		gmt_cart_to_geo (GMT, &yy[k], &xx[k], X, true);	/* Get lon/lat of this point along arc */
	}
	*xp = xx;	*yp = yy;
	return (n);
}

GMT_LOCAL uint64_t gmtplot_small_circle_arc (struct GMT_CTRL *GMT, double *A, double step, double P[], double rot, double **xp, double **yp) {
	/* Given vectors A and B, return small circle path sampled every step. */
	/* Use small circle pole P and its opening rot */
	uint64_t k, n;
	double X[3], R[3][3], R0[3][3], w, *xx = NULL, *yy = NULL;

	if (step == 0 || gmt_M_is_zero (step)) step = GMT->current.map.path_step;	/* Use default map-step if given as 0 */
	n = lrint (ceil (fabs (rot) / step)) + 1;	/* Number of segments needed for smooth curve from A to B inclusive */
	step = D2R * rot / (n - 1);	/* Adjust step for exact fit, convert to radians */
	gmt_M_malloc2 (GMT, xx, yy, n, NULL, double);	/* Allocate space for arrays */
	gmtlib_init_rot_matrix (R0, P);			/* Get partial rotation matrix since no actual angle is applied yet */
	for (k = 0; k < n; k++) {	/* March along the arc */
		w = k * step;					/* Opening angle from A to this point X */
		gmt_M_memcpy (R, R0, 9, double);			/* Get a copy of the "0-angle" rotation matrix */
		gmtlib_load_rot_matrix (w, R, P);			/* Build the actual rotation matrix for this angle */
		gmt_matrix_vect_mult (GMT, 3U, R, A, X);			/* Rotate point A towards B and get X */
		gmt_cart_to_geo (GMT, &yy[k], &xx[k], X, true);	/* Get lon/lat of this point along arc */
	}
	*xp = xx;	*yp = yy;
	return (n);
}

GMT_LOCAL double gmtplot_get_local_scale (struct GMT_CTRL *GMT, double lon0, double lat0, double length, double azimuth) {
	/* Determine the local scale at lon0,lat in the direction azimuth using a test distance length in degrees.
	 * The scale returned can be used to convert a map distance in inch to great circle degrees.
	 * This is approximate only. */

	double tlon, tlat, x0, y0, x1, y1;
	gmtlib_get_point_from_r_az (GMT, lon0, lat0, length, azimuth, &tlon, &tlat);	/* Arbitrary 2nd point near (lon0,lat0) in the azimuth direction */
	gmt_geo_to_xy (GMT, lon0, lat0, &x0, &y0);	/* Get map position in inches for (lon0,lat0) */
	gmt_geo_to_xy (GMT, tlon, tlat, &x1, &y1);	/* Get map position in inches for close point */
	return (length / hypot (x1 - x0, y1 - y0));	/* This scales a length in inches to degrees, approximately */
}

GMT_LOCAL void gmtplot_circle_pen_poly (struct GMT_CTRL *GMT, double *A, double *B, bool longway, double rot, struct GMT_PEN *pen, struct GMT_SYMBOL *S, struct GMT_CIRCLE *C, double scale) {
	/* Given vectors A and B, return a small circle polygon path sampled every step that approximates a pen of given width
	 * drawn on the map.  Use small circle pole and its opening rot */
	uint64_t k, n, n2;
	double Ai[3], Ao[3], Px[3], X[3], R[3][3], R0[3][3], w, step;
	struct GMT_DATASEGMENT *L = NULL;
	gmt_M_unused(B);

	gmt_cross3v (GMT, A, C->P, Px);	/* Px is Pole to plane through A and P  */
	gmt_normalize3v (GMT, Px);			/* Rotation pole unit vector */
	/* Rotate A back/fore by rotation angle of +/- pen halfwidth about Px */
	w = 0.5 * scale * pen->width * GMT->session.u2u[GMT_PT][GMT_INCH] * S->v.scale;		/* Half-width of pen in degrees */
	gmt_make_rot_matrix2 (GMT, Px, +w, R);		/* Rotation of rot_v degrees about pole P */
	gmt_matrix_vect_mult (GMT, 3U, R, A, Ai);	/* Get Ai = R * A */
	gmt_make_rot_matrix2 (GMT, Px, -w, R);		/* Rotation of rot_v degrees about pole P */
	gmt_matrix_vect_mult (GMT, 3U, R, A, Ao);	/* Get Ao = R * A */
	if (longway) rot = 360.0 - rot;

	step = GMT->current.map.path_step;		/* Use default map-step if given as 0 */
	n = lrint (ceil (fabs (rot) / step)) + 1;	/* Number of segments needed for smooth curve from A to B inclusive */
	step = D2R * rot / (n - 1);			/* Adjust step for exact fit, convert to radians */
	L = GMT_Alloc_Segment (GMT->parent, GMT_NO_STRINGS, 2*n+1, 2, NULL, NULL);	/* Allocate polygon to draw filled path */
	n2 = 2*n-1;
	gmtlib_init_rot_matrix (R0, C->P);			/* Get partial rotation matrix since no actual angle is applied yet */
	for (k = 0; k < n; k++) {	/* March along the arc */
		w = k * step;					/* Opening angle from A to this point X */
		gmt_M_memcpy (R, R0, 9U, double);			/* Get a copy of the "0-angle" rotation matrix */
		gmtlib_load_rot_matrix (w, R, C->P);			/* Build the actual rotation matrix for this angle */
		gmt_matrix_vect_mult (GMT, 3U, R, Ai, X);		/* Rotate point Ai towards B and get X */
		gmt_cart_to_geo (GMT, &L->data[GMT_Y][k], &L->data[GMT_X][k], X, true);	/* Get lon/lat of this point along arc */
		gmt_matrix_vect_mult (GMT, 3U, R, Ao, X);		/* Rotate point Ai towards B and get X */
		gmt_cart_to_geo (GMT, &L->data[GMT_Y][n2-k], &L->data[GMT_X][n2-k], X, true);	/* Get lon/lat of this point along arc */
	}
	L->data[GMT_X][2*n] = L->data[GMT_X][0];	/* Explicitly close the polygon */
	L->data[GMT_Y][2*n] = L->data[GMT_Y][0];

	/* Plot pen as a closed filled polygon without outline */
	PSL_command (GMT->PSL, "V\n");
	gmt_setpen (GMT, pen);		/* Set pen width just so later setpen's will work */
	PSL_setfill (GMT->PSL, pen->rgb, 0);	/* Fill, no outline */
	gmt_geo_polygons (GMT, L);	/* "Draw" the line */
	PSL_command (GMT->PSL, "U\n");

	gmt_free_segment (GMT, &L);
}

GMT_LOCAL void gmtplot_gcircle_sub (struct GMT_CTRL *GMT, double lon0, double lat0, double azimuth, double length, struct GMT_SYMBOL *S, struct GMT_CIRCLE *C) {
	/* We must determine points A and B, whose great-circle connector is the arc we seek to draw */

	int justify = PSL_vec_justify (S->v.status);	/* Return justification as 0-3 */
	double x, y;
	gmt_M_memset (C, 1, struct GMT_CIRCLE);	/* Set all to zero */

	if (S->v.status & PSL_VEC_JUST_S) { /* Was given coordinates of A and B */
		justify = 3;
	}
	switch (justify) {	/* A and B depends on chosen justification */
		case 0: /* Was given coordinates of A; determine B */
			C->lon[0] = lon0;	C->lat[0] = lat0;
			gmt_geo_to_cart (GMT, C->lat[0], C->lon[0], C->A, true);
			C->r0 = C->r = length / GMT->current.proj.DIST_KM_PR_DEG;	/* Arch length in spherical degrees */
			if (C->r > 180.0) {C->longway = true; C->r -= 180.0;}	/* Temporarily adjust if arcs > 180 degrees are chosen */
			gmtlib_get_point_from_r_az (GMT, C->lon[0], C->lat[0], C->r, azimuth, &C->lon[1], &C->lat[1]);
			if (C->longway) C->lon[1] += 180.0, C->lat[1] = -C->lat[1];	/* Undo adjustment */
			gmt_geo_to_cart (GMT, C->lat[1], C->lon[1], C->B, true);	/* Get B */
			break;
		case 1: /* Was given coordinates of halfway point; determine A and B */
			C->r0 = C->r = length / GMT->current.proj.DIST_KM_PR_DEG;	/* Arch length in spherical degrees */
			if (C->r > 180.0) C->longway = true;	/* Temporarily adjust if arcs > 180 degrees are chosen */
			gmtlib_get_point_from_r_az (GMT, lon0, lat0, 0.5*C->r, azimuth, &C->lon[1], &C->lat[1]);
			gmt_geo_to_cart (GMT, C->lat[1], C->lon[1], C->B, true);	/* Get B */
			gmtlib_get_point_from_r_az (GMT, lon0, lat0, 0.5*C->r, azimuth+180.0, &x, &y);
			C->lon[0] = x;	C->lat[0] = y;	/* Replace the original A point */
			gmt_geo_to_cart (GMT, C->lat[0], C->lon[0], C->A, true);	/* Get A */
			break;
		case 2: /* Was given coordinates of B point; determine A */
			C->lon[0] = lon0;	C->lat[0] = lat0;
			gmt_geo_to_cart (GMT, C->lat[0], C->lon[0], C->B, true);
			C->r0 = C->r = length / GMT->current.proj.DIST_KM_PR_DEG;	/* Arch length in spherical degrees */
			if (C->r > 180.0) {C->longway = true; C->r -= 180.0;}	/* Temporarily adjust if arcs > 180 degrees are chosen */
			gmtlib_get_point_from_r_az (GMT, C->lon[0], C->lat[0], C->r, azimuth+180.0, &C->lon[1], &C->lat[1]);
			if (C->longway) C->lon[1] += 180.0, C->lat[1] = -C->lat[1];	/* Undo adjustment */
			gmt_geo_to_cart (GMT, C->lat[1], C->lon[1], C->A, true);	/* Get A */
			gmt_M_double_swap (C->lon[0], C->lon[1]);	gmt_M_double_swap (C->lat[0], C->lat[1]);	/* Now A is first and B is second */
			break;
		case 3: /* Was given coordinates of B instead of azimuth and length; can never be longway */
			C->lon[0] = lon0;	C->lat[0] = lat0;
			C->lat[1] = length;	C->lon[1] = azimuth;
			gmt_geo_to_cart (GMT, C->lat[0], C->lon[0], C->A, true);	/* Get A */
			gmt_geo_to_cart (GMT, C->lat[1], C->lon[1], C->B, true);	/* Get B */
			C->r0 = C->r = d_acosd (gmt_dot3v (GMT, C->A, C->B));		/* Arc length in degrees */
			break;
	}
	gmt_cross3v (GMT, C->A, C->B, C->P);	/* Rotation pole */
	gmt_normalize3v (GMT, C->P);		/* Rotation pole unit vector */

	if (C->longway) {	/* Want to go the long way */
		C->P[0] = -C->P[0], C->P[1] = -C->P[1], C->P[2] = -C->P[2];
	}
}

GMT_LOCAL void gmtplot_scircle_sub (struct GMT_CTRL *GMT, double lon0, double lat0, double angle_1, double angle_2, struct GMT_SYMBOL *S, struct GMT_CIRCLE *C) {
	/* We must determine points A and B, whose small-circle connector about pole P is the arc we seek to draw */

	int justify = PSL_vec_justify (S->v.status);	/* Return justification as 0-3 */
	double R[3][3], M[3];
	gmt_M_memset (C, 1, struct GMT_CIRCLE);	/* Set all to zero */
	/* Requires the rotation matrix for pole S->v.pole */

	/* Here angle_1, angle_2 are not necessarily that, depending on S->v.status:
	 * S->v.pole & PSL_VEC_ANGLES : angle_1 is opening angle1 and angle_2 is opening angle2 about the pole.
	 * Otherwise:	angle_2 is the length of the arc in km */
	gmt_geo_to_cart (GMT, lat0, lon0, M, true);	/* Given input point */
	gmt_geo_to_cart (GMT, S->v.pole[GMT_Y], S->v.pole[GMT_X], C->P, true);
	C->colat = d_acosd (gmt_dot3v (GMT, M, C->P));	/* Colatitude of input point relative to pole, in degrees */

	if (S->v.status & PSL_VEC_ANGLES) {
		/* Was given the two opening angles; compute A and B accordingly */
		gmt_make_rot_matrix (GMT, S->v.pole[GMT_X], S->v.pole[GMT_Y], angle_1, R);
		gmt_matrix_vect_mult (GMT, 3U, R, M, C->A);	/* Get A */
		gmt_cart_to_geo (GMT, &C->lat[0], &C->lon[0], C->A, true);
		gmt_make_rot_matrix (GMT, S->v.pole[GMT_X], S->v.pole[GMT_Y], angle_2, R);
		gmt_matrix_vect_mult (GMT, 3U, R, M, C->B);	/* Get B */
		gmt_cart_to_geo (GMT, &C->lat[1], &C->lon[1], C->B, true);
		C->rot = C->r0 = C->r = angle_2 - angle_1;
	}
	else {
		/* Here A, B, or midpoint was given, + the arc length via angle_2 */
		/* Determine co-latitude for this point */
		C->rot = C->r0 = C->r = (angle_1 / GMT->current.proj.DIST_KM_PR_DEG) / sind (C->colat);	/* Opening angle in spherical degrees */
		switch (justify) {	/* A and B depends on chosen justification */
			case 0: /* Was given coordinates of A; determine B */
				gmt_M_memcpy (C->A, M, 3, double);
				C->lon[0] = lon0;	C->lat[0] = lat0;
				if (C->r > 180.0) {C->longway = true; C->r -= 180.0;}	/* Temporarily adjust if arcs > 180 degrees are chosen */
				/* Rotate A by C->r0 degrees about P to get B */
				gmt_make_rot_matrix (GMT, S->v.pole[GMT_X], S->v.pole[GMT_Y], C->r0, R);
				gmt_matrix_vect_mult (GMT, 3U, R, C->A, C->B);	/* Get B */
				gmt_cart_to_geo (GMT, &C->lat[1], &C->lon[1], C->B, true);
				break;
			case 1: /* Was given coordinates of halfway point; determine A and B */
				if (C->r > 180.0) C->longway = true;	/* Temporarily adjust if arcs > 180 degrees are chosen */
				/* Rotate M by -C->r0/2 degrees about P to get A */
				gmt_make_rot_matrix (GMT, S->v.pole[GMT_X], S->v.pole[GMT_Y], -0.5 * C->r0, R);
				gmt_matrix_vect_mult (GMT, 3U, R, M, C->A);	/* Get A */
				gmt_cart_to_geo (GMT, &C->lat[0], &C->lon[0], C->A, true);
				/* Rotate M by +C->r0/2 degrees about P to get B */
				gmt_make_rot_matrix (GMT, S->v.pole[GMT_X], S->v.pole[GMT_Y], +0.5 * C->r0, R);
				gmt_matrix_vect_mult (GMT, 3U, R, M, C->B);	/* Get B */
				gmt_cart_to_geo (GMT, &C->lat[1], &C->lon[1], C->B, true);
				break;
			case 2: /* Was given coordinates of B point; determine A */
				gmt_M_memcpy (C->B, M, 3, double);
				C->lon[1] = lon0;	C->lat[1] = lat0;
				if (C->r > 180.0) {C->longway = true; C->r -= 180.0;}	/* Temporarily adjust if arcs > 180 degrees are chosen */
				/* Rotate B by -C->r0 degrees about P to get A */
				gmt_make_rot_matrix (GMT, S->v.pole[GMT_X], S->v.pole[GMT_Y], -C->r0, R);
				gmt_matrix_vect_mult (GMT, 3U, R, C->B, C->A);	/* Get A */
				gmt_cart_to_geo (GMT, &C->lat[0], &C->lon[0], C->A, true);
				break;
		}
	}
}

GMT_LOCAL double gmtplot_smallcircle_az (struct GMT_CTRL *GMT, double P[], struct GMT_SYMBOL *S) {
	/* Compute the azimuth at P along small circle given pole */
	double R[3][3], X[3], xlon1, xlat1, xlon2, xlat2, az;
	/* Make rotation matrix for a +0.005 degree rotation */
	gmt_make_rot_matrix (GMT, S->v.pole[GMT_X], S->v.pole[GMT_Y], -0.005, R);
	gmt_matrix_vect_mult (GMT, 3U, R, P, X);	/* Get point really close to P along small circle */
	gmt_cart_to_geo (GMT, &xlat1, &xlon1, X, true);	/* Get coordinates of X */
	gmt_make_rot_matrix (GMT, S->v.pole[GMT_X], S->v.pole[GMT_Y], +0.005, R);
	gmt_matrix_vect_mult (GMT, 3U, R, P, X);	/* Get point really close to P along small circle */
	gmt_cart_to_geo (GMT, &xlat2, &xlon2, X, true);	/* Get coordinates of X */
	az = gmt_az_backaz (GMT, xlon1, xlat1, xlon2, xlat2, false);	/* Compute the azimuth from P to X at A */
	return (az);
}

GMT_LOCAL void gmtplot_plot_vector_head (struct GMT_CTRL *GMT, double *xp, double *yp, uint64_t n, struct GMT_SYMBOL *S) {
	/* PW: Plots the polygon that makes up a vector head.  Because sometimes these head stick
	 * across a periodic boundary we must check if that is the case and plot the two parts separately.
	 * When that is the case we cannot draw the outline of the two new polygons since we wish to show
	 * the heads as "clipped" by the boundary; hence all the rigamorole below. */
	uint64_t start = 0, nin = 0;
	unsigned int n_use, *pin = NULL;	/* Copy of the pen moves */
	unsigned int cap = GMT->PSL->internal.line_cap;
	double *xin = NULL, *yin = NULL;	/* Temp vector with possibly clipped x,y line returned by gmt_geo_to_xy_line */
	if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, xp, yp, n)) == 0) return;	/* All outside, or use plot.x|y array */
	PSL_setlinecap (GMT->PSL, PSL_SQUARE_CAP);	/* In case there are clipped heads and we want to do the best we can with the lines */
	n = GMT->current.plot.n;	/* Possibly fewer points */
	if (PSL_vec_outline (S->v.status)) {
		bool close = (GMT->current.plot.n > 2 && gmt_polygon_is_open (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.n));
		nin = n;
		PSL_command (GMT->PSL, "O0\n");	/* Temporary turn off outline; must draw outline separately when head is split */
		if (close) nin++;
		xin = gmt_M_memory (GMT, NULL, nin, double);
		yin = gmt_M_memory (GMT, NULL, nin, double);
		pin = gmt_M_memory (GMT, NULL, nin, unsigned int);
		gmt_M_memcpy (xin, GMT->current.plot.x, n, double);
		gmt_M_memcpy (yin, GMT->current.plot.y, n, double);
		gmt_M_memcpy (pin, GMT->current.plot.pen, n, unsigned int);
		if (close) {	/* Explicitly close the polygon outline */
			xin[n] = xin[0];
			yin[n] = yin[0];
		}
	}

	if ((*GMT->current.map.will_it_wrap) (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.n, &start)) {	/* Polygon does indeed wrap */
		double *xtmp = NULL, *ytmp = NULL;	/* Temp vector for map truncating */
		GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Vector head polygon will wrap at periodic boundary and will be split into two sections\n");
		xtmp = gmt_M_memory (GMT, NULL, n, double);
		ytmp = gmt_M_memory (GMT, NULL, n, double);
		gmt_M_memcpy (xtmp, GMT->current.plot.x, n, double);
		gmt_M_memcpy (ytmp, GMT->current.plot.y, n, double);
		/* First truncate against left border */
		GMT->current.plot.n = gmt_map_truncate (GMT, xtmp, ytmp, n, start, -1);
		n_use = (unsigned int)gmt_compact_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.n, false, 0);
		PSL_beginclipping (GMT->PSL, GMT->current.plot.x, GMT->current.plot.y, (int)n_use, GMT->session.no_rgb, 3);
		PSL_plotpolygon (GMT->PSL, GMT->current.plot.x, GMT->current.plot.y, (int)n_use);
		if (PSL_vec_outline (S->v.status)) gmt_plot_line (GMT, xin, yin, pin, nin, PSL_LINEAR);
		PSL_endclipping (GMT->PSL, 1);
		/* Then truncate against right border */
		GMT->current.plot.n = gmt_map_truncate (GMT, xtmp, ytmp, n, start, +1);
		n_use = (unsigned int)gmt_compact_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.n, false, 0);
		PSL_beginclipping (GMT->PSL, GMT->current.plot.x, GMT->current.plot.y, (int)n_use, GMT->session.no_rgb, 3);
		PSL_plotpolygon (GMT->PSL, GMT->current.plot.x, GMT->current.plot.y, (int)n_use);
		if (PSL_vec_outline (S->v.status)) gmt_plot_line (GMT, xin, yin, pin, nin, PSL_LINEAR);
		PSL_endclipping (GMT->PSL, 1);
		gmt_M_free (GMT, xtmp);		gmt_M_free (GMT, ytmp);
	}
	else {	/* No wrapping but may be clipped */
		PSL_beginclipping (GMT->PSL, GMT->current.plot.x, GMT->current.plot.y, (int)GMT->current.plot.n, GMT->session.no_rgb, 3);
		PSL_plotpolygon (GMT->PSL, GMT->current.plot.x, GMT->current.plot.y, (int)GMT->current.plot.n);
		if (PSL_vec_outline (S->v.status)) gmt_plot_line (GMT, xin, yin, pin, nin, PSL_LINEAR);
		PSL_endclipping (GMT->PSL, 1);
	}
	PSL_setlinecap (GMT->PSL, cap);
	if (PSL_vec_outline (S->v.status)) {	/* Turn on outline again and free temp memory */
		PSL_command (GMT->PSL, "O1\n");
		gmt_M_free (GMT, xin);
		gmt_M_free (GMT, yin);
		gmt_M_free (GMT, pin);
	}
}

GMT_LOCAL unsigned int gmtplot_geo_vector_smallcircle (struct GMT_CTRL *GMT, double lon0, double lat0, double azimuth, double length, struct GMT_PEN *ppen, struct GMT_SYMBOL *S) {
	/* Draws a small-circle vector with or without heads, etc. There are some complications to consider:
	 * When there are no heads it is simple.  If +n is active we may shrink the line thickness.
	 * With heads there are these cases:
	 * Head length is longer than 90% of the vector length.  We then skip the head and return 1
	 * +n<norm> is in effect.  We shrink vector pen and head length.  Still, the shrunk head
	 * may be longer than 90% of the vector length.  We then shrink head (not pen) further and return 2
	*/

	uint64_t n1, n2, n, add;
	int heads, side[2], outline = 0;
	unsigned int warn = 0;
	size_t n_alloc;
	bool perspective, pure;
	double P[3], Pa[3], Ax[3], Bx[3], Ax2[3], Bx2[3], R[3][3], h_length_limit, max_length;
	double dr[2] = {0.0, 0.0}, az[2] = {0.0, 0.0}, oaz[2] = {0.0, 0.0}, scl[2];
	double da = 0.0, dshift[2] = {0.0, 0.0}, s = 1.0, s1 = 1.0, s2 = 1.0, olon[2], olat[2], head_length, arc_width, n_az, arc;
	double rot[2] = {0.0, 0.0}, rot_v[2] = {0.0, 0.0};
	double *xp = NULL, *yp = NULL, *xp2 = NULL, *yp2 = NULL;
	double *rgb = S->v.fill.rgb;
	struct GMT_CIRCLE C;

#if 0
	/* We must determine points A and B, whose great-circle connector is the arc we seek to draw */
	justify = PSL_vec_justify (S->v.status);	/* Return justification as 0-3 */
#endif

	gmtplot_scircle_sub (GMT, lon0, lat0, azimuth, length, S, &C);
	perspective = gmt_M_is_perspective (GMT);

	/* Here we have the endpoints A and B of the great (or small) circle arc */

	/* If shrink-option (+n) is active we may have to scale down head attributes and pen width */
	/* If shrink-option (+n) is active we may have to scale down head attributes and pen width */
	heads = PSL_vec_head (S->v.status);	/* Return head selection as 0-3 */
	pure = (S->v.v_norm == -1.0f);	/* True if no shrinking has been specified */
	h_length_limit = (1.0 - S->v.v_stem) * C.r0;	/* Max length of arrow in degrees to ensure the stem is still showing */
	if (heads == 3) h_length_limit *= 0.5;		/* Split this length between the two heads */
	if (heads && !pure) {	/* Need to determine head length in degrees */
		az[0] = gmtplot_smallcircle_az (GMT, C.A, S);	/* Compute the azimuth from A to B at A along small circle */
		scl[0] = (perspective) ? S->v.scale : gmtplot_get_local_scale (GMT, C.lon[0], C.lat[0], 0.001 * C.r, az[0]);	/* Get local deg/inch scale at A in az[0] direction */
		dr[0] = scl[0] * S->size_x;	/* This is arrow head length in degrees, approximately */
		az[1] = -gmtplot_smallcircle_az (GMT, C.B, S);	/* Compute the azimuth from B to A at B along small circle */
		scl[1] = (perspective) ? S->v.scale : gmtplot_get_local_scale (GMT, C.lon[1], C.lat[1], 0.01 * C.r, az[1]);	/* Get local deg/inch scale */
		dr[1] = scl[1] * S->size_x;	/* This is arrow head length in degrees, approximately, adjusted for ~pen thickness to ensure no gap between head and line */
		max_length = MAX (dr[0], dr[1]);
		if (max_length > h_length_limit) {
			s2 = h_length_limit / max_length;
			warn = 2;
		}
	}

	/* Might have to shrink things */
	if (C.r0 < S->v.v_norm)
		s1 = C.r0 / S->v.v_norm;
	if (s1 < s2) {	/* Pick the smallest scale required for head shrinking */
		s = s1;
		warn = 0;	/* When requesting +n there is no warning unless s2 is the smaller scale */
	}
	else
		s = s2;
	head_length = s * S->size_x;
	arc_width   = s1 * S->v.v_width;	/* Use scale s1 for pen shrinking */

	gmt_M_memcpy (olon, C.lon, 2, double);	gmt_M_memcpy (olat, C.lat, 2, double);	/* Keep copy of original coordinates */

	/* When only one side of a vector head is requested (side = -1/+1) there are complications that leads to some
	 * extra work: Since we are clipping the head polygon, the head outline pen is effectively half that of the
	 * vector.  Thus, there will be an offset of 1/2 penwidth at the end of the vector line and the start of the
	 * back-end of the half vector head.  We adjust this changing the colatitude from the pole by the equivalent
	 * distance of 1/2 pen width away from the side with the half arrowhead.  This makes the outline of the head
	 * align with the vector line. */
	for (n = 0; n < 2; n++) {
		side[n] = PSL_vec_side (S->v.status, n);	/* Return side selection as 0,-1,+1 for this head */
		dshift[n] = (side[n]) ? 0.5 * arc_width : 0.0;	/* Half-width of arc thickness if side != 0 */
	}
	if (heads & 1) {	/* Placing head at A means we must shorten the arc and use Ax instead of A */
		az[0] = gmtplot_smallcircle_az (GMT, C.A, S);	/* Compute the azimuth from A to B at A along small circle */
		scl[0] = (perspective) ? S->v.scale : gmtplot_get_local_scale (GMT, C.lon[0], C.lat[0], 0.001 * C.r, az[0]);	/* Get local deg/inch scale at A in az[0] direction */
		dr[0] = scl[0] * (head_length - 1.1*dshift[0]);	/* This is arrow head length in degrees, approximately, adjusted for ~pen thickness to ensure no gap between head and line (the 1.1 slop) */
		if (pure && dr[0] > h_length_limit) {	/* Head length too long, refuse to plot it */
			gmt_M_memcpy (Ax, C.A, 3, double);	/* No need to shorten arc at beginning */
			heads -= 1;	/* Not purse this head any further */
			warn = 1;
		}
		else {
			dr[0] /= sind (C.colat);	/* Scale dr[0] to opening angle degrees given colatitude */
			/* Determine mid-back point of arrow head by rotating A back by chosen fraction rot of dr[0] */
			rot_v[0] = 0.5 * dr[0] * (1.95 - S->v.v_shape);	/* 1.95 instead of 2 to allow for slop */
			gmt_make_rot_matrix2 (GMT, C.P, rot_v[0], R);	/* Rotation of rot_v[0] degrees about pole P */
			gmt_matrix_vect_mult (GMT, 3U, R, C.A, Ax);	/* Get Ax = R * A*/
			C.rot -= rot_v[0];	/* Shorten full arc by the same amount */
			dr[0] = scl[0] * head_length;	/* This is arrow head length in great-circle degrees, approximately, without any pen-width compensation */
		}
	}
	else
		gmt_M_memcpy (Ax, C.A, 3, double);	/* No need to shorten arc at beginning */

	if (heads & 2) { /* Place arrow head at B */
		az[1] = -gmtplot_smallcircle_az (GMT, C.B, S);	/* Compute the azimuth from B to A at B along small circle */
		scl[1] = (perspective) ? S->v.scale : gmtplot_get_local_scale (GMT, C.lon[1], C.lat[1], 0.01 * C.r, az[1]);	/* Get local deg/inch scale */
		dr[1] = S->v.scale * (head_length - 1.1*dshift[1]);	/* This is arrow head length in degrees, approximately, adjusted for ~pen thickness to ensure no gap between head and line */
		if (pure && dr[1] > h_length_limit) {	/* Head length too long, refuse to plot it */
			gmt_M_memcpy (Bx, C.B, 3, double);	/* No need to shorten arc at end */
			heads -= 1;	/* Not purse this head any further */
			warn = 1;
		}
		else {
			dr[1] /= sind (C.colat);	/* Scale dr[1] to opening angle degrees given colatitude */
			/* Determine mid-back point of arrow head by rotating B back by chosen fraction rot of dr[1] */
			rot_v[1] = 0.5 * dr[1] * (1.95 - S->v.v_shape);	/* 1.95 instead of 2 to allow for slop */
			gmt_make_rot_matrix2 (GMT, C.P, -rot_v[1], R);	/* Rotation of -rot_v[1] degrees about pole P */
			gmt_matrix_vect_mult (GMT, 3U, R, C.B, Bx);	/* Get Bx = R * B*/
			C.rot -= rot_v[1];	/* Shorten full arc by the same amount */
			dr[1] = scl[1] * head_length;	/* This is arrow head length in great-circle degrees, approximately, without any pen-width compensation */
		}
	}
	else
		gmt_M_memcpy (Bx, C.B, 3, double);	/* No need to shorten arc at end */

	gmt_M_memcpy (oaz, az, 2, double);	/* Keep copy of original azimuths */

	/* Get array of lon,lat points that defines the arc */

	if (side[0] || side[1]) {	/* Must adjust the distance from pole to A, B by 1/2 the pen width */
		double Scl, Off, xlon, xlat, tlon, tlat;
		gmt_cart_to_geo (GMT, &xlat, &xlon, Ax, true);
		n_az = gmt_az_backaz (GMT, xlon, xlat, S->v.pole[GMT_X], S->v.pole[GMT_Y], false);	/* Compute the azimuth from Ax to P at Ax along great circle */
		if (side[0] == +1) n_az += 180.0;	/* Might be for side == +1, check */
		Scl = (perspective) ? S->v.scale : gmtplot_get_local_scale (GMT, xlon, xlat, 0.01, n_az);	/* Get deg/inch scale at A perpendicular to arc */
		Off = Scl * dshift[0];	/* Offset in degrees due to 1/2 pen thickness */
		gmtlib_get_point_from_r_az (GMT, xlon, xlat, Off, n_az, &tlon, &tlat);	/* Adjusted Ax */
		gmt_geo_to_cart (GMT, tlat, tlon, Ax2, true);
		gmt_make_rot_matrix2 (GMT, C.P, C.rot, R);		/* Rotation of C->rot degrees about pole P */
		gmt_matrix_vect_mult (GMT, 3U, R, Ax2, Bx2);		/* Get revised Bx = R * Ax */
	}
	else {
		gmt_M_memcpy (Ax2, Ax, 3, double);	/* No need to shorten arc at end */
		gmt_M_memcpy (Bx2, Bx, 3, double);	/* No need to shorten arc at end */
	}

	gmtplot_circle_pen_poly (GMT, Ax2, Bx2, false, C.rot, ppen, S, &C, s1);

	if (!heads) return (warn);	/* All done */

	PSL_command (GMT->PSL, "V\n");
	/* Get half-angle at head and possibly change pen */
	da = 0.5 * S->v.v_angle;	/* Half-opening angle at arrow head */
	if ((S->v.status & PSL_VEC_OUTLINE) == 0)
		PSL_command (GMT->PSL, "O0\n");	/* Turn off outline */
	else {
		gmt_setpen (GMT, &S->v.pen);
		outline = 1;
	}
	if ((S->v.status & PSL_VEC_FILL) == 0)
		PSL_command (GMT->PSL, "FQ\n");	/* Turn off vector head fill */
	else
		PSL_setfill (GMT->PSL, rgb, outline);

	if (heads & 1) { /* Place arrow head at A */
		if (C.longway) az[0] += 180.0;
		rot[0] = dr[0] / sind (C.colat);	/* Small circle rotation angle to draw arrow outline */
		if (side[0] != +1) {	/* Want to draw left side of arrow */
			gmt_make_rot_matrix2 (GMT, C.A, da, R);		/* Rotation of da degrees about A */
			gmt_matrix_vect_mult (GMT, 3U, R, C.P, Pa);	/* Rotate pole C.P to inner arc pole location Pa */
			gmt_make_rot_matrix2 (GMT, Pa, rot[0], R);	/* Rotation of rot[0] degrees about Pa */
			gmt_matrix_vect_mult (GMT, 3U, R, C.A, P);	/* Rotate A to inner arc end point P */
			arc = rot[0];
		}
		else {
			gmt_M_memcpy (P, Ax, 3, double);		/* Start from (possibly adjusted) mid point instead... */
			gmt_M_memcpy (Pa, C.P, 3, double);	/* ...and use circle pole */
			arc = rot_v[0];
		}
		n1 = (int)gmtplot_small_circle_arc (GMT, P, 0.0, Pa, -arc, &xp, &yp);	/* Compute small circle arc from P to A */
		if (side[0] != -1) {	/* Want to draw right side of arrow */
			gmt_make_rot_matrix2 (GMT, C.A, -da, R);	/* Rotation of -da degrees about A */
			gmt_matrix_vect_mult (GMT, 3U, R, C.P, Pa);	/* Rotate pole C.P to outer arc pole location Pa */
			gmt_make_rot_matrix2 (GMT, Pa, rot[0], R);	/* Rotation of rot[0] degrees about Pa */
			gmt_matrix_vect_mult (GMT, 3U, R, C.A, P);	/* Rotate A to outer arc end point P */
			arc = rot[0];
		}
		else {
			gmt_M_memcpy (P, Ax, 3, double);		/* Start from (adjusted) mid point instead... */
			gmt_M_memcpy (Pa, C.P, 3, double);	/* ...and use circle pole */
			arc = rot_v[0];
		}
		n2 = (unsigned int)gmtplot_small_circle_arc (GMT, C.A, 0.0, Pa, arc, &xp2, &yp2);	/* Compute great circle arc from A to P */
		add = (side[0] == 0) ? 1 : 0;	/* Need to add mid point explicitly */
		n_alloc = n = n1 + n2 + add;
		gmt_M_malloc2 (GMT, xp, yp, 0U, &n_alloc, double);	/* Allocate space for total path */
		gmt_M_memcpy (&xp[n1], xp2, n2, double);
		gmt_M_memcpy (&yp[n1], yp2, n2, double);
		if (add) {	/* Mid point of arrow */
			gmt_cart_to_geo (GMT, &yp[n-1], &xp[n-1], Ax, true);	/* Add geo coordinates of this new back mid point for arc */
		}
		gmtplot_plot_vector_head (GMT, xp, yp, n, S);
		gmt_M_free (GMT, xp);	gmt_M_free (GMT, yp);
		gmt_M_free (GMT, xp2);	gmt_M_free (GMT, yp2);
	}
	if (heads & 2) { /* Place arrow head at B */
		if (C.longway) az[1] += 180.0;
		rot[1] = dr[1] / sind (C.colat);	/* Small circle rotation angle to draw arrow outline */
		if (side[1] != +1) {	/* Want to draw left side of arrow */
			gmt_make_rot_matrix2 (GMT, C.B, -da, R);	/* Rotation of -da degrees about B */
			gmt_matrix_vect_mult (GMT, 3U, R, C.P, Pa);	/* Rotate pole C.P to inner arc pole location Pa */
			gmt_make_rot_matrix2 (GMT, Pa, -rot[1], R);	/* Rotation of -rot[1] degrees about Pa */
			gmt_matrix_vect_mult (GMT, 3U, R, C.B, P);	/* Rotate B to inner arc end point P */
			arc = rot[1];
		}
		else {
			gmt_M_memcpy (P, Bx, 3, double);		/* Start from (adjusted) mid point instead... */
			gmt_M_memcpy (Pa, C.P, 3, double);	/* ...and use circle pole */
			arc = rot_v[1];
		}
		n1 = gmtplot_small_circle_arc (GMT, P, 0.0, Pa, arc, &xp, &yp);	/* Compute small circle arc from P to B */
		if (side[1] != -1) {	/* Want to draw right side of arrow */
			gmt_make_rot_matrix2 (GMT, C.B, da, R);		/* Rotation of da degrees about B */
			gmt_matrix_vect_mult (GMT, 3U, R, C.P, Pa);	/* Rotate pole C.P to outer arc pole location Pa */
			gmt_make_rot_matrix2 (GMT, Pa, -rot[1], R);	/* Rotation of -rot[1] degrees about Pa */
			gmt_matrix_vect_mult (GMT, 3U, R, C.B, P);	/* Rotate B to outer arc end point P */
			arc = rot[1];
		}
		else {
			gmt_M_memcpy (P, Bx, 3, double);		/* Start from (adjusted) mid point instead */
			gmt_M_memcpy (Pa, C.P, 3, double);	/* ...and use circle pole */
			arc = rot_v[1];
		}
		n2 = (unsigned int)gmtplot_small_circle_arc (GMT, C.B, 0.0, Pa, -arc, &xp2, &yp2);	/* Compute small circle arc from B to P */
		add = (side[1] == 0) ? 1 : 0;	/* Need to add mid point explicitly */
		n_alloc = n = n1 + n2 + add;
		gmt_M_malloc2 (GMT, xp, yp, 0U, &n_alloc, double);	/* Allocate space for total path */
		gmt_M_memcpy (&xp[n1], xp2, n2, double);
		gmt_M_memcpy (&yp[n1], yp2, n2, double);
		if (add) {	/* Mid point of arrow */
			gmt_cart_to_geo (GMT, &yp[n-1], &xp[n-1], Bx, true);	/* Add geo coordinates of this new back-mid point for arc */
		}
		gmtplot_plot_vector_head (GMT, xp, yp, n, S);
		gmt_M_free (GMT, xp);	gmt_M_free (GMT, yp);
		gmt_M_free (GMT, xp2);	gmt_M_free (GMT, yp2);
	}
	PSL_command (GMT->PSL, "U\n");
	return (warn);
}

GMT_LOCAL unsigned int gmtplot_geo_vector_greatcircle (struct GMT_CTRL *GMT, double lon0, double lat0, double azimuth, double length, struct GMT_PEN *ppen, struct GMT_SYMBOL *S) {
	/* Draws a great-circle vector with our without heads, etc. There are some complications to consider:
	 * When there are no heads it is simple.  If +n is on we may shrink the line thickness.
	 * With heads there are these cases:
	 * Head length is longer than 90% of the vector length.  We then skip the head and return 1
	 * +n<norm> is in effect.  We shrink vector pen and head length.  Still, the shrunk head
	 * may be longer than 90% of the vector length.  We then shrink head (not pen) further and return 2
	*/

	uint64_t n1, n2, n, add;
	int heads, side[2], outline = 0;
	unsigned int warn = 0;
	size_t n_alloc;
	bool perspective, pure;
	double tlon, tlat, mlon, mlat, P[3], Ax[3], Bx[3], h_length_limit, max_length;
	double dr[2] = {0.0, 0.0}, az[2] = {0.0, 0.0}, oaz[2] = {0.0, 0.0}, off[2] = {0.0, 0.0};
	double da = 0.0, dshift[2] = {0.0, 0.0}, s = 1.0, s1 = 1.0, s2 = 1.0, olon[2], olat[2], head_length, arc_width, rot, scl[2];
	double *xp = NULL, *yp = NULL, *xp2 = NULL, *yp2 = NULL;
	double *rgb = S->v.fill.rgb;
	struct GMT_CIRCLE C;

	/* We must determine points A and B, whose great-circle connector is the arc we seek to draw */
	gmtplot_gcircle_sub (GMT, lon0, lat0, azimuth, length, S, &C);
	perspective = gmt_M_is_perspective (GMT);

	/* Here we have the endpoints A and B of the great (or small) circle arc */

	/* If shrink-option (+n) is active we may have to scale down head attributes and pen width */
	pure = (S->v.v_norm == -1.0f);	/* True if no shrinking has been specified */
	heads = PSL_vec_head (S->v.status);	/* Return head selection as 0-3 */
	h_length_limit = (1.0 - S->v.v_stem) * C.r0;	/* Max length of arrow in degrees to ensure the stem is still showing */
	if (heads == 3) h_length_limit *= 0.5;		/* Split this length between the two heads */
	if (heads && !pure) {	/* Need to determine head length in degrees */
		az[0]  = gmt_az_backaz (GMT, C.lon[0], C.lat[0], C.lon[1], C.lat[1], false);	/* Compute the azimuth from A to B at A along great circle */
		scl[0] = (perspective) ? S->v.scale : gmtplot_get_local_scale (GMT, C.lon[0], C.lat[0], 0.001 * C.r, az[0]);	/* Get local deg/inch scale at A in az[0] direction */
		dr[0] = scl[0] * S->size_x;	/* This is arrow head length in degrees, approximately */
		az[1]  = gmt_az_backaz (GMT, C.lon[1], C.lat[1], C.lon[0], C.lat[0], false);	/* Compute the azimuth from B to A at B along great circle */
		scl[1] = (perspective) ? S->v.scale : gmtplot_get_local_scale (GMT, C.lon[1], C.lat[1], 0.01 * C.r, az[1]);	/* Get local deg/inch scale */
		dr[1] = scl[1] * S->size_x;	/* This is arrow head length in degrees, approximately, adjusted for ~pen thickness to ensure no gap between head and line */
		max_length = MAX (dr[0], dr[1]);
		if (max_length > h_length_limit) {
			s2 = h_length_limit / max_length;
			warn = 2;
		}
	}

	/* Might have to shrink things */
	if (C.r0 < S->v.v_norm)
		s1 = C.r0 / S->v.v_norm;
	if (s1 < s2) {	/* Pick the smallest scale required for head shrinking */
		s = s1;
		warn = 0;	/* When requesting +n there is no warning unless s2 is the smaller scale */
	}
	else
		s = s2;
	head_length = s * S->size_x;
	arc_width   = s1 * S->v.v_width;	/* Use scale s1 for pen shrinking */
	gmt_M_memcpy (olon, C.lon, 2, double);	gmt_M_memcpy (olat, C.lat, 2, double);	/* Keep copy of original coordinates */

	/* When only one side of a vector head is requested (side = -1/+1) there are complications that leads to some
	 * extra work: Since we are clipping the head polygon, the head outline pen is effectively half that of the
	 * vector.  Thus, there will be an offset of 1/2 penwidth at the end of the vector line and the start of the
	 * back-end of the half vector head.  We adjust this (similar to straight and curved vectors in postscriptlight) by moving
	 * the vector tip and the mid-vector back point the equivalent distance of 1/2 pen width away from the side with
	 * the half arrowhead.  This makes the outline of the head align with the vector line. */

	for (n = 0; n < 2; n++) {
		side[n] = PSL_vec_side (S->v.status, n);	/* Return side selection as 0,-1,+1 for this head */
		dshift[n] = (side[n]) ? 0.5 * arc_width : 0.0;	/* Half-width of arc thickness if side != 0 */
	}
	side[0] = -side[0];	/* Since implmenented backwards */
	if (heads & 1) {	/* Placing head at A means we must shorten the arc and use Ax instead of A */
		az[0] = gmt_az_backaz (GMT, C.lon[0], C.lat[0], C.lon[1], C.lat[1], false);	/* Compute the azimuth from A to B at A along great circle */
		scl[0] = (perspective) ? S->v.scale : gmtplot_get_local_scale (GMT, C.lon[0], C.lat[0], 0.001 * C.r, az[0]);	/* Get local deg/inch scale at A in az[0] direction */
		dr[0] = scl[0] * (head_length - 1.1*dshift[0]);	/* This is arrow head length in degrees, approximately, adjusted for ~pen thickness to ensure no gap between head and line */
		if (pure && dr[0] > h_length_limit) {	/* Head length too long, refuse to plot it */
			gmt_M_memcpy (Ax, C.A, 3, double);	/* No need to shorten arc at beginning */
			heads -= 1;	/* Not purse this head any further */
			warn = 1;
		}
		else {
			gmtlib_get_point_from_r_az (GMT, C.lon[0], C.lat[0], 0.5*dr[0]*(1.95 - S->v.v_shape), az[0], &tlon, &tlat);	/* Back mid-point of arrow */
			gmt_geo_to_cart (GMT, tlat, tlon, Ax, true);	/* Get Cartesian coordinates of this new start point for arc */
			dr[0] = scl[0] * head_length;	/* This is arrow head length in degrees, approximately, without any pen-width compensation */
		}
	}
	else
		gmt_M_memcpy (Ax, C.A, 3, double);	/* No need to shorten arc at beginning */

	if (heads & 2) { /* Place arrow head at B */
		az[1] = gmt_az_backaz (GMT, C.lon[1], C.lat[1], C.lon[0], C.lat[0], false);	/* Compute the azimuth from B to A at B along great circle */
		scl[1] = (perspective) ? S->v.scale : gmtplot_get_local_scale (GMT, C.lon[1], C.lat[1], 0.01 * C.r, az[1]);	/* Get local deg/inch scale */
		dr[1] = scl[1] * (head_length - 1.1*dshift[1]);	/* This is arrow head length in degrees, approximately, adjusted for ~pen thickness to ensure no gap between head and line */
		if (pure && dr[1] > h_length_limit) {	/* Head length too long, refuse to plot it */
			gmt_M_memcpy (Bx, C.B, 3, double);	/* No need to shorten arc at end */
			heads -= 2;	/* Not purse this head any further */
			warn = 1;
		}
		else {
			gmtlib_get_point_from_r_az (GMT, C.lon[1], C.lat[1], 0.5*dr[1]*(1.95 - S->v.v_shape), az[1], &tlon, &tlat);	/* Back mid-point of arrow */
			gmt_geo_to_cart (GMT, tlat, tlon, Bx, true);	/* Get Cartesian coordinates of this new end point for arc */
			dr[1] = scl[1] * head_length;	/* This is arrow head length in degrees, approximately, without any pen-width compensation */
		}
	}
	else
		gmt_M_memcpy (Bx, C.B, 3, double);	/* No need to shorten arc at end */

	gmt_M_memcpy (oaz, az, 2, double);	/* Keep copy of original azimuths */

	rot = d_acosd (gmt_dot3v (GMT, Ax, Bx));	/* opening angle in degrees */
	gmtplot_circle_pen_poly (GMT, Ax, Bx, C.longway, rot, ppen, S, &C, s1);

	if (!heads) return (warn);	/* All done */

	/* Get half-angle at head and possibly change pen */
	da = 0.5 * S->v.v_angle;	/* Half-opening angle at arrow head */
	PSL_command (GMT->PSL, "V\n");
	if ((S->v.status & PSL_VEC_OUTLINE) == 0)
		PSL_command (GMT->PSL, "O0\n");	/* Turn off outline */
	else {
		gmt_setpen (GMT, &S->v.pen);
		outline = 1;
	}
	if ((S->v.status & PSL_VEC_FILL) == 0)
		PSL_command (GMT->PSL, "FQ\n");	/* Turn off vector head fill */
	else
		PSL_setfill (GMT->PSL, rgb, outline);

	if (heads & 1) { /* Place arrow head at A */
		if (C.longway) az[0] += 180.0;
		gmtlib_get_point_from_r_az (GMT, C.lon[0], C.lat[0], 0.5*dr[0]*(2.0 - S->v.v_shape), az[0], &mlon, &mlat);	/* Back mid-point of arrow  */
		if (side[0]) {	/* Must adjust the back mid- and end point by 1/2 the pen width */
			az[0] = gmt_az_backaz (GMT, mlon, mlat, C.lon[1], C.lat[1], false);	/* Compute the azimuth from M to B at M */
			scl[0] = (perspective) ? S->v.scale : gmtplot_get_local_scale (GMT, mlon, mlat, tand (da) * dr[0], az[0]+side[0]*90.0);	/* Get deg/inch scale at M perpendicular to arc */
			off[0] = scl[0] * dshift[0];	/* Offset in degrees due to 1/2 pen thickness */
			gmtlib_get_point_from_r_az (GMT, mlon, mlat, off[0], az[0]+side[0]*90.0, &tlon, &tlat);	/* Adjusted back mid-point of arrow head */
			mlon = tlon;	mlat = tlat;	/* Update shifted mid-point */
			gmtlib_get_point_from_r_az (GMT, C.lon[0], C.lat[0], off[0], oaz[0]+side[0]*90.0, &tlon, &tlat);	/* Adjusted tip of arrow head A */
			C.lon[0] = tlon;	C.lat[0] = tlat;	/* Update shifted A location */
			gmt_geo_to_cart (GMT, tlat, tlon, C.A, true);	/* New A vector */
		}
		if (side[0] != +1) {	/* Want to draw left side of arrow */
			gmtlib_get_point_from_r_az (GMT, olon[0], olat[0], dr[0]+off[0], oaz[0]+da, &tlon, &tlat);	/* Start point of arrow on left side */
			gmt_geo_to_cart (GMT, tlat, tlon, P, true);
		}
		else
			gmt_geo_to_cart (GMT, mlat, mlon, P, true);	/* Start from (adjusted) mid point instead */
		n1 = gmtplot_great_circle_arc (GMT, P, C.A, 0.0, false, &xp, &yp);	/* Compute great circle arc from P to A */
		if (side[0] != -1) {	/* Want to draw right side of arrow */
			gmtlib_get_point_from_r_az (GMT, olon[0], olat[0], dr[0]+off[0], oaz[0]-da, &tlon, &tlat);	/* End point of arrow on right side */
			gmt_geo_to_cart (GMT, tlat, tlon, P, true);
		}
		else
			gmt_geo_to_cart (GMT, mlat, mlon, P, true);	/* End at (adjusted) mid point instead */
		n2 = gmtplot_great_circle_arc (GMT, C.A, P, 0.0, false, &xp2, &yp2);	/* Compute great circle arc from A to P */
		add = (side[0] == 0) ? 1 : 0;	/* Need to add mid point explicitly */
		n_alloc = n = n1 + n2 + add;
		gmt_M_malloc2 (GMT, xp, yp, 0U, &n_alloc, double);	/* Allocate space for total path */
		gmt_M_memcpy (&xp[n1], xp2, n2, double);
		gmt_M_memcpy (&yp[n1], yp2, n2, double);
		if (add) {	/* Mid point of arrow */
			xp[n-1] = mlon;	yp[n-1] = mlat;
		}
		gmtplot_plot_vector_head (GMT, xp, yp, n, S);
		gmt_M_free (GMT, xp);	gmt_M_free (GMT, yp);
		gmt_M_free (GMT, xp2);	gmt_M_free (GMT, yp2);
	}
	if (heads & 2) { /* Place arrow head at B */
		if (C.longway) az[1] += 180.0;
		gmtlib_get_point_from_r_az (GMT, C.lon[1], C.lat[1], 0.5*dr[1]*(2.0 - S->v.v_shape), az[1], &mlon, &mlat);	/* Mid point of arrow */
		if (side[1]) {	/* Must adjust the mid-point and end point by 1/2 the pen width */
			az[1] = gmt_az_backaz (GMT, C.lon[1], C.lat[1], C.lon[0], C.lat[0], false);	/* Compute the azimuth from M to A at M */
			scl[1] = (perspective) ? S->v.scale : gmtplot_get_local_scale (GMT, mlon, mlat, tand (da) * dr[1], az[1]+side[1]*90.0);	/* Get deg/inch scale at M perpendicular to arc */
			off[1] = scl[1] * dshift[1];	/* Offset in degrees due to 1/2 pen thickness */
			gmtlib_get_point_from_r_az (GMT, mlon, mlat, off[1], az[1]+side[1]*90.0, &tlon, &tlat);	/* Adjusted back mid-point of arrow head  */
			mlon = tlon;	mlat = tlat;	/* Update shifted mid-point */
			gmtlib_get_point_from_r_az (GMT, C.lon[1], C.lat[1], off[1], oaz[1]+side[1]*90.0, &tlon, &tlat);	/* Adjusted tip of arrow head */
			C.lon[1] = tlon;	C.lat[1] = tlat;	/* Update shifted B location */
			gmt_geo_to_cart (GMT, tlat, tlon, C.B, true);	/* New B vector */
		}
		if (side[1] != +1) {	/* Want to draw left side of arrow */
			gmtlib_get_point_from_r_az (GMT, olon[1], olat[1], dr[1]+off[1], oaz[1]+da, &tlon, &tlat);	/* Start point of arrow on left side */
			gmt_geo_to_cart (GMT, tlat, tlon, P, true);
		}
		else
			gmt_geo_to_cart (GMT, mlat, mlon, P, true);	/* Start from (adjusted)mid point instead */
		n1 = gmtplot_great_circle_arc (GMT, P, C.B, 0.0, false, &xp, &yp);	/* Compute great circle arc from P to B */
		if (side[1] != -1) {	/* Want to draw right side of arrow */
			gmtlib_get_point_from_r_az (GMT, olon[1], olat[1], dr[1]+off[1], oaz[1]-da, &tlon, &tlat);	/* Start point of arrow on other side */
			gmt_geo_to_cart (GMT, tlat, tlon, P, true);
		}
		else
			gmt_geo_to_cart (GMT, mlat, mlon, P, true);	/* End at (adjusted) mid point instead */
		n2 = gmtplot_great_circle_arc (GMT, C.B, P, 0.0, false,  &xp2, &yp2);	/* Compute great circle arc from B to P */
		add = (side[1] == 0) ? 1 : 0;	/* Need to add mid point explicitly */
		n_alloc = n = n1 + n2 + add;
		gmt_M_malloc2 (GMT, xp, yp, 0U, &n_alloc, double);	/* Allocate space for total path */
		gmt_M_memcpy (&xp[n1], xp2, n2, double);
		gmt_M_memcpy (&yp[n1], yp2, n2, double);
		if (add) {	/* Mid point of arrow */
			xp[n-1] = mlon;	yp[n-1] = mlat;
		}
		gmtplot_plot_vector_head (GMT, xp, yp, n, S);
		gmt_M_free (GMT, xp);	gmt_M_free (GMT, yp);
		gmt_M_free (GMT, xp2);	gmt_M_free (GMT, yp2);
	}
	PSL_command (GMT->PSL, "U\n");
	return (warn);
}

/*----------------------------------------------------------|
 * Public functions that are part of the GMT Devel library  |
 *----------------------------------------------------------|
 */

void gmt_map_text (struct GMT_CTRL *GMT, double x, double y, struct GMT_FONT *font, char *label, double angle, int just, unsigned int form) {
	/* Function to plot single-line text in pstext.c */
	struct PSL_CTRL *PSL= GMT->PSL;

	if (gmt_text_is_latex (GMT, label)) {	/* Detected LaTeX commands, i.e., "....@[LaTeX...@[ ..." or  "....<math>LaTeX...</math> ..." */
		double w, h;
		unsigned char *eps = NULL;
		struct imageinfo header;

		if ((eps = gmtplot_latex_eps (GMT, font, label, &header)) == NULL) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "gmt_map_text: Conversion of LaTeX to EPS failed\n");
			return;	/* Done */
		}
		/* Scale up EPS dimensions by the ratio of label font size to LaTeX default size of 10p */
		w = (header.width  / 72.0) * (font->size / 10.0);
		h = (header.height / 72.0) * (font->size / 10.0);
		/* Place EPS file as label, then free eps */
		GMT_Report (GMT->parent, GMT_MSG_DEBUG, "gmt_map_text: Conversion of LaTeX gave dimensions %g x %g\n", w, h);
		PSL_command (PSL, "V\n");	/* Keep the relative changes inside a save/restore block */
		PSL_setorigin (PSL, x, y, angle, PSL_FWD);		/* Move to desired point and possibly rotate to angle */
		PSL_plotlatexeps (PSL, 0.0, 0.0, w, h, just, eps, font->fill.rgb, &header);	/* Place the EPS plot */
		PSL_command (PSL, "U\n");	/* Close up the block */
		PSL_free (eps);
	}
	else	/* Regular text label */
		PSL_plottext (PSL, x, y, font->size, label, angle, just, form);
}

void gmt_xy_axis (struct GMT_CTRL *GMT, double x0, double y0, double length, double val0, double val1, struct GMT_PLOT_AXIS *A, bool below, unsigned side) {
	unsigned int k, i, nx, nx1, np = 0;/* Misc. variables */
	unsigned int annot_pos;	/* Either 0 for upper annotation or 1 for lower annotation */
	unsigned int primary;		/* Axis item number of annotation with largest interval/unit */
	unsigned int axis = A->id;	/* Axis id (GMT_X, GMT_Y, GMT_Z) */
	unsigned int justify;
	unsigned int lx_just = PSL_BC, ly_just = PSL_BC;
	bool horizontal;		/* true if axis is horizontal */
	bool annotate = ((side & GMT_AXIS_ANNOT) > 0);
	bool neg = below;		/* true if annotations are to the left of or below the axis */
	bool faro;			/* true if the anchor point of annotations is on the far side of the axis */
	bool first = true;
	bool is_interval;		/* true when the annotation is interval annotation and not tick annotation */
	bool do_annot;		/* true unless we are dealing with Gregorian weeks */
	bool do_tick;		/* true unless we are dealing with bits of weeks */
	bool form;			/* true for outline font */
	bool ortho = false;		/* true if annotations are orthogonal to axes */
	bool flip = false;		/* true if annotations are inside axes */
	bool MM_set = false;		/* true after we define the MM PS macro for label offsets */
	bool angled = false;		/* True if user used +angle to select a slanted annotation */
	bool skip = false, just_set = false;
	bool save_pi = GMT->current.plot.substitute_pi;
	double *knots = NULL, *knots_p = NULL;	/* Array pointers with tick/annotation knots, the latter for primary annotations */
	double x, t_use, text_angle, cos_a = 0.0, sin_a = 0.0, delta;	/* Misc. variables */
	double x_angle_add = 0.0, y_angle_add = 0.0;	/* Used when dealing with perspectives */
	struct GMT_FONT font;			/* Annotation font (FONT_ANNOT_PRIMARY or FONT_ANNOT_SECONDARY) */
	struct GMT_PLOT_AXIS_ITEM *T = NULL;	/* Pointer to the current axis item */
	char string[GMT_LEN256] = {""};	/* Annotation string */
	char format[GMT_LEN256] = {""};	/* format used for non-time annotations */
	char *axis_chr[3] = {"ns", "ew", "zz"};	/* Characters corresponding to axes */
	char **label_c = NULL;
	double (*xyz_fwd) (struct GMT_CTRL *, double) = NULL;
	struct PSL_CTRL *PSL= GMT->PSL;

	/* Initialize parameters for this axis */

	if (GMT->current.proj.three_D && axis != GMT_Z) {
		/* Because perspective is now done in PostScript we must sometimes reverse the label orientations
		 * so they don't appear upside down when viewed from certain quadrants */
		switch (GMT->current.proj.z_project.quadrant) {
			case 1: x_angle_add = 180.0; lx_just = PSL_TC; break;
			case 3: y_angle_add = 180.0; ly_just = PSL_TC; break;
			case 4: x_angle_add = y_angle_add = 180.0; lx_just = ly_just = PSL_TC; break;
			default: x_angle_add = y_angle_add = 0.0; break;
		}
	}
	horizontal = (axis == GMT_X);	/* This is a horizontal axis */
	xyz_fwd = ((axis == GMT_X) ? &gmt_x_to_xx : (axis == GMT_Y) ? &gmt_y_to_yy : &gmt_z_to_zz);
	primary = gmtplot_get_primary_annot (A);			/* Find primary axis items */
	if (A->use_angle) {	/* Must honor the +a modifier */
		if (axis != GMT_X && doubleAlmostEqualZero (A->angle, 90.0)) ortho = false, just_set = true;	/* Y/Z-Annotations are parallel */
		else if (axis != GMT_X && doubleAlmostEqualZero (A->angle, 0.0)) ortho = true, just_set = true;	/* Y/Z-Annotations are normal */
		if (doubleAlmostEqualZero (A->angle, 0.0)) skip = true;	/* Annotations are parallel so do nothing */
	}
	else if (strchr (GMT->current.setting.map_annot_ortho, axis_chr[axis][below]))
		ortho = true;	/* Annotations are orthogonal for this axis */
	if (GMT->current.setting.map_frame_type & GMT_IS_INSIDE) neg = !neg;	/* Annotations go either below or above the axis */
	faro = (neg == (horizontal && !ortho));			/* Current point is at the far side of the tickmark? */
	if (A->type != GMT_TIME) {						/* Set the annotation format template */
		int n_dec = gmt_get_format (GMT, gmtlib_get_map_interval (GMT, A->type, &A->item[GMT_ANNOT_UPPER]), A->unit, A->prefix, format);
		if (!doubleAlmostEqualZero (A->phase, 0.0)) {	/* Phase is nonzero, must see if adding it gives more decimals */
			char format_w_phase[GMT_LEN256] = {""};	/* Format used if phase is added */
			int ndec_p = gmt_get_format (GMT, gmtlib_get_map_interval (GMT, A->type, &A->item[GMT_ANNOT_UPPER]) + A->phase, A->unit, A->prefix, format_w_phase);
			if (ndec_p > n_dec)	/* Use the phase format which has more decimals */
				strncpy (format, format_w_phase, GMT_LEN256);
		}
	}
	text_angle = (ortho == horizontal) ? 90.0 : 0.0;
	justify = (ortho) ? PSL_MR : PSL_BC;
	if (A->use_angle && !skip) {	/* User override annotation angle */
		text_angle = A->angle;
		angled = true;
		if (!just_set) {
			if (text_angle > 0.0)
				justify = (below) ? PSL_MR : PSL_ML;
			else
				justify = (below) ? PSL_ML : PSL_MR;
		}
		if (axis == GMT_Y) {
			ortho = true;
			if (just_set && !below) justify = PSL_TC;
		}
		if (axis == GMT_X) {
			cos_a = 0.5 * cosd (text_angle);	/* Half-height of text at an angle */
			sin_a = fabs (sind (text_angle));	/* Fraction of y offset due to slanted annotation */
		}
		else {
			cos_a = cosd (text_angle);	/* Full-height of text at an angle */
			sin_a = fabs (sind (text_angle));	/* The y offset due to slanted annotation */
		}
	}
	flip = (GMT->current.setting.map_frame_type & GMT_IS_INSIDE);	/* Inside annotation */
	if (axis != GMT_Z && GMT->current.proj.three_D && GMT->current.proj.z_project.cos_az > 0) {	/* Rotate x/y-annotations when seen "from North" */
		if (!flip) justify = gmt_flip_justify (GMT, justify);
		text_angle += 180.0;
	}
//	else if (flip)
//		justify = gmt_flip_justify (GMT, justify);

	/* Ready to draw axis */
	if (axis == GMT_X)
		PSL_comment (PSL, below ? "Start of lower x-axis\n" : "Start of upper x-axis\n");
	else if (axis == GMT_Y)
		PSL_comment (PSL, below ? "Start of left y-axis\n" : "Start of right y-axis\n");
	else
		PSL_comment (PSL, below ? "Start of front z-axis\n" : "Start of back z-axis\n");
	PSL_setorigin (PSL, x0, y0, 0.0, PSL_FWD);

	PSL_comment (PSL, "Axis tick marks and annotations\n");
	gmt_setpen (GMT, &GMT->current.setting.map_frame_pen);
	if (horizontal)
		PSL_plotsegment (PSL, 0.0, 0.0, length, 0.0);
	else
		PSL_plotsegment (PSL, 0.0, length, 0.0, 0.0);
	if (GMT->current.setting.map_frame_type & GMT_IS_GRAPH) {	/* Extend axis with an arrow */
		struct GMT_FILL arrow;
		double vector_width, dim[PSL_MAX_DIMS], g_scale_begin = 0.0, g_scale_end = 0.0, g_ext = 0.0;
		gmt_init_fill (GMT, &arrow, GMT->current.setting.map_frame_pen.rgb[0], GMT->current.setting.map_frame_pen.rgb[1], GMT->current.setting.map_frame_pen.rgb[2]);
		gmt_setfill (GMT, &arrow, 0);
		gmt_M_memset (dim, PSL_MAX_DIMS, double);
		vector_width = rint (PSL_DOTS_PER_INCH * GMT->current.setting.map_frame_pen.width / PSL_POINTS_PER_INCH) / PSL_DOTS_PER_INCH;	/* Round off vector width same way as pen width */
		dim[PSL_VEC_TAIL_WIDTH]    = vector_width;
		dim[PSL_VEC_HEAD_LENGTH]   = 10.0 * vector_width;
		dim[PSL_VEC_HEAD_WIDTH]    = 5.0 * vector_width;
		dim[PSL_VEC_HEAD_SHAPE]    = GMT->current.setting.map_vector_shape;
		dim[PSL_VEC_STATUS]        = PSL_VEC_END | PSL_VEC_FILL;
		dim[PSL_VEC_HEAD_PENWIDTH] = 0.5 * GMT->current.setting.map_frame_pen.width;
		PSL_defpen (PSL, "PSL_vecheadpen", GMT->current.setting.map_frame_pen.width, GMT->current.setting.map_frame_pen.style, GMT->current.setting.map_frame_pen.offset, GMT->current.setting.map_frame_pen.rgb);

		/* The start of the vector at the end of a positive axis is computed as x = g_scale_end * length + g_ext;
		 * The start of the vector at the start of a negative axis is computed as x = g_scale_begin * length - g_ext;
		 * The required constants are set below */
		if (GMT->current.setting.map_graph_extension_unit == GMT_GRAPH_EXTENSION_UNIT) {	/* Gave scaling in percent */
			/* Here, the computed length is to the tip of the vector */
			g_scale_end = 1.0 + 0.01 * GMT->current.setting.map_graph_extension;
			g_scale_begin = -0.01 * GMT->current.setting.map_graph_extension;
		}
		else {	/* Gave actual length extension (now in inches) */
			/* Here, the computed length is to the base of the vector */
			g_scale_end = 1.0;
			g_ext = GMT->current.setting.map_graph_extension + dim[3];
		}

		if (horizontal) {
			double x = 0.0;
			if (GMT->current.proj.xyz_pos[axis]) {
				x = length;
				dim[0] = g_scale_end * length + g_ext;
			}
			else
				dim[0] = g_scale_begin * length - g_ext;
			PSL_plotsymbol (PSL, x, 0.0, dim, PSL_VECTOR);
		}
		else {
			double y = 0.0;
			if (GMT->current.proj.xyz_pos[axis]) {
				y = length;
				dim[1] = g_scale_end * length + g_ext;
			}
			else
				dim[1] = g_scale_begin * length - g_ext;
			PSL_plotsymbol (PSL, 0.0, y, dim, PSL_VECTOR);
		}
	}

	if (side == GMT_AXIS_DRAW) {	/* Just drawing the axis for this one */
		PSL_setorigin (PSL, -x0, -y0, 0.0, PSL_INV);
		return;
	}

	GMT->current.plot.substitute_pi = A->substitute_pi;
	np = gmtlib_coordinate_array (GMT, val0, val1, &A->item[primary], &knots_p, NULL);	/* Get all the primary tick annotation knots */

	/* Axis items are in order: GMT_ANNOT_UPPER, GMT_ANNOT_LOWER, GMT_TICK_UPPER, GMT_TICK_LOWER, GMT_GRID_UPPER, GMT_GRID_LOWER */

	for (k = 0; k < 2; k++)	/* Define annotation offsets */
		PSL_command (PSL, "/PSL_A%d_y %d def\n", k, A->item[k].active || A->item[k+2].active ? PSL_IZ (PSL, MAX(0.0,GMT->current.setting.map_tick_length[k])) : 0);	/* Length of primary/secondary tickmark */

	for (k = 0; k < GMT_GRID_UPPER; k++) {	/* For each one of the 6 axis items (gridlines are done separately) */

		T = &A->item[k];		/* Get pointer to this item */
		if (!T->active) continue;	/* Do not want this item plotted - go to next item */

		gmt_setpen (GMT, &GMT->current.setting.map_tick_pen[GMT_PRIMARY]);

		is_interval = (T->type == 'i' || T->type == 'I');	/* Interval or tick mark annotation? */
		nx = gmtlib_coordinate_array (GMT, val0, val1, &A->item[k], &knots, &label_c);	/* Get all the annotation tick knots */
		delta = (label_c) ? GMT->session.d_NaN : A->item[k].interval;
		do_annot = (nx && k < GMT_TICK_UPPER && annotate && !gmt_M_axis_is_geo_strict (GMT, axis) && T->unit != 'r');	/* Cannot annotate a Gregorian week */
		do_tick = !((T->unit == 'K' || T->unit == 'k') && T->interval > 1 && fmod (T->interval, 7.0) > 0.0);	/* Do we want tick marks? */
		nx1 = (nx > 0 && is_interval) ? nx - 1 : nx;
		/* First plot all the tick marks */

		if (do_tick) {
			for (i = 0; i < nx; i++) {
				x = (*xyz_fwd) (GMT, knots[i]);	/* Convert to inches on the page */
				if (x < (- GMT_CONV4_LIMIT) || x > (length + GMT_CONV4_LIMIT)) continue;	/* Outside the range */
				if (gmtplot_skip_second_annot (k, knots[i], knots_p, np, primary)) continue;	/* Minor tick marks skipped when coinciding with major */
				if (horizontal)
					PSL_plotsegment (PSL, x, 0.0, x, ((neg) ? -1.0 : 1.0) * GMT->current.setting.map_tick_length[k]);
				else
					PSL_plotsegment (PSL, 0.0, x, ((neg) ? -1.0 : 1.0) * GMT->current.setting.map_tick_length[k], x);
			}
		}

		/* Then do annotations too - here just set text height/width parameters in PostScript */

		if (do_annot) {
			annot_pos = (T->type == 'A' || T->type == 'I') ? 1 : 0;	/* 1 means lower annotation, 0 means upper (close to axis) */
			font = GMT->current.setting.font_annot[annot_pos];	/* Set the font to use */
			form = gmt_setfont (GMT, &font);
			if (A->type != GMT_TIME) PSL_settextmode (PSL, PSL_TXTMODE_MINUS);	/* Replace hyphens with minus */
			if (first) {	/* Change up/down (neg) and/or flip coordinates (exch) */
				PSL_command (PSL, "/MM {%s%sM} def\n", neg ? "neg " : "", (axis != GMT_X) ? "exch " : "");
				first = false;
				MM_set = true;
			}
			PSL_command (PSL, "/PSL_AH%d 0\n", annot_pos);	/* Largest annotation width (or height) so far */

			for (i = 0; i < nx1; i++) {
				if (gmtlib_annot_pos (GMT, val0, val1, T, &knots[i], &t_use)) continue;			/* Outside range */
				if (axis == GMT_Z && fabs (knots[i] - GMT->current.proj.z_level) < GMT_CONV8_LIMIT) continue;	/* Skip z annotation coinciding with z-level plane */
				x = (*xyz_fwd) (GMT, knots[i]);	/* Convert to inches on the page */
				if (gmtplot_skip_end_annotation (GMT, axis, x, length)) continue;	/* Don't want annotations exactly at one or both ends of the axis */
				if (!is_interval && gmtplot_skip_second_annot (k, knots[i], knots_p, np, primary)) continue;	/* Secondary annotation skipped when coinciding with primary annotation */
				if (label_c && label_c[i] && label_c[i][0]) {
					strncpy (string, label_c[i], GMT_LEN256-1);
				}
				else
					gmtlib_get_coordinate_label (GMT, string, &GMT->current.plot.calclock, format, T, knots[i], delta);	/* Get annotation string */
				PSL_deftextdim (PSL, ortho ? "-w" : "-h", font.size, string);
				PSL_command (PSL, "mx\n");		/* Update the longest annotation stored in PSL_AH? */
			}
			PSL_command (PSL, "def\n");	/* Finalize the definition of longest (y-axis) or tallest (x-annotation) found */
			if (angled && axis == GMT_X) {	/* Also need annotation width since a component of length is projected in the y-direction so label and title must be displaced by this amount */
				/* Note: PSL_slant_y will also be used when placing a Cartesian frame title */
				PSL_command (PSL, "/PSL_slant_y ");
				if (label_c && label_c[nx1-1] && label_c[nx1-1][0])	/* We only use the string for the max annotation value at index nx1-1 */
					strncpy (string, label_c[nx1-1], GMT_LEN256-1);
				else
					gmtlib_get_coordinate_label (GMT, string, &GMT->current.plot.calclock, format, T, knots[nx1-1], delta);	/* Get annotation string */
				PSL_deftextdim (PSL, "-w", font.size, string);	/* Compute the width */
				PSL_command (PSL, " %.12g mul def\n", sin_a);	/* Multiply this width by sine of the angle to get the y-component */
				//PSL_command (PSL, " %.12g mul %g add def\n", sin_a, PSL_IZ (PSL, 72.0 * GMT->current.setting.font_annot[GMT_PRIMARY].size));	/* Multiply this width by sine of the angle to get the y-component */
			}
			else if (angled && axis == GMT_Y) {	/* Need a slant in x and reset PSL_AH? to 0 */
				//PSL_command (PSL, "/PSL_slant_x PSL_AH%d %.12g mul %g add def\n", annot_pos, cos_a, PSL_IZ (PSL, 72.0 * GMT->current.setting.font_annot[GMT_PRIMARY].size));	/* Largest annotation width (or height) so far */
				PSL_command (PSL, "/PSL_slant_x PSL_AH%d %.12g mul  def\n", annot_pos, cos_a);	/* Largest annotation width (or height) so far */
				PSL_command (PSL, "/PSL_AH%d 0 def\n", annot_pos);	/* Largest annotation width (or height) so far */
			}

			if (annot_pos == 0)
				PSL_command (PSL, "/PSL_A0_y PSL_A0_y %d add ", PSL_IZ (PSL, GMT->current.setting.map_annot_offset[annot_pos]));
			else
				PSL_command (PSL, "/PSL_A1_y PSL_A0_y PSL_A1_y mx %d add ", PSL_IZ (PSL, GMT->current.setting.map_annot_offset[annot_pos]));
			if (faro) PSL_command (PSL, "PSL_AH%d add ", annot_pos);
			PSL_command (PSL, "def\n");

			for (i = 0; i < nx1; i++) {
				if (gmtlib_annot_pos (GMT, val0, val1, T, &knots[i], &t_use)) continue;			/* Outside range */
				if (axis == GMT_Z && fabs (knots[i] - GMT->current.proj.z_level) < GMT_CONV8_LIMIT) continue;	/* Skip z annotation coinciding with z-level plane */
				x = (*xyz_fwd) (GMT, t_use);	/* Convert to inches on the page */
				if (gmtplot_skip_end_annotation (GMT, axis, x, length)) continue;	/* Don't want annotations exactly at one or both ends of the axis */
				if (!is_interval && gmtplot_skip_second_annot (k, knots[i], knots_p, np, primary)) continue;	/* Secondary annotation skipped when coinciding with primary annotation */
				/* Move to new anchor point */
				PSL_command (PSL, "%d PSL_A%d_y MM\n", PSL_IZ (PSL, x), annot_pos);
				if (angled) {	/* Must compensate for rotated textbox sticking too close to axis */
					if (below) /* S axis */
						PSL_command (PSL, "0 PSL_AH%d 1 %.12g sub mul G\n", annot_pos, cos_a);
					else /* N axis */
						PSL_command (PSL, "0 PSL_AH%d %.12g mul G\n", annot_pos, cos_a);
				}
				if (label_c && label_c[i] && label_c[i][0])
					strncpy (string, label_c[i], GMT_LEN256-1);
				else
					gmtlib_get_coordinate_label (GMT, string, &GMT->current.plot.calclock, format, T, knots[i], delta);	/* Get annotation string */
				PSL_plottext (PSL, 0.0, 0.0, -font.size, string, text_angle, justify, form);
			}
			if (!faro) PSL_command (PSL, "/PSL_A%d_y PSL_A%d_y PSL_AH%d add def\n", annot_pos, annot_pos, annot_pos);
		}

		if (nx) gmt_M_free (GMT, knots);
		if (label_c) {
			for (i = 0; i < nx; i++) gmt_M_str_free (label_c[i]);
			gmt_M_free (GMT, label_c);
		}
		PSL_settextmode (PSL, PSL_TXTMODE_HYPHEN);	/* Replace hyphens with minus */
	}
	if (np) gmt_M_free (GMT, knots_p);

	/* Finally do axis label */

	if (A->label[0] && annotate && !gmt_M_axis_is_geo_strict (GMT, axis)) {
		unsigned int far_ = !below, l_just, L_axis = (axis == GMT_Z) ? GMT_X : axis;
		double label_angle;
		char *this_label = (far_ && A->secondary_label[0]) ? A->secondary_label : A->label;	/* Get primary or secondary axis label */

		if (axis == GMT_Y && A->label_mode) {
			l_just = (below) ? PSL_MR : PSL_ML;
			label_angle = 0.0 + y_angle_add;
		}
		else {
			double angle_add = (axis == GMT_X) ? x_angle_add : y_angle_add;
			l_just = (axis == GMT_X) ? lx_just : ly_just;
			label_angle = (horizontal ? 0.0 : 90.0) + angle_add;
		}
		if (!MM_set) PSL_command (PSL, "/MM {%s%sM} def\n", neg ? "neg " : "", (axis != GMT_X) ? "exch " : "");
		form = gmt_setfont (GMT, &GMT->current.setting.font_label);
		PSL_command (PSL, "/PSL_LH ");	/* PSL_LH is the height of the label text based on height of letter M */
		PSL_deftextdim (PSL, "-h", GMT->current.setting.font_label.size, "M");
		PSL_command (PSL, "def\n");
		if (GMT->current.setting.map_label_mode[L_axis] == GMT_LABEL_AXIS) /* Base of label is a fixed distance from axis */
			PSL_command (PSL, "/PSL_L_y %d %sdef\n", PSL_IZ (PSL, fabs (GMT->current.setting.map_label_offset[L_axis])), (neg == horizontal) ? "PSL_LH add " : "");
		else
			PSL_command (PSL, "/PSL_L_y PSL_A0_y PSL_A1_y mx %d add %sdef\n", PSL_IZ (PSL, GMT->current.setting.map_label_offset[L_axis]), (neg == horizontal) ? "PSL_LH add " : "");
		/* Move to new anchor point for label */
		if (angled && axis == GMT_X)	/* Add offset due to angled x-annotations */
			PSL_command (PSL, "%d PSL_L_y PSL_slant_y add MM\n", PSL_IZ (PSL, 0.5 * length));
		else if (angled && axis == GMT_Y)
			PSL_command (PSL, "%d PSL_L_y PSL_slant_x add MM\n", PSL_IZ (PSL, 0.5 * length));
		else
			PSL_command (PSL, "%d PSL_L_y MM\n", PSL_IZ (PSL, 0.5 * length));
		gmtplot_map_label (GMT, 0.0, 0.0, this_label, label_angle, l_just, axis, below);
	}
	else
		PSL_command (PSL, "/PSL_LH 0 def /PSL_L_y PSL_A0_y PSL_A1_y mx def\n");
	if (axis == GMT_X)
		PSL_comment (PSL, below ? "End of lower x-axis\n" : "End of upper x-axis\n");
	else if (axis == GMT_Y)
		PSL_comment (PSL, below ? "End of left y-axis\n" : "End of right y-axis\n");
	else
		PSL_comment (PSL, below ? "End of front z-axis\n" : "End of back z-axis\n");
	PSL_setorigin (PSL, -x0, -y0, 0.0, PSL_INV);

	GMT->current.plot.substitute_pi = save_pi;
}

GMT_LOCAL unsigned int gmtplot_the_side_we_cut (struct GMT_CTRL *GMT, double x, double y) {
	/* We know the (x,y) sits on the map boundary.  Depending on what map scheme
	 * we have we use different strategies to find which side the cut is on. */
	if (gmt_M_is_azimuthal (GMT) && !GMT->current.proj.polar) {	/* Not implemented yet */
		/* Here the map boundaries are things like circles and there is no "side" */
		return 9;
	}
	else if (GMT->common.R.oblique) {	/* Rectangular box, work with plot units */
		/* We first check if the x value equal either the left or right x-coordinate
		 * on the border for this y.  If that fails then we know it sits on the south
		 * or north border and we can just pick the closest one. */
		double width = gmt_half_map_width (GMT, y);	/* Width of map at current latitude (not all projections have straight w/e boundaries */
		double xx = GMT->current.map.half_width - width;	/* x-coordinate on west boundary for this y */
		if (doubleAlmostEqualZero (x, xx)) return 3;
		xx = GMT->current.map.half_width + width;		/* x-coordinate on east boundary for this y */
		if (doubleAlmostEqualZero (x, xx)) return 1;
		if (y < GMT->current.map.half_height) return 0;		/* Sits on south border */
		return 2;	/* OK, so it was north */
	}
	else {	/* Boundaries are meridians and parallels so just detect which extrema */
		double lon, lat;
		gmt_xy_to_geo (GMT, &lon, &lat, x, y);
		if (doubleAlmostEqualZero (lon, GMT->common.R.wesn[XHI])) return 1;
		else if (doubleAlmostEqualZero (lon, GMT->common.R.wesn[XLO])) return 3;
		else if (doubleAlmostEqualZero (lat, GMT->common.R.wesn[YHI])) return 2;
		else if (doubleAlmostEqualZero (lat, GMT->common.R.wesn[YLO])) return 0;
		else return 9;	/* SHould not happen but if it does a 9 means "do not extend" */
	}
}

GMT_LOCAL double gmtplot_get_border_angle (struct GMT_CTRL *GMT, double xc, double yc, unsigned int side) {
	/* Return the angle of the tangent to the map border at the point x, y.  This does
	 * not work for circular boundaries yet */
	double lon, lat, lon0, lat0, lon1, lat1, del, angle, x0, y0, x1, y1;
	gmt_xy_to_geo (GMT, &lon, &lat, xc, yc);
	switch (side) {
		case 0:	case 2:	/* Pick two points along the parallel */
			if (GMT->common.R.oblique) return (0.0);	/* We know it is horizontal */
			del = (GMT->common.R.wesn[XHI] - GMT->common.R.wesn[XLO]) / 360.0;
			//if (del > 0.1) del = 0.1;
			lon0 = lon - del;	lon1 = lon + del;	lat0 = lat1 = lat;
			break;
		default:	/* Pick two points along the meridian */
			if (GMT->common.R.oblique) return (90.0);	/* We know it is vertical */
			del = (GMT->common.R.wesn[YHI] - GMT->common.R.wesn[YLO]) / 180.0;
			//if (del > 0.1) del = 0.1;
			/* Make sure lats do not exceed the limits */
			lat0 = MAX (-90.0, lat - del);
			lat1 = MIN (+90.0, lat + del);	lon0 = lon1 = lon;
			break;
	}
	gmt_geo_to_xy (GMT, lon0, lat0, &x0, &y0);
	gmt_geo_to_xy (GMT, lon1, lat1, &x1, &y1);
	angle = atan2 (y1 - y0, x1 - x0) * R2D;	/* Get angle of line in direction from "on" point to "in" point */

	return angle;
}

#if 0
/* Saving this for now in case something blows */
GMT_LOCAL void gmtplot_get_outside_point_extension_orig (struct GMT_CTRL *GMT, double x_on, double y_on, double x_in, double y_in, double *x_off, double *y_off) {
	/* The point (x_on, y_on) is known to sit on a rectangular map border, and the point (x_in, y_in) is either inside or is on another border point.
	 * A line will be drawn between the two points, but given the finite pen width we must adjust the starting and|or end point so that
	 * the pen is not clipped, resulting in a gap between the end and the border. */
	double W = 0.5 * GMT->current.setting.ps_penwidth * GMT->session.u2u[GMT_PT][GMT_INCH];	/* Half the current pen width in inches */
	double L, angle, dx, dy, tan_angle;
	static char side[] = "SENW";
	unsigned int k = gmtplot_the_side_we_cut (GMT, x_on, y_on);	/* Which border side is this point on? */
	dx = x_in - x_on;	dy = y_in - y_on;		/* Get coordinate increments from the "on" point to the "in" point */
	angle = atan2 (dy, dx) * R2D;				/* Get angle of line in direction from "on" point to "in" point */
	tan_angle = tand (angle);				/* Compute the tangent of that line direction (-180 to +180)  */
	/* Need to modify the stuff below which assumes border_angle = 0|90 */
	if (k%2 == 0) {	/* Cutting the S or N border */
		if (doubleAlmostEqualZero (dx, 0.0))		/* Line is almost vertical; no need to extend it when at a horizontal border */
			L = 0.0;
		else if (doubleAlmostEqualZero (dy, 0.0))	/* Line is almost horizontal; would get infinity so truncate to width of map */
			L = GMT->current.map.width;
		else	/* Safe to get width (but still truncate to width of map) */
			L = MIN (fabs (W / tan_angle), GMT->current.map.half_width);	/* L is positive */
	}
	else {	/* Cutting the E or W border */
		if (doubleAlmostEqualZero (dx, 0.0))		/* Line is almost vertical; would get infinity so clip truncate height of map */
			L = GMT->current.map.height;
		else if (doubleAlmostEqualZero (dy, 0.0))	/* Line is almost horizontal; no need to extend it when at a vertical border */
			L = 0.0;
		else	/* Safe to get width (but still truncate to height of map) */
			L = MIN (fabs (W * tan_angle), GMT->current.map.height);	/* L is positive */
	}
	/* Compute the coordinates of the point at a distance L away from clip point in the direction of angle */
	*x_off = x_on - L * cosd (angle);
	*y_off = y_on - L * sind (angle);
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Extend from (%g, %g) to (%g, %g) in direction %g by %g\" on %c side\n", *x_off, *y_off, x_on, y_on, angle, L, side[k]);
}
#endif

GMT_LOCAL void gmtplot_get_outside_point_extension (struct GMT_CTRL *GMT, double x_on, double y_on, double x_in, double y_in, double *x_off, double *y_off) {
	/* The point (x_on, y_on) is known to sit on a rectangular map border, and the point (x_in, y_in) is either inside or is on another border point.
	 * A line will be drawn between the two points, but given the finite pen width we must adjust the starting and|or end point so that
	 * the pen is not clipped, resulting in a gap between the end and the border. */
	double W = 0.5 * GMT->current.setting.ps_penwidth * GMT->session.u2u[GMT_PT][GMT_INCH];	/* Half the current pen width in inches */
	double L, angle, dx, dy, tan_angle, border_angle, d_angle;
	static char side[] = "SENW";
	unsigned int k = gmtplot_the_side_we_cut (GMT, x_on, y_on);	/* Which border side is this point on? */
	if (k == 9) {	/* No can do yet */
		*x_off = x_on; *y_off = y_on; return;
	}
	dx = x_in - x_on;	dy = y_in - y_on;		/* Get coordinate increments from the "on" point to the "in" point */
	angle = atan2 (dy, dx) * R2D;				/* Get angle of line in direction from "on" point to "in" point */
	tan_angle = tand (angle);				/* Compute the tangent of that line direction (-180 to +180)  */
	border_angle = gmtplot_get_border_angle (GMT, x_on, y_on, k);
	d_angle = angle - border_angle;		/* Compute the tangent of that line direction (-180 to +180)  */
	tan_angle = tand (d_angle);		/* Compute the tangent of that line direction (-180 to +180)  */
	/* The is_dnan catches bad cases yielding a NaN angle */
	if (gmt_M_is_dnan (d_angle) || doubleAlmostEqualZero (fabs (d_angle), 90.0))		/* Line is orthogonal to border; no need to extend it when at a horizontal border */
		L = 0.0;
	else if (doubleAlmostEqualZero (d_angle, 0.0))	/* Line is parallel to border; would get infinity so truncate to width of map */
		//L = GMT->current.map.width; /* What we want but still picking up stray horizontal lines, hence L = 0 for now */
		L = 0.0;
	else	/* Safe to get width (but still truncate to width of map) */
		L = MIN (fabs (W / tan_angle), GMT->current.map.half_width);	/* L is positive */
	/* Compute the coordinates of the point at a distance L away from clip point in the direction of angle */
	*x_off = x_on - L * cosd (angle);
	*y_off = y_on - L * sind (angle);
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Extend from (%g, %g) to crossing point (%g, %g) in direction %g by %g\" on %c side with border angle %g and delta angle %g\n", *x_off, *y_off, x_on, y_on, angle, L, side[k], border_angle, d_angle);
}

GMT_LOCAL bool gmtplot_these_are_duplicates (double x0, double y0, double x1, double y1) {
	/* Return true if (x0,y0) and (x1,y1) is the same point within bitnoise */
	if (!doubleAlmostEqualZero (x0, x1)) return false;
	if (!doubleAlmostEqualZero (y0, y1)) return false;
	return true;
}

void gmt_plot_line (struct GMT_CTRL *GMT, double *x, double *y, unsigned int *pen, uint64_t n, unsigned int mode) {
	/* Mode = PSL_LINEAR [0] or PSL_BEZIER [1] */
	/* Drawing lines of finite thickness that cross map boundaries is tricky because the exit or entry point
	 * refers to the center line of a finite thickness line but the wider part of the line do not share those
	 * coordinates.  To avoid ugly line gaps we must extend such lines a bit beyond the border and clip so that
	 * things look pretty.  For now we only do this for a few situations:
	 * 1) The map border must be rectangular
	 * 2) The line must be a linear spline and not the Bezier spline
	 * Hopefully we can work to avoid the first case soon */

	uint64_t i, j, im1, ip1, k;
	int way;
	bool close, stop;
	double x_cross[2], y_cross[2];
	double x_ext, y_ext;	/* Coordinates of the external helper point */
	unsigned int start_pen = PSL_MOVE, end_pen = PSL_STROKE;	/* Default pen code of given line endings */
	struct PSL_CTRL *PSL= GMT->PSL;

	if (n < 2) return;

	gmtplot_NaN_pen_up (x, y, pen, n);	/* Ensure we don't have NaNs in the coordinates */

	/* First skip any repeating PSL_MOVE in the beginning since only the last one will matter */
	i = 0;
	while (i < (n-1) && (pen[i+1] & PSL_MOVE)) i++;	/* Skip repeating pen & PSL_MOVE in beginning */
	if ((n-i) < 2) return;	/* Less than 2 points is not a line */
	j = i;	/* j is now the first valid point with a PSL_MOVE */
	/* Then we want to prevent duplicate points since the first point may have a different pen (e.g., include PSL_CLIP) */
	while (i < (n-1) && gmtplot_these_are_duplicates (x[j], y[j], x[i+1], y[i+1])) i++;	/* Skip duplicate points in beginning */
	if ((n-i) < 2) return;	/* Less than 2 points is not a line */
	if (i > j) pen[i] = pen[j];	/* Skipped initial duplicates but must maintain the initial pen code */
	/* Likewise, skip any repeating PSL_MOVE at the end */
	while (n > 1 && (pen[n-1] & PSL_MOVE)) n--;	/* Cut off repeating pen & PSL_MOVE at end */
	if ((n-i) < 2) return;	/* Less than 2 points is not a line */
	/* Then we want to prevent duplicate points at end since the last valid point may have a different pen (e.g., include PSL_CLIP) */
	j = n - 1;	/* j is now last point after initial skip test */
	while (n > 1 && gmtplot_these_are_duplicates (x[j], y[j], x[n-2], y[n-2])) n--;	/* Skip duplicate points at end */
	if ((n-i) < 2) return;	/* Less than 2 points is not a line */
	if (n <= j) {	/* Got a duplicate at the end */
		if ((pen[n-1] & PSL_CLIP) == 0) pen[n-1] = pen[j];	/* Skipped trailing duplicates but must maintain initial pen code of last valid point unless has clip information */
		while (n > 1 && (pen[n-1] & PSL_MOVE)) n--;	/* Cut off repeating pen & PSL_MOVE at end that might result after removing that duplicate */
	}
	for (j = i + 1; j < n && !(pen[j] & PSL_MOVE); j++);	/* j == n means no PSL_MOVEs present */
	close = (j == n) ? (hypot (x[n-1] - x[i], y[n-1] - y[i]) < GMT_CONV4_LIMIT) : false;

	/* First see if we can use the PSL_plotline call directly to save points */

	for (j = i + 1, stop = false; !stop && j < n; j++) stop = (pen[j] & PSL_MOVE || (*GMT->current.map.jump) (GMT, x[j-1], y[j-1], x[j], y[j]));
	if (!stop) {	/* Safe to draw single line directly since no jumping around */
		if (mode == PSL_BEZIER)
			PSL_plotcurve (PSL, &x[i], &y[i], (int)(n - i), PSL_MOVE + PSL_STROKE + close * PSL_CLOSE);
		else {	/* Pay more attention here since exiting lines of finite thickness may need some corrections */
			//bool rect = gmt_M_is_rect_graticule (GMT) || GMT->common.R.oblique;	/* We have a rectangular map boundary */
			bool rect = !(gmt_M_is_azimuthal (GMT) && !GMT->current.proj.polar);
			for (j = i + 1, stop = false; !stop && j < n; j++) stop = (pen[j] & PSL_CLIP);	/* Look for clipped end points */
			if (!rect || !stop)	/* No troubled clip points or it is not a rectangular border */
			//if (!stop)	/* No troubled clip points or it is not a rectangular border */
				PSL_plotline (PSL, &x[i], &y[i], (int)(n - i), PSL_MOVE + PSL_STROKE + close * PSL_CLOSE);
			else {	/* At least one clip point on a rectangular border will need a line extension  */
				if (pen[i] & PSL_CLIP) { 	/* Must add extra point before the first point */
					gmtplot_get_outside_point_extension (GMT, x[i], y[i], x[i+1], y[i+1], &x_ext, &y_ext);
					start_pen = PSL_DRAW;	/* Cannot have a move when drawing line after this point */
					PSL_plotpoint (PSL, x_ext, y_ext, PSL_MOVE);	/* Lay down new start point */
				}
				j = n - 1;	/* Last point index */
				if (pen[j] & PSL_CLIP) {	/* Must draw to extra point after the last point, so set end_pen to draw since not done yet */
					gmtplot_get_outside_point_extension (GMT, x[j], y[j], x[j-1], y[j-1], &x_ext, &y_ext);
					end_pen = PSL_DRAW;
				}
				PSL_plotline (PSL, &x[i], &y[i], (int)(n - i), start_pen|end_pen);	/* Lay down the original line, possibly with adjusted pen code */
				if (end_pen == PSL_DRAW) PSL_plotpoint (PSL, x_ext, y_ext, PSL_STROKE);	/* Add extension from last point to outside helper point */
			}
		}
		return;
	}

	/* Here we must check for jumps, pen changes etc */

	 if (GMT->common.R.oblique) {
		  double www = GMT->current.proj.rect[XHI] - GMT->current.proj.rect[XLO];
		  for (k = i+1; k < n; k++) {
				if (doubleAlmostEqual (fabs (x[k-1] - x[k]), www) && pen[k] & PSL_CLIP) {
					 if (pen[k] == PSL_CLIP) /* Fix a mistake made earlier */
						  pen[k] = PSL_MOVE;
				}
		  }
	 }
	start_pen = pen[i];
	if (pen[i] & PSL_CLIP) { 	/* Must add extra extension point before the first point */
		gmtplot_get_outside_point_extension (GMT, x[i], y[i], x[i+1], y[i+1], &x_ext, &y_ext);
		start_pen = PSL_DRAW;	/* Cannot have a move when drawing line after this point */
		PSL_plotpoint (PSL, x_ext, y_ext, PSL_MOVE);	/* Lay down new start point */
	}
	PSL_plotpoint (PSL, x[i], y[i], start_pen);

	i++;
	while (i < n) {
		im1 = i - 1;
		ip1 = i + 1;
		if (pen[im1] & PSL_MOVE && pen[i] == PSL_DRAW && (ip1 == n || pen[ip1] & PSL_MOVE) && GMT->current.proj.projection_GMT == GMT_OBLIQUE_MERC && fabs (y[i]-y[im1]) < 0.001) {
			double mw = 2.0 * gmt_half_map_width (GMT, y[i]);	/* Get map width */
			/* We have a single 2-point ~horizontal line segment with move, draw, move.  Check if the draw is across the entire oblique Mercator map */
			if (doubleAlmostEqual (fabs (x[i]-x[im1]), mw)) {	/* Yes, so skip such stray lines across map */
				/* This fix was implemented in response to the problem first illustrated in test/psbasemap/oblique.sh.
				 * Ideally, we should fix this upstream but not that easy to follow the logic. */
				i++;
				GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Skipping a stray line across map.\n");
				continue;
			}
		}
		if (pen[i] == pen[im1] && (way = (*GMT->current.map.jump) (GMT, x[im1], y[im1], x[i], y[i]))) {	/* Jumped across the map */
			(*GMT->current.map.get_crossings) (GMT, x_cross, y_cross, x[im1], y[im1], x[i], y[i]);
			if (way == -1) {	/* Add left border point */
				PSL_plotpoint (PSL, x_cross[0], y_cross[0], PSL_DRAW);	/* Draw to left boundary... */
				gmtplot_get_outside_point_extension (GMT, x[im1], y[im1], x_cross[0], y_cross[0], &x_ext, &y_ext);	/* And extend further to get proper clipping */
				PSL_plotpoint (PSL, x_ext, y_ext, PSL_DRAW);	/* Draw to outside point */
				gmtplot_get_outside_point_extension (GMT, x_cross[1], y_cross[1], x[i], y[i], &x_ext, &y_ext);	/* Find outside point on the right */
				PSL_plotpoint (PSL, x_ext, y_ext, PSL_MOVE);	/* Lay down new outside point */
				PSL_plotpoint (PSL, x_cross[1], y_cross[1], PSL_DRAW);	/* ...then jump to the right boundary */
			}
			else {
				PSL_plotpoint (PSL, x_cross[1], y_cross[1], PSL_DRAW);	/* Draw to right boundary... */
				gmtplot_get_outside_point_extension (GMT, x[im1], y[im1], x_cross[1], y_cross[1], &x_ext, &y_ext);	/* And extend further to get proper clipping */
				PSL_plotpoint (PSL, x_ext, y_ext, PSL_DRAW);	/* Draw to outside point */
				gmtplot_get_outside_point_extension (GMT, x[i], y[i], x_cross[0], y_cross[0], &x_ext, &y_ext);	/* Find outside point on the left */
				PSL_plotpoint (PSL, x_ext, y_ext, PSL_MOVE);	/* Lay down new outside point */
				PSL_plotpoint (PSL, x_cross[0], y_cross[0], PSL_DRAW);	/* ...then jump to the left boundary */
			}
			close = false;
		}
		if (pen[i] & PSL_CLIP) { 	/* Must add extra point before/after the next point */
			if (pen[i] & PSL_MOVE) {	/* Start of new line, so find external point and lay down before clipped point */
				gmtplot_get_outside_point_extension (GMT, x[i], y[i], x[i+1], y[i+1], &x_ext, &y_ext);
				PSL_plotpoint (PSL, x_ext, y_ext, PSL_MOVE);	/* Lay down new start point */
				PSL_plotpoint (PSL, x[i], y[i], PSL_DRAW);
			}
			else {	/* End of a line, add extension after placing the last point of the line  */
				PSL_plotpoint (PSL, x[i], y[i], PSL_DRAW);
				gmtplot_get_outside_point_extension (GMT, x[i], y[i], x[i-1], y[i-1], &x_ext, &y_ext);
				PSL_plotpoint (PSL, x_ext, y_ext, PSL_DRAW);	/* Lay down new start point */
			}
		}
		else	/* Just another point in the middle of the line */
			PSL_plotpoint (PSL, x[i], y[i], pen[i]);
		i++;
	}
	PSL_command (PSL, close ? "P S\n" : "S\n");
}

void gmt_xy_axis2 (struct GMT_CTRL *GMT, double x0, double y0, double length, double val0, double val1, struct GMT_PLOT_AXIS *A, bool below, bool annotate, unsigned side) {
    /* Only used in psscale.c.  Because gmt_xy_axis does not do gridlines (done in gmt_map_basemap at a higher level),
     * we do that separately in psscale.c */
	if (annotate) side |= GMT_AXIS_BARB;
	gmt_xy_axis (GMT, x0, y0, length, val0, val1, A, below, side);
}

void gmt_linearx_grid (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double w, double e, double s, double n, double dval) {
	double *x = NULL, ys, yn, p_cap = 0.0, cap_start[2] = {0.0, 0.0}, cap_stop[2] = {0.0, 0.0};
	unsigned int idup = 0, i, nx;
	bool cap = false;
	char *type = (gmt_M_x_is_lon (GMT, GMT_IN)) ? "meridian" : "x";

	/* Do we have duplicate e and w boundaries ? */
	idup = (gmt_M_is_azimuthal(GMT) && doubleAlmostEqual(e-w, 360.0)) ? 1 : 0;

	if (gmt_M_pole_is_point(GMT)) {	/* Might have two separate domains of gridlines */
		if (GMT->current.proj.projection_GMT == GMT_POLAR) {	/* Different for polar graphs since "lat" = 0 is at the center */
			ys = cap_stop[0] = cap_stop[1] = p_cap = 90.0 - GMT->current.setting.map_polar_cap[0];
			yn = n;
			cap_start[0] = cap_start[1] = 0.0;
		}
		else {
			p_cap = GMT->current.setting.map_polar_cap[0];
			ys = MAX (s, -p_cap);
			yn = MIN (n, p_cap);
			cap_start[0] = s;
			cap_stop[0]  = -p_cap ;
			cap_start[1] = p_cap;
			cap_stop[1]  = n;
		}
		cap = !doubleAlmostEqual (GMT->current.setting.map_polar_cap[0], 90.0);
	}
	else {
		ys = s;
		yn = n;
	}
	nx = gmtlib_linear_array (GMT, w, e, dval, GMT->current.map.frame.axis[GMT_X].phase, &x);
	if (idup && !gmt_M_360_range(x[0],x[nx-1])) idup = 0;	/* Probably due to phase we don't need to remove any duplicate */
	for (i = 0; i < nx - idup; i++)  {
		GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Draw %s = %g from %g to %g\n", type, x[i], ys, yn);
		gmtplot_map_lonline (GMT, PSL, x[i], ys, yn);
	}
	if (nx) gmt_M_free (GMT, x);

	if (cap) {	/* Also draw the polar cap(s) */
		nx = 0;
		if (s < -GMT->current.setting.map_polar_cap[0]) {	/* Must draw some or all of the S polar cap */
			nx = gmtlib_linear_array (GMT, w, e, GMT->current.setting.map_polar_cap[1], GMT->current.map.frame.axis[GMT_X].phase, &x);
			for (i = 0; i < nx - idup; i++) {
				GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Draw S polar cap %s = %g from %g to %g\n", type, x[i], ys, yn);
				gmtplot_map_lonline (GMT, PSL, x[i], cap_start[0], cap_stop[0]);
			}
			gmtplot_map_latline (GMT, PSL, -p_cap, w, e);
		}
		if (n > GMT->current.setting.map_polar_cap[0]) {	/* Must draw some or all of the N polar cap */
			if (nx == 0) nx = gmtlib_linear_array (GMT, w, e, GMT->current.setting.map_polar_cap[1], GMT->current.map.frame.axis[GMT_X].phase, &x);
			for (i = 0; i < nx - idup; i++) {
				GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Draw N polar cap %s = %g from %g to %g\n", type, x[i], ys, yn);
				gmtplot_map_lonline (GMT, PSL, x[i], cap_start[1], cap_stop[1]);
			}
			gmtplot_map_latline (GMT, PSL, p_cap, w, e);
		}
		if (nx) gmt_M_free (GMT, x);
	}
}

GMT_LOCAL void gmtplot_check_primary_secondary (struct GMT_CTRL *GMT) {
	/* We requires the primary interval to be finer than the secondary interval (when given).
	 * However, it is easy to get this mixed up. Here, we make the check and if we find
	 * the primary interval is larger than the secondary we warn. */
	struct GMT_PLOT_AXIS *A = NULL;
	struct GMT_PLOT_AXIS_ITEM *P = NULL, *S = NULL;
	static char *kind[3] = {"annotation", "tick", "grid-line"};
	static char *axis[2][3] = { {"x", "y", "z"}, {"longitude", "latitude", "z"}};
	unsigned int no, k, type;
	double dP, dS;

	type = gmt_M_is_geographic (GMT, GMT_IN);
	for (no = 0; no <= GMT_Z; no++) {
		if (no < GMT_Z && !(GMT->current.map.frame.side[no] && GMT->current.map.frame.side[no+2])) continue;	/* That axis will not be annotated */
		A = &GMT->current.map.frame.axis[no];
		if (A->type == GMT_TIME) continue;	/* We assume those are set correctly */
		for (k = 0; k < 3; k++) {	/* For each axis */
			P = &(A->item[2*k]);	/* Primary item */
			S = &(A->item[2*k+1]);	/* Secondary item */
			if ((P->active + S->active) < 2) continue;	/* Primary and secondary not both set */
			if (P->special || S->special) continue;		/* Primary and/or secondary are custom so no fixed interval to compare with */
			/* Here they are both set, check the intervals */
			dP = gmtlib_get_map_interval (GMT, A->type, P);
			dS = gmtlib_get_map_interval (GMT, A->type, S);
			if (dP > dS) {	/* Must warn since primary should be the finer-grained interval */
				GMT_Report (GMT->parent, GMT_MSG_WARNING, "Your primary %s %s interval exceeds the secondary interval.\n", axis[type][no], kind[k]);
				GMT_Report (GMT->parent, GMT_MSG_WARNING, "GMT expects it to be the other way around (primary annotations are closest to axis, secondary are further away)\n");
				GMT_Report (GMT->parent, GMT_MSG_WARNING, "Consider correcting your command - proceeding with your selections\n");
			}
		}
	}
}


GMT_LOCAL void gmtplot_map_griditems (struct GMT_CTRL *GMT) {
	double w, e, s, n;
	struct PSL_CTRL *PSL= GMT->PSL;

	if (GMT->current.map.frame.gridline_plotted) return;	/* Already plotted */
	if (GMT->current.map.frame.order == GMT_BASEMAP_BEFORE && GMT->current.map.frame.basemap_flag & GMT_BASEMAP_GRID_AFTER) return;	/* Wrong order */
	if (GMT->current.map.frame.order == GMT_BASEMAP_AFTER  && !(GMT->current.map.frame.basemap_flag & GMT_BASEMAP_GRID_AFTER)) return;	/* Wrong order */

	if (GMT->common.B.active[GMT_PRIMARY] && GMT->common.B.active[GMT_SECONDARY]) {
		/* Make sure primary intervals are < than secondary intervals, otherwise we swap them */
		gmtplot_check_primary_secondary (GMT);
	}

	w = GMT->common.R.wesn[XLO], e = GMT->common.R.wesn[XHI], s = GMT->common.R.wesn[YLO], n = GMT->common.R.wesn[YHI];

	PSL_comment (PSL, "Start of gridlines - if any\n");

	PSL_setdash (PSL, NULL, 0);	/* To ensure no dashed pens are set prior */

	if (GMT->current.proj.got_azimuths) gmt_M_uint_swap (GMT->current.map.frame.side[E_SIDE], GMT->current.map.frame.side[W_SIDE]);	/* Temporary swap to trick justify machinery */

	gmtplot_map_gridlines (GMT, PSL, w, e, s, n);	/* At most only one of these three would kick in */
	gmtplot_map_gridcross (GMT, PSL, w, e, s, n);
	gmtplot_map_gridticks (GMT, PSL, w, e, s, n);

	if (GMT->current.proj.got_azimuths) gmt_M_uint_swap (GMT->current.map.frame.side[E_SIDE], GMT->current.map.frame.side[W_SIDE]);	/* Undo temporary swap */

	GMT->current.map.frame.gridline_plotted = true;	/* Since gmtplot_map_gridlines is called in gmt_map_basemap we flag if we already have done this step */
}

GMT_LOCAL void gmtplot_map_tick_marks (struct GMT_CTRL *GMT) {
	double w, e, s, n;
	struct PSL_CTRL *PSL= GMT->PSL;

	if (GMT->current.map.frame.order == GMT_BASEMAP_BEFORE && GMT->current.map.frame.basemap_flag & GMT_BASEMAP_ANNOT_AFTER) return;	/* Wrong order */
	if (GMT->current.map.frame.order == GMT_BASEMAP_AFTER  && !(GMT->current.map.frame.basemap_flag & GMT_BASEMAP_ANNOT_AFTER)) return;	/* Wrong order */

	w = GMT->common.R.wesn[XLO], e = GMT->common.R.wesn[XHI], s = GMT->common.R.wesn[YLO], n = GMT->common.R.wesn[YHI];

	gmtplot_map_tickmarks (GMT, PSL, w, e, s, n);
}

GMT_LOCAL void gmtplot_map_annotations (struct GMT_CTRL *GMT) {
	double w, e, s, n;
	struct PSL_CTRL *PSL= GMT->PSL;

	if (GMT->current.map.frame.order == GMT_BASEMAP_BEFORE && GMT->current.map.frame.basemap_flag & GMT_BASEMAP_ANNOT_AFTER) return;	/* Wrong order */
	if (GMT->current.map.frame.order == GMT_BASEMAP_AFTER  && !(GMT->current.map.frame.basemap_flag & GMT_BASEMAP_ANNOT_AFTER)) return;	/* Wrong order */

	w = GMT->common.R.wesn[XLO], e = GMT->common.R.wesn[XHI], s = GMT->common.R.wesn[YLO], n = GMT->common.R.wesn[YHI];

	gmtplot_map_annotate (GMT, PSL, w, e, s, n);
}

void gmtplot_title_breaks_decode (struct GMT_CTRL *GMT, const char *in_string, char *out_string) {
	/* Deal with long-form @^ or <break> strings in title and subtitle and replace with GMT_ASCII_GS */
	unsigned int i, o, kl[2] = {2, 7}, id;
	char *kw[2] = {"@^", "<break>"};
	gmt_M_unused (GMT);
	if (in_string[0] == '\0') return;	/* Got nothing */
	if (strstr (in_string, kw[1]))
		id = 1;
	else if (strstr (in_string, kw[0]))
		id = 0;
	else {	/* No markers given */
		strncpy (out_string, in_string, GMT_LEN256);
		return;
	}
	/* Here we must replace kl[id] with GMT_ASCII_GS */
	for (i = o = 0; i < strlen (in_string); i++) {
		if (!strncmp (&in_string[i], kw[id], kl[id]))
			out_string[o++] = GMT_ASCII_GS, i += (kl[id] - 1);	/* Skip one less, allowing for i++ */
		else
			out_string[o++] = in_string[i];
	}
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Converted %s to %s\n", in_string, out_string);
}

GMT_LOCAL double gmtplot_place_latex_eps (struct GMT_CTRL *GMT, double x, double y, struct GMT_FONT *F, const char *string) {
	/* Detected LaTeX commands, i.e., "....@[LaTeX...@[ ..." or  "....<math>LaTeX...</math> ..." */
	bool pos_set = (gmt_M_is_zero (x) && gmt_M_is_zero (y));
	double w, h;
	unsigned char *eps = NULL;
	struct imageinfo header;

	if ((eps = gmtplot_latex_eps (GMT, F, string, &header)) == NULL) {
		GMT_Report (GMT->parent, GMT_MSG_ERROR, "Conversion of LaTeX string to EPS failed\n");
		return 0;	/* Done */
	}
	/* Scale up EPS dimensions by the ratio of title font size to LaTeX default size of 10p */
	w = (header.width  / 72.0) * (F->size / 10.0);
	h = (header.height / 72.0) * (F->size / 10.0);
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "gmtplot_place_latex_eps: Conversion of LaTeX string gave dimensions %g x %g\n", w, h);
	/* Place EPS file instead of text */
	PSL_command (GMT->PSL, "V\n");	/* Keep the relative changes inside a save/restore block */
	if (pos_set)
		PSL_command (GMT->PSL, "currentpoint T\n");	/* Translate to currentpoint since already set by calling function */
	PSL_plotlatexeps (GMT->PSL, x, y, w, h, PSL_BC, eps, F->fill.rgb, &header);
	PSL_command (GMT->PSL, "U\n");
	PSL_free (eps);
	return h;
}

void gmt_map_title (struct GMT_CTRL *GMT, double x, double y) {
	/* Place plot title and optionally subtitle, which both may be a single line or multiple lines.
	 * Note, when x = y = 0 it means current point has already been selected so we must store that and keep
	 * moving up for each line in the multi-line title.
	 */
	bool pos_set = (gmt_M_is_zero (x) && gmt_M_is_zero (y)), head_latex = false;
	double sign = (pos_set) ? -1.0 : +1.0, y_next = 0.0, line_spacing;
	unsigned int n_breaks_T = 0, n_breaks_S, k, form;
	char *word = NULL, title[GMT_LEN256] = {""}, subtitle[GMT_LEN256] = {""}, sep[2] = {""};
	struct PSL_CTRL *PSL= GMT->PSL;

	if (!GMT->current.map.frame.header[0]) return;	/* No title given, so cannot have subtitle either */

	if ((head_latex = gmt_text_is_latex (GMT, GMT->current.map.frame.header)) && gmtplot_has_title_breaks (GMT, GMT->current.map.frame.header)) {
		GMT_Report (GMT->parent, GMT_MSG_ERROR, "LaTeX expressions are only allowed in single-line titles\n");
		return;	/* Done */
	}
	if (GMT->current.map.frame.sub_header[0] && gmt_text_is_latex (GMT, GMT->current.map.frame.sub_header) && gmtplot_has_title_breaks (GMT, GMT->current.map.frame.sub_header)) {
		GMT_Report (GMT->parent, GMT_MSG_ERROR, "LaTeX expressions are only allowed in single-line subtitles\n");
		return;	/* Done */
	}

	/* OK, here we know that if there is LaTeX it is only a single-line string, else it could be multi-line titles and subtitles */

	gmtplot_title_breaks_decode (GMT, GMT->current.map.frame.header, title);
	gmtplot_title_breaks_decode (GMT, GMT->current.map.frame.sub_header, subtitle);

	n_breaks_T = gmt_char_count (title, GMT_ASCII_GS);		/* Is there a title spilling over several lines */
	n_breaks_S = gmt_char_count (subtitle, GMT_ASCII_GS);	/* Is there a subtitle spilling over several lines */

	if (!(n_breaks_T || n_breaks_S || subtitle[0])){	/* Just a single title string on one line */
		if (gmt_text_is_latex (GMT, title)) {
			/* Detected LaTeX commands, i.e., "....@[LaTeX...@[ ..." or  "....<math>LaTeX...</math> ..." */
			(void)gmtplot_place_latex_eps (GMT, x, y, &GMT->current.setting.font_title, title);
		}
		else {
			form = gmt_setfont (GMT, &GMT->current.setting.font_title);
			PSL_plottext (PSL, x, y, sign * GMT->current.setting.font_title.size, GMT->current.map.frame.header, 0.0, -PSL_BC, form);
			GMT->current.map.frame.plotted_header = true;
		}
		return;	/* Done */
	}

	sep[0] = GMT_ASCII_GS;
	/* Must put everything inside a gsave/grestore block */
	PSL_command (PSL, "V\n");	/* Keep the relative changes inside a save/restore block */
	if (pos_set) PSL_command (PSL, "currentpoint /PSL_text_y edef /PSL_text_x edef\n");	/* Remember the position set for the text */

	/* If there are subtitles then these must be placed first since we start from base and move up/backwards */

	if (subtitle[0]) {	/* Plot a subtitle over one or more several lines */
		if (gmt_text_is_latex (GMT, subtitle)) {
			/* Detected LaTeX commands, i.e., "....@[LaTeX...@[ ..." or  "....<math>LaTeX...</math> ..." */
			double h = gmtplot_place_latex_eps (GMT, x, y, &GMT->current.setting.font_subtitle, subtitle);
			if (!head_latex) h += GMT->current.setting.map_title_offset;
			if (pos_set) {	/* Must move up based on height above initial current point for every line since title will come later */
				y_next += h;
				PSL_command (PSL, "PSL_text_x PSL_text_y M 0 %d G\n", PSL->internal.y0 + (int)lrint (y_next * PSL->internal.y2iy));
			}
			else
				y += h;
		}
		else {
			form = gmt_setfont (GMT, &GMT->current.setting.font_subtitle);
			line_spacing = 1.1 * GMT->current.setting.font_subtitle.size / PSL_POINTS_PER_INCH;
			for (k = 0; k <= n_breaks_S; k++) {
				word = gmt_get_word (subtitle, sep, n_breaks_S - k);	/* Pick from the end going forward */
				PSL_plottext (PSL, x, y, sign * GMT->current.setting.font_subtitle.size, word, 0.0, -PSL_BC, form);
				gmt_M_str_free (word);
				if (pos_set) {	/* Must move up based on height above initial current point for every line since title will come later */
					y_next += line_spacing;
					PSL_command (PSL, "PSL_text_x PSL_text_y M 0 %d G\n", PSL->internal.y0 + (int)lrint (y_next * PSL->internal.y2iy));
				}
				else	/* Just increment the change in y location */
					y += line_spacing;
			}
		}
	}

	if (gmt_text_is_latex (GMT, GMT->current.map.frame.header)) {
		/* Detected LaTeX commands, i.e., "....@[LaTeX...@[ ..." or  "....<math>LaTeX...</math> ..." */
		(void)gmtplot_place_latex_eps (GMT, x, y, &GMT->current.setting.font_title, GMT->current.map.frame.header);
	}
	else {
		/* Now plot the title string(s), i.e., more than one line or LaTeX-free string */
		form = gmt_setfont (GMT, &GMT->current.setting.font_title);
		line_spacing = 1.1 * GMT->current.setting.font_title.size / PSL_POINTS_PER_INCH;
		for (k = 0; k <= n_breaks_T; k++) {
			word = gmt_get_word (title, sep, n_breaks_T - k);	/* Pick from the end going forward */
			PSL_plottext (PSL, x, y, sign * GMT->current.setting.font_title.size, word, 0.0, -PSL_BC, form);
			gmt_M_str_free (word);
			if (k < n_breaks_T) {	/* If there are more lines above this one */
				if (pos_set) {
					y_next += line_spacing;
					PSL_command (PSL, "PSL_text_x PSL_text_y M 0 %d G\n", PSL->internal.y0 + (int)lrint (y_next * PSL->internal.y2iy));
				}
				else	/* Just increment the change in y location */
					y += line_spacing;
			}
		}
	}
	/* Since we are calling gmt_setfont inside a gsave/grestore block we must reset the font and stroke color to "not set" so it will be
	 * correctly executed by the next gmt_setfont or gmt_setcolor calls elsewhere */
	gmt_M_memcpy (PSL->current.rgb[PSL_IS_FONT], GMT->session.no_rgb, 3, double);	/* Reset to -1,-1,-1 since text setting must set the color desired */
	gmt_M_memcpy (PSL->current.rgb[PSL_IS_STROKE], GMT->session.no_rgb, 3, double);	/* Reset to -1,-1,-1 since text setting must set the color desired */
	PSL_command (PSL, "U\n");

	GMT->current.map.frame.plotted_header = true;
}

void gmt_map_basemap (struct GMT_CTRL *GMT) {
	/* This function is usually called twice by modules: Once before data-plotting starts and
	 * once after all data-plotting has ended.  This is because different modules have different
	 * needs related to visibility and order of items. Each call will consult a bitflag that determines
	 * if a particular aspect should be plotted this time.  The three items are:
	 *
	 * 1. The base frame (fancy or plain). Flag 0 means plot before data, 1 means after data
	 * 2. The grid lines or grid ticks. Flag 0 means plot before data, 2 means after data
	 * 3. Frame annotations and tick marks. Flag 0 means plot before data, 4 means after data
	 *
	 * In addition, your -B selections may not actually include all of those choices, of course.
	 */

	unsigned int side;
	bool clip_on = false, too_crowded = false;
	char *order[2] = {"before", "after"};
	struct PSL_CTRL *PSL= GMT->PSL;

	/* 0. Determine if we need to be here and set a few parameters */

	/* If a user only gave -Bx<stuff> OR -By<stuff> then we override the -BWESN settings to turn off the unspecified axes */
	if (!GMT->current.map.frame.set_frame[GMT_PRIMARY]) {
		if (GMT->current.map.frame.set[GMT_Y] && !GMT->current.map.frame.set[GMT_X]) GMT->current.map.frame.side[S_SIDE] = GMT->current.map.frame.side[N_SIDE] = GMT_AXIS_NONE;
		if (GMT->current.map.frame.set[GMT_X] && !GMT->current.map.frame.set[GMT_Y]) GMT->current.map.frame.side[W_SIDE] = GMT->current.map.frame.side[E_SIDE] = GMT_AXIS_NONE;
	}

	if (!GMT->common.B.active[GMT_PRIMARY] && !GMT->common.B.active[GMT_SECONDARY]) return;	/* No frame annotation/ticks/gridlines specified */

	if (GMT->common.B.active[GMT_PRIMARY] && GMT->common.B.active[GMT_SECONDARY]) {
		/* Make sure primary intervals are < than secondary intervals, otherwise we swap them */
		gmtplot_check_primary_secondary (GMT);
	}

	if (GMT->current.setting.map_annot_oblique & GMT_OBL_ANNOT_LON_HORIZONTAL) GMT->current.map.frame.horizontal = 2;
	if (GMT->current.map.frame.horizontal == 2) GMT->current.setting.map_annot_oblique |= GMT_OBL_ANNOT_LON_HORIZONTAL;
	if (GMT->current.setting.map_frame_type & GMT_IS_GRAPH && gmt_M_is_geographic (GMT, GMT_IN)) GMT->current.setting.map_frame_type = GMT_IS_PLAIN;
	if (GMT->current.setting.map_frame_type & GMT_IS_FANCY && !gmtplot_is_fancy_boundary(GMT)) GMT->current.setting.map_frame_type = GMT_IS_PLAIN;

	PSL_comment (PSL, "Start of basemap (placed %s the plot contents)\n", order[GMT->current.map.frame.order]);

	PSL_setdash (PSL, NULL, 0);	/* To ensure no dashed pens are set prior */

	/* These three commands resets the memory of PSL regarding pen width, color, and outline */
	gmt_M_memcpy (PSL->current.rgb[PSL_IS_STROKE], GMT->session.no_rgb, 3, double);	/* Reset to -1,-1,-1 so it can be reset below */
	PSL->current.linewidth = GMT_NOTSET;	/* For a reset of internal setting in PSL */
	PSL->current.outline = GMT_NOTSET;		/* Will now be changed by first PSL_setfill */

	if (GMT->current.proj.three_D && GMT->current.map.frame.drawz) GMT->current.map.frame.plotted_header = true;	/* Just so it is not plotted by gmtplot_map_boundary first */

	if (GMT->current.proj.got_azimuths) gmt_M_uint_swap (GMT->current.map.frame.side[E_SIDE], GMT->current.map.frame.side[W_SIDE]);	/* Temporary swap to trick justify machinery */

	if (GMT->current.setting.map_frame_type & GMT_IS_INSIDE && !GMT->current.map.frame.header[0]) {
		gmt_map_clip_on (GMT, GMT->session.no_rgb, 3);	/* Must clip to ensure things are inside */
		clip_on = true;
	}

	/* 1. Lowest plot level is to place grid lines, ticks, crosses, if they are requested and expected at the current order.
	 *    For some modules, gridlines are placed on top of the data (polygons, land-masses). */

	gmtplot_map_griditems (GMT);

	/* 2. Next is tick marks */

	gmtplot_map_tick_marks (GMT);

	/* 3. Next is map frame */

	gmtplot_map_boundary (GMT);	/* This sets frame.side[] = true|false so MUST come before gmtplot_map_annotations */

	/* 4. End with annotations, if requested in the current order */

	gmtplot_map_annotations (GMT);

	if (clip_on) gmt_map_clip_off (GMT);

	PSL_setdash (PSL, NULL, 0);

	/* 4. Undo various temporary changes */

	if (GMT->current.proj.got_azimuths) gmt_M_uint_swap (GMT->current.map.frame.side[E_SIDE], GMT->current.map.frame.side[W_SIDE]);	/* Undo temporary swap */

	if (GMT->current.proj.three_D && GMT->current.map.frame.drawz) GMT->current.map.frame.plotted_header = false;	/* Now we can plot the title [if selected via -B+t] */

	if ((GMT->current.map.frame.order == GMT_BASEMAP_BEFORE && (GMT->current.plot.mode_3D & 1)) || (GMT->current.map.frame.order == GMT_BASEMAP_AFTER && (GMT->current.plot.mode_3D & 2)))
		gmt_vertical_axis (GMT, GMT->current.plot.mode_3D);

	PSL_comment (PSL, "End of basemap (placed %s the plot contents)\n", order[GMT->current.map.frame.order]);

	for (side = 0; side < 4; side++) {	/* Reset annotation crowdedness arrays */
		if (GMT_n_annotations[side]) {
			gmt_M_free (GMT, GMT_x_annotation[side]);
			gmt_M_free (GMT, GMT_y_annotation[side]);
			GMT_n_annotations[side] = GMT_alloc_annotations[side] = 0;
		}
		if (GMT_n_annotations_skip[side]) {
			static char *name[4] = {"bottom", "right", "top", "left"};
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "%d annotations along the %s border were skipped due to crowding.\n", GMT_n_annotations_skip[side], name[side]);
			too_crowded = true;
			GMT_n_annotations_skip[side] = 0;
		}
	}
	if (too_crowded) {
		GMT_Report (GMT->parent, GMT_MSG_WARNING, "Crowding decisions is controlled by MAP_ANNOT_MIN_SPACING, currently set to %s.\n", GMT->current.setting.map_annot_min_spacing_txt);
		GMT_Report (GMT->parent, GMT_MSG_WARNING, "Decrease or increase MAP_ANNOT_MIN_SPACING to see more or fewer annotations, with 0 showing all annotations.\n");
	}
	PSL_setcolor (PSL, GMT->current.setting.map_default_pen.rgb, PSL_IS_STROKE);

	if (GMT->current.map.frame.order == GMT_BASEMAP_AFTER) {	/* Undo at end in case of multi processes */
		GMT->current.map.frame.gridline_plotted = false;
		GMT->current.map.frame.basemap_flag = GMT_BASEMAP_BEFORE;
	}
	else
		GMT->current.map.frame.order = GMT_BASEMAP_AFTER;	/* Move to next order */
}

void gmt_set_basemap_orders (struct GMT_CTRL *GMT, unsigned int frame, unsigned int grid, unsigned int annot) {
	/* Helper function to fill out the basemap flags for the calling module based on its peculiarities and initial settings */
	/* First apply some general over-ruling depending on 3-D and inside annotations */
	static char *place[2] = {"below", "above"};
	if (GMT->current.proj.three_D && GMT->common.J.zactive) {
		frame = GMT_BASEMAP_FRAME_BEFORE;	/* In true 3-D plots we must lay down the x-y frame first regardless of desire to place it at the end */
		annot = GMT_BASEMAP_ANNOT_BEFORE;
		grid  = GMT_BASEMAP_GRID_BEFORE;
	}
	else if (GMT->current.setting.map_frame_type == GMT_IS_INSIDE)	/* Must do annotations and ticks at end since inside the map */
		annot = GMT_BASEMAP_ANNOT_AFTER;

		/* Finally, since ticks overprint on top of the frame we make sure annot/ticks are never done after the frame */
	if (annot == GMT_BASEMAP_ANNOT_AFTER && frame == GMT_BASEMAP_FRAME_BEFORE)
		annot = GMT_BASEMAP_ANNOT_BEFORE;
	GMT->current.map.frame.basemap_flag = frame + grid + annot;
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Basemap order: Frame = %s  Grid = %s  Tick/Annot = %s\n", place[frame], place[grid/2], place[annot/4]);
}

GMT_LOCAL bool gmtplot_z_axis_side (struct GMT_CTRL *GMT, unsigned int axis, unsigned int quadrant) {
	bool below;
	gmt_M_unused(GMT);
	axis++;	/* 1-4 */
	switch (axis) {
		case 1: below = (quadrant == 2);	break;
		case 2: below = (quadrant != 3);	break;
		case 3: below = (quadrant != 2);	break;
		case 4: below = (quadrant == 3);	break;
		default: below = true;	/* Just to avoid Coverity issues */
	}
	return below;
}

void gmt_vertical_axis (struct GMT_CTRL *GMT, unsigned int mode) {
	/* Mode means: 1 = background walls and title, 2 = foreground walls and axis, 3 = all */
	unsigned int fore, back, old_plane;
	double nesw[4], old_level, xx, yy, az;
	struct PSL_CTRL *PSL= GMT->PSL;

	if (!GMT->current.proj.three_D) return;
	if (!GMT->current.map.frame.drawz) return;

	if (GMT->current.map.frame.axis[GMT_Z].item[GMT_GRID_UPPER].active || GMT->current.map.frame.axis[GMT_Z].item[GMT_GRID_LOWER].active) {
		GMT->current.map.frame.draw_box |= GMT_3D_WALL;
		GMT->current.map.frame.draw_wall = true;
	}

	nesw[0] = GMT->current.proj.rect[YHI], nesw[1] = GMT->current.proj.rect[XHI], nesw[2] = GMT->current.proj.rect[YLO], nesw[3] = GMT->current.proj.rect[XLO];

	fore = mode & 2, back = mode & 1;

	/* Since this routine messes with the perspective, we save the state first, then restore later */
	old_plane = GMT->current.proj.z_project.plane;
	old_level = GMT->current.proj.z_project.level;

	/* Vertical walls */

	if (GMT->current.map.frame.draw_box) {
		PSL_setfill (PSL, GMT->session.no_rgb, 1);
		gmt_setpen (GMT, &GMT->current.setting.map_grid_pen[GMT_PRIMARY]);
		if (back) {	/* Draw back walls and cube boxes */
			gmtplot_vertical_wall (GMT, PSL, GMT->current.proj.z_project.quadrant + 1, nesw, true, GMT->current.map.frame.draw_box);
			gmtplot_vertical_wall (GMT, PSL, GMT->current.proj.z_project.quadrant + 2, nesw, true, GMT->current.map.frame.draw_box);
			if (GMT->current.map.frame.draw_box & GMT_3D_BOX) {
				gmtplot_cube_box (GMT, PSL, GMT->current.proj.z_project.quadrant + 1, nesw);
				gmtplot_cube_box (GMT, PSL, GMT->current.proj.z_project.quadrant + 2, nesw);
			}
		}
		if (fore) {	/* Draw front walls and cube boxes */
			gmtplot_vertical_wall (GMT, PSL, GMT->current.proj.z_project.quadrant + 3, nesw, false, GMT->current.map.frame.draw_box);
			gmtplot_vertical_wall (GMT, PSL, GMT->current.proj.z_project.quadrant    , nesw, false, GMT->current.map.frame.draw_box);
			if (GMT->current.map.frame.draw_box & GMT_3D_BOX) {
				gmtplot_cube_box (GMT, PSL, GMT->current.proj.z_project.quadrant + 3, nesw);
				gmtplot_cube_box (GMT, PSL, GMT->current.proj.z_project.quadrant    , nesw);
			}
		}
	}

	/* Vertical axis */

	if (fore && GMT->current.map.frame.side[Z_SIDE]) {
		unsigned int k, n_z, quadrant, corner_to_quadrant[5] = {0, 2, 1, 4, 3}, z_axis[4];	/* Given corner ID 1-4, return quadrant, or vice versa (0 is unused) */
		bool below;
		gmt_M_memcpy (z_axis, GMT->current.map.frame.z_axis, 4, unsigned int);
		for (k = n_z = 0; k < 4; k++)	/* Count # of vertical axes specified; if 0 then we do an auto-select */
			if (z_axis[k]) n_z++;
		if (n_z == 0) z_axis[corner_to_quadrant[GMT->current.proj.z_project.quadrant]-1] = 1;	/* Set the default corner given the quadrant */
		gmt_plane_perspective (GMT, -1, 0.0);
		for (k = 0; k < 4; k++) {
			if (z_axis[k] == 0) continue;	/* Not drawing this vertical axis */
			quadrant = corner_to_quadrant[k+1];	/* Given corner (k+1), return quadrant */
			below = gmtplot_z_axis_side (GMT, k, GMT->current.proj.z_project.quadrant);
			gmt_xyz_to_xy (GMT, nesw[(quadrant/2*2+1)%4], nesw[((quadrant+1)/2*2)%4], GMT->common.R.wesn[ZLO], &xx, &yy);
			/* Restrict reduced azimuth to -45 to 45 range */
			az = GMT->current.proj.z_project.view_azimuth - 90.0 - floor ((GMT->current.proj.z_project.view_azimuth - 45.0) / 90.0) * 90.0;
			PSL_command (PSL, "/PSL_GPP matrix currentmatrix def [%.12g %.12g %.12g %.12g %.12g %.12g] concat\n",
				cosd(az), sind(az) * GMT->current.proj.z_project.sin_el, 0.0, GMT->current.proj.z_project.cos_el, xx * PSL->internal.x2ix, yy * PSL->internal.y2iy);
			gmt_xy_axis (GMT, 0.0, -GMT->common.R.wesn[ZLO], GMT->current.proj.zmax - GMT->current.proj.zmin, GMT->common.R.wesn[ZLO],
				GMT->common.R.wesn[ZHI], &GMT->current.map.frame.axis[GMT_Z], below, GMT->current.map.frame.side[Z_SIDE]);
			PSL_command (PSL, "PSL_GPP setmatrix\n");
		}
	}

	/* Title */

	if (back && GMT->current.map.frame.header[0] && !GMT->current.map.frame.plotted_header) {	/* No header today */
		gmt_plane_perspective (GMT, -1, 0.0);
		gmt_map_title (GMT, 0.5 * (GMT->current.proj.z_project.xmin + GMT->current.proj.z_project.xmax), GMT->current.proj.z_project.ymax + GMT->current.setting.map_title_offset);
	}

	gmt_plane_perspective (GMT, old_plane, old_level);
}

void gmt_map_clip_on (struct GMT_CTRL *GMT, double rgb[], unsigned int flag) {
	/* This function sets up a clip path so that only plotting
	 * inside the map area will be drawn on paper. map_setup
	 * must have been called first.  If r >= 0, the map area will
	 * first be painted in the r,g,b colors specified.  flag can
	 * be 0-3, as described in PSL_beginclipping():
	 * flag : 0 = continue adding pieces to the current clipping path
	 *        1 = start new clipping path (more must follows)
	 *        2 = end clipping path (this is the last segment added)
	 *        3 = this is the complete clipping path (start to end)
	 * 	  Add 4 to select even-odd clipping [nonzero-winding rule].
	 */

	uint64_t np;
	bool donut;
	double *work_x = NULL, *work_y = NULL;
	struct PSL_CTRL *PSL= GMT->PSL;

	np = gmt_map_clip_path (GMT, &work_x, &work_y, &donut);

	PSL_comment (PSL, "Activate Map clip path\n");
	if (donut) {
		PSL_beginclipping (PSL, work_x, work_y, (int)np, rgb, 1);
		PSL_beginclipping (PSL, &work_x[np], &work_y[np], (int)np, rgb, 2);
	}
	else
		PSL_beginclipping (PSL, work_x, work_y, (int)np, rgb, flag);

	gmt_M_free (GMT, work_x);
	gmt_M_free (GMT, work_y);
}

void gmt_map_clip_off (struct GMT_CTRL *GMT) {
	/* Restores the original clipping path for the plot */

	PSL_comment (GMT->PSL, "Deactivate Map clip path\n");
	PSL_endclipping (GMT->PSL, 1);		/* Reduce polygon clipping by one level */
}

void gmt_BB_clip_on (struct GMT_CTRL *GMT, double rgb[], unsigned int flag) {
	/* This function sets up a clip path so that only plotting
	 * inside the bounding box rectangular area will be drawn on paper. map_setup
	 * must have been called first.  If r >= 0, the map area will
	 * first be painted in the r,g,b colors specified.  flag can
	 * be 0-3, as described in PSL_beginclipping():
	 * flag : 0 = continue adding pieces to the current clipping path
	 *        1 = start new clipping path (more must follows)
	 *        2 = end clipping path (this is the last segment added)
	 *        3 = this is the complete clipping path (start to end)
	 * 	  Add 4 to select even-odd clipping [nonzero-winding rule].
	 */

	double work_x[5], work_y[5];
	struct PSL_CTRL *PSL= GMT->PSL;

	work_x[0] = work_x[3] = work_x[4] = GMT->current.proj.rect[XLO];
	work_x[1] = work_x[2] = GMT->current.proj.rect[XHI];
	work_y[0] = work_y[1] = work_y[4] = GMT->current.proj.rect[YLO];
	work_y[2] = work_y[3] = GMT->current.proj.rect[YHI];

	PSL_comment (PSL, "Activate BoundingBox Map clip path\n");
	PSL_beginclipping (PSL, work_x, work_y, 5, rgb, flag);
}

void gmt_setfill (struct GMT_CTRL *GMT, struct GMT_FILL *fill, int outline) {
	struct PSL_CTRL *PSL= GMT->PSL;
	if (!fill) /* NO fill pointer = no fill */
		PSL_setfill (PSL, GMT->session.no_rgb, outline);
	else if (fill->use_pattern) {
		/* Fill with a pattern */
		double rgb[4] = {-3.0, -3.0, -3.0, 0.0};
		rgb[1] = (double)PSL_setimage (PSL, fill->pattern_no, fill->pattern, fill->image, fill->dpi, fill->dim, fill->f_rgb, fill->b_rgb);
		if (rgb[1] < 0.0) {	/* Error in PSL_setimage */
			gmt_M_memset (rgb, 4, double);
			PSL_comment (GMT->PSL, "PSL_setimage failed: Setting fill to black\n");
		}
		PSL_setfill (PSL, rgb, outline);
	}
	else	/* Fill with a color */
		PSL_setfill (PSL, fill->rgb, outline);
}

unsigned int gmt_setfont (struct GMT_CTRL *GMT, struct GMT_FONT *F) {
	/* Set all attributes of the selected font in the PS */
	unsigned int outline;

	PSL_setfont (GMT->PSL, F->id);	/* Set the current font ID */
	if (F->form & 2) {	/* Outline font requested; set pen and fill (rgb[0] == -1 means no fill) */
		gmt_setpen (GMT, &F->pen);		/* Stroke the text outline with this pen */
		gmt_setfill (GMT, &F->fill, 1);	/* Use this color or pattern (if any) for the text fill */
		outline = 1;	/* Indicates outline font is needed and should be stroked */
		if (F->form & 8) outline = 3;	/* Indicates that the outline is allowed to overlap the font */
	}
	else if (F->form & 4) {	/* Want to use a pattern to fill the text but do not draw outline */
		gmt_setfill (GMT, &F->fill, 0);
		outline = 2;	/* Indicates outline font is needed for filling but will not be stroked */
	}
	else {	/* Regular, solid text fill is set via stroke color */
		PSL_setcolor (GMT->PSL, F->fill.rgb, PSL_IS_FONT);
		outline = 0;	/* Indicates we will fill text using "show" which takes current color (i.e., stroke color) */
	}
	return (outline);
}

void gmt_draw_map_inset (struct GMT_CTRL *GMT, struct GMT_MAP_INSET *B, bool clip) {
	/* Place a rectangle on the map, as defined by center point and dimensions or w/e/s/n in geo or projected coordinates */
	unsigned int k;
	double rect[4], dim[3], s;
	struct GMT_MAP_PANEL *panel = B->panel;

	if (B->refpoint) gmt_set_refpoint (GMT, B->refpoint);	/* Finalize reference point plot coordinates, if needed */

	/* First convert the information we have into the center and dimensions of a rectangle */

	if (B->refpoint || B->unit || B->oblique) {	/* Dealing with projected coordinates and dimensions or got oblique box */
		if (B->unit) {
			if (gmt_init_distaz (GMT, B->unit, GMT_GREATCIRCLE, GMT_MAP_DIST) == GMT_NOT_A_VALID_TYPE)	/* Get scales for this unit */
				return;
		}
		if (B->refpoint) {	/* Got a geographic center point and width/height for a rectangular box */
			gmt_M_memcpy (dim, B->dim, 2, double);		/* Duplicate the width/height of rectangle */
			if (B->unit) {	/* Gave dimensioned box */
				for (k = 0; k < 2; k++) {
					dim[k] /= GMT->current.map.dist[GMT_MAP_DIST].scale;	/* Convert units to meters */
					dim[k] *= GMT->current.proj.scale[k];		/* Turns meters into inches on map */
				}
			}
			/* Now in inches */
			gmt_adjust_refpoint (GMT, B->refpoint, dim, B->off, B->justify, PSL_BL);	/* Adjust to bottom left corner */
			rect[XLO] = B->refpoint->x;	rect[XHI] = B->refpoint->x + dim[GMT_X];	/* Get the min/max map coordinates of the rectangle */
			rect[YLO] = B->refpoint->y;	rect[YHI] = B->refpoint->y + dim[GMT_Y];
		}
		else if (B->oblique) {	/* Got lower left and upper right coordinates of rectangular box */
			gmt_geo_to_xy (GMT, B->wesn[XLO], B->wesn[YLO], &rect[XLO], &rect[YLO]);	/* Lower left corner in inches */
			gmt_geo_to_xy (GMT, B->wesn[XHI], B->wesn[YHI], &rect[XHI], &rect[YHI]);	/* Lower left corner in inches */
		}
		else {	/* Got 4 sides in projected units or a rectangular box*/
			gmt_M_memcpy (rect, B->wesn, 4, double);
			for (k = 0; k < 4; k++) rect[k] /= GMT->current.map.dist[GMT_MAP_DIST].scale;	/* Turns units to meters */
			/* Turns meters into inches on map */
			rect[XLO] = rect[XLO] * GMT->current.proj.scale[GMT_X] + GMT->current.proj.origin[GMT_X];
			rect[XHI] = rect[XHI] * GMT->current.proj.scale[GMT_X] + GMT->current.proj.origin[GMT_X];
			rect[YLO] = rect[YLO] * GMT->current.proj.scale[GMT_Y] + GMT->current.proj.origin[GMT_Y];
			rect[YHI] = rect[YHI] * GMT->current.proj.scale[GMT_Y] + GMT->current.proj.origin[GMT_Y];
		}
	}
	else {	/* Got 4 geographic coordinates */
		if (gmt_M_is_rect_graticule (GMT)) {	/* Will give rectangle */
			gmt_geo_to_xy (GMT, B->wesn[XLO], B->wesn[YLO], &rect[XLO], &rect[YLO]);	/* Lower left corner in inches */
			gmt_geo_to_xy (GMT, B->wesn[XHI], B->wesn[YHI], &rect[XHI], &rect[YHI]);	/* Lower left corner in inches */
		}
		else {	/* Curved map inset */
			uint64_t np;
			int outline = ((panel->mode & GMT_PANEL_OUTLINE) == GMT_PANEL_OUTLINE) ? 1 : 0;	/* Does the panel have an outline? */
			double *lon = NULL, *lat = NULL;
			struct GMT_DATASEGMENT *S = GMT_Alloc_Segment (GMT->parent, GMT_NO_STRINGS, 0, 2, NULL, NULL);	/* Just get empty array pointers */
			np = gmt_graticule_path (GMT, &lon, &lat, 1, false, B->wesn[XLO], B->wesn[XHI], B->wesn[YLO], B->wesn[YHI]);
			S->data[GMT_X] = lon;	S->data[GMT_Y] = lat;
			S->n_rows = np;
			if ((panel->mode & 1)) gmt_setfill (GMT, &panel->fill, outline);
			if ((panel->mode & 2)) gmt_setpen (GMT, &panel->pen1);
			gmt_geo_polygons (GMT, S);
			gmt_free_segment (GMT, &S);
			return;	/* Done here */
		}
	}

	/* Deal with rectangular inset */
	/* Determine panel dimensions */

	dim[GMT_X] = rect[XHI] - rect[XLO];	dim[GMT_Y] = rect[YHI] - rect[YLO];
	if (B->refpoint == NULL) {	/* Need to set the BL refpoint since needed in inset to set the temporary origin */
		char fake[GMT_LEN64] = {""};
		snprintf (fake, GMT_LEN64, "x%.16lgi/%.16lgi+jBL", rect[XLO], rect[YLO]);
		if ((B->refpoint = gmt_get_refpoint (GMT, fake, 'D')) == NULL)
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Unable to parse arg %s\n", fake);
	}

	/* Report position and dimensions */
	s = GMT->session.u2u[GMT_INCH][GMT->current.setting.proj_length_unit];
	GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "Map inset lower left corner and dimensions (in %s): %g %g %g %g\n",
		GMT->session.unit_name[GMT->current.setting.proj_length_unit], rect[XLO]*s, rect[YLO]*s, dim[GMT_X]*s, dim[GMT_Y]*s);
	if (B->file) {	/* Save x0 y0 w h to file */
		FILE *fp = fopen (B->file, "w");
		if (fp) {
			fprintf (fp, "%.12g %.12g %.12g %.12g\n", rect[XLO]*s, rect[YLO]*s, dim[GMT_X]*s, dim[GMT_Y]*s);
			fclose (fp);
		}
		else
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Unable to create file %s\n", B->file);
		gmt_M_str_free (B->file);
	}
	if (panel) {	/* Requested to draw a panel */
		panel->width = dim[GMT_X];	panel->height = dim[GMT_Y];
		if (!panel->clearance) gmt_M_memset (panel->padding, 4, double);	/* No clearance is default for map insets unless actually specified */
		gmt_draw_map_panel (GMT, 0.5 * (rect[XHI] + rect[XLO]), 0.5 * (rect[YHI] + rect[YLO]), 3U, panel);
	}
	if (clip) {	/* Set up clip path for this inset */
		double xc[4], yc[4];
		/* Adjust for the padding so that clipping matches the panel rectangle which may be larger than inset */
		xc[0] = xc[3] = rect[XLO];	xc[1] = xc[2] = rect[XHI];
		yc[0] = yc[1] = rect[YLO];	yc[2] = yc[3] = rect[YHI];
		if (panel) {	/* Adjust for clearance, if any */
			xc[0] -= panel->padding[XLO]; xc[3] -= panel->padding[XLO];	xc[1] += panel->padding[XHI]; xc[2] += panel->padding[XHI];
			yc[0] -= panel->padding[YLO]; yc[1] -= panel->padding[YLO];	yc[2] += panel->padding[YHI]; yc[3] += panel->padding[YHI];
		}
		PSL_comment (GMT->PSL, "Start of inset clip path\n");
		PSL_command (GMT->PSL, "clipsave\n");
		PSL_plotline (GMT->PSL, xc, yc, 4, PSL_MOVE | PSL_CLOSE);	/* Must not close path since first point not given ! */
		PSL_command (GMT->PSL, "clip N\n");
		PSL_command (GMT->PSL, "/PSL_inset_clip 1 def\n");	/* Remember that inset clipping is on */
	}
	else	/* No inset clipping set */
		PSL_command (GMT->PSL, "/PSL_inset_clip 0 def\n");

	if (B->translate)	/* Translate the plot origin */
		PSL_setorigin (GMT->PSL, rect[XLO], rect[YLO], 0.0, PSL_FWD);
}

int gmt_draw_map_scale (struct GMT_CTRL *GMT, struct GMT_MAP_SCALE *ms) {
	/* unit_width is array of approximate widths of one space plus the various unit abbreviations in multiples of GMT_LET_WIDTH for Helvetica */
	/* name_width is array of approximate widths of the various labels in multiples of GMT_LET_WIDTH for Helvetica */
	unsigned int form, j;
	double unit_width[GMT_N_UNITS] = {1.1, 1.5, 1.5, 1.5, 2.0, 2.0, 2.0, 1.1, 2.25};
	double name_width[GMT_N_UNITS] = {1.0, 2.3, 4.02, 10.7, 3.2, 2.0, 2.0, 3.05, 8.6};
	enum gmt_enum_units unit = GMT_IS_METER;
	double x1, x2, y1, y2, tx, ty, dist_to_annot, scale_height, x_left, x_right, bar_length_km, dim[4];
	double XL, YL, XR, YR, dist, scl, bar_width, dx, x0_scl, y0_scl;
	char txt[GMT_LEN256] = {""}, format[GMT_LEN64] = {""}, *this_label = NULL;
	/* inch, cm, pt is not used here but in the array since gmtlib_get_unit_number uses this sequence */
	char *label[GMT_N_UNITS] = {"m", "km", "miles", "nautical miles", "inch", "cm", "pt", "feet", "survey feet"};
	char *units[GMT_N_UNITS] = {"m", "km", "mi", "nm", "in", "cm", "pt", "ft", "usft"}, measure;
	struct PSL_CTRL *PSL= GMT->PSL;
	struct GMT_MAP_PANEL *panel = ms->panel;

	if (!ms->plot) return GMT_OK;

	gmt_set_refpoint (GMT, ms->refpoint);	/* Finalize reference point plot coordinates, if needed */

	if (gmt_M_is_cartesian (GMT, GMT_IN))
		bar_length_km = ms->length;	/* Just as is */
	else {
		measure = (ms->measure == 0) ? 'k' : ms->measure;	/* Km is default distance unit */
		if ((unit = gmtlib_get_unit_number (GMT, measure)) == GMT_IS_NOUNIT) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Bad distance unit %c\n", measure);
			return GMT_PARSE_ERROR;
		}
		bar_length_km = 0.001 * GMT->current.proj.m_per_unit[unit] * ms->length;	/* Now in km */
		if (ms->origin_mode == GMT_SCALE_ORIGIN_PLACE) {	/* Pick the lon/lat of the scale placement as map scale origin */
			gmt_xy_to_geo (GMT, &ms->origin[GMT_X], &ms->origin[GMT_Y], ms->refpoint->x, ms->refpoint->y);
			GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "Map scale origin selected to be %g/%g\n", ms->origin[GMT_X], ms->origin[GMT_Y]);
		}
		else if (ms->origin_mode == GMT_SCALE_ORIGIN_MIDDLE) {	/* Pick middle of map as map scale origin */
			gmt_xy_to_geo (GMT, &ms->origin[GMT_X], &ms->origin[GMT_Y], 0.5 * GMT->current.map.width, 0.5 * GMT->current.map.height);
			GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "Map scale origin selected to be %g/%g\n", ms->origin[GMT_X], ms->origin[GMT_Y]);
		}
	}

	/* We will determine the scale at the point ms->origin[GMT_X], ms->origin[GMT_Y] */
	if (gmt_M_is_dnan (ms->origin[GMT_X])) ms->origin[GMT_X] = GMT->current.proj.central_meridian;
	gmt_geo_to_xy (GMT, ms->origin[GMT_X], ms->origin[GMT_Y], &x0_scl, &y0_scl);
	/* 1. Pick a reasonably small dx, e.g., 5% of map width */
	dx = 0.05 * GMT->current.map.width;
	/* 2. Compute test x1, x2 to either side of ms->origin[GMT_X], ms->origin[GMT_Y] */
	x1 = x0_scl - dx;	x2 = x0_scl + dx;
	/* 3. Convert (x1,y0), (x2,y0) to lon,lat coordinates */
	gmt_xy_to_geo (GMT, &XL, &YL, x1, y0_scl);
	gmt_xy_to_geo (GMT, &XR, &YR, x2, y0_scl);
	/* 4. Get distances between (XL,YL) and (XR,YR) */
	if (gmt_M_is_cartesian (GMT, GMT_IN))
		dist = hypot (XL - XR, YL - YR);
	else	/* in km */
		dist = 0.001 * gmt_great_circle_dist_meter (GMT, XL, YL, XR, YR);
	/* Get local scale of desired length to this reference length */
	scl = bar_length_km / dist;
	/* Revise the selection of dx, x1, x2 using this scale.  Use center y as coordinates for a horizontal bar */
	dx *= scl; x1 = ms->refpoint->x - dx;	x2 = ms->refpoint->x + dx;
	y1 = y2 = ms->refpoint->y;

	dim[GMT_X] = bar_width = hypot (x2 - x1, y2 - y1);		/* Width of scale bar in inches */
	dim[GMT_Y] = scale_height = fabs (GMT->current.setting.map_scale_height);	/* Nominal scale bar height */
	dist_to_annot = scale_height + 0.75 * GMT->current.setting.map_annot_offset[GMT_PRIMARY];	/* Dist from top of scalebar to top of annotations */

	gmt_adjust_refpoint (GMT, ms->refpoint, dim, ms->off, ms->justify, PSL_TC);	/* Adjust refpoint to top center */

	x_left  = ms->refpoint->x - 0.5 * bar_width;	/* x-coordinate of leftmost  scalebar point */
	x_right = ms->refpoint->x + 0.5 * bar_width;	/* x-coordinate of rightmost scalebar point */

	if (ms->fancy && gmt_M_is_geographic (GMT, GMT_IN)) {	/* Fancy scale */
		unsigned int i, justify, form;
		unsigned int n_f_ticks[10] = {5, 4, 6, 4, 5, 6, 7, 4, 3, 5};
		unsigned int n_a_ticks[10] = {1, 2, 3, 2, 1, 3, 1, 2, 1, 1};
		int js;
		double base, d_base, dx_f, dx_a, bar_height, bar_tick_len;
		PSL_comment (PSL, "Draw fancy map scale\n");
		/* Based on magnitue of length, determine reasonable annotation spacings */
		js = irint (floor (d_log10 (GMT, ms->length / 0.95)));
		base = pow (10.0, (double)js);
		i = urint (ms->length / base) - 1;
		d_base = ms->length / n_a_ticks[i];
		dx_f = bar_width / n_f_ticks[i];
		dx_a = bar_width / n_a_ticks[i];
		bar_height = 0.5 * fabs (GMT->current.setting.map_scale_height);	/* Height of the black/white checkered fancy bar */
		bar_tick_len = 0.75 * fabs (GMT->current.setting.map_scale_height);	/* Length of tickmarks */
		if (panel && panel->mode) {	/* Place rectangle behind the map scale */
			double x_center, y_center, l_width = 0.0, l_height = 0.0, l_shift = 0.0;

			/* Adjustment for size of largest annotation is ~half its length (= js+1) times approximate dimensions: */
			dim[XLO] = dim[XHI] = 0.5 * (js+1) * GMT_DEC_WIDTH * GMT->current.setting.font_annot[GMT_PRIMARY].size / PSL_POINTS_PER_INCH;
			if (ms->unit) {	/* Adjust for units at end of each annotation */
				dim[XLO] += unit_width[unit] * GMT_LET_WIDTH * GMT->current.setting.font_annot[GMT_PRIMARY].size / PSL_POINTS_PER_INCH;
				dim[XHI] += unit_width[unit] * GMT_LET_WIDTH * GMT->current.setting.font_annot[GMT_PRIMARY].size / PSL_POINTS_PER_INCH;
			}
			if (ms->do_label) {
				/* We estimate the width of the label similarly; if no label given then we use the name_widths from above rather than strlen */
				l_width = bar_tick_len + (ms->label[0] ? strlen (ms->label) : name_width[unit]) * GMT_LET_WIDTH * GMT->current.setting.font_label.size * GMT->session.u2u[GMT_PT][GMT_INCH];
				/* When label is placed to left|right it sticks up by l_shift and we must ensure dim[YHI] is at least that large */
				l_height = GMT_LET_HEIGHT * GMT->current.setting.font_label.size / PSL_POINTS_PER_INCH;	/* Approximate height of label */
				l_shift  = l_height - scale_height;	/* Adjust for the shift in y-coordinate */
			}
			if (ms->alignment == 'l' && l_width > dim[XLO]) dim[XLO] = l_width;	/* Extend rectangle on the left to accommodate the label */
			else if (ms->alignment == 'r' && l_width > dim[XHI]) dim[XHI] = l_width;	/* Extend rectangle on the right to accommodate the label */
			/* Estimate approximate distance from anchor point down to base of annotations */
			dim[YLO] = dist_to_annot + GMT_LET_HEIGHT * GMT->current.setting.font_annot[GMT_PRIMARY].size / PSL_POINTS_PER_INCH;
			dim[YHI] = 0.0;	/* Normally nothing above the scale bar */
			/* If label is above or below bar, add label offset and approximate label height to the space dim */
			if (ms->do_label && ms->alignment == 'b') dim[YLO] += fabs (GMT->current.setting.map_label_offset[GMT_Y]) + l_height;
			else if (ms->do_label && ms->alignment == 't') dim[YHI] += fabs (GMT->current.setting.map_label_offset[GMT_Y]) + l_height;
			else if (ms->do_label && l_shift > dim[YHI]) dim[YHI] = l_shift;
			/* Determine center of gravity for panel */
			x_center = ms->refpoint->x + 0.5 * (dim[XHI] - dim[XLO]);
			y_center = ms->refpoint->y + 0.5 * (dim[YHI] - dim[YLO]);
			/* Determine panel dimensions */
			panel->width = bar_width + dim[XHI] + dim[XLO];	panel->height = dim[YHI] + dim[YLO];
			gmt_draw_map_panel (GMT, x_center, y_center, 3U, panel);
		}
		/* Draw bar scale ticks using tick pen as well as checkerboard using map_default_pen color [black] and page color [white] */
		gmt_setpen (GMT, &GMT->current.setting.map_tick_pen[GMT_PRIMARY]);
		PSL_plotsegment (PSL, x_left, ms->refpoint->y - bar_tick_len, x_left, ms->refpoint->y);
		for (j = 0; j < n_f_ticks[i]; j++) {
			PSL_setfill (PSL, (j%2) ? GMT->PSL->init.page_rgb : GMT->current.setting.map_default_pen.rgb, 1);
			PSL_plotbox (PSL, x_left + j * dx_f, ms->refpoint->y, x_left + (j+1) * dx_f, ms->refpoint->y - bar_height);
			PSL_plotsegment (PSL, x_left + (j+1) * dx_f, ms->refpoint->y - bar_tick_len, x_left + (j+1) * dx_f, ms->refpoint->y);
		}
		/* Place annotations */
		ty = ms->refpoint->y - dist_to_annot;	/* The y-coordinate at the top of the annotations */
		form = gmt_setfont (GMT, &GMT->current.setting.font_annot[GMT_PRIMARY]);
		for (j = 0; j <= n_a_ticks[i]; j++) {
			gmt_sprintf_float (GMT, format, GMT->current.setting.format_float_map, j * d_base);
			if (ms->unit) /* Must append unit */
				snprintf (txt, GMT_LEN256, "%s %s", format, units[unit]);
			else
				snprintf (txt, GMT_LEN256, "%s", format);
			tx = x_left + j * dx_a;	/* center x-coordinate for this annotation */
			PSL_plotsegment (PSL, tx, ms->refpoint->y - scale_height, tx, ms->refpoint->y);	/* Draw the tick mark */
			PSL_plottext (PSL, tx, ty, GMT->current.setting.font_annot[GMT_PRIMARY].size, txt, 0.0, PSL_TC, form);	/* Place annotation */
		}
		if (ms->do_label) {	/* Label was requested */
			/* Determine placement of the label */
			switch (ms->alignment) {
				case 'l':	/* Left */
					tx = x_left - bar_tick_len;
					ty = ms->refpoint->y - scale_height;
					justify = PSL_BR;	/* Left side annotation are right-justified, etc. */
					break;
				case 'r':	/* Right */
					tx = x_right + bar_tick_len;
					ty = ms->refpoint->y - scale_height;
					justify = PSL_BL;
					break;
				case 't':	/* Top */
					tx = ms->refpoint->x;
					ty = ms->refpoint->y + fabs (GMT->current.setting.map_label_offset[GMT_Y]);
					justify = PSL_BC;
					break;
				default:	/* Bottom */
					tx = ms->refpoint->x;
					ty = ms->refpoint->y - dist_to_annot - GMT_LET_HEIGHT * GMT->current.setting.font_annot[GMT_PRIMARY].size / PSL_POINTS_PER_INCH - fabs (GMT->current.setting.map_label_offset[GMT_Y]);
					justify = PSL_TC;
					break;
			}
			this_label = (ms->label[0]) ? ms->label : label[unit];	/* Use user label or default to the unit names */
			form = gmt_setfont (GMT, &GMT->current.setting.font_label);
			PSL_plottext (PSL, tx, ty, GMT->current.setting.font_label.size, this_label, 0.0, justify, form);
		}
	}
	else {	/* Simple scale has no annotation and just the length and unit centered beneath */
		double xp[4], yp[4];	/* Line for simple scale */
		PSL_comment (PSL, "Draw plain map scale\n");
		if (panel && panel->mode) {	/* Place rectangle behind the map scale */
			double x_center, y_center, dim[4];

			/* Adjustment for size of largest annotation is half the length times dimensions: */
			dim[XLO] = dim[XHI] = fabs (GMT->current.setting.map_annot_offset[GMT_PRIMARY]);
			dim[YLO] = dist_to_annot + GMT_LET_HEIGHT * GMT->current.setting.font_annot[GMT_PRIMARY].size / PSL_POINTS_PER_INCH;
			dim[YHI] = 0.0;
			/* Determine center of gravity for panel */
			x_center = ms->refpoint->x + 0.5 * (dim[XHI] - dim[XLO]);
			y_center = ms->refpoint->y + 0.5 * (dim[YHI] - dim[YLO]);
			/* Determine panel dimensions */
			panel->width = bar_width + dim[XHI] + dim[XLO];	panel->height = dim[YHI] + dim[YLO];
			gmt_draw_map_panel (GMT, x_center, y_center, 3U, panel);
		}
		/* Draw the simple scale bar using tick pen */
		gmt_setpen (GMT, &GMT->current.setting.map_tick_pen[GMT_PRIMARY]);
		xp[0] = xp[1] = x_left;	xp[2] = xp[3] = x_right;
		yp[0] = yp[3] = ms->refpoint->y - scale_height;	yp[1] = yp[2] = ms->refpoint->y;
		PSL_plotline (PSL, xp, yp, 4, PSL_MOVE|PSL_STROKE);
		/* Make a basic label using the length and chosen unit and place below the scale */
		gmt_sprintf_float (GMT, format, GMT->current.setting.format_float_map, ms->length);
		if (gmt_M_is_geographic (GMT, GMT_IN))
			snprintf (txt, GMT_LEN256, "%s %s", format, (ms->unit) ? units[unit] : label[unit]);
		else {	/* No unit for Cartesian scales, only label */
			if (ms->do_label)
				snprintf (txt, GMT_LEN256, "%s %s", format, ms->label);
			else
				strcpy (txt, format);
		}
		form = gmt_setfont (GMT, &GMT->current.setting.font_annot[GMT_PRIMARY]);
		PSL_plottext (PSL, ms->refpoint->x, ms->refpoint->y - dist_to_annot, GMT->current.setting.font_annot[GMT_PRIMARY].size, txt, 0.0, PSL_TC, form);
	}
	return GMT_OK;
}

void gmt_draw_vertical_scale (struct GMT_CTRL *GMT, struct GMT_MAP_SCALE *ms) {
	/* Draws a basic vertical scale bar at given refpoint and labels it as specified , */
	double half_scale_length, scale, x0, y0, xx[4], yy[4], off, dim[2], sign = 1.0;
	char txt[GMT_LEN256] = {""};
	int form, just = PSL_ML;

	gmt_set_refpoint (GMT, ms->refpoint);	/* Finalize reference point plot coordinates, if needed */
	/* The ms->origin_mode checks only kick in for geographic plots */
	if (ms->origin_mode == GMT_SCALE_ORIGIN_PLACE) {	/* Pick the lon/lat of the scale placement as map scale origin */
		gmt_xy_to_geo (GMT, &ms->origin[GMT_X], &ms->origin[GMT_Y], ms->refpoint->x, ms->refpoint->y);
		GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "Vertical map scale origin selected to be %g/%g\n", ms->origin[GMT_X], ms->origin[GMT_Y]);
	}
	else if (ms->origin_mode == GMT_SCALE_ORIGIN_MIDDLE) {	/* Pick middle of map as map scale origin */
		gmt_xy_to_geo (GMT, &ms->origin[GMT_X], &ms->origin[GMT_Y], 0.5 * GMT->current.map.width, 0.5 * GMT->current.map.height);
		GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "Vertical map scale origin selected to be %g/%g\n", ms->origin[GMT_X], ms->origin[GMT_Y]);
	}
	if (ms->label[0]) /* Append data unit to the scale length */
		snprintf (txt, GMT_LEN256, "%g %s", ms->length, ms->label);
	else
		snprintf (txt, GMT_LEN256, "%g", ms->length);

	if (ms->zdata) /* Append data unit to the scale length */
		scale = ms->z_scale;
	else
		scale = GMT->current.proj.scale[GMT_Y];
	half_scale_length = 0.5 * ms->length * scale;

	off = ((GMT->current.setting.map_scale_height > 0.0) ? GMT->current.setting.map_tick_length[GMT_PRIMARY] : 0.0) + GMT->current.setting.map_annot_offset[GMT_PRIMARY];
	dim[GMT_X] = strlen (txt) * GMT_DEC_WIDTH * GMT->current.setting.font_annot[GMT_PRIMARY].size / PSL_POINTS_PER_INCH + off;
	dim[GMT_Y] = 2.0 * half_scale_length;

	gmt_adjust_refpoint (GMT, ms->refpoint, dim, ms->off, ms->justify, PSL_ML);	/* Adjust refpoint to ML */

	x0 = ms->refpoint->x;	y0 = ms->refpoint->y;
	if (ms->panel && ms->panel->mode) {	/* Place rectangle behind the map scale */
		double x_center, y_center;
		x_center = x0 + 0.5 * dim[GMT_X]; y_center = y0;
		ms->panel->width = dim[GMT_X];	ms->panel->height = dim[GMT_Y];
		gmt_draw_map_panel (GMT, x_center, y_center, 3U, ms->panel);
	}
	/* Compute the 4 coordinates on the scale line */
	if (ms->alignment == 'l') {	/* Flip scale with label on left */
		x0 += dim[GMT_X];
		just = PSL_MR;
		sign = -1;
	}
	gmt_xyz_to_xy (GMT, x0 + sign * GMT->current.setting.map_scale_height, y0 - half_scale_length, 0.0, &xx[0], &yy[0]);
	gmt_xyz_to_xy (GMT, x0, y0 - half_scale_length, 0.0, &xx[1], &yy[1]);
	gmt_xyz_to_xy (GMT, x0, y0 + half_scale_length, 0.0, &xx[2], &yy[2]);
	gmt_xyz_to_xy (GMT, x0 + sign * GMT->current.setting.map_scale_height, y0 + half_scale_length, 0.0, &xx[3], &yy[3]);
	gmt_setpen (GMT, &GMT->current.setting.map_tick_pen[GMT_PRIMARY]);
	PSL_plotline (GMT->PSL, xx, yy, 4, PSL_MOVE|PSL_STROKE);
	form = gmt_setfont (GMT, &GMT->current.setting.font_annot[GMT_PRIMARY]);
	PSL_plottext (GMT->PSL, x0 + sign * off, y0, GMT->current.setting.font_annot[GMT_PRIMARY].size, txt, 0.0, just, form);
}

void gmt_draw_vertical_scale_old (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL, double x0, double y0, double length, double zscale, int gave_xy, char *units) {
	/* Draws a basic vertical scale bar at (x0,y0) and labels it as specified , struct GMT_MAP_SCALE *ms*/
	int form;
	double dy, off, xx[4], yy[4];
	char txt[GMT_LEN256] = {""};

	gmt_setpen (GMT, &GMT->current.setting.map_tick_pen[GMT_PRIMARY]);

	if (!gave_xy) {	/* Project lon,lat to get position of scale */
		gmt_geo_to_xy (GMT, x0, y0, &xx[0], &yy[0]);
		x0 = xx[0];	y0 = yy[0];
	}

	if (units) /* Append data unit to the scale length */
		snprintf (txt, GMT_LEN256, "%g %s", length, units);
	else
		snprintf (txt, GMT_LEN256, "%g", length);
	dy = 0.5 * length * zscale;
	/* Compute the 4 coordinates on the scale line */
	gmt_xyz_to_xy (GMT, x0 + GMT->current.setting.map_scale_height, y0 - dy, 0.0, &xx[0], &yy[0]);
	gmt_xyz_to_xy (GMT, x0, y0 - dy, 0.0, &xx[1], &yy[1]);
	gmt_xyz_to_xy (GMT, x0, y0 + dy, 0.0, &xx[2], &yy[2]);
	gmt_xyz_to_xy (GMT, x0 + GMT->current.setting.map_scale_height, y0 + dy, 0.0, &xx[3], &yy[3]);
	PSL_plotline (PSL, xx, yy, 4, PSL_MOVE|PSL_STROKE);
	off = ((GMT->current.setting.map_scale_height > 0.0) ? GMT->current.setting.map_tick_length[GMT_PRIMARY] : 0.0) + GMT->current.setting.map_annot_offset[GMT_PRIMARY];
	form = gmt_setfont (GMT, &GMT->current.setting.font_annot[GMT_PRIMARY]);
	PSL_plottext (PSL, x0 + off, y0, GMT->current.setting.font_annot[GMT_PRIMARY].size, txt, 0.0, PSL_ML, form);
}

void gmt_draw_map_rose (struct GMT_CTRL *GMT, struct GMT_MAP_ROSE *mr) {
	int tmp_join, tmp_limit;
	struct PSL_CTRL *PSL= GMT->PSL;
	struct GMT_MAP_PANEL *panel = mr->panel;
	double dim[2];
	if (!mr->plot) return;

	dim[GMT_X] = dim[GMT_Y] = mr->size;
	gmt_set_refpoint (GMT, mr->refpoint);	/* Finalize reference point plot coordinates, if needed */
	gmt_adjust_refpoint (GMT, mr->refpoint, dim, mr->off, mr->justify, PSL_MC);	/* Adjust refpoint to MC */

	if (panel && panel->mode) {	/* Place rectangle behind the map rose */
		/* Determine panel dimensions */
		panel->width = panel->height = mr->size;
		gmt_draw_map_panel (GMT, mr->refpoint->x, mr->refpoint->y, 3U, panel);
	}

	/* Temporarily use miter to get sharp points to compass rose */
	tmp_join = PSL->internal.line_join;	PSL_setlinejoin (PSL, 0);
	tmp_limit = PSL->internal.miter_limit;	PSL_setmiterlimit (PSL, 0);

	if (mr->type == GMT_ROSE_MAG)	/* Do magnetic compass rose */
		gmtplot_draw_mag_rose (GMT, PSL, mr);
	else
		gmtplot_draw_dir_rose (GMT, PSL, mr);

	/* Switch line join style back */
	PSL_setlinejoin (PSL, tmp_join);
	PSL_setmiterlimit (PSL, tmp_limit);
}

void gmt_draw_map_panel (struct GMT_CTRL *GMT, double x, double y, unsigned int mode, struct GMT_MAP_PANEL *P) {
	/* Draw a rectangular backpanel behind things like logos, scales, legends, images.
	 * Here, (x,y) is the center-point of the panel.
	 * mode is a bit flag that can be 1,2, or 3:
	 * mode = 1.  Lay down fills for background (if any fills) but no outlines
	 * mode = 2.  Just draw any outlines requested
	 * mode = 3.  Do both at the same time. */
	double dim[3] = {0.0, 0.0, 0.0};
	int outline;
	struct GMT_FILL *fill = NULL;	/* Default is no fill */
	if (!P) return;	/* No panel given */
	outline = ((P->mode & GMT_PANEL_OUTLINE) == GMT_PANEL_OUTLINE) ? 1 : 0;	/* Does the panel have an outline? */
	GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "Place rectangular back panel\n");
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Clearance: %g/%g/%g/%g\n", P->padding[XLO], P->padding[XLO], P->padding[YLO], P->padding[YHI]);
	dim[GMT_X] = P->width  + P->padding[XLO] + P->padding[XHI];	/* Rectangle width */
	dim[GMT_Y] = P->height + P->padding[YLO] + P->padding[YHI];	/* Rectangle height */
	dim[GMT_Z] = P->radius;	/* Corner radius, or zero */
	/* In case clearances are not symmetric we need to shift the symbol center accordingly */
	x += 0.5 * (P->padding[XHI] - P->padding[XLO]);
	y += 0.5 * (P->padding[YHI] - P->padding[YLO]);
	if (mode == 1) outline = 0;	/* Do not draw outlines (even if requested) at this time since mode == 1*/
	if ((mode & 1) && (P->mode & GMT_PANEL_SHADOW)) {	/* Draw offset background shadow first */
		gmt_setfill (GMT, &P->sfill, 0);	/* The shadow has no outline */
		PSL_plotsymbol (GMT->PSL, x + P->off[GMT_X], y + P->off[GMT_Y], dim, (P->mode & GMT_PANEL_ROUNDED) ? PSL_RNDRECT : PSL_RECT);
	}
	if ((mode & 2) && outline) gmt_setpen (GMT, &P->pen1);	/* Need to set frame outline pen */
	if (mode & 1) fill = &P->fill;		/* Select fill (which may be NULL) unless we are just doing outlines */
	if (fill || outline) gmt_setfill (GMT, fill, outline);	/* Activate frame fill, with optional outline */
	if ((mode&1 && (P->mode & GMT_PANEL_FILL)) || (mode&2 && (P->mode & GMT_PANEL_OUTLINE))) PSL_plotsymbol (GMT->PSL, x, y, dim, (P->mode & GMT_PANEL_ROUNDED) ? PSL_RNDRECT : PSL_RECT);
	if ((mode & 2) && (P->mode & GMT_PANEL_INNER)) {	/* Also draw secondary frame on the inside */
		dim[GMT_X] -= 2.0 * P->gap;	/* Shrink dimension of panel by the uniform gap on all sides */
		dim[GMT_Y] -= 2.0 * P->gap;
		gmt_setpen (GMT, &P->pen2);	/* Set inner border pen */
		gmt_setfill (GMT, NULL, 1);	/* Never fill for inner frame */
		PSL_plotsymbol (GMT->PSL, x, y, dim, (P->mode & GMT_PANEL_ROUNDED) ? PSL_RNDRECT : PSL_RECT);
	}
	/* Reset color */
	PSL_setcolor (GMT->PSL, GMT->current.setting.map_frame_pen.rgb, PSL_IS_STROKE);
}

void gmt_setpen (struct GMT_CTRL *GMT, struct GMT_PEN *pen) {
	/* gmt_setpen issues PostScript code to set the specified pen. */

	if (!pen) return;
	GMT->current.setting.ps_penwidth = pen->width;	/* Remember current pen width */
	PSL_setlinewidth (GMT->PSL, pen->width);
	PSL_setdash (GMT->PSL, pen->style, pen->offset);
	PSL_setcolor (GMT->PSL, pen->rgb, PSL_IS_STROKE);
}

#define GMT_N_COND_LEVELS	10	/* Number of max nesting level for conditionals */

GMT_LOCAL void gmtplot_get_the_pen (struct GMT_PEN *p, struct GMT_CUSTOM_SYMBOL_ITEM *s, struct GMT_PEN *cp, struct GMT_FILL *f) {
	/* Returns the pointer to the pen we should use.  If this is an action pen then
	 * we set restore to true so that we can reset it later.  If var_pen contains a -GMT_USE_FILL_RGB
	 * then we overwrite the pen color with the current fill color */
	if (s->pen) {
		gmt_M_memcpy (p, s->pen, 1, struct GMT_PEN);
	}
	else if (cp) {
		gmt_M_memcpy (p, cp, 1, struct GMT_PEN);
	}
	if (s->var_pen < 0 && (abs(s->var_pen) & GMT_USE_FILL_RGB) && f) gmt_M_rgb_copy (p->rgb, f->rgb);
}

GMT_LOCAL void gmtplot_get_the_fill (struct GMT_FILL *f, struct GMT_CUSTOM_SYMBOL_ITEM *s, struct GMT_PEN *cp, struct GMT_FILL *cf) {
	/* Returns pointer to the chosen fill: The one specified by -G inside the macro or -G on command line.
	 * If var_pen contains a -GMT_USE_PEN_RGB then we copy in the color of the current pen instead */
	if (s->fill)
		gmt_M_memcpy (f, s->fill, 1, struct GMT_FILL);
	else if (cf)
		gmt_M_memcpy (f, cf, 1, struct GMT_FILL);
	else
		f->rgb[0] = -1;	/* No fill */
	if (f->rgb[0] >= 0.0 && s->var_pen < 0 && (abs(s->var_pen) & GMT_USE_PEN_RGB) && cp) gmt_M_rgb_copy (f->rgb, cp->rgb);
}

GMT_LOCAL void gmtplot_draw_eps_symbol (struct GMT_CTRL *GMT, double x0, double y0, double size[], struct GMT_CUSTOM_SYMBOL_ITEM *symbol) {
	/* Special Encapsulated PostScript-only symbol */
	struct GMT_CUSTOM_SYMBOL_EPS *E = symbol->eps;	/* SHorthand */
	double off = 0.5*size[0], fy = (E->BB[3] - E->BB[2]) / (E->BB[1] - E->BB[0]);
	struct PSL_CTRL *PSL= GMT->PSL;
	if (!E->placed) {	/* First time we must dump the PS code definition */
		double scl = 72.0 / (E->BB[1] - E->BB[0]);
		PSL_comment (PSL, "Start of symbol %s\n", E->name);
		PSL_command (PSL, "/Sk_%s {\nPSL_eps_begin\n", E->name);
		/* We use the symbol's bounding box and scale its width to 1 inch since PSL uses inches */
		PSL_command (PSL, "%.8f dup scale\n", scl);
		if (!E->GMT_made)	/* Non-GMT-produced EPS macro - must scale points to GMT's unit */
			PSL_command (PSL, "1200 72 div dup scale\n");
		PSL_command (PSL, "%%%%BeginDocument: %s.eps\n", E->name);
		PSL_copy (PSL, E->macro);	/* Since it may be quite large */
		PSL_command (PSL, "%%%%EndDocument\n");
		PSL_command (PSL, "PSL_eps_end } def\n");
		PSL_comment (PSL, "End of symbol %s\n", E->name);
		E->placed = true;	/* Flag now says we have dumped the EPS code */
	}
	PSL_command (PSL, "V ");
	PSL_setorigin (PSL, x0-off, y0-fy*off, 0.0, PSL_FWD);
	PSL_command (PSL, "%.12g dup scale ", size[0]);
	PSL_command (PSL, "Sk_%s U\n", E->name);
}

GMT_LOCAL bool gmtplot_is_azimuth (struct GMT_CUSTOM_SYMBOL *symbol, struct GMT_CUSTOM_SYMBOL_ITEM *s, unsigned int var_no) {
	/* We either know a particular column has variable azimuths or we were given a constant angle flagged by a trailing 'a' */
    if (s->angular == GMT_IS_AZIMUTH) return true; /* Must convert a constant azimuth to Cartesian angle */
    if (s->is_var[var_no] && symbol->type[s->var[var_no]-1] == GMT_IS_AZIMUTH ) return true; /* Must convert variable azimuth to Cartesian angle */
	return false;	/* Got Cartesian angle */
}

int gmt_draw_custom_symbol (struct GMT_CTRL *GMT, double x0, double y0, double size[], char *tr_text, struct GMT_CUSTOM_SYMBOL *symbol, struct GMT_PEN *pen, struct GMT_FILL *fill, unsigned int outline) {
	int action, ifs, this_outline = 0, error = GMT_NOERROR;
	unsigned int na, i, id = 0, level;
	bool flush = false, skip[GMT_N_COND_LEVELS+1], done[GMT_N_COND_LEVELS+1];
	uint64_t n = 0;
	size_t n_alloc = 0;
	double x, y, lon, lat, az, angle1, angle2, p_width = 0.0, *xx = NULL, *yy = NULL, *xp = NULL, *yp = NULL, dim[PSL_MAX_DIMS];
	char user_text[GMT_LEN256] = {""};
	struct GMT_CUSTOM_SYMBOL_ITEM *s = NULL;
	struct GMT_FILL f, *current_fill = fill;
	struct GMT_PEN save_pen, p, *current_pen = pen;
	struct GMT_FONT font = GMT->current.setting.font_annot[GMT_PRIMARY];
	struct PSL_CTRL *PSL= GMT->PSL;

	/* Regular macro symbol */

	/* Remember current settings as we wish to restore at the end */
	gmtplot_savepen (GMT, &save_pen);
	gmt_M_memset (&f, 1, struct GMT_FILL);
	gmt_M_memset (&p, 1, struct GMT_PEN);
	gmt_M_memset (skip, GMT_N_COND_LEVELS+1, bool);
	gmt_M_memset (done, GMT_N_COND_LEVELS+1, bool);

	if (symbol->text) {	/* This symbol places text, so we must set macros for fonts and fontsizes outside the gsave/grestore around each symbol */
		symbol->text = 0;	/* Only do this formatting once */
		s = symbol->first;	/* Start at first item */
		while (s) {		/* Examine all items for possible text */
			if (s->action == GMT_SYMBOL_TEXT || s->action == GMT_SYMBOL_VARTEXT) {	/* Text item found */
				gmtplot_format_symbol_string (GMT, s, size, user_text);
				if (s->p[0] < 0.0)	/* Fixed point size for text */
					s->font.size = -s->p[0];
				else	/* Fractional size that depends on symbol size */
					s->font.size = s->p[0] * size[0] * PSL_POINTS_PER_INCH;
				/* Set PS macro for fetching this font and size */
				gmtplot_encodefont (PSL, s->font.id, symbol->name, id++);
			}
			s = s->next;
		}
	}

	/* We encapsulate symbol with gsave and translate origin to (x0, y0) first */
	PSL_command (PSL, "V ");
	PSL_setorigin (PSL, x0, y0, 0.0, PSL_FWD);
	gmt_set_meminc (GMT, GMT_SMALL_CHUNK);

	s = symbol->first;
	id = level = 0;
	while (s) {
		if (s->conditional > GMT_BEGIN_SINGLE_IF) {	/* Process if/elseif/else and endif by updating level and skip array, then go to next item */
			/* We keep track of all the nested levels of tests via the skip array.  If a higher level test fails, then we will skip anything inside
			 * it (e.g., lower-level nested test) since those tests don't matter since the upper test failed.  Hence skip is set to true for all deeper
			 * tests (regardless of their actual test result) since we will not get there anyway if the earlier test failed. Finally, when we have
			 * a series if if,elseif,else at the same level then we consult the done[level] array.  This is set to true if we pass a test and actually
			 * draw something and once that is done none of the other tests at the same level can pass. */
			if (s->conditional == GMT_BEGIN_BLOCK_IF) {	/* Beginning of a new if branch. If we are inside an earlier branch whose test was false then all inside shall be false */
				if (level > 0 && skip[level])
					skip[level+1] = true;
				else {
					if ((error = gmtplot_custum_failed_bool_test (GMT, s, size, tr_text, &skip[level+1]))) return (error);
				}
				level++;
				if (level == GMT_N_COND_LEVELS) {
					GMT_Report (GMT->parent, GMT_MSG_ERROR, "Symbol macro (%s) logical nesting too deep [> %d]\n", symbol->name, GMT_N_COND_LEVELS);
					return GMT_DIM_TOO_LARGE;
				}
				done[level] = false;	/* Have not yet taken any action at this level */
			}
			else if (s->conditional == GMT_END_IF) {	/* Simply reduce indent */
				if (level == 0) {
					GMT_Report (GMT->parent, GMT_MSG_ERROR, "Symbol macro (%s) logical nesting error\n", symbol->name);
					return GMT_RUNTIME_ERROR;
				}
				level--;
			}
			else if (s->conditional == GMT_END_IF_ELSE)	/* else branch */
				skip[level] = (!done[level] && skip[level] && !skip[level-1]) ? false : true;
			else if (s->conditional == GMT_BEGIN_ELSEIF) {	/* Skip if prior if/elseif was true, otherwise evaluate test at this level */
				if (!done[level] && skip[level] && !skip[level-1]) {
					if ((error = gmtplot_custum_failed_bool_test (GMT, s, size, tr_text, &skip[level]))) return (error);
				}
				else
					skip[level] = true;
			}
			s = s->next;
			continue;
		}
		if (level && skip[level]) {	/* We are inside an if-block but the block test was false, so we skip */
			s = s->next;
			continue;
		}
		/* Finally, check for 1-line if tests */
		if (s->conditional == GMT_BEGIN_SINGLE_IF) {
			bool retval;
			if ((error = gmtplot_custum_failed_bool_test (GMT, s, size, tr_text, &retval))) return (error);
			if (retval) {
				s = s->next;	/* Done here, move to next item */
				continue;
			}
		}
		done[level] = true;	/* Here we will actually draw something */

		/* Scale coordinates and size parameters by the scale in size[0] */

		x = s->x * size[0];
		y = s->y * size[0];
		gmt_M_memset (dim, PSL_MAX_DIMS, double);
		dim[0] = s->p[0] * size[0];
		dim[1] = s->p[1] * size[0];
		dim[2] = s->p[2] * size[0];
		if (s->pen) {	/* This action has a pen setting */
			p_width = s->pen->width;	/* Remember what it was before messing with it below */
			if (s->pen->width < 0.0)	/* Convert the normalized pen width to points given current size */
				s->pen->width = fabs (s->pen->width * size[0] * GMT->session.u2u[GMT_INCH][GMT_PT]);
			else if (s->var_pen > 0)	/* Convert the specified variable to points */
				s->pen->width = size[s->var_pen] * GMT->session.u2u[GMT_INCH][GMT_PT];
		}
		if (s->action == '?')	/* Reset to what is last in the input record */
			action = GMT->current.io.curr_trailing_text[strlen(GMT->current.io.curr_trailing_text)-1];
		else
			action = s->action;

		switch (action) {
			case GMT_SYMBOL_MOVE:	/* Flush existing polygon and start a new path */
				if (flush) gmtplot_flush_symbol_piece (GMT, PSL, xx, yy, &n, &p, &f, this_outline, &flush);
				n = 0;
				if (n >= n_alloc) gmt_M_malloc2 (GMT, xx, yy, n, &n_alloc, double);
				xx[n] = x, yy[n] = y, n++;
				gmtplot_get_the_pen (&p, s, current_pen, current_fill);
				gmtplot_get_the_fill (&f, s, current_pen, current_fill);
				this_outline = (p.rgb[0] == -1) ? false : outline;
				break;

			case GMT_SYMBOL_STROKE:	/* To force the drawing of a line (outline == 2), not a closed polygon */
				if (flush) gmtplot_flush_symbol_piece (GMT, PSL, xx, yy, &n, &p, &f, 2, &flush);
				n = 0;
				break;

			case GMT_SYMBOL_DRAW:	/* Append another point to the path */
				flush = true;
				if (n >= n_alloc) gmt_M_malloc2 (GMT, xx, yy, n, &n_alloc, double);
				xx[n] = x, yy[n] = y, n++;
				break;

			case GMT_SYMBOL_ARC:	/* Append a circular arc to the path */
				flush = true;
				angle1 = (s->is_var[1]) ? size[s->var[1]] : s->p[1];
				angle2 = (s->is_var[2]) ? size[s->var[2]] : s->p[2];

				if (gmtplot_is_azimuth (symbol, s, 1) || gmtplot_is_azimuth (symbol, s, 2)) {
					angle1 = gmt_azim_to_angle (GMT, lon, lat, 0.1, angle1);
					angle2 = gmt_azim_to_angle (GMT, lon, lat, 0.1, angle2);
				}
				na = gmtlib_get_arc (GMT, x, y, 0.5 * s->p[0] * size[0], angle1, angle2, &xp, &yp);
				for (i = 0; i < na; i++) {
					if (n >= n_alloc) gmt_M_malloc2 (GMT, xx, yy, n, &n_alloc, double);
					xx[n] = xp[i], yy[n] = yp[i], n++;
				}
				gmt_M_free (GMT, xp);
				gmt_M_free (GMT, yp);
				break;

			case GMT_SYMBOL_ROTATE:		/* Rotate the symbol coordinate system by a fixed amount */
				if (flush) gmtplot_flush_symbol_piece (GMT, PSL, xx, yy, &n, &p, &f, this_outline, &flush);
				PSL_setorigin (PSL, 0.0, 0.0, s->p[0], PSL_FWD);
				break;

			case GMT_SYMBOL_AZIMROTATE:	/* Rotate the symbol y-axis to the a fixed azimuth */
				if (flush) gmtplot_flush_symbol_piece (GMT, PSL, xx, yy, &n, &p, &f, this_outline, &flush);
				/* Need to recover actual lon,lat location of symbol first */
				gmt_xy_to_geo (GMT, &lon, &lat, x0, y0);
				angle1 = gmt_azim_to_angle (GMT, lon, lat, 0.1, 90.0 - s->p[0]);
				PSL_setorigin (PSL, 0.0, 0.0, angle1, PSL_FWD);
				break;

			case GMT_SYMBOL_VARROTATE:	/* Rotate the symbol coordinate system by a variable amount */
				if (flush) gmtplot_flush_symbol_piece (GMT, PSL, xx, yy, &n, &p, &f, this_outline, &flush);
				if (symbol->type[s->var[0]-1] == GMT_IS_AZIMUTH) {	/* Must convert azimuth to map angle then  rotatate by -pseudoazimuth */
					gmt_xy_to_geo (GMT, &lon, &lat, x0, y0);
					az = 90.0 - gmt_azim_to_angle (GMT, lon, lat, 0.1, s->var_sign[0] * size[s->var[0]]);
					PSL_setorigin (PSL, 0.0, 0.0, -az, PSL_FWD);
				}
				else	/* Use rotation angle as given */
					PSL_setorigin (PSL, 0.0, 0.0, s->var_sign[0] * size[s->var[0]], PSL_FWD);
				break;

			case GMT_SYMBOL_TEXTURE:	/* Change the current pen/fill settings */
				if (s->fill) current_fill = s->fill;
				if (s->pen) current_pen = s->pen;
				break;

			case (int)'C':
				if (gmt_M_compat_check (GMT, 4)) {	/* Warn and purposefully fall through to assign the rest of the statements */
					GMT_Report (GMT->parent, GMT_MSG_COMPAT, "Circle macro symbol C is deprecated; use c instead\n");
					action = s->action = PSL_CIRCLE;	/* Backwards compatibility, circles are now 'c' */
				}
				else {
					GMT_Report (GMT->parent, GMT_MSG_ERROR,
									"Unrecognized symbol code (%d = '%c') passed to gmt_draw_custom_symbol\n", action, (char)action);
					return GMT_PARSE_ERROR;
					break;
				}
				/* Intentionally fall through - to assign the rest of the statements */
			case PSL_CROSS:
			case PSL_CIRCLE:
			case PSL_SQUARE:
			case PSL_TRIANGLE:
			case PSL_DIAMOND:
			case PSL_STAR:
			case PSL_DOT:
			case PSL_HEXAGON:
			case PSL_OCTAGON:
			case PSL_PENTAGON:
			case PSL_PLUS:
			case PSL_INVTRIANGLE:
			case PSL_RECT:
			case PSL_RNDRECT:
			case PSL_XDASH:
			case PSL_YDASH:
				if (flush) gmtplot_flush_symbol_piece (GMT, PSL, xx, yy, &n, &p, &f, this_outline, &flush);
				gmtplot_get_the_fill (&f, s, current_pen, current_fill);
				gmtplot_get_the_pen (&p, s, current_pen, current_fill);
				this_outline = (p.rgb[0] == -1) ? 0 : outline;
				if (this_outline) gmt_setpen (GMT, &p);
				gmt_setfill (GMT, &f, this_outline);
				PSL_plotsymbol (PSL, x, y, dim, action);
				break;

			case PSL_ELLIPSE:
			case PSL_ROTRECT:
				if (flush) gmtplot_flush_symbol_piece (GMT, PSL, xx, yy, &n, &p, &f, this_outline, &flush);
				gmtplot_get_the_fill (&f, s, current_pen, current_fill);
				gmtplot_get_the_pen (&p, s, current_pen, current_fill);
				this_outline = (p.rgb[0] == -1) ? 0 : outline;
				if (this_outline) gmt_setpen (GMT, &p);
				gmt_setfill (GMT, &f, this_outline);
				angle1 = (s->is_var[0]) ? size[s->var[0]] : s->p[0];
				if (gmtplot_is_azimuth (symbol, s, 0))
					angle1 = gmt_azim_to_angle (GMT, lon, lat, 0.1, angle1);
				dim[0] = angle1;
				if (s->is_var[1]) dim[1] = size[s->var[1]];
				PSL_plotsymbol (PSL, x, y, dim, PSL_ELLIPSE);
				break;

			case PSL_MARC:
				if (flush) gmtplot_flush_symbol_piece (GMT, PSL, xx, yy, &n, &p, &f, this_outline, &flush);
				gmtplot_get_the_fill (&f, s, current_pen, current_fill);
				gmtplot_get_the_pen (&p, s, current_pen, current_fill);
				this_outline = (p.rgb[0] == -1) ? 0 : outline;
				if (this_outline) gmt_setpen (GMT, &p);
				gmt_setfill (GMT, &f, this_outline);
				angle1 = (s->is_var[1]) ? size[s->var[1]] : s->p[1];
				angle2 = (s->is_var[2]) ? size[s->var[2]] : s->p[2];
				if (gmtplot_is_azimuth (symbol, s, 1) || gmtplot_is_azimuth (symbol, s, 2)) {
					angle1 = gmt_azim_to_angle (GMT, lon, lat, 0.1, angle1);
					angle2 = gmt_azim_to_angle (GMT, lon, lat, 0.1, angle2);
				}
				dim[0] *= 0.5;	/* Give diameter */
				dim[1] = angle1;
				dim[2] = angle2;
				dim[5] = p.width * GMT->session.u2u[GMT_PT][GMT_INCH];
				PSL_plotsymbol (PSL, x, y, dim, PSL_MARC);
				break;

			case PSL_WEDGE:
				if (flush) gmtplot_flush_symbol_piece (GMT, PSL, xx, yy, &n, &p, &f, this_outline, &flush);
				gmtplot_get_the_fill (&f, s, current_pen, current_fill);
				gmtplot_get_the_pen (&p, s, current_pen, current_fill);
				this_outline = (p.rgb[0] == -1) ? 0 : outline;
				if (this_outline) gmt_setpen (GMT, &p);
				gmt_setfill (GMT, &f, this_outline);
				angle1 = (s->is_var[1]) ? size[s->var[1]] : s->p[1];
				angle2 = (s->is_var[2]) ? size[s->var[2]] : s->p[2];
				if (gmtplot_is_azimuth (symbol, s, 1) || gmtplot_is_azimuth (symbol, s, 2)) {
					angle1 = gmt_azim_to_angle (GMT, lon, lat, 0.1, angle1);
					angle2 = gmt_azim_to_angle (GMT, lon, lat, 0.1, angle2);
				}
				dim[0] *= 0.5;	/* Give diameter */
				dim[1] = angle1;
				dim[2] = angle2;
				dim[3] = 0;
				dim[7] = 3;
				PSL_plotsymbol (PSL, x, y, dim, PSL_WEDGE);
				break;

			case PSL_VECTOR:
				/* Sets sensible vector head size derived from current stem pen thickness (-W). Head outline is always drawn, fill is optional */
				gmtplot_flush_symbol_piece (GMT, PSL, xx, yy, &n, &p, &f, this_outline, &flush);
				gmtplot_get_the_fill (&f, s, current_pen, current_fill);
				gmtplot_get_the_pen (&p, s, current_pen, current_fill);
				this_outline = (p.rgb[0] == -1) ? 0 : outline;
				if (this_outline) gmt_setpen (GMT, &p);
				gmt_setfill (GMT, &f, this_outline);
				gmt_xy_to_geo (GMT, &lon, &lat, x0, y0);
				angle1 = (s->is_var[0]) ? size[s->var[0]] : s->p[0];
				if (gmtplot_is_azimuth (symbol, s, 0))
					angle1 = gmt_azim_to_angle (GMT, lon, lat, 0.1, angle1);
				PSL_setorigin (PSL, 0.0, 0.0, angle1, PSL_FWD);
				dim[PSL_VEC_XTIP]        = 2.0 * s->p[1] * size[0];	/* Must undo a 0.5 in psxy */
				dim[PSL_VEC_YTIP]        = 0.0;
				dim[PSL_VEC_TAIL_WIDTH]  = current_pen->width / PSL_POINTS_PER_INCH;
				dim[PSL_VEC_HEAD_LENGTH] = 8.0 * current_pen->width / PSL_POINTS_PER_INCH;
				dim[PSL_VEC_HEAD_WIDTH]  = dim[PSL_VEC_HEAD_LENGTH] * 0.5358983848621;	/* 2*tan (15) as default apex is 30 */
				dim[PSL_VEC_HEAD_SHAPE]  = 0.5;
				dim[PSL_VEC_HEAD_PENWIDTH]   = 0.5 * current_pen->width;
				dim[PSL_VEC_STATUS]      = PSL_VEC_END | PSL_VEC_FILL | PSL_VEC_OUTLINE;
				PSL_defpen (PSL, "PSL_vecheadpen", current_pen->width, current_pen->style, current_pen->offset, current_pen->rgb);
				PSL_plotsymbol (PSL, x, y, dim, PSL_VECTOR);
				PSL_setorigin (PSL, 0.0, 0.0, -angle1, PSL_INV);
				break;

			case GMT_SYMBOL_TEXT:
			case GMT_SYMBOL_VARTEXT:
				if (flush) gmtplot_flush_symbol_piece (GMT, PSL, xx, yy, &n, &p, &f, this_outline, &flush);
				gmtplot_get_the_fill (&f, s, current_pen, current_fill);
				gmtplot_get_the_pen (&p, s, current_pen, current_fill);
				this_outline = (p.rgb[0] == -1) ? 0 : outline;
				if (this_outline) gmt_setpen (GMT, &p);
				gmtplot_format_symbol_string (GMT, s, size, user_text);
				if (s->p[0] < 0.0)	/* Fixed point size */
					font.size = -s->p[0];
				else	/* Fractional size */
					font.size = s->p[0] * size[0] * PSL_POINTS_PER_INCH;
				ifs = (int)lrint (font.size * PSL->internal.dpp);	/* Basically psl_ip */
				gmt_setfont (GMT, &s->font);
				if (f.rgb[0] >= 0.0 && this_outline)
					gmt_setfill (GMT, &f, this_outline);
				else if (f.rgb[0] >= 0.0)
					PSL_setcolor (PSL, f.rgb, PSL_IS_FONT);
				else
					PSL_setfill (PSL, GMT->session.no_rgb, this_outline);
				PSL_command (PSL, "%d PSL_symbol_%s_setfont_%d\n", ifs, symbol->name, id++);
				PSL_plottext (PSL, x, y, font.size, user_text, 0.0, s->justify, this_outline);
				break;

			case GMT_SYMBOL_EPS:
				/* Special Encapsulated PostScript-only symbol */
				gmtplot_draw_eps_symbol (GMT, x, y, dim, s);
				break;

			default:
				GMT_Report (GMT->parent, GMT_MSG_ERROR,
								"Unrecognized symbol code (%d = '%c') passed to gmt_draw_custom_symbol\n", action, (char)action);
				return GMT_PARSE_ERROR;
				break;
		}
		if (s->pen) s->pen->width = p_width;	/* Reset to what it was before scaling */
		s = s->next;
	}
	if (flush) gmtplot_flush_symbol_piece (GMT, PSL, xx, yy, &n, &p, &f, this_outline, &flush);
	PSL_command (PSL, "U\n");
	PSL_comment (PSL, "End of symbol %s\n", symbol->name);
	gmt_reset_meminc (GMT);

	/* Restore settings */
	gmt_setpen (GMT, &save_pen);
	if (xx) gmt_M_free (GMT, xx);
	if (yy) gmt_M_free (GMT, yy);

	return (GMT_OK);
}

/* Plotting functions related to contours */

void gmt_add_label_record (struct GMT_CTRL *GMT, struct GMT_DATASET *T, double x, double y, double angle, char *label) {
	/* Add one record to the output file */
	double geo[2];
	uint64_t col, rec = T->table[0]->segment[0]->n_rows;	/* Current record */
	struct GMT_DATASEGMENT_HIDDEN *SH = gmt_get_DS_hidden (T->table[0]->segment[0]);
	/* Convert lon/lat and save x/y */
	gmt_xy_to_geo (GMT, &geo[GMT_X], &geo[GMT_Y], x, y);
	T->table[0]->segment[0]->data[GMT_X][rec] = geo[GMT_X];
	T->table[0]->segment[0]->data[GMT_Y][rec] = geo[GMT_Y];
	/* Also output the label angle */
	T->table[0]->segment[0]->data[GMT_Z][rec] = angle;
	/* Then save the label */
	T->table[0]->segment[0]->text[rec] = strdup (label);
	/* Increment counter */
	T->table[0]->segment[0]->n_rows++;
	if (T->table[0]->segment[0]->n_rows == SH->n_alloc) {
		SH->n_alloc <<= 1;
		T->table[0]->segment[0]->text = gmt_M_memory (GMT, T->table[0]->segment[0]->text, SH->n_alloc, char *);
		for (col = 0; col < T->n_columns; col++)
			T->table[0]->segment[0]->data[col] = gmt_M_memory (GMT, T->table[0]->segment[0]->data[col], SH->n_alloc, double);
	}
}

int gmt_contlabel_save_begin (struct GMT_CTRL *GMT, struct GMT_CONTOUR *G) {
	int kind;
	uint64_t k, seg, dim[4] = {1, 1, GMT_SMALL_CHUNK, 3};	/* Single table with 1 segment and unknown rows, each with 3 columns and trailing text */
	char record[GMT_BUFSIZ] = {""};
	char *xname[2] = {"x", "lon"}, *yname[2] = {"y", "lat"};
	double angle = 0.0;
	struct GMT_CONTOUR_LINE *L = NULL;

	/* Save the lon, lat, angle, text for each annotation to specified file */

	kind = gmt_M_is_geographic (GMT, GMT_IN);
	if ((G->Out = GMT_Create_Data (GMT->parent, GMT_IS_DATASET, GMT_IS_POINT, GMT_WITH_STRINGS, dim, NULL, NULL, 0, 0, NULL)) == NULL) {
		GMT_Report (GMT->parent, GMT_MSG_ERROR, "Unable to create a dataset\n");
		return (GMT_MEMORY_ERROR);	/* Establishes data output */
	}
	/* Write lon, lat, angle, label record */
	snprintf (record, GMT_BUFSIZ, "# %s%s%s%sangle%slabel", xname[kind], GMT->current.setting.io_col_separator, yname[kind],
		GMT->current.setting.io_col_separator, GMT->current.setting.io_col_separator);
	GMT_Set_Comment (GMT->parent, GMT_IS_DATASET, GMT_COMMENT_IS_TEXT | GMT_COMMENT_IS_COMMAND, record, G->Out);
	G->Out->table[0]->segment[0]->n_rows = 0;
	for (seg = 0; seg < G->n_segments; seg++) {
		L = G->segment[seg];	/* Pointer to current segment */
		if (!L->annot || L->n_labels == 0) continue;
		for (k = 0; k < L->n_labels; k++) {
			angle = fmod (2.0 * (L->L[k].angle + 360.0), 360.0) / 2.0;		/* Get text line in 0-180 range */
			gmt_add_label_record (GMT, G->Out, L->L[k].x, L->L[k].y, angle, L->L[k].label);	/* Store text record */
		}
	}
	gmtlib_finalize_dataset (GMT, G->Out);
	return (GMT_NOERROR);
	/* To finish and close the file, call gmt_contlabel_save_end */
}

int gmt_contlabel_save_end (struct GMT_CTRL *GMT, struct GMT_CONTOUR *G) {
	/* Finalize this dataset and write it out */
	gmt_set_dataset_minmax (GMT, G->Out);
	if (GMT_Write_Data (GMT->parent, GMT_IS_DATASET, GMT_IS_FILE, GMT_IS_POINT, GMT_WRITE_SET, NULL, G->label_file, G->Out) != GMT_NOERROR) {
		GMT_Report (GMT->parent, GMT_MSG_ERROR, "Unable to create/write to file %s\n", G->label_file);
		return (GMT_ERROR_ON_FOPEN);	/* Establishes data output */
	}
	GMT_Destroy_Data (GMT->parent, &(G->Out));
	return (GMT_NOERROR);
}

void gmt_textpath_init (struct GMT_CTRL *GMT, struct GMT_PEN *BP, double Brgb[]) {
	PSL_comment (GMT->PSL, "Pen and fill for text boxes (if enabled):\n");
	PSL_defpen (GMT->PSL, "PSL_setboxpen", BP->width, BP->style, BP->offset, BP->rgb);
	PSL_defcolor (GMT->PSL, "PSL_setboxrgb", Brgb);
}

void gmt_contlabel_plot (struct GMT_CTRL *GMT, struct GMT_CONTOUR *G) {
	unsigned int i, mode;
	bool no_labels;
	struct PSL_CTRL *PSL= GMT->PSL;

	if (!G->n_segments) return;	/* Northing to do here */

	/* See if there are labels at all */
	for (i = 0, no_labels = true; i < G->n_segments && no_labels; i++)
		if (G->segment[i]->n_labels) no_labels = false;

	if (!G->delay) PSL_command (GMT->PSL, "V\n");	/* Plotting something, so protect current graphics state */

	if (G->debug) gmtplot_contlabel_debug (GMT, PSL, G);		/* Debugging lines and points */

	if (no_labels) {	/* No labels, just draw lines; no clipping required */
		gmtplot_contlabel_drawlines (GMT, PSL, G, 0);
		PSL_command (GMT->PSL, "U\n");	/* Restore to where we were */
		return;
	}

	PSL_settextmode (PSL, PSL_TXTMODE_MINUS);	/* Replace hyphens with minus signs */
	gmt_setfont (GMT, &G->font_label);

	if (G->must_clip) {		/* Transparent boxes means we must set up plot text, then set up clip paths, then draw lines, then deactivate clipping */
		/* Place PSL variables, plot labels, set up clip paths, draw lines */
		mode = PSL_TXT_INIT | PSL_TXT_SHOW | PSL_TXT_CLIP_ON;
		if (G->draw) mode |= PSL_TXT_DRAW;
		if (!G->delay) mode |= PSL_TXT_CLIP_OFF;	/* Also turn off clip path when done */
		gmtplot_contlabel_plotlabels (GMT, PSL, G, mode);	/* Take the above actions */
	}
	else {	/* Opaque text boxes */
		mode = PSL_TXT_INIT;	/* Plot text */
		if (G->draw) mode |= PSL_TXT_DRAW;	/* Draw lines */
		gmtplot_contlabel_plotlabels (GMT, PSL, G, mode);	/* Place PSL variables and draw lines */
		mode = PSL_TXT_SHOW;				/* Plot text */
		if (G->delay) mode |= PSL_TXT_CLIP_ON;		/* Also turn on clip path after done */
		gmtplot_contlabel_plotlabels (GMT, PSL, G, mode);	/* Plot labels and possibly turn on clipping if delay */
	}
	PSL_command (GMT->PSL, "[] 0 B\n");	/* Ensure no pen textures remain in effect */
	PSL_settextmode (PSL, PSL_TXTMODE_HYPHEN);	/* Back to leave as is */

	if (!G->delay) PSL_command (GMT->PSL, "U\n");	/* Restore to where we were */
}

#if 0        	// We have no current need for this anymore
GMT_LOCAL void gmtplot_wipe_substr(char *str1, char *str2) {
	/* Set the substring str2 of str1 to blanks */
	return;
	char *pch;
	if ((pch = strstr(str1, str2)) != NULL) {
		size_t k;
		for (k = 0; k < strlen(str2); k++)
			pch[k] = ' ';
	}
}
#endif

#ifdef HAVE_GDAL
char *gmt_importproj4 (struct GMT_CTRL *GMT, char *pStr, int *scale_pos) {
	/* Take a PROJ.4 projection string or EPSG code and try to find the equivalent -J syntax
		scale_pos is position on the return string where starts the scale sub-string.
	*/
	unsigned int pos = 0;
	//bool got_lonlat = false;
	char opt_J[GMT_LEN256] = {""}, szProj4[GMT_LEN256] = {""}, prjcode[16] = {""};
	char token[GMT_LEN256] = {""}, scale_c[GMT_LEN32] = {""}, *pch = NULL, *pStrOut = NULL;
	char lon_0[32] = {""}, lat_0[32] = {""}, lat_1[32] = {""}, lat_2[32] = {""}, lat_ts[32] = {""};

	snprintf(szProj4, GMT_LEN256-1, "%s", pStr);

	if ((pch = strchr(szProj4, '/')) != NULL) {		/* Get the scale factor and chop it out of proj string */
		strncpy(scale_c, &pch[1], GMT_LEN32-1);
		pch[0] = '\0';
	}
	else {
		GMT->current.ps.active ? sprintf(scale_c, "14c") : sprintf(scale_c, "1:1");
	}

	if (isdigit(szProj4[1])) {		/* A EPSG code. By looking at 2nd char instead of 1st both +epsg and epsg work */
		int   EPSGID;
		char *pszResult = NULL;
		OGRSpatialReferenceH  hSRS;
		OGRErr eErr = OGRERR_NONE;
#if 0
		bool found = false;
		char buffer [GMT_LEN256] = {""};
		FILE *fp = NULL;
#endif

		if ((pch = strstr(szProj4, "+width=")) != NULL || (pch = strstr(szProj4, "+scale=")) != NULL) {
			snprintf (scale_c, GMT_LEN32-1, "%s", pch);
			pch[0] = '\0';			/* Strip it */
		}

		if (szProj4[0] == '+')		/* Let it both work: -J+epsg or -Jepsg */
			EPSGID = atoi(&szProj4[1]);
		else
			EPSGID = atoi(szProj4);

		/* Rely on GDAL to tell us the proj4 string of this EPSG code */
		hSRS = OSRNewSpatialReference(NULL);
		if ((eErr = OSRImportFromEPSG(hSRS, EPSGID)) != OGRERR_NONE) {
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "Did not get the SRS from input EPSG  %d\n", EPSGID);
			return (pStrOut);
		}
		if ((eErr = OSRExportToProj4(hSRS, &pszResult)) != OGRERR_NONE) {
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "Failed to convert the SRS to proj4 syntax\n");
			return (pStrOut);
		}
		snprintf(szProj4, GMT_LEN256-1, "%s", pszResult);
		if (scale_c[0] != '\0') strcat (szProj4, scale_c);		/* Add the width/scale found above */
		CPLFree(pszResult);
		OSRDestroySpatialReference(hSRS);

#if 0
		gmt_getsharepath (GMT, "", "epsg", ".txt", buffer, R_OK);
		if ((fp = fopen (buffer, "r")) == NULL) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Cannot open %s file\n", buffer);
			return (pStrOut);
		}
		while (fgets (buffer, GMT_LEN256, fp)) {
			if (buffer[0] == '#') continue;
			if ((pch = strstr(buffer, "+proj")) != NULL) {
				pch[0] = '\0';			/* Break the line before the +proj=... */
				if (EPSGID == atoi(buffer)) {
					pch[0] = '+';		/* Undo previous break */
					snprintf(szProj4, GMT_LEN256-1, "%s", pch);
					found = true;
					break;
				}
			}
		}
		fclose (fp);
		if (!found) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Could not find the EPGS code %d in database\n", EPSGID);
			return (pStrOut);
		}
#endif
	}

	if (gmt_strtok(szProj4, " \t+", &pos, token)) {
		snprintf(prjcode, 16, "%s",(token[0] == '+' ? &token[6] : &token[5]));	/* PROJ4 projection code. */
		//gmtplot_wipe_substr(szProj4, token);	/* Consumed, clear it from list */
	}

	if (!strcmp(prjcode, "longlat") || !strcmp(prjcode, "latlong")) {
		strcat (opt_J, "X");
		GMT->current.proj.projection_GMT = GMT_LINEAR;
		//got_lonlat = true;			/* At the end we need to append a 'd' to the scale */
	}
	/* Cylindrical projections */
	else if (!strcmp(prjcode, "cea") || !strcmp(prjcode, "eqc") || !strcmp(prjcode, "tmerc") ||
			  !strcmp(prjcode, "mill") || !strcmp(prjcode, "merc") || !strcmp(prjcode, "cass")) {
		if (!strcmp(prjcode, "tmerc"))     strcat (opt_J, "T");
		else if (!strcmp(prjcode, "cea"))  strcat (opt_J, "Y");
		else if (!strcmp(prjcode, "eqc"))  strcat (opt_J, "Q");
		else if (!strcmp(prjcode, "merc")) strcat (opt_J, "M");
		else if (!strcmp(prjcode, "mill")) strcat (opt_J, "J");
		else strcat (opt_J, "C");
		while (gmt_strtok (szProj4, " \t+", &pos, token)) {
			if ((pch = strstr(token, "lon_0=")) != NULL) {
				strncat(lon_0, &token[6], 31);
				//gmtplot_wipe_substr(szProj4, token);
			}
			else if ((pch = strstr(token, "lat_0=")) != NULL) {
				strncat(lat_0, &token[6], 31);
				//gmtplot_wipe_substr(szProj4, token);
			}
		}
		if (!strcmp(prjcode, "cea") || !strcmp(prjcode, "eqc") || !strcmp(prjcode, "cass") ||
			 !strcmp(prjcode, "tmerc") || !strcmp(prjcode, "merc")) {
			if (!lon_0[0]) strcat(lon_0, "0");
			if (!lat_0[0]) strcat(lat_0, "0");
			strcat(opt_J, lon_0);	strcat (opt_J, "/");
			strcat(opt_J, lat_0);	strcat (opt_J, "/");
		}
		else {		// "mill"
			if (lon_0[0]) strcat(opt_J, lon_0),	strcat (opt_J, "/");
		}
	}
	else if (!strcmp(prjcode, "omerc")) {
		char lon_1[32] = {""}, lon_2[32] = {""}, lonc[32] = {""}, alpha[32] = {""};
		strcat (opt_J, "OC");
		while (gmt_strtok (szProj4, " \t+", &pos, token)) {
			if ((pch = strstr(token, "lon_1=")) != NULL)
				strncat(lon_1, &token[6], 31);
			else if ((pch = strstr(token, "lat_1=")) != NULL)
				strncat(lat_1, &token[6], 31);
			else if ((pch = strstr(token, "lon_2=")) != NULL)
				strncat(lon_2, &token[6], 31);
			else if ((pch = strstr(token, "lat_2=")) != NULL)
				strncat(lat_2, &token[6], 31);
			else if ((pch = strstr(token, "lat_0=")) != NULL)
				strncat(lat_0, &token[6], 31);
			else if ((pch = strstr(token, "lonc=")) != NULL)
				strncat(lonc, &token[5], 31);
			else if ((pch = strstr(token, "alpha=")) != NULL)
				strncat(alpha, &token[6], 31);
		}
		if (lon_2[0] && lat_2[0]) {
			if (!lon_1[0]) strcat(lon_1, "0");
			if (!lat_1[0]) strcat(lat_1, "0");
			strcat(opt_J, lon_1);	strcat (opt_J, "/");	strcat(opt_J, lat_1);	strcat (opt_J, "/");
			strcat(opt_J, lon_2);	strcat (opt_J, "/");	strcat(opt_J, lat_2);	strcat (opt_J, "/");
		}
		else if (lonc[0] && alpha[0]) {
			if (!lat_0[0]) strcat(lat_0, "0");
			strcat(opt_J, lonc);	strcat (opt_J, "/");	strcat(opt_J, lat_0);	strcat (opt_J, "/");
			strcat(opt_J, alpha);	strcat (opt_J, "/");
			opt_J[1] = 'A';
		}
		else {
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "In projection %s proj parameters\n", prjcode);
			return (pStrOut);
		}
	}
	else if (!strcmp(prjcode, "utm")) {
		int zone = 100;
		strcat (opt_J, "U");
		while (gmt_strtok (szProj4, " \t+", &pos, token)) {
			if ((pch = strstr(token, "zone=")) != NULL) {
				zone = atoi(&token[5]);
				//gmtplot_wipe_substr(szProj4, token);
			}
			else if ((pch = strstr(token, "south=")) != NULL) {
				zone *= -1;
				//gmtplot_wipe_substr(szProj4, token);
			}
			else if ((pch = strstr(token, "lon_0=")) != NULL) {
				double x;
				gmt_scanf (GMT, &token[6], GMT_IS_GEO, &x);
				if (x > 180) x -= 360;
				zone = (int)(x + 180) / 6 + 1;
				//gmtplot_wipe_substr(szProj4, token);
			}
		}
		if (zone == 100) {
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "UTM proj selected but no info about utm zone.\n");
			return (pStrOut);
		}
		else if (zone > 64) {
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "UTM proj The lon_0 argument was not correctly parsed.\n");
		}
		else {
			char t[4];
			sprintf(t, "%d/", zone);
			strcat(opt_J, t);
		}
	}

	/* Conic projections */
	else if (!strcmp(prjcode, "aea") || !strcmp(prjcode, "eqdc") || !strcmp(prjcode, "lcc") || !strcmp(prjcode, "poly")) {
		if (!strcmp(prjcode, "aea")) strcat (opt_J, "B");
		else if (!strcmp(prjcode, "eqdc")) strcat (opt_J, "D");
		else if (!strcmp(prjcode, "lcc")) strcat (opt_J, "L");
		else strcat (opt_J, "Poly/");
		while (gmt_strtok (szProj4, " \t+", &pos, token)) {
			if ((pch = strstr(token, "lon_0=")) != NULL) {
				strncat(lon_0, &token[6], 31);
				//gmtplot_wipe_substr(szProj4, token);
			}
			else if ((pch = strstr(token, "lat_0=")) != NULL) {
				strncat(lat_0, &token[6], 31);
				//gmtplot_wipe_substr(szProj4, token);
			}
			else if ((pch = strstr(token, "lat_1=")) != NULL) {
				strncat(lat_1, &token[6], 31);
				//gmtplot_wipe_substr(szProj4, token);
			}
			else if ((pch = strstr(token, "lat_2=")) != NULL) {
				strncat(lat_2, &token[6], 31);
				//gmtplot_wipe_substr(szProj4, token);
			}
		}
		/* Check if user errors */
		if (!lon_0[0]) strcat(lon_0, "0");
		if (!lat_0[0]) strcat(lat_0, "0");
		if (strcmp(prjcode, "poly")) {			/* i.e. if NOT Poly */
			if (!lat_1[0] || !lat_2[0]) {
				GMT_Report (GMT->parent, GMT_MSG_WARNING, "Projection %s needs the lat_1 & lat_2 proj parameters\n", prjcode);
				return (pStrOut);
			}
			strcat(opt_J, lon_0);	strcat (opt_J, "/");	strcat(opt_J, lat_0);	strcat (opt_J, "/");
			strcat(opt_J, lat_1);	strcat (opt_J, "/");	strcat(opt_J, lat_2);	strcat (opt_J, "/");
		}
		else {
			strcat(opt_J, lon_0);	strcat (opt_J, "/");	strcat(opt_J, lat_0);	strcat (opt_J, "/");
		}
	}

	/* Azimuthal projections */
	else if (!strcmp(prjcode, "stere") || !strcmp(prjcode, "laea") || !strcmp(prjcode, "aeqd") || !strcmp(prjcode, "gnom") ||
				!strcmp(prjcode, "sterea") || !strcmp(prjcode, "stere")) {
		if (!strcmp(prjcode, "stere")) strcat (opt_J, "S");
		else if (!strcmp(prjcode, "laea")) strcat (opt_J, "A");
		else if (!strcmp(prjcode, "aeqd")) strcat (opt_J, "E");
		else if (!strcmp(prjcode, "gnom")) strcat (opt_J, "F");
		else strcat (opt_J, "S");
		while (gmt_strtok (szProj4, " \t+", &pos, token)) {
			if ((pch = strstr(token, "lon_0=")) != NULL) {
				strncat(lon_0, &token[6], 31);
				//gmtplot_wipe_substr(szProj4, token);
			}
			else if ((pch = strstr(token, "lat_0=")) != NULL) {
				strncat(lat_0, &token[6], 31);
				//gmtplot_wipe_substr(szProj4, token);
			}
			else if ((pch = strstr(token, "lat_ts=")) != NULL) {
				strncat(lat_ts, &token[7], 31);
				//gmtplot_wipe_substr(szProj4, token);
			}
		}
		//if (!strcmp(prjcode, "ups")) lon_0[0] = lat_0[0] = '\0';
		if (!strcmp(prjcode, "stere")) {
			if (!lat_0[0]) strcat(lat_0, "90");		/* ptoj4 says lat_0 = 90 but what if in southerm hemisphere? */
		}
		if (!lon_0[0]) strcat(lon_0, "0");
		if (!lat_0[0]) strcat(lat_0, "0");
		strcat(opt_J, lon_0);	strcat (opt_J, "/");	strcat(opt_J, lat_0);	strcat (opt_J, "/");
		if (!strcmp(prjcode, "stere")) strcat(opt_J, "75/");
	}

	/* Misc projections */
	else if (!strcmp(prjcode, "moll") || !strcmp(prjcode, "sinu") || !strcmp(prjcode, "vandg") ||
				!strcmp(prjcode, "robin") || !strcmp(prjcode, "wintri") || !strcmp(prjcode, "hammer") ||
				!strcmp(prjcode, "eck4") || !strcmp(prjcode, "eck6")) {
		if (!strcmp(prjcode, "moll")) strcat (opt_J, "W");
		else if (!strcmp(prjcode, "sinu"))   strcat (opt_J, "I");
		else if (!strcmp(prjcode, "vandg"))  strcat (opt_J, "V");
		else if (!strcmp(prjcode, "robin"))  strcat (opt_J, "N");
		else if (!strcmp(prjcode, "wintri")) strcat (opt_J, "R");
		else if (!strcmp(prjcode, "hammer")) strcat (opt_J, "H");
		else if (!strcmp(prjcode, "eck4"))   strcat (opt_J, "Kf");
		else strcat (opt_J, "Ks");
		while (gmt_strtok (szProj4, " \t+", &pos, token)) {
			if ((pch = strstr(token, "lon_0=")) != NULL) {
				strncat(opt_J, &token[6], GMT_LEN256-1);	strcat (opt_J, "/");
				//gmtplot_wipe_substr(szProj4, token);
			}
			else if ((pch = strstr(token, "lun_0=")) != NULL) {
				strncat(opt_J, &token[6], GMT_LEN256-1);	strcat (opt_J, "/");
			}
		}
		//if (opt_J[strlen(opt_J)-1] != '/')		/* Not strictly needed by GMT but needed in gmt_parse_common_options() */
			//strcat(opt_J, "0/");
	}

	else	/* We don't return yet because we may have a +width/+scale to parse */
		GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "This projection '%s' is not natively supported in GMT\n", prjcode);

#if 0
	bool got_a = false, got_b = false;
	char ename[GMT_LEN16] = {""};
	/* Work on common flags */
	if ((pch = strstr(szProj4, "+a=")) != NULL) {	/* Check for major axis +a=xxxx */
		pos = 0;	gmt_strtok (pch, " \t+", &pos, token);
		GMT->current.setting.ref_ellipsoid[GMT->current.setting.proj_ellipsoid].eq_radius = atof(&token[2]);
		got_a = true;
		gmtplot_wipe_substr(szProj4, token);	/* Set the token part in szProj4 to blanks */
	}
	if ((pch = strstr(szProj4, "+b=")) != NULL) {	/* Check for minor axis +b=xxxx */
		double f;
		if (!got_a) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "+b= in proj4 string without the companion +a\n");
			return (pStrOut);
		}
		pos = 0;	gmt_strtok (pch, " \t+", &pos, token);
		/* f =  1 - b / a; */
		f =  1 - atof(&token[2]) / GMT->current.setting.ref_ellipsoid[GMT->current.setting.proj_ellipsoid].eq_radius;
		GMT->current.setting.ref_ellipsoid[GMT->current.setting.proj_ellipsoid].flattening = f;
		got_b = true;
		gmtplot_wipe_substr(szProj4, token);
	}
	if ((pch = strstr(szProj4, "+k_0=")) != NULL) {		/* Check for scale factor */
		pos = 0;	gmt_strtok (pch, " \t+", &pos, token);
		GMT->current.setting.proj_scale_factor = atof (&token[4]);
		gmtplot_wipe_substr(szProj4, token);
	}
	if ((pch = strstr(szProj4, "+k=")) != NULL) {		/* Check for scale factor */
		pos = 0;	gmt_strtok (pch, " \t+", &pos, token);
		GMT->current.setting.proj_scale_factor = atof (&token[2]);
		gmtplot_wipe_substr(szProj4, token);
	}

	if ((pch = strstr(szProj4, "+ellps=")) != NULL) {	/* Check for ellipsoids */
		pos = 0;	gmt_strtok (pch, " \t+", &pos, token);
		gmtplot_ellipsoid_name_convert2(&token[6], ename);
		if (ename[0] != '\0')
			;//sprintf(GMT->current.setting.ref_ellipsoid[GMT->current.setting.proj_ellipsoid].name, "%s", ename);
		else {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "could not translate the ellipsoid name %s\n", &token[6]);
			return (pStrOut);
		}
		gmtplot_wipe_substr(szProj4, token);
	}

	if ((pch = strstr(szProj4, "+datum=")) != NULL) {	/* Check for datum */
		char t[128] = {""};
		pos = 0;	gmt_strtok (pch, " \t+", &pos, token);

		if (!strcmp(&token[6], "WGS84"))
			gmtplot_ellipsoid_name_convert2("WGS84", ename);
		else if (!strcmp(&token[6], "GGRS87")) {
			gmtplot_ellipsoid_name_convert2("GRS80", ename);
			sprintf(t, " +towgs84=-199.87,74.79,246.62");
		}
		else if (!strcmp(&token[6], "NAD83"))
			gmtplot_ellipsoid_name_convert2("WGS84", ename);
		else if (!strcmp(&token[6], "potsdam")) {
			gmtplot_ellipsoid_name_convert2("bessel", ename);
			sprintf(t, " +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7");
		}
		else {
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "this datum %s is not yet ported to GMT\n", &token[6]);
			return (pStrOut);
		}

		if (ename[0] != '\0')
			sprintf(GMT->current.setting.ref_ellipsoid[GMT->current.setting.proj_ellipsoid].name, "%s", ename);
		else {
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "could not translate the ellipsoid name %s\n", &token[6]);
			return (pStrOut);
		}
		strcat(szProj4, t);			/* Append to the proj4 string so that the +towgs84 case below will handle this */
		gmtplot_wipe_substr(szProj4, token);
	}

	if ((pch = strstr(szProj4, "+towgs84=")) != NULL) {
		char *txt, t[128] = {""};
		struct GMT_DATUM to, from;

		memset(&to, 0, sizeof(struct GMT_DATUM));
		memset(&from, 0, sizeof(struct GMT_DATUM));

		if (ename[0] == '\0' && (!got_a || !got_b)) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR,
							" Cannot convert to WGS84 if you don't tell me the ellipsoid of origin\n"
						 "(miss +ellips=xxxx OR +a=xxxx +b=xxxx)\n");
			return (pStrOut);
		}
		pos = 0;	gmt_strtok (pch, " \t+", &pos, token);
		txt = strdup (&token[8]);
		if (ename[0] != '\0')                   /* We know the ellipsoid name */
			sprintf(t, "%s:%s", ename, txt);    /* Create an ellip:dx,dy,dz string */
		else                                    /* Need to pass a,1/f */
			sprintf(t, "%12g,%.10f:%s",
					  GMT->current.setting.ref_ellipsoid[GMT->current.setting.proj_ellipsoid].eq_radius,
					1 / GMT->current.setting.ref_ellipsoid[GMT->current.setting.proj_ellipsoid].flattening,
					txt);

		// Need to create and fill 2 GMT_DATUM structs (from & to). See gmt_set_datum() & gmt_datum_init()
		gmt_set_datum (GMT, "-", &from);		/* 'to' is WG84 */
		gmt_set_datum (GMT, t, &to);
		gmt_datum_init(GMT, &from, &to, false);

		for (k = 0; k < 7; k++)
			GMT->current.proj.datum.bursa[k] = GMT->current.proj.datum.to.xyz[k];

		free(txt);
		gmtplot_wipe_substr(szProj4, token);
	}

	if ((pch = strstr(szProj4, "+x_0=")) != NULL) {
		pos = 0;	gmt_strtok (pch, " \t+", &pos, token);
		GMT->current.proj.proj4_x0 = atof(&token[4]);
		gmtplot_wipe_substr(szProj4, token);
	}
	if ((pch = strstr(szProj4, "+y_0=")) != NULL) {
		pos = 0;	gmt_strtok (pch, " \t+", &pos, token);
		GMT->current.proj.proj4_y0 = atof(&token[4]);
		gmtplot_wipe_substr(szProj4, token);
	}

	if ((pch = strstr(szProj4, "+units=")) != NULL) {
		pos = 0;	gmt_strtok (pch, " \t+", &pos, token);
		if (token[6] != 'm') {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "using units other than meters is not yet implemented\n");
		}
		gmtplot_wipe_substr(szProj4, token);
	}

	/* Remove known unused fields in proj4 string */
	if ((pch = strstr(szProj4, "+no_defs")) != NULL) {
		pos = 0;	gmt_strtok (pch, " \t+", &pos, token);
		gmtplot_wipe_substr(szProj4, token);
	}
	if ((pch = strstr(szProj4, "+wktext")) != NULL) {	/* +wktext is used by GDAL to capture the proj.4 string of unsupported projs */
		pos = 0;	gmt_strtok (pch, " \t+", &pos, token);
		gmtplot_wipe_substr(szProj4, token);
	}
#endif

	/* Override the /1:xxx scale set at the beginning of this function if a +scale=scale is found */
	if ((pch = strstr(szProj4, "+scale=")) != NULL) {
		pos = 0;
		if (gmt_strtok (pch, " \t+", &pos, token)) sprintf(scale_c, "%s", &token[6]);
		//gmtplot_wipe_substr(szProj4, token);
		GMT->current.proj.gave_map_width = 0;
		GMT->current.proj.proj4_scl = gmt_M_to_inch (GMT, scale_c);
	}
	/* If a +width=xx is given, append it a 'W' so that we identify this in gmt_parse_common_options() and act */
	if ((pch = strstr(szProj4, "+width=")) != NULL) {
		pos = 0;
		if (gmt_strtok (pch, " \t+", &pos, token)) sprintf(scale_c, "%sW", &token[6]);
		//gmtplot_wipe_substr(szProj4, token);
		GMT->current.proj.gave_map_width = 1;
		GMT->current.proj.proj4_scl = gmt_M_to_inch (GMT, &token[6]);
	}

	*scale_pos = (int)strlen(opt_J);		/* The position at which the scale string will be appended */

	if (!opt_J[0])		/* No corresponding GMT proj found but we need the scale separated with a slash */
		sprintf(opt_J, "/%s", scale_c);
	else
		strcat(opt_J, scale_c);	/* Append the scale */

	if (strchr(scale_c, ':')) {	/* If we have a scale in the 1:xxxx form use lower case codes */
		opt_J[0] = (char)tolower(opt_J[0]);
		if (!strcmp(prjcode, "omerc")) opt_J[1] = (char)tolower(opt_J[1]);
	}

#if 0
	/* This part is now not particularly useful but let it in case we may need it in future */
	k = 0;
	while (szProj4[k] && (szProj4[k] == ' ' || szProj4[k] == '+')) k++;
	if (k < strlen(szProj4))
		GMT_Report (GMT->parent, GMT_MSG_WARNING, "Original proj4 string was not all consumend. Remaining options:\n\t%s\n", szProj4);
#endif

	GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "Converted to -J syntax = -J%s\n", opt_J);
	GMT->current.proj.is_proj4 = true;		/* Used so far in map|grdproject to set local -C */

	pStrOut = strdup(opt_J);
	return pStrOut;
}
#endif

char *gmt_export2proj4 (struct GMT_CTRL *GMT) {
	char *pStrOut = NULL;
	char szProj4[GMT_LEN512], proj4_ename[GMT_LEN16];
	double scale_factor, false_easting = 0.0, false_northing = 0.0, a, b, f;

	/* If we already have a proj.4 string, obviously use it and return right away */
	if (GMT->common.J.proj4string[0] == '+') {
		pStrOut = strdup (GMT->common.J.proj4string);
		return pStrOut;
	}

	scale_factor = GMT->current.setting.proj_scale_factor;
	if (scale_factor < 0) scale_factor = 1;
	szProj4[0] = 0;

	switch (GMT->current.proj.projection_GMT) {
	/* Cylindrical projections */
	case GMT_UTM:
		snprintf (szProj4, GMT_LEN512, "+proj=utm +zone=%d +units=m", (int)GMT->current.proj.pars[0]);
		if (GMT->current.proj.utm_hemisphere < 0) strcat (szProj4, " +south");
		break;
	case GMT_MERCATOR:
		snprintf (szProj4, GMT_LEN512, "+proj=merc +lon_0=%.16g +k=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[0] >= -360 ? GMT->current.proj.pars[0] : 0, scale_factor, false_easting, false_northing);
		break;
	case GMT_CYL_EQ:
		snprintf (szProj4, GMT_LEN512, "+proj=cea +lon_0=%.16g +lat_ts=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[1], GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_CYL_EQDIST:
		snprintf (szProj4, GMT_LEN512, "+proj=eqc +lat_ts=%.16g +lat_0=%.16g +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[1], 0.0, GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_CYL_STEREO:
		break;
	case GMT_MILLER:
		snprintf (szProj4, GMT_LEN512, "+proj=mill +lat_0=%.16g +lon_0=%.16g +x_0=%.16g +y_0=%.16g +R_A +units=m",
			GMT->current.proj.pars[1], GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_TM:
		snprintf (szProj4, GMT_LEN512, "+proj=tmerc +lat_0=%.16g +lon_0=%.16g +k=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[1], GMT->current.proj.pars[0], scale_factor, false_easting, false_northing);
		break;
	case GMT_CASSINI:
		snprintf (szProj4, GMT_LEN512, "+proj=cass +lat_0=%.16g +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[1], GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_OBLIQUE_MERC:
		snprintf (szProj4, GMT_LEN512, "+unavailable");
		/*snprintf (szProj4, GMT_LEN512, "+proj=omerc +lat_0=%.16g +lonc=%.16g +alpha=%.16g +k=%.16g +x_0=%.16g +y_0=%.16g",
		0.0,0.0,0.0,0.0,0.0,0.0 );*/
		break;
	case GMT_OBLIQUE_MERC_POLE:
		sprintf (szProj4, "+unavailable");
		break;

	/* Conic projections */
	case GMT_ALBERS:
		snprintf (szProj4, GMT_LEN512, "+proj=aea +lat_1=%.16g +lat_2=%.16g +lat_0=%.16g +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[2], GMT->current.proj.pars[3], GMT->current.proj.pars[1], GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_ECONIC:
		snprintf (szProj4, GMT_LEN512, "+proj=eqdc +lat_1=%.16g +lat_2=%.16g +lat_0=%.16g +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[2], GMT->current.proj.pars[3], GMT->current.proj.pars[1], GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_LAMBERT:
		snprintf (szProj4, GMT_LEN512, "+proj=lcc +lat_1=%.16g +lat_2=%.16g +lat_0=%.16g +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[2], GMT->current.proj.pars[3], GMT->current.proj.pars[1], GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_POLYCONIC:
		snprintf (szProj4, GMT_LEN512, "+proj=poly +lat_0=%.16g +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[1], GMT->current.proj.pars[0], false_easting, false_northing);
		break;

	/* Azimuthal projections */
	case GMT_STEREO:
		snprintf (szProj4, GMT_LEN512, "+proj=stere +lat_0=%.16g +lon_0=%.16g +k=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[1], GMT->current.proj.pars[0], scale_factor, false_easting, false_northing);
		break;
	case GMT_LAMB_AZ_EQ:
		snprintf (szProj4, GMT_LEN512, "+proj=laea +lat_0=%.16g +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[1], GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_ORTHO:
		sprintf (szProj4, "+unavailable");
		break;
	case GMT_AZ_EQDIST:
		snprintf (szProj4, GMT_LEN512, "+proj=aeqd +lat_0=%.16g +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[1], GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_GNOMONIC:
		snprintf (szProj4, GMT_LEN512, "+proj=gnom +lat_0=%.16g +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[1], GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_GENPER:
		sprintf (szProj4, "+unavailable");
		break;
	case GMT_POLAR:
		sprintf (szProj4, "+unavailable");
		break;

	/* Misc projections */
	case GMT_MOLLWEIDE:
		snprintf (szProj4, GMT_LEN512, "+proj=moll +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_HAMMER:
		sprintf (szProj4, "+unavailable");
		break;
	case GMT_SINUSOIDAL:
		snprintf (szProj4, GMT_LEN512, "+proj=sinu +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_VANGRINTEN:
		snprintf (szProj4, GMT_LEN512, "+proj=vandg +lon_0=%.16g +x_0=%.16g +y_0=%.16g +R_A +units=m",
			GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_ROBINSON:
		snprintf (szProj4, GMT_LEN512, "+proj=robin +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_ECKERT4:
		snprintf (szProj4, GMT_LEN512, "+proj=eck4 +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_ECKERT6:
		snprintf (szProj4, GMT_LEN512, "+proj=eck6 +lon_0=%.16g +x_0=%.16g +y_0=%.16g +units=m",
			GMT->current.proj.pars[0], false_easting, false_northing);
		break;
	case GMT_WINKEL:
		 printf (szProj4, "+unavailable");
		break;
	default:
		if (gmt_M_is_geographic (GMT, GMT_IN))
			sprintf (szProj4, "+proj=latlong");
		else
			sprintf (szProj4, "+xy");	/* Probably useless as an info, but put there something */
	}

	if (strcmp(szProj4, "+xy")) {
		size_t len = strlen (szProj4);
		a = GMT->current.setting.ref_ellipsoid[GMT->current.setting.proj_ellipsoid].eq_radius;
		f = GMT->current.setting.ref_ellipsoid[GMT->current.setting.proj_ellipsoid].flattening;
		b = a * (1 - f);
		snprintf (szProj4+len, GMT_LEN512-len, " +a=%.3f +b=%.3f", a, b);
		len = strlen (szProj4);
		if (fabs(a - b) > 1) {		/* WGS84 is not spherical */
			gmtlib_ellipsoid_name_convert (GMT->current.setting.ref_ellipsoid[GMT->current.setting.proj_ellipsoid].name, proj4_ename);
			snprintf(szProj4+len, GMT_LEN512-len, " +ellps=%s", proj4_ename);
			len = strlen (szProj4);
			if (!strcmp(proj4_ename, "WGS84"))
				snprintf(szProj4+strlen(szProj4), GMT_LEN512-len, " +datum=WGS84");
		}
		len = strlen (szProj4);
		snprintf(szProj4+len, GMT_LEN512-len, " +units=m +no_defs");
	}

	pStrOut = strdup(szProj4);
	return (pStrOut);
}

GMT_LOCAL void gmtplot_just_f_xy (int justify, double *fx, double *fy) {
	/* Factors to multiply feature width/height and add to x,y that updates its mid-point location */
	*fy = -0.5 * ((justify / 4) - 1);
	*fx = -0.5 * ((justify % 4) - 2);
}

/* Here are 6 (a-f) different progress indicator plotting functions.  See the movie documentation for how they look */

GMT_LOCAL void gmtplot_prog_indicator_A (struct GMT_CTRL *GMT, double x, double y, double t, double w, int justify, char *F1, char *F2) {
	/* Place default progress indicator pie */
	double dim[PSL_MAX_DIMS], fx, fy;
	struct GMT_FILL fill;
	gmt_M_memset (dim, PSL_MAX_DIMS, double);
	gmtplot_just_f_xy (justify, &fx, &fy);
	x += fx * w;	y += fy * w;	/* Move to center of circle */
	dim[0] = w;
	if (t < 1.0) {	/* Need the background full circle since partly visible */
		gmt_getfill (GMT, F2, &fill);	/* Want to paint background full circle*/
		PSL_setfill (GMT->PSL, fill.rgb, 0);	/* Full circle color */
		PSL_plotsymbol (GMT->PSL, x, y, dim, PSL_CIRCLE);	/* Plot full circle */
	}
	gmt_getfill (GMT, F1, &fill);	/* Want to paint inside of tag box */
	PSL_setfill (GMT->PSL, fill.rgb, 0);	/* Wedge color */
	if (doubleAlmostEqual (t, 1.0)) {
		PSL_plotsymbol (GMT->PSL, x, y, dim, PSL_CIRCLE);	/* Plot full circle */
	}
	else {
		dim[0] = 0.5 * w;	/* Apparently we take radius for wedge */
		dim[1] = 90.0 - 360 * t;	/* Go clockwise. Convert t to 0-360 degrees */
		dim[2] = 90.0;	/* Start is 12 'oclock */
		dim[7] = 1;	/* Lay down filled wedge */
		PSL_command (GMT->PSL, "/PSL_spiderpen {} def\n");	/* So wedge won't fuss about not being set (like in psxy) */
		PSL_plotsymbol (GMT->PSL, x, y, dim, PSL_WEDGE);	/* Plot wedge */
	}
}

GMT_LOCAL void gmtplot_prog_indicator_B (struct GMT_CTRL *GMT, double x, double y, double t, double w, int justify, char *P1, char *P2, char *label, char kind, double fsize) {
	/* Place growing ring */
	double fx, fy, dr2;
	struct GMT_PEN pen;

	gmtplot_just_f_xy (justify, &fx, &fy);
	gmt_M_memset (&pen, 1, struct GMT_PEN);
	gmt_getpen (GMT, P2, &pen);	/* Want to draw full circle */
	dr2 = pen.width / PSL_POINTS_PER_INCH;	/* Half pen width */
	x += fx * (w+dr2);	y += fy * (w+dr2);	/* Move to center of circle */
	PSL_command (GMT->PSL, "FQ %% Force turn off any prior fill\n");
	if (kind == 'B') {
		if (fsize == 0.0) fsize = 0.3 * w * PSL_POINTS_PER_INCH;
		PSL_plottext (GMT->PSL, x, y, fsize, label, 0.0, PSL_MC, 0);
	}
	if (t < 1.0) {	/* Need to plot the background full circle since partly visible */
		gmt_setpen (GMT, &pen);	/* Full circle pen */
		PSL_setfill (GMT->PSL, GMT->session.no_rgb, 1);
		PSL_plotsymbol (GMT->PSL, x, y, &w, PSL_CIRCLE);	/* Plot full circle */
	}
	gmt_M_memset (&pen, 1, struct GMT_PEN);
	gmt_getpen (GMT, P1, &pen);	/* Always draw foreground circle */
	gmt_setpen (GMT, &pen);	/* Full circle pen */
	PSL_setfill (GMT->PSL, GMT->session.no_rgb, 1);
	if (doubleAlmostEqual (t, 1.0))
		PSL_plotsymbol (GMT->PSL, x, y, &w, PSL_CIRCLE);	/* Plot full circle */
	else
		PSL_plotarc (GMT->PSL, x, y, 0.5*w, 90.0 - 360 * t, 90.0, PSL_MOVE | PSL_STROKE);	/* Draw the arc */
}

GMT_LOCAL void gmtplot_prog_indicator_C (struct GMT_CTRL *GMT, double x, double y, double t, double w, int justify, char *P1, char *P2, char *label, char kind, double fsize) {
	/* Place growing math arrow */
	double fx, fy, dim[PSL_MAX_DIMS];
	struct GMT_PEN pen;

	gmt_M_memset (dim, PSL_MAX_DIMS, double);
	gmtplot_just_f_xy (justify, &fx, &fy);
	gmt_M_memset (&pen, 1, struct GMT_PEN);
	gmt_getpen (GMT, P2, &pen);	/* Want to draw full circle */
	x += fx * 1.15 * w;	y += fy * 1.15 * w;	/* Move to center of circle, add 15% to get more space so arrow head won't clip at canvas edge */
	PSL_command (GMT->PSL, "FQ %% Force turn off any prior fill\n");
	if (kind == 'C') {	/* Place center label */
		if (fsize == 0.0) fsize = 0.3 * w * PSL_POINTS_PER_INCH;
		PSL_plottext (GMT->PSL, x, y, fsize, label, 0.0, PSL_MC, 0);
	}
	if (t < 1.0) {	/* Need to plot the background full circle */
		gmt_setpen (GMT, &pen);	/* Full circle pen */
		PSL_setfill (GMT->PSL, GMT->session.no_rgb, 1);
		PSL_plotsymbol (GMT->PSL, x, y, &w, PSL_CIRCLE);	/* Plot full circle */
	}
	gmt_M_memset (&pen, 1, struct GMT_PEN);
	gmt_getpen (GMT, P1, &pen);	/* Always draw foreground circle */
	gmt_setpen (GMT, &pen);	/* Front circle pen */
	PSL_setfill (GMT->PSL, pen.rgb, 1);
	dim[PSL_MATHARC_RADIUS]          = 0.5 * w;	/* Apparently we take radius for wedge */
	dim[PSL_MATHARC_ANGLE_BEGIN]     = 90.0 - 360 * t;	/* Go clockwise. Convert t to 0-360 degrees */
	dim[PSL_MATHARC_ANGLE_END]       = 90.0;	/* Start is 12 'oclock */
	dim[PSL_MATHARC_HEAD_LENGTH]     = 0.2 * w;
	dim[PSL_MATHARC_HEAD_WIDTH]      = 0.2 * w;
	dim[PSL_MATHARC_ARC_PENWIDTH]    = pen.width / PSL_POINTS_PER_INCH;
	dim[PSL_MATHARC_HEAD_SHAPE]      = 0.75;	/* Fixed shape setting */
	dim[PSL_MATHARC_STATUS]          = (double)(PSL_VEC_BEGIN | PSL_VEC_FILL);
	dim[PSL_MATHARC_HEAD_TYPE_BEGIN] = (double)PSL_VEC_ARROW;
	dim[PSL_MATHARC_HEAD_TYPE_END]   = (double)PSL_VEC_ARROW;
	dim[PSL_MATHARC_HEAD_PENWIDTH]   = 0.5 * pen.width;
	PSL_command (GMT->PSL, "/PSL_vecheadpen {} def\n");	/* So wedge won't fuss about not being set (like in psxy) */
	PSL_plotsymbol (GMT->PSL, x, y, dim, PSL_MARC);
}

GMT_LOCAL void gmtplot_prog_indicator_D (struct GMT_CTRL *GMT, double x, double y, double t, double w, int justify, char *P1, char *P2, char *label, char kind, double fsize) {
	/* Place rounded line indicator with a single crossmark line and optional label.
	 * If label is in percent then we strip off the percent and place it at end of line insteadl */
	int text_justify, just_p;
	char *p = NULL;
	double fx, fy, dy, dy2, xt, del_x = 0.0, del_y = 0.0, angle = 0.0, xp;
	struct GMT_PEN pen;

	gmtplot_just_f_xy (justify, &fx, &fy);
	gmt_M_memset (&pen, 1, struct GMT_PEN);
	gmt_getpen (GMT, P2, &pen);	/* Pen for the background fixed line */
	dy = pen.width / PSL_POINTS_PER_INCH;	/* Pen width in inches */
	dy2 = dy / 2.0;	/* Half pen-width */
	x += fx * (w + dy);	y += fy * dy;	/* Adjust center (x,y) depending on justification */
	xt = w * (t - 0.5);	/* Location of label and tick along x-axis */
	w /= 2;	/* From here on, w is half-length of fixed line */
	if (kind == 'D') {	/* Want a label at the current crossmark */
		if (fsize == 0.0) fsize = 2.0 * dy * PSL_POINTS_PER_INCH;	/* Set a scaled font size */
		if ((p = strchr (label, '%'))) p[0] = '\0';	/* Remove % here and add separately later */
		switch (justify) {	/* Deal with justification of text and possibly rotation */
			case PSL_BL: case PSL_BC: case PSL_BR:
				text_justify = PSL_BC; del_y = dy; break;
			case PSL_TL: case PSL_TC: case PSL_TR: case PSL_MC:
				text_justify = PSL_TC; del_y = -dy; break;
			case PSL_ML: text_justify = PSL_ML; angle = 90.0; x -= w; del_y = -dy; break;
			case PSL_MR: text_justify = PSL_MR; angle = 90.0; x += w; del_y = +dy; break;
			default:
				text_justify = PSL_BC; del_y = dy; break;
		}
		if (p) {	/* Got percentages; must plot it separately */
			if (justify == PSL_ML || justify == PSL_MR)
				just_p = PSL_BC, xp = w + dy;
			else {
				just_p = (justify % 4 == 3) ? PSL_MR : PSL_ML;
				xp = (just_p == PSL_MR) ? -dy-w : w+dy;
			}
		}
	}
	PSL_setorigin (GMT->PSL, x, y, angle, PSL_FWD);	/* Origin is now mid-point of rounded line */
	if (kind == 'D') {
		PSL_plottext (GMT->PSL, xt+del_x, del_y, fsize, label, -angle, text_justify, 0);
		if (p) {
			PSL_plottext (GMT->PSL, xp, 0.0, fsize, "%", -angle, just_p, 0);
			p[0] = '%';	/* Restore it */
		}
	}
	PSL_command (GMT->PSL, "V\n");	/* Do a save/restore around the line cap change */
	gmt_setpen (GMT, &pen);
	PSL_setlinecap (GMT->PSL, PSL_ROUND_CAP);	/* Want the rounded line effect */
	PSL_plotsegment (GMT->PSL, -w, 0.0, +w, 0.0);	/* Draw thick rounded line */
	PSL_command (GMT->PSL, "U\n");	/* grestore to restore to whatever line cap was set to before */
	gmt_M_memset (&pen, 1, struct GMT_PEN);	/* Wipe pen */
	gmt_getpen (GMT, P1, &pen);	/* Get pen for crossmark */
	gmt_setpen (GMT, &pen);
	PSL_plotsegment (GMT->PSL, xt, -dy2, xt, dy2);	/* Draw foreground crossmark at label time */
	PSL_setorigin (GMT->PSL, -x, -y, -angle, PSL_INV);	/* Undo coordinate transformation */
	PSL_setlinecap (GMT->PSL, PSL_BUTT_CAP);	/* Want the rounded line effect */
}

GMT_LOCAL void gmtplot_prog_indicator_E (struct GMT_CTRL *GMT, double x, double y, double t, double w, int justify, char *P1, char *P2, char *label, char kind, double fsize) {
	/* Place solid line indicator with optional start/end labels */
	int text_justify, just_p;
	double fx, fy, dy, dy2, xt, xp, del_x = 0.0, del_y = 0.0, angle = 0.0;
	struct GMT_PEN pen;

	gmtplot_just_f_xy (justify, &fx, &fy);
	gmt_M_memset (&pen, 1, struct GMT_PEN);
	gmt_getpen (GMT, P2, &pen);	/* Pen for the background fixed line */
	dy = pen.width / PSL_POINTS_PER_INCH;	/* Pen width in inches */
	dy2 = dy / 2.0;	/* Half pen-width */
	x += fx * (w + dy2);	y += fy * dy;	/* Adjust center (x,y) depending on justification */
	xt = w * (t - 0.5);	/* Point of current time on axis */
	w /= 2;	/* From here on, w is half-length of fixed line length */
	if (kind == 'E') {	/* Want a label at the current crossmark */
		if (fsize == 0.0) fsize = 2.0 * dy * PSL_POINTS_PER_INCH;	/* Set a scaled font size */
		switch (justify) {	/* Deal with justification of text and possibly rotation */
			case PSL_BL: case PSL_BC: case PSL_BR:
				text_justify = PSL_BC; del_y = dy; break;
			case PSL_TL: case PSL_TC: case PSL_TR: case PSL_MC:
				text_justify = PSL_TC; del_y = -dy; break;
			case PSL_ML: text_justify = PSL_ML; angle = 90.0; x -= w; del_y = -dy; break;
			case PSL_MR: text_justify = PSL_MR; angle = 90.0; x += w; del_y = +dy; break;
			default:
				text_justify = PSL_BC; del_y = dy; break;
		}
	}
	PSL_setorigin (GMT->PSL, x, y, angle, PSL_FWD);	/* Origin is now mid-point of rounded line and we have possibly rotated 90 degrees */
	if (kind == 'E') {					/* Place the two labels */
		char *c = strchr (label, ';'), *p = NULL;
		if (c) c[0] = '\0';				/* Chop of label at t == 1 and place the first at t == 0 */
		if ((p = strchr (label, '%')))	/* Watch for percentage sign; if found place 0 and add % to axis label */
			PSL_plottext (GMT->PSL, del_x-w, del_y, fsize, "0", -angle, text_justify, 0);
		else
			PSL_plottext (GMT->PSL, del_x-w, del_y, fsize, label, -angle, text_justify, 0);
		if ((p = strchr (&c[1], '%')))	/* Watch for percentage sign; if found place 100 and add % to axis label */
			PSL_plottext (GMT->PSL, del_x+w, del_y, fsize, "100", -angle, text_justify, 0);
		else
			PSL_plottext (GMT->PSL, del_x+w, del_y, fsize, &c[1], -angle, text_justify, 0);
		if (p) {	/* Since we found percentage we add it as an axis label */
			if (justify == PSL_ML || justify == PSL_MR)
				just_p = PSL_BC, xp = w + dy;
			else {
				just_p = (justify % 4 == 3) ? PSL_MR : PSL_ML;
				xp = (just_p == PSL_MR) ? -dy-w : w+dy;
			}
			PSL_plottext (GMT->PSL, xp, 0.0, fsize, "%", -angle, just_p, 0);
		}
		c[0] = ';';	/* Restore the semi-colon separator */
	}
	gmt_setpen (GMT, &pen);
	PSL_plotsegment (GMT->PSL, -w, 0.0, +w, 0.0);	/* Draw thick rounded line */
	gmt_M_memset (&pen, 1, struct GMT_PEN);
	gmt_getpen (GMT, P1, &pen);	/* Pen for time-variable line */
	gmt_setpen (GMT, &pen);
	PSL_plotsegment (GMT->PSL, w, 0.0, xt, 0.0);
	PSL_setorigin (GMT->PSL, -x, -y, -angle, PSL_INV);	/* Undo coordinate transformation */
}

GMT_LOCAL void gmtplot_prog_indicator_F (struct GMT_CTRL *GMT, double x, double y, double t, double w, int justify, char *P1, char *F1, char *label, char kind, double width, double fsize, struct GMT_FONT *F) {
	/* Place time axis indicator via call to basemap plus adding triangle here */
	int symbol = PSL_INVTRIANGLE;	/* Time marker when labels are below the line */
	char cmd[GMT_LEN256] = {""}, region[GMT_LEN64] = {""}, unit[4] = {""}, axis = 0;
	bool was = GMT->current.map.frame.init;
	double fx, fy, dy, dy2, xt, s = 1.0, angle = 0.0, h;
	struct GMT_PEN pen;
	struct GMT_FILL fill;

	gmtplot_just_f_xy (justify, &fx, &fy);
	gmt_M_memset (&pen, 1, struct GMT_PEN);
	gmt_getpen (GMT, P1, &pen);	/* Pen for the time axis */
	dy = pen.width / PSL_POINTS_PER_INCH;	/* Need to get pen width to determine triangle size */
	dy2 = dy / 2.0;	/* Half pen width */
	x += fx * w;	y += fy * dy;	/* Adjust the (x,y) to be mid-point of axis */
	xt = w * t;	/* Relative fraction of w along the axis from left point */
	/* Must allow for the space needed for tick-mark, annot-offset, and possibly annotations */
	h = GMT->current.setting.map_tick_length[GMT_PRIMARY] + MAX (0.0, GMT->current.setting.map_annot_offset[GMT_PRIMARY]);
	if (kind == 'F') h += GMT_LET_HEIGHT * GMT->current.setting.font_annot[GMT_PRIMARY].size / PSL_POINTS_PER_INCH;
	x -= 0.5 * w;	/* Let x be start position of axis rather than middle */
	switch (justify) {	/* Deal with justification and adjust x,y for starting point of axis */
		case PSL_BL: case PSL_BC: case PSL_BR:	/* At bottom we need to move up to give space for annotations below axis */
			axis = 'S'; y += h; break;
		case PSL_TL: case PSL_TC: case PSL_TR: case PSL_MC:	/* At top (and middle) we likewise need space above, plus switch to TRIANGLE */
			axis = 'N'; y-= h; s = -1; symbol = PSL_TRIANGLE; break;
		case PSL_ML: angle = 90.0; y -= 0.5*w; axis = 'N'; x += h; s = -1; symbol = PSL_TRIANGLE; break;	/* Left is like TC but rotated 90 */
		case PSL_MR: angle = 90.0; y -= 0.5*w; axis = 'S'; x += w - h; break;	/* Right is like BC but rotated 90 */
	}
	if (strchr (label, '%')) {	/* Cannot have percentages in the -R string */
		sprintf (region, "-R0/100/0/1");	/* Always 0-100 % */
		strcpy (unit, "+u%");
	}
	else if (label[2] == ' ') { /* Eliminate leading spaces... */
		int k = 2; while (label[k] && label[k] == ' ') k++;	/* Find first non-space character */
		sprintf (region, "-R%s", &label[k]);
	}
	else	/* Use as is */
		strncpy (region, label, GMT_LEN64-1);
	PSL_setorigin (GMT->PSL, x, y, angle, PSL_FWD);	/* Origin (0,0) is now at left end-point of time axis */
		GMT->common.R.active[RSET] = GMT->common.J.active = false;
	if (fsize == 0.0)
			F->size = 3.0 * dy * PSL_POINTS_PER_INCH;	/* Set a scaled font size */
	if (kind == 'F')	/* Want annotations so add -Baf and set annotation font */
		sprintf (cmd, "%s -JX%gi/0.0001i -Baf%s -B%c -X0 -Y0 --MAP_FRAME_PEN=%s --FONT_ANNOT_PRIMARY=+%s --GMT_HISTORY=false", region, width, unit, axis, P1, gmt_putfont (GMT, F));
	else	/* Only axis with ticks so only -Bf is used */
		sprintf (cmd, "%s -JX%gi/0.0001i -Bf%s -B%c -X0 -Y0 --MAP_FRAME_PEN=%s --GMT_HISTORY=false", region, width, unit, axis, P1);
	GMT->current.map.frame.init = false;	/* To enable more -B parsing */
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Call basemap from gmtplot_prog_indicator_F with args %s\n", cmd);
	if (GMT_Call_Module (GMT->parent, "basemap", GMT_MODULE_CMD, cmd) != GMT_NOERROR) {
		GMT_Report (GMT->parent, GMT_MSG_ERROR, "Failed to call basemap with args %s ???\n", cmd);
	}
	GMT->current.map.frame.init = was;	/* Reset how we found it */
	if (gmt_getfill (GMT, F1, &fill))	/* Get and set color for the triangle; no outline */
		GMT_Report (GMT->parent, GMT_MSG_ERROR, "Bad fill argument %s\n", F1);
	PSL_setfill (GMT->PSL, fill.rgb, 0);
	w = dy * 2;	/* Set triangle size to 2 times the axis width */
	PSL_plotsymbol (GMT->PSL, xt, 1.2*s*dy2, &w, symbol);
}

GMT_LOCAL void gmtplot_reset_PSL (struct GMT_CTRL *GMT, struct PSL_CTRL *PSL) {
	/* Because PSL_*_completion procedures are called at the end of the module but crafted at the start, we must forget we used any fonts, fills, pens in creating them,
	 * otherwise the next PSL call requesting a fill or pen will not be properly set since PSL may think it is already the current choice. */
	if (PSL->current.font_no != GMT_NOTSET) PSL->internal.font[PSL->current.font_no].encoded = 0;	/* Forget about font encoding */
	PSL->current.font_no   = GMT_NOTSET;		/* Forget what the current font is */
	PSL->current.linewidth = GMT_NOTSET;		/* Forget we ever set a line width */
	PSL->current.outline   = GMT_NOTSET;		/* Forget we requested polygon outline */
	gmt_M_memcpy (PSL->current.rgb[PSL_IS_FONT],   GMT->session.no_rgb, 3, double);	/* Reset to -1,-1,-1 since text setting must set the color desired */
	gmt_M_memcpy (PSL->current.rgb[PSL_IS_FILL],   GMT->session.no_rgb, 3, double);	/* Reset to -1,-1,-1 since fill setting must set the color desired */
	gmt_M_memcpy (PSL->current.rgb[PSL_IS_STROKE], GMT->session.no_rgb, 3, double);	/* Reset to -1,-1,-1 since stroke setting must set the color desired */
}

struct PSL_CTRL *gmt_plotinit (struct GMT_CTRL *GMT, struct GMT_OPTION *options) {
	/* Shuffles parameters and calls PSL_beginplot, issues PS comments regarding the GMT options
	 * and places a time stamp, if selected */

	int k, id, fno[PSL_MAX_EPS_FONTS], n_fonts, last, form, justify;
	bool O_active = GMT->common.O.active, auto_media = false, do_paint = false;
	unsigned int this_proj, write_to_mem = 0, switch_charset = 0, n_movie_items[2] = {0, 0};
	char *mode[2] = {"w","a"}, *movie_item_arg[2][GMT_LEN32], not_used[GMT_LEN32] = {""};
	static char *ps_mode[2] = {"classic", "modern"}, *F_name[2] = {"label", "prog_indicator"};
	double media_size[2], plot_x, plot_y;
	FILE *fp = NULL;	/* Default which means stdout in PSL */
	struct GMT_FILL fill;
	struct GMT_OPTION *Out = NULL;
	struct PSL_CTRL *PSL = NULL;
	struct GMT_SUBPLOT *P = NULL;
	struct GMT_INSET *I = &(GMT->current.plot.inset);	/* I->active == 1 if an inset */

	PSL = GMT->PSL;	/* Shorthand */

	gmt_set_undefined_defaults (GMT, sqrt (GMT->current.map.width * GMT->current.map.height), false);	/* We must change any undefined defaults given max plot dimension */

	if (GMT->current.map.frame.paint[GMT_Z]) {	/* Must squirrel this away for now since we may call psbasemap during the movie-indicators below */
		do_paint = true;
		gmt_M_memcpy (&fill, &GMT->current.map.frame.fill[GMT_Z], 1U, struct GMT_FILL);
		GMT->current.map.frame.paint[GMT_Z] = false;	/* Turn off for now */
	}

	PSL->internal.verbose = GMT->current.setting.verbose;		/* Inherit verbosity level from GMT */
	if (gmt_M_compat_check (GMT, 4) && GMT->current.setting.ps_copies > 1) PSL->init.copies = GMT->current.setting.ps_copies;
	PSL_setdefaults (PSL, GMT->current.setting.ps_magnify, GMT->current.setting.ps_page_rgb, GMT->current.setting.ps_encoding.name);
	GMT->current.ps.memory = false;
	if (doubleAlmostEqual (GMT->current.setting.ps_page_size[GMT_X], GMT_PAPER_DIM) && doubleAlmostEqual (GMT->current.setting.ps_page_size[GMT_Y], GMT_PAPER_DIM))
		auto_media = true;
	gmt_M_memcpy (media_size, GMT->current.setting.ps_page_size, 2, double);

	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Running in PS mode %s\n", ps_mode[GMT->current.setting.run_mode]);
	if (GMT->current.setting.run_mode == GMT_MODERN) {	/* Write PS to hidden gmt_#.ps- file.  No -O -K allowed */
		char *verb[2] = {"Create", "Append to"};
		bool wants_PS;
		double paper_margin = GMT_PAPER_MARGIN_AUTO;

		if (gmtlib_fixed_paper_size (GMT->parent)) {	/* Must honor paper size and regular margin */
			paper_margin = GMT_PAPER_MARGIN_FIXED;
			gmt_M_memcpy (media_size, GMT->current.setting.ps_def_page_size, 2, double);
		}

		if ((k = gmt_set_psfilename (GMT)) == GMT_NOTSET) {	/* Get hidden file name for PS */
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "No workflow directory\n");
			GMT->parent->error = GMT_ERROR_ON_FOPEN;
			return NULL;
		}
		GMT_Report (GMT->parent, GMT_MSG_DEBUG, "%s hidden PS file %s\n", verb[k], GMT->current.ps.filename);
		if ((fp = PSL_fopen (PSL, GMT->current.ps.filename, mode[k])) == NULL) {	/* Must open inside PSL DLL */
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Cannot open %s with mode %s\n", GMT->current.ps.filename, mode[k]);
			GMT->parent->error = GMT_ERROR_ON_FOPEN;
			return NULL;
		}
		O_active = (k) ? true : false;	/* -O is determined by presence or absence of hidden PS file */
		/* Determine paper size */
		wants_PS = gmtlib_fig_is_ps (GMT);	/* True if we have requested a PostScript output format */
		if (wants_PS && !O_active) {	/* Requesting a new PostScript file in modern mode */
			if (auto_media) {	/* Cannot use "auto" if requesting a PostScript file */
				GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "You should specify a paper size when requesting a PostScript file.\n");
				if (GMT->current.setting.proj_length_unit == GMT_INCH) {	/* Use US settings */
					GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "Changing paper size to US Letter, but we cannot know if this is adequate for your plot; better to use PS_MEDIA explicitly.\n");
					media_size[GMT_X] = 612.0; media_size[GMT_Y] = 792.0;
				}
				else {	/* Use SI settings */
					GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "Changing paper size to A4, but we cannot know if this is adequate for your plot; better to use PS_MEDIA explicitly.\n");
					media_size[GMT_X] = 595.0; media_size[GMT_Y] = 842.0;
				}
			}
			if ((GMT->current.map.width > GMT->current.map.height) && (((GMT->current.map.width + GMT->current.setting.map_origin[GMT_X]) * 72) > media_size[GMT_X]) && GMT->current.setting.ps_orientation == PSL_PORTRAIT) {
				GMT->current.setting.ps_orientation = PSL_LANDSCAPE;
				GMT_Report (GMT->parent, GMT_MSG_WARNING, "Changing to PostScript landscape orientation based on your plot and paper dimensions, but we cannot be 100%% sure.\n");
				GMT_Report (GMT->parent, GMT_MSG_WARNING, "Use PS_MEDIA and/or PS_PAGE_ORIENTATION to specify correct paper dimensions and/or orientation if our guesses are inadequate.\n");
			}
		}
		else if (!O_active) {	/* Not desiring PS output so we can add safety margin of paper_margin inches for initial layer unless PS_MEDIA was set */
			if (!(GMT->common.X.active || GMT->common.Y.active) && auto_media)
				GMT->current.setting.map_origin[GMT_X] = GMT->current.setting.map_origin[GMT_Y] = paper_margin;
		}
		if (!O_active) {	/* See if special movie labeling file exists under modern mode */
			char file[PATH_MAX] = {""}, record[GMT_LEN256] = {""};
			FILE *fpl = NULL;
			for (k = 0; k < 2; k++) {
				snprintf (file, PATH_MAX, "%s/gmt.movie%ss", GMT->parent->gwf_dir, F_name[k]);
				if (!access (file, R_OK) && (fpl = fopen (file, "r"))) {	/* File exists and could be opened for reading */
					while (fgets (record, GMT_LEN256, fpl)) {
						if (record[0] == '#') continue;	/* Skip header */
						if (n_movie_items[k] == GMT_LEN32) {
							GMT_Report (GMT->parent, GMT_MSG_ERROR, "Number of movie %ss exceed capacity [%d] - skipped.\n", F_name[k], n_movie_items[k]);
							continue;
						}
						gmt_chop (record);		/* Chop off LF or CR/LF */
						GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Found MOVIE_LABEL_ARG%d = %s.\n", n_movie_items[k], record);
						movie_item_arg[k][n_movie_items[k]++] = strdup (record);
					}
					fclose (fpl);
					if (gmt_remove_file (GMT, file)) {	/* Now remove the file since done */
						GMT_Report (GMT->parent, GMT_MSG_ERROR, "Cannot delete file %s\n", file);
						GMT->parent->error = GMT_ERROR_ON_FOPEN;
						return NULL;
					}
				}
			}
		}
	}
	else if ((Out = GMT_Find_Option (GMT->parent, '>', options))) {	/* Want to use a specific output file */
		k = (Out->arg[0] == '>') ? 1 : 0;	/* Are we appending (k = 1) or starting a new file (k = 0) */
		if (O_active && k == 0 && !gmt_M_file_is_memory (Out->arg)) {
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "-O given but append-mode not selected for file %s\n", &(Out->arg[k]));
		}
		if (gmt_M_file_is_memory (&(Out->arg[k]))) {
			write_to_mem = 2;
			GMT->current.ps.memory = true;
			strncpy (GMT->current.ps.memname, &(Out->arg[k]), GMT_VF_LEN-1);
		}
		else {
			if ((fp = PSL_fopen (PSL, &(Out->arg[k]), mode[k])) == NULL) {	/* Must open inside PSL DLL */
				GMT_Report (GMT->parent, GMT_MSG_ERROR, "Cannot open %s with mode %s\n", &(Out->arg[k]), mode[k]);
				GMT->parent->error = GMT_ERROR_ON_FOPEN;
				return NULL;
			}
			GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Opened PS file %s\n", &(Out->arg[k]));
		}
	}

    if (GMT->common.P.active) GMT->current.setting.ps_orientation = PSL_PORTRAIT;

    if (GMT->current.setting.run_mode == GMT_CLASSIC && !O_active) {  /* Warn if plot dimensions are larger than current paper size */
		unsigned int X = (GMT->current.setting.ps_orientation == PSL_LANDSCAPE);
		unsigned int Y = (GMT->current.setting.ps_orientation == PSL_PORTRAIT);
		if ((((GMT->current.map.width + GMT->current.setting.map_origin[GMT_X]) * 72) > GMT->current.setting.ps_def_page_size[X]) ||
			(((GMT->current.map.height + GMT->current.setting.map_origin[GMT_Y]) * 72) > GMT->current.setting.ps_def_page_size[Y])) {
			double OX = GMT->current.setting.map_origin[GMT_X] * GMT->session.u2u[GMT_INCH][GMT->current.setting.proj_length_unit];
			double OY = GMT->current.setting.map_origin[GMT_Y] * GMT->session.u2u[GMT_INCH][GMT->current.setting.proj_length_unit];
			double W = GMT->current.map.width * GMT->session.u2u[GMT_INCH][GMT->current.setting.proj_length_unit];
			double H = GMT->current.map.height * GMT->session.u2u[GMT_INCH][GMT->current.setting.proj_length_unit];
			double PW = GMT->current.setting.ps_def_page_size[X] * GMT->session.u2u[GMT_PT][GMT->current.setting.proj_length_unit];
			double PH = GMT->current.setting.ps_def_page_size[Y] * GMT->session.u2u[GMT_PT][GMT->current.setting.proj_length_unit];
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "Your plot (WxH = %.2lf x %.2lf %s) placed at (%.2lf, %.2lf %s) may exceed your PS_MEDIA (WxH = %.2lf x %.2lf %s)\n",
				W, H, GMT->session.unit_name[GMT->current.setting.proj_length_unit], OX, OY, GMT->session.unit_name[GMT->current.setting.proj_length_unit], PW, PH, GMT->session.unit_name[GMT->current.setting.proj_length_unit]);
		}
	}

	/* Initialize the plot header and settings */

	/* Default for overlay plots is no shifting */

	if (!GMT->common.X.active && O_active) GMT->current.setting.map_origin[GMT_X] = 0.0;
	if (!GMT->common.Y.active && O_active) GMT->current.setting.map_origin[GMT_Y] = 0.0;

	/* Adjust offset when centering plot on center of page (PS does the rest) */

	if (GMT->current.ps.origin[GMT_X] == 'c') GMT->current.setting.map_origin[GMT_X] -= 0.5 * GMT->current.map.width;
	if (GMT->current.ps.origin[GMT_Y] == 'c') GMT->current.setting.map_origin[GMT_Y] -= 0.5 * GMT->current.map.height;

	/* Get font names used */

	id = 0;
	if (GMT->common.U.active) fno[id++] = GMT->current.setting.font_logo.id;	/* Add GMT logo font */
	/* Add title font if a title was used */
	if (GMT->current.map.frame.header[0]) fno[id++] = GMT->current.setting.font_title.id;
	/* Add the label font if labels were used */
	if (GMT->current.map.frame.axis[GMT_X].label[0] || GMT->current.map.frame.axis[GMT_Y].label[0] || GMT->current.map.frame.axis[GMT_Z].label[0])
		fno[id++] = GMT->current.setting.font_label.id;
	/* Always add annotation fonts */
	fno[id++] = GMT->current.setting.font_annot[GMT_PRIMARY].id;
	fno[id++] = GMT->current.setting.font_annot[GMT_SECONDARY].id;

	gmt_sort_array (GMT, fno, id, GMT_INT);

	last = GMT_NOTSET;
	for (k = n_fonts = 0; k < id; k++) {
		if (fno[k] != last) last = fno[n_fonts++] = fno[k]; /* To avoid duplicates */
	}
	for (k = n_fonts; k < PSL_MAX_EPS_FONTS; k++) fno[k] = GMT_NOTSET;	/* Terminate */

	P = &(GMT->current.plot.panel);

	if (P->active) {	/* Subplot panel mode is in effect */
		/* Consider offsets required to center the plot on the subplot panel [0/0] */
		GMT->current.setting.map_origin[GMT_X] += (P->dx + P->gap[XLO]);
		GMT->current.setting.map_origin[GMT_Y] += (P->dy + P->gap[YLO]);
		if (P->first && O_active)	/* Run completion script, if any */
			PSL_setexec (PSL, 1);
	}
	if (I->active && I->first) {
		GMT->current.setting.map_origin[GMT_X] += (I->dx);
		GMT->current.setting.map_origin[GMT_Y] += (I->dy);
		I->active = I->first = false;	/* For MATLAB mostly */
	}

	if (O_active && GMT->current.ps.switch_set) {	/* User used --PS_CHAR_ENCODING=<encoding> to change it */
		GMT->current.ps.switch_set = false;
		switch_charset = 4;
	}

	/* Get title */

	if (GMT->current.setting.run_mode == GMT_MODERN)
		snprintf (GMT->current.ps.title, GMT_LEN256, "GMT v%s Document from %s", GMT_VERSION, gmt_current_name (GMT->init.module_name, not_used));
	else
		snprintf (GMT->current.ps.title, GMT_LEN256, "GMT v%s Document from %s", GMT_VERSION, GMT->init.module_name);

	PSL_beginplot (PSL, fp, GMT->current.setting.ps_orientation|write_to_mem|switch_charset, O_active, GMT->current.setting.ps_color_mode,
		GMT->current.ps.origin, GMT->current.setting.map_origin, media_size, GMT->current.ps.title, fno);

	/* Issue the comments that allow us to trace down what command created this layer */

	gmtplot_echo_command (GMT, PSL, options);

	if (GMT->common.p.do_z_rotation) {	/* Need a rotation about z of the whole page */
		double x0 = 0.0, y0 = 0.0;	/* Default is to rotate around plot origin */
		if (GMT->current.proj.z_project.view_given) {	/* Rotation is about another z-axis than through the origin */
			x0 = GMT->current.proj.z_project.view_x;
			y0 = GMT->current.proj.z_project.view_y;
		}
		else if (GMT->current.proj.z_project.world_given)	/* Rotation is about another lon/lat pole than the origin */
			gmt_geo_to_xy (GMT, GMT->current.proj.z_project.world_x, GMT->current.proj.z_project.world_y, &x0, &y0);
		GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Transrot: Rotating plot by %g degrees about (%g, %g)\n", GMT->common.p.z_rotation, x0, y0);
		PSL_comment (GMT->PSL, "Possibly translate then rotate whole page\n");
		PSL_setorigin (PSL, x0, y0, GMT->common.p.z_rotation, PSL_FWD);
		PSL_setorigin (PSL, -x0, -y0, 0.0, PSL_FWD);
	}

	/* Create %%PROJ tag that psconvert can use to prepare a ESRI world file */

	this_proj = GMT->current.proj.projection_GMT;
	for (k = 0, id = GMT_NOTSET; id == GMT_NOTSET && k < GMT_N_PROJ4; k++)
		if (GMT->current.proj.proj4[k].id == this_proj) id = k;
	if (id >= 0) {			/* Valid projection for creating world file info */
		double Cartesian_m[4];	/* WESN equivalents in projected meters */
		char *pstr = NULL, proj4name[16] = {""};
		Cartesian_m[0] = (GMT->current.proj.rect[YLO] - GMT->current.proj.origin[GMT_Y]) * GMT->current.proj.i_scale[GMT_Y];
		Cartesian_m[1] = (GMT->current.proj.rect[XHI] - GMT->current.proj.origin[GMT_X]) * GMT->current.proj.i_scale[GMT_X];
		Cartesian_m[2] = (GMT->current.proj.rect[YHI] - GMT->current.proj.origin[GMT_Y]) * GMT->current.proj.i_scale[GMT_Y];
		Cartesian_m[3] = (GMT->current.proj.rect[XLO] - GMT->current.proj.origin[GMT_X]) * GMT->current.proj.i_scale[GMT_X];
		/* It woul be simpler if we had a cleaner way of telling when data is lon-lat */
		if (GMT->current.proj.projection_GMT == GMT_LINEAR && gmt_M_is_geographic (GMT, GMT_IN))
			strcpy (proj4name, "latlong");
		else
			strncpy (proj4name, GMT->current.proj.proj4[id].name, 15U);

		pstr = gmt_export2proj4 (GMT);
		PSL_command (PSL, "%%@PROJ: %s %.8f %.8f %.8f %.8f %.3f %.3f %.3f %.3f %s\n", proj4name,
			GMT->common.R.wesn[XLO], GMT->common.R.wesn[XHI], GMT->common.R.wesn[YLO], GMT->common.R.wesn[YHI],
			Cartesian_m[3], Cartesian_m[1], Cartesian_m[0], Cartesian_m[2], pstr);
		gmt_M_str_free (pstr);
	}
	if (!O_active) {	/* New plot, set GMT bounding box and reset layer counter */
		/* Add GMT bounding box comment which gives dimensions of inner frame */
		PSL_command (PSL, "%%GMTBoundingBox: %.12g %.12g %.12g %.12g\n",
			GMT->current.setting.map_origin[GMT_X] * 72.0, GMT->current.setting.map_origin[GMT_Y] * 72.0,
			GMT->current.map.width * 72.0, GMT->current.map.height * 72.0);
		GMT->current.ps.layer = 0;
	}

	PSL_beginlayer (GMT->PSL, ++GMT->current.ps.layer);
	/* Set layer transparency, if requested. Note that PSL_transp actually sets the opacity alpha, which is (1 - transparency) */
	if (GMT->common.t.active) {
		if (GMT->common.t.value[GMT_FILL_TRANSP] == 0.0 && GMT->common.t.value[GMT_PEN_TRANSP] == 0.0) {
			GMT_Report (GMT->parent, GMT_MSG_INFORMATION, "A transparency of 0/0 is the same as opaque. Skipped\n");
			GMT->common.t.active = false;
		}
		else 	/* Place both fill and stroke transparencies in 0-1 normalized range, plus the blend mode name */
			PSL_command (PSL, "%.12g %.12g /%s PSL_transp\n", 1.0 - 0.01 * GMT->common.t.value[GMT_FILL_TRANSP], 1.0 - 0.01 * GMT->common.t.value[GMT_PEN_TRANSP], GMT->current.setting.ps_transpmode);
	}
	/* If requested, place the timestamp */

	if (GMT->current.ps.logo_cmd) {
		char txt[4] = {' ', '-', 'X', 0}, not_used[GMT_LEN32] = {""};
		struct GMT_OPTION *opt = NULL;
		size_t len;
		/* -Uc was given as shorthand for "plot current command line" */
		if (GMT->current.setting.run_mode == GMT_MODERN)
			strncpy (GMT->current.ps.map_logo_label, gmt_current_name (GMT->init.module_name, not_used), GMT_LEN256-1);
		else
			strncpy (GMT->current.ps.map_logo_label, GMT->init.module_name, GMT_LEN256-1);
		len = strlen (GMT->current.ps.map_logo_label);
		for (opt = options; opt; opt = opt->next) {
			if (opt->option == GMT_OPT_INFILE || opt->option == GMT_OPT_OUTFILE) continue;	/* Skip file names */
			txt[2] = opt->option;
			strncat (GMT->current.ps.map_logo_label, txt, GMT_LEN256-len);
			len += 3;
			strncat (GMT->current.ps.map_logo_label, opt->arg, GMT_LEN256-len);
		}
		GMT->current.ps.logo_cmd = false;	/* Mission accomplished */
	}
	if (GMT->current.setting.map_logo)
		gmtplot_timestamp (GMT, PSL, GMT->current.setting.map_logo_pos[GMT_X], GMT->current.setting.map_logo_pos[GMT_Y], GMT->current.setting.map_logo_justify, GMT->current.ps.map_logo_label);

	PSL_settransparencymode (PSL, GMT->current.setting.ps_transpmode);	/* Set PDF transparency mode, if used */
	/* Enforce chosen line parameters */
	k = GMT->PSL->internal.line_cap;	GMT->PSL->internal.line_cap = GMT_NOTSET; PSL_setlinecap (PSL, k);
	k = GMT->PSL->internal.line_join;	GMT->PSL->internal.line_join = GMT_NOTSET; PSL_setlinejoin (PSL, k);
	k = GMT->PSL->internal.miter_limit;	GMT->PSL->internal.miter_limit = GMT_NOTSET; PSL_setmiterlimit (PSL, k);

	if (P->active && P->first) {
		/* Place the panel tag, once per panel (if requested), then update the gmt.panel file to say we have been there */
		if (strcmp (P->tag, "-")) {	/* Place the panel tag */
			int form, refpoint, justify;
			refpoint = gmt_just_decode (GMT, P->refpoint, PSL_NO_DEF);	/* Convert XX refpoint code to PSL number */
			gmtlib_refpoint_to_panel_xy (GMT, refpoint, P, &plot_x, &plot_y);	/* Convert just code to panel location */
			/* Undo any offsets above that was required to center the plot on the subplot panel */
			plot_x -= (P->dx);
			plot_y -= (P->dy);
			PSL_command (PSL, "/PSL_plot_completion {\nV\n");
			PSL_comment (PSL, "Start of panel tag for panel (%d,%d)\n", P->row, P->col);
			PSL_comment (PSL, "Will not execute until end of panel\n");
			PSL_setorigin (PSL, GMT->current.setting.map_origin[GMT_X] - P->gap[XLO], GMT->current.setting.map_origin[GMT_Y] - P->gap[YLO], 0.0, PSL_FWD);
			justify = gmt_just_decode (GMT, P->justify, PSL_NO_DEF);	/* Convert XX refpoint code to PSL number */
			gmt_smart_justify (GMT, justify, 0.0, P->off[GMT_X], P->off[GMT_Y], &plot_x, &plot_y, 1);	/* Shift as requested */
			form = gmt_setfont (GMT, &GMT->current.setting.font_tag);	/* Set the tag font */
			/* Because this is run later we must force the correct font color here */
			PSL_command (PSL, PSL_makecolor (PSL, GMT->current.setting.font_tag.fill.rgb));
			PSL_command (PSL, " ");
			PSL_setfont (PSL, GMT->current.setting.font_tag.id);
			if (P->pen[0] || P->fill[0] || P->shade[0]) {	/* Must deal with textbox fill/outline/shade */
				/* All fills and pens have already gone through a syntax check in subplot -A */
				int outline = 0;
				struct GMT_FILL fill;
				struct GMT_PEN pen;
				gmt_M_memset (&pen, 1, struct GMT_PEN);
				gmt_init_fill (GMT, &fill, -1.0, -1.0, -1.0);	/* No fill */
				PSL_command (PSL, "FQ O0\n");	/* Ensure fill/pen have been reset */
				if (P->shade[0] && !gmt_getfill (GMT, P->shade, &fill)) {  /* Want to paint an offset, shaded rectangle behind the tag box */
					PSL_setfill (PSL, fill.rgb, 0);  /* Shade color */
					PSL_plottextbox (PSL, plot_x+P->soff[GMT_X], plot_y+P->soff[GMT_Y], GMT->current.setting.font_tag.size, P->tag, 0.0, justify, P->clearance, 0);
				}
				if (P->pen[0] && !gmt_getpen (GMT, P->pen, &pen)) {	/* Want to draw the outline of the tag box */
					gmt_setpen (GMT, &pen);
					outline = 1;
				}
				if (P->fill[0] && !gmt_getfill (GMT, P->fill, &fill)) {  /* Want to paint inside of tag box */
					PSL_setfill (PSL, fill.rgb, outline);	/* Box color and possible outline */
					PSL_plottextbox (PSL, plot_x, plot_y, GMT->current.setting.font_tag.size, P->tag, 0.0, justify, P->clearance, 0);
				}
				form = gmt_setfont (GMT, &GMT->current.setting.font_tag);	/* Set the tag font */
				/* Finally place the tag text in the box */
				PSL_plottext (PSL, plot_x, plot_y, GMT->current.setting.font_tag.size, NULL, 0.0, justify, form);
			}
			else   /* Just place the tag text */
				PSL_plottext (PSL, plot_x, plot_y, GMT->current.setting.font_tag.size, P->tag, 0.0, justify, form);
			gmtplot_reset_PSL (GMT, PSL);	/* Reset anything we may have set in building the completion PS procedure */
			PSL_comment (PSL, "End of panel tag for panel (%d,%d)\n", P->row, P->col);
			PSL_command (PSL, "U\n}!\n");
		}
		/* Store first = 0 since we are done with -B and the optional tag */
		if (gmt_set_current_panel (GMT->parent, GMT->current.ps.figure, P->row, P->col, P->gap, P->tag, 0))
			return NULL;	/* Should never happen */
	}
	if (n_movie_items[MOVIE_ITEM_IS_LABEL]) {	/* Obtained movie frame labels, implement them via a completion PostScript procedure */
		/* Decode x/y/just/clearance_x/clearance_Y|offX|offY//pen/-/fill/-/font/txt in MOVIE_LABEL_ARG */
		double clearance[2] = {0.0, 0.0}, soff[2] = {0.0, 0.0};
		char FF[GMT_LEN64] = {""}, FS[GMT_LEN64] = {""}, PP[GMT_LEN64] = {""}, font[GMT_LEN64] = {""}, label[GMT_LEN256] = {""};
		int kk, nc, box;
		unsigned int T;
		struct GMT_FONT Tfont;
		k = MOVIE_ITEM_IS_LABEL;	/* Just a short hand to avoid repeating that long constant */
		/* Create special PSL_movie_label_completion PostScript procedure and define it here in the output PS.
		 * It will be called at the end of all plotting in PSL_endplot. It always exist but is NULL by default.
		 * If not NULL then it will plot 1-32 labels set via movie.c's -L option */

		PSL_command (PSL, "/PSL_movie_label_completion {\nV\n");
		PSL_comment (PSL, "Start of movie labels\n");
		PSL_comment (PSL, "Will not execute until end of plot\n");

		for (T = 0; T < n_movie_items[k]; T++) {
			GMT_Report (GMT->parent, GMT_MSG_DEBUG, "%d:  %s\n", T, movie_item_arg[k][T]);
			/* Parse -|x|y|-|-|just|clearance_x|clearance_Y|offX|offY|pen|-|fill|-|box|box_X|box_Y|sfill|font|txt in MOVIE_LABEL_ARG# strings (- means we don't parse but skip) */
			/* Replace the 19 leading slashes first with spaces */
			for (kk = nc = 0; movie_item_arg[k][T][kk] && nc < 19; kk++) if (movie_item_arg[k][T][kk] == '|') { movie_item_arg[k][T][kk] = ' '; nc++;}
			if ((kk = sscanf (movie_item_arg[k][T], "%*c %lg %lg %*s %*s %d %lg %lg %s %*s %s %*s %d %lg %lg %s %s %[^\n]", &plot_x, &plot_y, &justify, &clearance[GMT_X], &clearance[GMT_Y], PP, FF, &box, &soff[GMT_X], &soff[GMT_Y], FS, font, label)) != 13) {
				GMT_Report (GMT->parent, GMT_MSG_ERROR, "Unable to parse MOVIE_LABEL_ARG %s for 9 required items [only made %d of 13 conversions]\n", movie_item_arg[k][T], kk);
				return NULL;	/* Should never happen */
			}
			/* Because this PostScript procedure runs outside the main gsave/grestore block the origin is at (0,0) */
			if (font[0] == '+')	/* User provided a specific font/size which we must honor */
				gmt_getfont (GMT, &font[1], &Tfont);	/* We already parsed the font string in movie.c for correctness */
			else	/* Pick up the default font written in movie and scale it as needed */
				gmt_getfont (GMT, font, &Tfont);
			gmt_getfont (GMT, font, &Tfont);	/* We already parsed the font string in movie.c for correctness so no need to check here */
			form = gmt_setfont (GMT, &Tfont);	/* Obtain and set the tag font */
			PSL_setfont (PSL, Tfont.id);
			if (!(PP[0] == '-' && FF[0] == '-')) {	/* Requested textbox fill and/or outline */
				int outline = 0;	/* No outline */
				struct GMT_FILL fill;
				PSL_command (PSL, "FQ O0\n");	/* Ensure fill/pen have been reset */
				if (box == 1 || box == 5) {	/* Want a shaded box behind the text box */
					gmt_getfill (GMT, FS, &fill);	/* No need to check since we verified fill syntax in movie.c */
					PSL_setfill (PSL, fill.rgb, 0);	/* Box color (no outline for shade boxes) */
					PSL_plottextbox (PSL, plot_x+soff[GMT_X], plot_y+soff[GMT_Y], Tfont.size, label, 0.0, justify, clearance, box-1);
					box--;	/* Now 0 or 4 suitable for next call to PSL_plottextbox */
				}
				if (PP[0] != '-') {	/* Want to draw outline of tag box */
					struct GMT_PEN pen;
					gmt_M_memset (&pen, 1, struct GMT_PEN);
					gmt_getpen (GMT, PP, &pen);	/* No need to check since we verified pen syntax in movie.c */
					gmt_setpen (GMT, &pen);
					outline = 1;
				}
				gmt_init_fill (GMT, &fill, -1.0, -1.0, -1.0);	/* Initialize to no fill */
				if (FF[0] != '-')	/* Want to paint the background of tag box */
					gmt_getfill (GMT, FF, &fill);	/* No need to check since we verified fill syntax in movie.c */

				PSL_setfill (PSL, fill.rgb, outline);	/* Box color (or nothing) */
				PSL_plottextbox (PSL, plot_x, plot_y, Tfont.size, label, 0.0, justify, clearance, box);
				form = gmt_setfont (GMT, &Tfont);	/* Set the tag font (again?) */
				PSL_plottext (PSL, plot_x, plot_y, Tfont.size, NULL, 0.0, justify, form);
			}
			else
				PSL_plottext (PSL, plot_x, plot_y, Tfont.size, label, 0.0, justify, form);
			gmt_M_str_free (movie_item_arg[k][T]);	/* Done with this label */
			gmtplot_reset_PSL (GMT, PSL);	/* Reset anything we may have set in building the completion PS procedure */
		}
		PSL_comment (PSL, "End of movie labels\n");
		PSL_command (PSL, "U\n}!\n");
	}
	if (n_movie_items[MOVIE_ITEM_IS_PROG_INDICATOR]) {	/* Obtained movie frame progress indicators, implement them via a completion PostScript procedure */
		/* Decode kind|x|y|t|width|just|clearance_x|clearance_Y|offX|offY|pen|pen2|fill|fill2|-|-|-|-|font|txt in MOVIE_PROG_INDICATOR_ARG# strings */
		double clearance[2] = {0.0, 0.0}, width = 0.0, t, fsize = 0.0;
		char kind, F1[GMT_LEN64] = {""}, F2[GMT_LEN64] = {""}, P1[GMT_LEN64] = {""}, P2[GMT_LEN64] = {""}, font[GMT_LEN64] = {""}, label[GMT_LEN256] = {""};
		int kk, nc;
		unsigned int T;
		struct GMT_FONT Tfont;
		k = MOVIE_ITEM_IS_PROG_INDICATOR;	/* Just a short hand to avoid repeating that long constant */
		/* Create special PSL_movie_label_completion PostScript procedure and define it here in the output PS.
		 * It will be called at the end of all plotting in PSL_endplot. It always exist but is NULL by default.
		 * If not NULL then it will plot 1-32 progress indicators set via movie.c's -P option */

		PSL_command (PSL, "/PSL_movie_prog_indicator_completion {\n");
		PSL_comment (PSL, "Start of movie progress indicators\n");
		PSL_comment (PSL, "Will not execute until end of plot\n");

		for (T = 0; T < n_movie_items[k]; T++) {
			GMT_Report (GMT->parent, GMT_MSG_DEBUG, "%d:  %s\n", T, movie_item_arg[k][T]);
			PSL_command (PSL, "V\n");
			/* Replace the 19 leading bars first with spaces */
			for (kk = nc = 0; movie_item_arg[k][T][kk] && nc < 19; kk++) if (movie_item_arg[k][T][kk] == '|') { movie_item_arg[k][T][kk] = ' '; nc++;}
			if ((kk = sscanf (movie_item_arg[k][T], "%c %lg %lg %lg %lg %d %lg %lg %s %s %s %s %*d %*lg %*lg %*s %s %[^\n]", &kind, &plot_x, &plot_y, &t, &width, &justify, &clearance[GMT_X], &clearance[GMT_Y], P1, P2, F1, F2, font, label)) < 13) {
				GMT_Report (GMT->parent, GMT_MSG_ERROR, "Unable to parse MOVIE_PROG_INDICATOR_ARG %s for 14 required items [only made %d of 13-14 conversions]\n", movie_item_arg[k][T], kk);
				return NULL;	/* Should never happen */
			}
			/* Because this runs outside main gsave/grestore block the origin is (0,0) */
			if (isupper (kind) && kind != 'A') {	/* Requested text labels so initialize selected font */
				if (font[0] == '+') {	/* User provided a specific font/size which we must honor */
					gmt_getfont (GMT, &font[1], &Tfont);	/* We already parsed the font string in movie.c for correctness */
					fsize = Tfont.size;	/* Must honor the given size */
				}
				else {	/* Pick up the default font written in movie and scale it as needed */
					gmt_getfont (GMT, font, &Tfont);
					fsize = 0.0;	/* So we can scale at will for b-e */
				}
				form = gmt_setfont (GMT, &Tfont);	/* Set the font to be used */
				PSL_setfont (PSL, Tfont.id);
			}
			switch (tolower (kind)) {	/* Set the selected progress indicator a-f */
				case 'a': 	/* Default pie symbol (no label) */
					gmtplot_prog_indicator_A (GMT, plot_x, plot_y, t, width, justify, F1, F2);
					break;
				case 'b': 	/* Growing ring symbol with optional label */
					gmtplot_prog_indicator_B (GMT, plot_x, plot_y, t, width, justify, P1, P2, label, kind, fsize);
					break;
				case 'c': /* Growing circular arrow with optional label */
					gmtplot_prog_indicator_C (GMT, plot_x, plot_y, t, width, justify, P1, P2, label, kind, fsize);
					break;
				case 'd': /* Rounded time line with vertical bar and optional label */
					gmtplot_prog_indicator_D (GMT, plot_x, plot_y, t, width, justify, P1, P2, label, kind, fsize);
					break;
				case 'e': /* Straight line on top of fixed line with optional start/end labels */
					gmtplot_prog_indicator_E (GMT, plot_x, plot_y, t, width, justify, P1, P2, label, kind, fsize);
					break;
				case 'f': /* Moving triangle on basemap time-line */
					gmtplot_prog_indicator_F (GMT, plot_x, plot_y, t, width, justify, P2, F1, label, kind, width, fsize, &Tfont);
					break;
				default:	/* Just for Coverity */
					break;
			}
			gmt_M_str_free (movie_item_arg[k][T]);	/* Free since done with this string */
			gmtplot_reset_PSL (GMT, PSL);	/* Reset anything we may have set in building the completion PS procedure */
		PSL_command (PSL, "U\n");
		}
		PSL_comment (PSL, "End of movie progress indicators\n");
		PSL_command (PSL, "}!\n");
	}
	if (do_paint) {	/* Reset any canvas coloring here */
		GMT->current.map.frame.paint[GMT_Z] = true;
		gmt_M_memcpy (&GMT->current.map.frame.fill[GMT_Z], &fill, 1U, struct GMT_FILL);
	}

	return (PSL);
}

void gmt_plotcanvas (struct GMT_CTRL *GMT) {
	if (GMT->current.map.frame.paint[GMT_Z]) {	/* Paint the inside of the map (xy plane) with specified fill */
		double *x = NULL, *y = NULL;
		uint64_t np;
		bool donut;
		PSL_comment (GMT->PSL, "Fill the canvas %s\n", gmtlib_putfill (GMT, &GMT->current.map.frame.fill[GMT_Z]));
		np = gmt_map_clip_path (GMT, &x, &y, &donut);
		gmt_setfill (GMT, &GMT->current.map.frame.fill[GMT_Z], 0);
		PSL_plotpolygon (GMT->PSL, x, y, (int)((1 + donut) * np));
		gmt_M_free (GMT, x);
		gmt_M_free (GMT, y);
	}
}

int gmt_strip_layer (struct GMTAPI_CTRL *API, int nlayers) {
	/* Remove the last n layers from the current figure */
	char file[PATH_MAX] = {""}, buffer[GMT_LEN64] = {""};
	int k = 0, fig;
	size_t n_alloc = GMT_SMALL_CHUNK;
	FILE *fp = NULL;
	struct GMT_PSLAYER {	/* Used to hold the layer info locally */
		unsigned int id;
		size_t size;
	} *layer = NULL;

	fig = gmt_get_current_figure (API);
	/* Get the name of the corresponding gmt.layers.<fig> file */
	snprintf (file, PATH_MAX, "%s/gmt.layers.%d", API->gwf_dir, fig);
	if (nlayers == GMT_NOTSET) {	/* Reset to nothing, but still remain at current figure */
		if (gmt_remove_file (API->GMT, file))	/* Remove the layers file */
			GMT_Report (API, GMT_MSG_WARNING, "Failed to delete file: %s\n", file);
		snprintf (file, PATH_MAX, "%s/gmt_%d.ps-", API->gwf_dir, fig);
		if (gmt_remove_file (API->GMT, file))	/* Wipe the PS file */
			GMT_Report (API, GMT_MSG_WARNING, "Failed to delete file: %s\n", file);
		return GMT_NOERROR;
	}
	/* See if there is a layers file to read for this figure */
	if (access (file, R_OK)) {
		GMT_Report (API, GMT_MSG_ERROR, "No layers available [no file: %s]\n", file);
		return GMT_FILE_NOT_FOUND;
	}

	/* OK, the file exist, can we read it? */
	if ((fp = fopen (file, "r")) == NULL) {
		GMT_Report (API, GMT_MSG_ERROR, "Could not open file %s\n", file);
		return GMT_ERROR_ON_FOPEN;
	}
	/* Ingest the layer information */
	layer = gmt_M_memory (API->GMT, NULL, n_alloc, struct GMT_PSLAYER);
	while (fgets (buffer, GMT_LEN64, fp)) {
		sscanf (buffer, "%d %" PRIuS, &layer[k].id, &layer[k].size);
		if (++k == (int)n_alloc) {	/* Need more memory */
			n_alloc <<= 1;
			layer = gmt_M_memory (API->GMT, layer, n_alloc, struct GMT_PSLAYER);
		}
	}
	fclose (fp);
	if (nlayers >= k) {
		GMT_Report (API, GMT_MSG_ERROR, "Cannot revert %d layers as only %d layers found!\n", nlayers, k);
		gmt_M_free (API->GMT, layer);
		return GMT_RUNTIME_ERROR;
	}

	snprintf (file, PATH_MAX, "%s/gmt_%d.ps-", API->gwf_dir, fig);
	if (gmt_truncate_file (API, file, layer[k-nlayers-1].size)) {
		GMT_Report (API, GMT_MSG_ERROR, "Could not truncate file %s!\n", file);
		gmt_M_free (API->GMT, layer);
		return GMT_RUNTIME_ERROR;
	}
	/* Finally, rewrite the layers file to skip the reverted layers */
	snprintf (file, PATH_MAX, "%s/gmt.layers.%d", API->gwf_dir, fig);
	if ((fp = fopen (file, "w")) == NULL) {
		GMT_Report (API, GMT_MSG_ERROR, "Could not create new file %s\n", file);
		gmt_M_free (API->GMT, layer);
		return GMT_ERROR_ON_FOPEN;
	}
	nlayers = k - nlayers;	/* Number of surviving layers */
	for (k = 0; k < nlayers; k++)
		fprintf (fp, "%d %" PRIuS "\n", layer[k].id, layer[k].size);
	fclose (fp);
	gmt_M_free (API->GMT, layer);
	return GMT_NOERROR;
}

void gmt_plotend (struct GMT_CTRL *GMT) {
	unsigned int i;
	bool K_active = (GMT->current.setting.run_mode == GMT_MODERN) ? true : GMT->common.K.active;
	struct PSL_CTRL *PSL= GMT->PSL;
	PSL_endlayer (GMT->PSL);
	if (GMT->common.t.active) PSL_command (PSL, "1 1 /Normal PSL_transp\n"); /* Reset transparency to fully opaque, if required */

	if (GMT->common.p.do_z_rotation) {	/* Need a undo the rotation about z of the whole page */
		double x0 = 0.0, y0 = 0.0;	/* Default is to rotate around plot origin */
		if (GMT->current.proj.z_project.view_given) {	/* Rotation is about another z-axis than through the origin */
			x0 = GMT->current.proj.z_project.view_x;
			y0 = GMT->current.proj.z_project.view_y;
		}
		else if (GMT->current.proj.z_project.world_given)	/* Rotation is about another lon/lat pole than the origin */
			gmt_geo_to_xy (GMT, GMT->current.proj.z_project.world_x, GMT->current.proj.z_project.world_y, &x0, &y0);
		GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Transrot: Unrotating plot by %g degrees about (%g, %g)\n", -GMT->common.p.z_rotation, x0, y0);
		PSL_comment (GMT->PSL, "Possibly translate then unrotate rotate whole page\n");
		PSL_setorigin (PSL, x0, y0, -GMT->common.p.z_rotation, PSL_FWD);
		PSL_setorigin (PSL, -x0, -y0, 0.0, PSL_FWD);
	}

	/* Check expected change of clip level to achieved one. Update overall clip level. Check for pending clips. */

	if (abs (GMT->current.ps.nclip) == PSL_ALL_CLIP)	/* Special case where we reset all polygon clip levels */
		GMT->current.ps.clip_level = GMT->current.ps.nclip = PSL->current.nclip = 0;
	else
		GMT->current.ps.clip_level += GMT->current.ps.nclip;

	if (GMT->current.ps.nclip != PSL->current.nclip)
		GMT_Report (GMT->parent, GMT_MSG_INFORMATION,
						"Module was expected to change clip level by %d, but clip level changed by %d\n", GMT->current.ps.nclip, PSL->current.nclip);

	if (!K_active) {
		if (GMT->current.ps.clip_level > 0)
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "%d external clip operations were not terminated!\n", GMT->current.ps.clip_level);
		if (GMT->current.ps.clip_level < 0)
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "%d extra terminations of external clip operations!\n", -GMT->current.ps.clip_level);
		GMT->current.ps.clip_level = 0;	/* Reset to zero, so it will no longer show up in gmt.history */
	}
	for (i = 0; i < 3; i++) gmt_M_str_free (GMT->current.map.frame.axis[i].file_custom);
	PSL_endplot (PSL, !K_active);

	if (GMT->current.setting.run_mode == GMT_MODERN) {	/* Reset file pointer and name */
		struct stat buf;
		char file[PATH_MAX] = {""};
		FILE *fp = NULL;
		if (stat (GMT->current.ps.filename, &buf))
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "Could not determine size of file %s\n", GMT->current.ps.filename);
		else
			GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Current size of half-baked PS file %s = %" PRIuS ".\n", GMT->current.ps.filename, buf.st_size);
		GMT->current.ps.fp = NULL;
		GMT->current.ps.filename[0] = '\0';
		/* Write layer size to gmt.layers.<fig> in case of revert calls */
		snprintf (file, PATH_MAX, "%s/gmt.layers.%d", GMT->parent->gwf_dir, GMT->current.ps.figure);
		if ((fp = fopen (file, "a")) == NULL) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Could not open/create file %s\n", file);
			return;
		}
		fprintf (fp, "%d\t%" PRIuS "\n", GMT->current.ps.layer, (size_t)buf.st_size);
		fclose (fp);
	}
	else if (PSL->internal.memory) {    /* Time to write out buffer regardless of mode */
		struct GMT_POSTSCRIPT *P = gmt_get_postscript (GMT);
		struct GMT_POSTSCRIPT_HIDDEN *PH = gmt_get_P_hidden (P);
		if (GMT->current.ps.title[0]) {
			P->header = gmt_M_memory (GMT, NULL, 1, char *);
			P->header[0] = strdup (GMT->current.ps.title);
			P->n_headers = 1;
		}
		P->data = PSL_getplot (PSL);	/* Get the plot buffer */
		P->n_bytes = PSL->internal.n;   /* Length of plot buffer; note P->n_alloc = 0 since we did not allocate this string here */
		P->mode = PSL->internal.pmode;  /* Mode of plot (GMT_PS_{HEADER,TRAILER,COMPLETE}) */
		PH->alloc_mode = GMT_ALLOC_EXTERNALLY;	/* Since created in PSL */
		if (GMT_Write_Data (GMT->parent, GMT_IS_POSTSCRIPT, GMT_IS_REFERENCE, GMT_IS_TEXT, GMT_WRITE_NORMAL, NULL, GMT->current.ps.memname, P) != GMT_OK) {
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "Unable to write PS structure to file %s!\n", GMT->current.ps.memname);
		}
		/* coverity[leaked_storage] */	/* We can't free P because it was written into a 'memory file' */
	}
	GMT->current.ps.title[0] = '\0';	/* Reset title */
	if (GMT->current.ps.oneliner) GMT->current.ps.active = true;	/* Since we are plotting we reset this here in case other modules have turned it off */

	if (!K_active) GMT->current.plot.gridline_spacing[GMT_X] = GMT->current.plot.gridline_spacing[GMT_Y] = 0.0;	/* Done, if they ever were used */
#if 0
	if (GMT->current.setting.run_mode == GMT_CLASSIC) {	/* Remove any gridline file we may have made in /tmp */
		char file[PATH_MAX] = {""};
		snprintf (file, PATH_MAX, "%s/%s-gmt.gridlines", GMT->parent->tmp_dir, GMT->parent->session_name);
		gmt_remove_file (GMT, file);
	}
#endif
}

void gmt_geo_line (struct GMT_CTRL *GMT, double *lon, double *lat, uint64_t n) {
	/* When geographic lines are plotted, they may cross the boundaries, may need to be clipped,
	 * or may appear again on the other side of the map. This is all taken care of in this
	 * routine.
	 */
	if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, lon, lat, n)) == 0) return;	/* Nothing further to do */
	gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);	/* Separately plot the outline */
}

uint64_t gmt_geo_polarcap_segment (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S, double **lon, double **lat) {
	/* Special treatment for polar caps since they must add in parts of possibly curved periodic boundaries
	 * from the pole up or down to the intersection with the cap perimeter.  We handle this case separately here.
	 * This is in response to issue # 852. P. Wessel */
	struct GMT_DATASEGMENT_HIDDEN *SH = gmt_get_DS_hidden (S);
	uint64_t k0, perim_n, n_new, m, n = S->n_rows, k;
	double start_lon, stop_lon, yc = 0.0, dx, pole_lat = 90.0 * SH->pole;
	double *x_perim = NULL, *y_perim = NULL, *plon = NULL, *plat = NULL;
	static char *pole = "S N";
	int type;
#if 0
	FILE *fp;
	bool aap = gmtplot_at_pole (S->data[GMT_Y], S->n_rows);
#endif
	if (GMT->common.R.oblique) return n;	/* Algorithm assumes meridian boundaries */
	/* We want this code to be used for the misc. global projections but also global cylindrical or linear(if degrees) maps */
	if (!(gmt_M_is_misc(GMT) || (GMT->current.map.is_world  && (gmt_M_is_cylindrical(GMT) || (gmt_M_is_linear(GMT) && gmt_M_is_geographic(GMT,GMT_IN)))))) return n;	/* We are only concerned with the global misc projections here */
	/* Global projection need to handle pole path properly */
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Try to include %c pole in polar cap path\n", pole[SH->pole+1]);
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "West longitude = %g.  East longitude = %g\n", GMT->common.R.wesn[XLO], GMT->common.R.wesn[XHI]);
	type = gmtlib_determine_pole (GMT, S->data[GMT_X], S->data[GMT_Y], n);
	if (abs(type) == 2)	/* The algorithm only works for clockwise polygon so anything CCW we simply reverse... */
		gmtplot_reverse_polygon (GMT, S);
	start_lon = GMT->common.R.wesn[XHI];
	stop_lon  = GMT->common.R.wesn[XLO];

	for (k = 0; k < n; k++) 	/* Make negative longitudes only */
		if (S->data[GMT_X][k] >= 180.0) S->data[GMT_X][k] -= 360.0;
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "First longitude = %g.  Last longitude = %g\n", S->data[GMT_X][0], S->data[GMT_X][n-1]);
	for (k = 1, k0 = 0; k0 == 0 && k < n; k++) {	/* Determine where the perimeter crossing with the west boundary occurs */
		/* Must assume longitudes randomly adds/subtracts 360 and hence carefully check for the crossing */
		//double R0, L0;
		double R, L, D;
		//R0 = GMT->common.R.wesn[XLO]-S->data[GMT_X][k];
		//L0 = GMT->common.R.wesn[XLO]-S->data[GMT_X][k-1];
		gmt_M_set_delta_lon (S->data[GMT_X][k],   GMT->common.R.wesn[XLO], R);	/* Handles the 360 jump cases */
		gmt_M_set_delta_lon (S->data[GMT_X][k-1], GMT->common.R.wesn[XLO], L);	/* Handles the 360 jump cases */
		D = R - L;
		//fprintf (stderr, "k = %d L0 R0 = %g %g  L R = %g %g  D = %g\n", (int)k, L0, R0, L, R, D);
		if (D >= 0.0 && D < 180.0 && L <= 0.0 && R >= 0.0) k0 = k;
	}
	/* Determine the latitude of that crossing */
	if (k0) {	/* Occurred somewhere along the perimeter between points k0 and k0-1 */
		//double x_dist = S->data[GMT_X][k0-1] - GMT->common.R.wesn[XLO];
		double x_dist;
		gmt_M_set_delta_lon (S->data[GMT_X][k0-1], S->data[GMT_X][k0], dx);	/* Handles the 360 jump cases */
		gmt_M_set_delta_lon (S->data[GMT_X][k0-1], GMT->common.R.wesn[XLO], x_dist);	/* Handles the 360 jump cases */
		yc = S->data[GMT_Y][k0-1] + (S->data[GMT_Y][k0] - S->data[GMT_Y][k0-1]) * x_dist / dx;
	}
	else	/* Very first point is at the right longitude */
		yc = S->data[GMT_Y][k0];
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Crossing at %g,%g\n", GMT->common.R.wesn[XLO], yc);
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "k at point closest to lon %g is = %d [n = %d]\n", GMT->common.R.wesn[XLO], (int)k0, (int)n);
	/* Then, add path from pole to start longitude, then copy perimeter path, then add path from stop longitude back to pole */
	/* With cylindrical projections we may not go all the way to the pole, so we adjust pole_lat based on -R: */
	if (pole_lat < GMT->common.R.wesn[YLO]) pole_lat = GMT->common.R.wesn[YLO];
	if (pole_lat >GMT->common.R.wesn[YHI])  pole_lat = GMT->common.R.wesn[YHI];
	/* 1. Calculate the path from yc to pole: */
	perim_n = gmtlib_lonpath (GMT, start_lon, pole_lat, yc, &x_perim, &y_perim);
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Created path from %g/%g to %g/%g [%d points]\n", start_lon, pole_lat, start_lon, yc, perim_n);
	/* 2. Allocate enough space for new polar cap polygon */
	n_new = 2 * perim_n + n;
	plon = gmt_M_memory (GMT, NULL, n_new, double);
	plat = gmt_M_memory (GMT, NULL, n_new, double);
	/* Start off with the path from the pole to the crossing */
	if (perim_n) {
		gmt_M_memcpy (plon, x_perim, perim_n, double);
		gmt_M_memcpy (plat, y_perim, perim_n, double);
	}
	/* Now walk from k0 to the end of polygon, wrapping around if needed */
	m = perim_n;	/* Index of next output point */
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Add perimeter data from k0->n [%d->%d], then 0->k0 [%d]\n", k0, n, k0);
	for (k = k0; k < n; k++, m++) {
		plon[m] = S->data[GMT_X][k];
		plat[m] = S->data[GMT_Y][k];
	}
	for (k = 0; k < k0; k++, m++) {
		plon[m] = S->data[GMT_X][k];
		plat[m] = S->data[GMT_Y][k];
	}
	/* Then add the opposite path to the pole, switching the longitude to stop_lon */
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Add path from %g/%g to %g/%g [%d points]\n", stop_lon, yc, stop_lon, pole_lat, perim_n);
	if (perim_n) {
		for (k = perim_n-1; k > 0; k--, m++) {
			plon[m] = stop_lon;
			plat[m] = y_perim[k];
		}
	}
	/* Finally add the duplicate pole at the end of polygon so it is closed */
	plon[m] = stop_lon;
	plat[m++] = pole_lat;
	if (perim_n) {
		gmt_M_free (GMT, x_perim);	gmt_M_free (GMT, y_perim);	/* No longer needed */
	}
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "New path has %d points, we allocated %d points\n", m, n_new);
#if 0
	fp = fopen ("shit.txt", "w");
	for (k = 0; k < m; k++) {
		fprintf (fp, "%lg\t%lg\n", plon[k], plat[k]);
	}
	fclose (fp);
#endif
	*lon = plon;	*lat = plat;
	return (m);	/* Number of points in new polygon */
}

void gmt_geo_polygons (struct GMT_CTRL *GMT, struct GMT_DATASEGMENT *S) {
	/* Deal with plotting of one or more polygons that may wrap across the map.
	 * Multi-polygons occur if composed of a perimeter and one or more holes.
	 * This is marked by S->next being set to point to the next hole.
	 * Also, if the perimeter is a polar cap we must add a helping line that
	 * connects to the pole but this line should not be drawn.  This is why
	 * we must lay down path twice (first for fill; then for line) since the
	 * two paths are not the same.  If no fill is requested then we just draw lines.
	 */
	struct GMT_DATASEGMENT *S2 = NULL;
	struct GMT_DATASEGMENT_HIDDEN *SH = gmt_get_DS_hidden (S);
	bool add_pole, separate;
	int outline = 0, P_handedness = GMT_NOTSET, H_handedness;
	int geo = gmt_M_is_geographic (GMT, GMT_IN);
	uint64_t used = 0;
	char *type[2] = {"Perimeter", "Polar cap perimeter"};
	char *use[2] = {"fill only", "fill and outline"};
	char comment[GMT_LEN64] = {""};
	struct PSL_CTRL *PSL= GMT->PSL;

	/* CASE 1: NO FILL REQUESTED -- JUST DRAW OUTLINE */

	if (gmt_M_eq (PSL->current.rgb[PSL_IS_FILL][0], -1.0)) {
		if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, S->data[GMT_X], S->data[GMT_Y], S->n_rows)) == 0) return;	/* Nothing further to do */
		PSL_comment (PSL, "Perimeter polygon for outline only\n");
		gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);	/* Separately plot the outline */
		for (S2 = gmt_get_next_S (S); S2; S2 = gmt_get_next_S (S2)) {
			if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, S2->data[GMT_X], S2->data[GMT_Y], S2->n_rows)) == 0) continue;	/* Nothing for this hole */
			PSL_comment (PSL, "Hole polygon for outline only\n");
			gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);	/* Separately plot the outline */
		}
		return;	/* Done with the simple task of drawing lines */
	}

	/* CASE 2: FILL REQUESTED -- WITH OR WITHOUT OUTLINE */

	add_pole = (abs (SH->pole) == 1);	/* true if a polar cap */
	separate = ((add_pole || gmt_get_next_S (S)) && PSL->current.outline);	/* Multi-polygon (or polar cap) fill with outline handled by doing fill and outline separately */
	if (separate) {				/* Do fill and outline separately */
		outline = PSL->current.outline;	/* Keep a copy of what we wanted */
		PSL->current.outline = false;	/* Turns off outline for now (if set) */
		PSL_command (PSL, "O0\n");	/* Temporarily switch off outline in the PS */
	}

	/* Here we must lay down the perimeter and then the holes.  */

	if (PSL->internal.comments) snprintf (comment, GMT_LEN64, "%s polygon for %s\n", type[add_pole], use[PSL->current.outline]);
	used = gmtplot_geo_polygon_segment (GMT, S, add_pole, true, comment);	/* First lay down perimeter */
	for (S2 = gmt_get_next_S (S); S2; S2 = gmt_get_next_S (S2)) {	/* Process all holes [none processed if there aren't any holes] */
		if (P_handedness == GMT_NOTSET)
			P_handedness = gmt_polygon_orientation (GMT, S->data[GMT_X], S->data[GMT_Y], S->n_rows, geo);
		H_handedness = gmt_polygon_orientation (GMT, S2->data[GMT_X], S2->data[GMT_Y], S2->n_rows, geo);

		if (PSL->internal.comments) snprintf (comment, GMT_LEN64, "Hole polygon for %s\n", use[PSL->current.outline]);
		if (H_handedness == P_handedness) {
			uint64_t row_f, row_l, col;
			for (row_f = 0, row_l = S2->n_rows - 1; row_f < S2->n_rows/2; row_f++, row_l--) {
				for (col = 0; col < S2->n_columns; col++) gmt_M_double_swap (S2->data[col][row_f], S2->data[col][row_l]);
			}
		}
		used += gmtplot_geo_polygon_segment (GMT, S2, false, false, comment);	/* Add this hole to the path */
	}
	if (used) {
		PSL_comment (PSL, "Reset FO and fill the path\n");
		PSL_command (PSL, "/FO {fs os}!\nFO\n");	/* Reset FO to its original settings, then force the fill */
	}
	if (separate) {	/* Must draw outline separately */
		if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, S->data[GMT_X], S->data[GMT_Y], S->n_rows)) == 0) return;	/* Nothing further to do */
		PSL_command (PSL, "O1\n");	/* Switch on outline again */
		PSL_comment (PSL, "%s polygon for outline only\n", type[add_pole]);
		PSL->current.outline = outline;	/* Reset outline to what it was originally */
		gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);	/* Separately plot the outline */
		for (S2 = gmt_get_next_S (S); S2; S2 = gmt_get_next_S (S2)) {
			if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, S2->data[GMT_X], S2->data[GMT_Y], S2->n_rows)) == 0) continue;	/* Nothing for this hole */
			PSL_comment (PSL, "Hole polygon for outline only\n");
			gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);	/* Separately plot the outline */
		}
	}
}

void gmt_plot_geo_ellipse (struct GMT_CTRL *GMT, double lon, double lat, double major, double minor, double azimuth) {
	/* gmt_plot_geo_ellipse takes the location, axes (in km), and azimuth of an ellipse
		and draws an approximate ellipse using N-sided polygon and the chosen map projection.
		We determine N by looking at the spacing between successive points for a trial N and
		then scale N up or down to satisfy the minimum point spacing criteria. */

	struct GMT_DATASEGMENT *S = gmt_get_geo_ellipse (GMT, lon, lat, major, minor, azimuth, 0);

	gmt_geo_polygons (GMT, S);

	gmt_free_segment (GMT, &S);
}

GMT_LOCAL void gmtplot_get_far_point (struct GMT_CTRL *GMT, double lon, double lat, double radius, double *P) {
	/* Get a point P that is radius degrees away along the meridian through our point X */
	double plat, plon;
	plat = lat + radius;
	plon = lon;
	if (plat > 90.0) {	/* Went over the pole */
		plon += 180.0;
		plat = 180 - plat;
	}
	gmt_geo_to_cart (GMT, plat, plon, P, true);	/* Vector <radius> degrees away from (lon,lat) along meridian */
}

GMT_LOCAL void gmtplot_geo_spider (struct GMT_CTRL *GMT, double xlon, double xlat, double radius_i, double radius_o, double dr, double az_start, double az_stop, double da, unsigned int wmode) {
	/* gmtplot_geo_spider takes the location, radius_i (in unit), radius_o, dr, and start/stop azimuths of an geo-wedge and da
		and draws an approximate circular spider-web using N-sided polygon.
		mode & 1: Just draw arc(s)
		mode & 2: Just draw radii(s)
	*/

	int n_arc, k, ko, n_seg;
	uint64_t n_new;
	bool windshield = (radius_i > 0.0);	/* Windshield web */
	enum GMT_enum_wedgetype mode = wmode;
	double qlon, qlat, az, rot_start, E[3], P_o[3], P_i[3], Q[3], R[3][3];
	struct GMT_DATASEGMENT *S = NULL;

	radius_i = gmtlib_conv_distance (GMT, radius_i, 'k', 'd');	/* Convert inner radius from km to degrees */
	radius_o = gmtlib_conv_distance (GMT, radius_o, 'k', 'd');	/* Convert outer radius from km to degrees */
	dr = gmtlib_conv_distance (GMT, dr, 'k', 'd');			/* Convert radial increments from km to degrees */

	gmt_geo_to_cart (GMT, xlat, xlon, E, true);	/* Euler rotation pole at point location */

	/* Get a point P_o that is radius_o degrees away along the meridian through our point X */
	gmtplot_get_far_point (GMT, xlon, xlat, radius_o, P_o);
	if (windshield)	/* Get a point P_i that is radius_i degrees away along the meridian through our point X */
		gmtplot_get_far_point (GMT, xlon, xlat, radius_i, P_i);

	PSL_comment (GMT->PSL, "Placing Spider Web\n");
	PSL_command (GMT->PSL, "V PSL_spiderpen\n");	/* Place spider under gsave/grestore and change to spiderpen */
	if (mode == GMT_WEDGE_RADII || mode == GMT_WEDGE_SPIDER) {	/* Need to draw radial start and end points */
		double *azim = NULL;
		if (da > 0.0)	/* Must draw several radial lines */
			n_seg = gmtlib_linear_array (GMT, az_start, az_stop, da, 0.0, &azim);
		else {	/* Two radial lines only */
			azim = gmt_M_memory (GMT, NULL, n_seg = 2, double);
			azim[0] = az_start;	azim[1] = az_stop;
		}
		PSL_comment (GMT->PSL, "Drawing Spider Radials\n");
		for (k = 0; k < n_seg; k++) {	/* For each radial line */
			S = GMT_Alloc_Segment (GMT->parent, GMT_NO_STRINGS, 2, 2, NULL, NULL);	/* Segment with two rows */
			gmt_make_rot_matrix2 (GMT, E, -azim[k], R);	/* Since we have a right-handed rotation but gave azimuths */
			if (windshield) {	/* Start at inner arc point */
				gmt_matrix_vect_mult (GMT, 3U, R, P_i, Q);
				gmt_cart_to_geo (GMT, &S->data[GMT_Y][0], &S->data[GMT_X][0], Q, true);	/* Rotated inner point */
			}
			else {	/* Start at origin */
				S->data[GMT_X][0] = xlon;	S->data[GMT_Y][0] = xlat;
			}
			gmt_matrix_vect_mult (GMT, 3U, R, P_o, Q);
			gmt_cart_to_geo (GMT, &S->data[GMT_Y][1], &S->data[GMT_X][1], Q, true);	/* Rotated outer point */
			if ((n_new = gmt_fix_up_path (GMT, &S->data[GMT_X], &S->data[GMT_Y], S->n_rows, 0.0, 0)) == 0) {
				gmt_free_segment (GMT, &S);
				continue;
			}
			S->n_rows = n_new;
			if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, S->data[GMT_X], S->data[GMT_Y], S->n_rows)) == 0) continue;
			gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);
			gmt_free_segment (GMT, &S);
		}
		gmt_M_free (GMT, azim);
	}
	if (mode == GMT_WEDGE_ARCS || mode == GMT_WEDGE_SPIDER) {	/* Need to draw one or more arcs */
		double d_az, px, py, qx, qy, L, plon_o, plat_o;
		double *r = NULL;
		plat_o = xlat + radius_o;
		plon_o = xlon;
		if (plat_o > 90.0) {	/* Went over the pole */
			plon_o += 180.0;
			plat_o = 180 - plat_o;
		}
		/* Compute distance between P_o and Q and compare to map_line_step to determine azimuthal sampling */
		gmt_make_rot_matrix2 (GMT, E, 1.0, R);		/* Test point 1 degree away */
		gmt_matrix_vect_mult (GMT, 3U, R, P_o, Q);	/* Get Q = R * P_o */
		gmt_cart_to_geo (GMT, &qlat, &qlon, Q, true);	/* Coordinates of rotated point */
		gmt_geo_to_xy (GMT, plon_o, plat_o, &px, &py);	/* P_o projected on map */
		gmt_geo_to_xy (GMT, qlon, qlat, &qx, &qy);	/* Q projected on map */
		L = hypot (px - qx, py - qy);			/* Distance in inches along arc for 1 degree of azimuth change */
		/* Estimate how many intermediate points by dividing the estinate arc length by line_step */
		n_arc = MAX (2, irint (fabs (az_stop - az_start) * L /GMT->current.setting.map_line_step));
		d_az = (az_stop - az_start) / (n_arc - 1);	/* Azimuthal sampling rate */
		GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Arcs will be approximated by %d-point lines\n", n_arc);
		PSL_comment (GMT->PSL, "Drawing Spider Arcs\n");

		if (dr > 0.0)	/* Must draw several arcs */
			n_seg = gmtlib_linear_array (GMT, (radius_i > 0.0) ? radius_i : dr, radius_o, dr, 0.0, &r);
		else {	/* One or two arcs only */
			if (radius_i > 0.0) {
				r = gmt_M_memory (GMT, NULL, n_seg = 2, double);
				r[0] = radius_i;	r[1] = radius_o;
			}
			else {
				r = gmt_M_memory (GMT, NULL, n_seg = 1, double);
				r[0] = radius_o;
			}
		}

		d_az = -d_az;	/* Since we have a right-handed rotation but gave azimuths */
		for (k = 0; k < n_seg; k++) {
			S = GMT_Alloc_Segment (GMT->parent, GMT_NO_STRINGS, n_arc, 2, NULL, NULL);
			gmtplot_get_far_point (GMT, xlon, xlat, r[k], P_o);
			rot_start = -az_start;	/* Since we have a right-handed rotation but gave azimuths */
			for (ko = 0; ko < n_arc; ko++) {
				az = rot_start + ko * d_az;
				gmt_make_rot_matrix2 (GMT, E, az, R);
				gmt_matrix_vect_mult (GMT, 3U, R, P_o, Q);
				gmt_cart_to_geo (GMT, &S->data[GMT_Y][ko], &S->data[GMT_X][ko], Q, true);	/* Rotated outer point */
			}
			if ((n_new = gmt_fix_up_path (GMT, &S->data[GMT_X], &S->data[GMT_Y], S->n_rows, 0.0, 0)) == 0) {
				gmt_free_segment (GMT, &S);
				continue;
			}
			S->n_rows = n_new;
			if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, S->data[GMT_X], S->data[GMT_Y], S->n_rows)) == 0) continue;
			gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);
			gmt_free_segment (GMT, &S);
		}
		gmt_M_free (GMT, r);
	}
	PSL_command (GMT->PSL, "U\n");	/* End block where spiderpen is active; this resets previous pen (for outline) */
}

GMT_LOCAL void gmtplot_geo_wedge_fill (struct GMT_CTRL *GMT, double xlon, double xlat, double radius_i, double radius_o, double az_start, double az_stop, bool fill, bool outline) {
	/* gmt_geo_wedge takes the location, radius_i (in unit), radius_o (in unit), and start/stop azimuths of an geo-wedge
		and draws an approximate circluar wedge using N-sided polygon.
		mode = 3: Draw the entire closed wedge.
		mode = 1: Just draw arc [no fill possible]
		mode = 2: Just draw jaw [no fill possible]
	*/

	int n_arc, k, ki, ko, n_path;
	uint64_t n_new;
	bool windshield = (radius_i > 0.0);	/* Windshield plot */
	double d_az, px, py, qx, qy, L, plon_o, plat_o, qlon, qlat, az, rot_start, E[3], P_o[3], P_i[3], Q[3], R[3][3];
	struct GMT_DATASEGMENT *S = NULL;

	radius_i = gmtlib_conv_distance (GMT, radius_i, 'k', 'd');	/* Convert inner radius from km to degrees */
	radius_o = gmtlib_conv_distance (GMT, radius_o, 'k', 'd');	/* Convert outer radius from km to degrees */

	gmt_geo_to_cart (GMT, xlat, xlon, E, true);	/* Euler rotation pole at point location */

	/* Get a point P_o that is radius_o degrees away along the meridian through our point X */
	gmtplot_get_far_point (GMT, xlon, xlat, radius_o, P_o);
	if (windshield)	/* Get a point P_i that is radius_i degrees away along the meridian through our point X */
		gmtplot_get_far_point (GMT, xlon, xlat, radius_i, P_i);
	/* Compute distance between P_o and Q and compare to map_line_step to determine azimuthal sampling */
	plat_o = xlat + radius_o;
	plon_o = xlon;
	if (plat_o > 90.0) {	/* Went over the pole */
		plon_o += 180.0;
		plat_o = 180 - plat_o;
	}
	gmt_make_rot_matrix2 (GMT, E, 1.0, R);		/* Test point 1 degree away */
	gmt_matrix_vect_mult (GMT, 3U, R, P_o, Q);	/* Get Q = R * P_o */
	gmt_cart_to_geo (GMT, &qlat, &qlon, Q, true);	/* Coordinates of rotated point */
	gmt_geo_to_xy (GMT, plon_o, plat_o, &px, &py);	/* P_o projected on map */
	gmt_geo_to_xy (GMT, qlon, qlat, &qx, &qy);	/* Q projected on map */
	L = hypot (px - qx, py - qy);			/* Distance in inches along arc for 1 degree of azimuth change */
	/* Estimate how many intermediate points by dividing the estinate arc length by line_step */
	n_arc = MAX (2, irint (fabs (az_stop - az_start) * L /GMT->current.setting.map_line_step));
	d_az = (az_stop - az_start) / (n_arc - 1);	/* Azimuthal sampling rate */
	n_path = n_arc;			/* Total number of points along outer arc */
	if (windshield)	/* Need points along inner arc */
		n_path *= 2;
	else	/* Just place apex */
		n_path++;
	n_path++;	/* Closed the polygon */
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Wedge will be approximated by %d-sided polygon\n", n_path);
	S = GMT_Alloc_Segment (GMT->parent, GMT_NO_STRINGS, n_path, 2, NULL, NULL);	/* Polygon array */
	rot_start = -az_start;	/* Since we have a right-handed rotation but gave azimuths */
	d_az = -d_az;		/* Same reason */
	for (k = ko = 0; k < n_arc; k++, ko++) {
		az = rot_start + k * d_az;
		gmt_make_rot_matrix2 (GMT, E, az, R);
		gmt_matrix_vect_mult (GMT, 3U, R, P_o, Q);
		gmt_cart_to_geo (GMT, &S->data[GMT_Y][ko], &S->data[GMT_X][ko], Q, true);	/* Rotated outer point */
		if (windshield) {	/* Get the corresponding point on the inner arc but offset its index */
			gmt_matrix_vect_mult (GMT, 3U, R, P_i, Q);
			ki = n_path - 2 - ko;	/* Point index for this angle on the return arc */
			gmt_cart_to_geo (GMT, &S->data[GMT_Y][ki], &S->data[GMT_X][ki], Q, true);	/* Rotated inner point */
		}
	}
	if (!windshield)	/* Add apex of wedge */
		S->data[GMT_X][ko] = xlon, S->data[GMT_Y][ko++] = xlat;
	/* Close polygon */
	S->data[GMT_X][n_path-1] = S->data[GMT_X][0];	S->data[GMT_Y][n_path-1] = S->data[GMT_Y][0];
	/* Resample along great circles */
	if ((n_new = gmt_fix_up_path (GMT, &S->data[GMT_X], &S->data[GMT_Y], S->n_rows, 0.0, 0)) == 0) {
		gmt_free_segment (GMT, &S);
		return;
	}
	S->n_rows = n_new;
	gmt_set_seg_minmax (GMT, GMT_IS_POLY, 2, S);	/* Update min/max of x/y only */

	if (fill && outline) {
		PSL_comment (GMT->PSL, "Drawing Wedge fill and outline\n");
		gmt_geo_polygons (GMT, S);
	}
	else if (fill) {
		PSL_comment (GMT->PSL, "Drawing Wedge fill only\n");
		PSL_command (GMT->PSL, "V O0\n");	/* Temporarily disable outlines */
		gmt_geo_polygons (GMT, S);
		PSL_command (GMT->PSL, "U\n");		/* Undo */
	}
	else if (outline) {
		PSL_comment (GMT->PSL, "Drawing Wedge outline only\n");
		if ((GMT->current.plot.n = gmt_geo_to_xy_line (GMT, S->data[GMT_X], S->data[GMT_Y], S->n_rows)))
			gmt_plot_line (GMT, GMT->current.plot.x, GMT->current.plot.y, GMT->current.plot.pen, GMT->current.plot.n, PSL_LINEAR);
	}

	gmt_free_segment (GMT, &S);
}

void gmt_geo_wedge (struct GMT_CTRL *GMT, double xlon, double xlat, double radius_i, double radius_o, double dr, double az_start, double az_stop, double da, unsigned int wmode, bool fill, bool outline) {
	/* gmt_geo_wedge takes the location, radius_i (in km), radius_o (in km), dr (in km) and start/stop/da azimuths of an geo-wedge
		and draws an approximate circular wedge or windshield using N-sided polygon.
		If fill is true then we paint this polygon.
		If mode is not GMT_WEDGE_NORMAL then we draw the spider web on top.
		The form of that depends on dr > 0 and da > 0, otherwise we just do the outline.
	*/

	enum GMT_enum_wedgetype mode = wmode;

	if (mode == GMT_WEDGE_NORMAL) {	/* Just a regular filled/outlined geowedge */
		gmtplot_geo_wedge_fill (GMT, xlon, xlat, radius_i, radius_o, az_start, az_stop, fill, outline);
		return;
	}

	/* Here we need to lay down fill first (if active), then spider web, then wedge outline (if active) */

	if (fill) gmtplot_geo_wedge_fill (GMT, xlon, xlat, radius_i, radius_o, az_start, az_stop, true, false);	/* Just fill (if active) */
	gmtplot_geo_spider (GMT, xlon, xlat, radius_i, radius_o, dr, az_start, az_stop, da, wmode);			/* Spiderweb */
	if (outline) gmtplot_geo_wedge_fill (GMT, xlon, xlat, radius_i, radius_o, az_start, az_stop, false, true);	/* Just outline (if active) */
}

unsigned int gmt_geo_vector (struct GMT_CTRL *GMT, double lon0, double lat0, double azimuth, double length, struct GMT_PEN *pen, struct GMT_SYMBOL *S) {
	/* gmt_geo_vector takes the location lon0, lat0, azimuth of the vector at that point, and the
		length (in km), and and draws the vector using the chosen map projection.  If arrow heads
		have been requested we compute an arc length in degrees that is equivalent to the chosen
		symbol size.  With arrow heads we also shorten the vector arc so that unfilled vector heads
		are possible. If a small-circle vector is chosen then azimuth, length may be opening angles
		1 and 2 if PSL_VEC_ANGLES is set as well. */
	unsigned int warn;
	if ((S->v.status & PSL_VEC_SCALE) == 0) {	/* Must determine the best inch to degree scale for this map */
		if (gmt_M_is_perspective (GMT)) {
			double clon, clat;
			gmt_xy_to_geo (GMT, &clon, &clat, GMT->current.map.half_width, GMT->current.map.half_height);	/* Geographic coordinates of middle map point */
			S->v.scale = gmtplot_inch_to_degree_scale (GMT, clon, clat, azimuth);
			S->v.status |= PSL_VEC_SCALE;
			GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Vector stem scale is %g degrees/inch at (%g, %g) for az = %g\n", S->v.scale, clon, clat, azimuth);
		}
		else {	/* Set scale each time locally */
			S->v.scale = gmtplot_inch_to_degree_scale (GMT, lon0, lat0, azimuth);
			GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Vector stem scale is %g degrees/inch at (%g, %g) for az = %g\n", S->v.scale, lon0, lat0, azimuth);
		}
	}

	if (S->v.status & PSL_VEC_POLE)
		warn = gmtplot_geo_vector_smallcircle (GMT, lon0, lat0, azimuth, length, pen, S);
	else
		warn = gmtplot_geo_vector_greatcircle (GMT, lon0, lat0, azimuth, length, pen, S);
	return (warn);
}

void gmt_geo_rectangle (struct GMT_CTRL *GMT, double lon, double lat, double width, double height, double azimuth) {
	/* gmt_geo_rectangle takes the location, axes (in km), and azimuth of a rectangle
		and draws the rectangle using the chosen map projection */

	int jump;
	double sin_azimuth, cos_azimuth, sinp, cosp, x, y, x_prime, y_prime, rho, c, dim[3];
	double sin_c, cos_c, center, lon_w, lat_w, lon_h, lat_h, xp, yp, xw, yw, xh, yh;
	struct PSL_CTRL *PSL= GMT->PSL;

	azimuth = gmt_azim_to_angle (GMT, lon, lat, 0.1, azimuth);
	gmt_geo_to_xy (GMT, lon, lat, &xp, &yp);		/* Center of rectangle */

	width *= 500.0, height *= 500.0;	/* Convert to meters and get half the size */
	dim[0] = azimuth;
	sincosd (azimuth, &sin_azimuth, &cos_azimuth);
	sincosd (lat, &sinp, &cosp);		/* Set up azimuthal equidistant projection */

	center = (GMT->current.proj.central_meridian < GMT->common.R.wesn[XLO] || GMT->current.proj.central_meridian > GMT->common.R.wesn[XHI]) ? 0.5 * (GMT->common.R.wesn[XLO] + GMT->common.R.wesn[XHI]) : GMT->current.proj.central_meridian;

	/* Get first point width away from center */
	sincos (0.0, &y, &x);
	x *= width;
	y *= height;
	/* Get rotated coordinates in m */
	x_prime = x * cos_azimuth - y * sin_azimuth;
	y_prime = x * sin_azimuth + y * cos_azimuth;
	/* Convert m back to lon lat */
	rho = hypot (x_prime, y_prime);
	c = rho / GMT->current.proj.EQ_RAD;
	sincos (c, &sin_c, &cos_c);
	lat_w = d_asind (cos_c * sinp + (y_prime * sin_c * cosp / rho));
	if ((lat - 90.0) > -GMT_CONV8_LIMIT)	/* origin in Northern hemisphere */
		lon_w = lon + d_atan2d (x_prime, -y_prime);
	else if ((lat + 90.0) < GMT_CONV8_LIMIT)	/* origin in Southern hemisphere */
		lon_w = lon + d_atan2d (x_prime, y_prime);
	else
		lon_w = lon + d_atan2d (x_prime * sin_c, (rho * cosp * cos_c - y_prime * sinp * sin_c));
	while ((lon_w - center) < -180.0) lon_w += 360.0;
	while ((lon_w - center) > +180.0) lon_w -= 360.0;
	gmt_geo_to_xy (GMT, lon_w, lat_w, &xw, &yw);	/* Get projected x,y coordinates */
	if ((jump = (*GMT->current.map.jump) (GMT, xp, yp, xw, yw)))	/* Adjust for map jumps */
		xw += jump * 2.0 * gmt_half_map_width (GMT, yp);
	dim[1] = 2.0 * hypot (xp - xw, yp - yw);	/* Estimate of rectangle width in plot units (inch) */
	/* Get 2nd point height away from center */
	sincos (M_PI_2, &y, &x);
	x *= width;
	y *= height;
	/* Get rotated coordinates in m */
	x_prime = x * cos_azimuth - y * sin_azimuth;
	y_prime = x * sin_azimuth + y * cos_azimuth;
	/* Convert m back to lon lat */
	rho = hypot (x_prime, y_prime);
	c = rho / GMT->current.proj.EQ_RAD;
	sincos (c, &sin_c, &cos_c);
	lat_h = d_asind (cos_c * sinp + (y_prime * sin_c * cosp / rho));
	if ((lat - 90.0) > -GMT_CONV8_LIMIT)	/* origin in Northern hemisphere */
		lon_h = lon + d_atan2d (x_prime, -y_prime);
	else if ((lat + 90.0) < GMT_CONV8_LIMIT)	/* origin in Southern hemisphere */
		lon_h = lon + d_atan2d (x_prime, y_prime);
	else
		lon_h = lon + d_atan2d (x_prime * sin_c, (rho * cosp * cos_c - y_prime * sinp * sin_c));
	while ((lon_h - center) < -180.0) lon_h += 360.0;
	while ((lon_h - center) > +180.0) lon_h -= 360.0;
	gmt_geo_to_xy (GMT, lon_h, lat_h, &xh, &yh);
	if ((jump = (*GMT->current.map.jump) (GMT, xp, yp, xh, yh)))	/* Adjust for map jumps */
		xh += jump * 2.0 * gmt_half_map_width (GMT, yp);
	dim[2] = 2.0 * hypot (xp - xh, yp - yh);	/* Estimate of rectangle width in plot units (inch) */
	PSL_plotsymbol (PSL, xp, yp, dim, PSL_ROTRECT);
}

#define M_PI_3          1.047197551196598	/* 60 degrees in radians */

void gmt_draw_front (struct GMT_CTRL *GMT, double x[], double y[], uint64_t n, struct GMT_FRONTLINE *f) {
	int ngap, tmp_join = 0, tmp_limit = 0;
	bool skip;
	uint64_t i;
	double *s = NULL, xx[6], yy[6], dist = 0.0, w, frac, dx, dy, angle, dir1, dir2;
	double gap, x0, y0, xp, yp, len2, len3, len4, cosa, sina, sa, ca, offx, offy, dim[PSL_MAX_DIMS];
	double a1, a2, sina1, sina2, cosa1, cosa2;
	struct PSL_CTRL *PSL= GMT->PSL;

	if (n < 2) return;

	gmt_M_memset (dim, PSL_MAX_DIMS, double);
	s = gmt_M_memory (GMT, NULL, n, double);
	for (i = 1, s[0] = 0.0; i < n; i++) {
		/* Watch out for longitude wraps */
		dx = x[i] - x[i-1];
		w = gmt_half_map_width (GMT, y[i]);
		if (GMT->current.map.is_world && fabs (dx) > w) dx = copysign (2.0 * w - fabs (dx), -dx);
		s[i] = s[i-1] + hypot (dx, y[i] - y[i-1]);
	}

	if (f->f_gap > 0.0) {	/* Gave positive interval */
		if (f->f_exact) {	/* Use gap exactly as given, not worry about ending line with tick */
			ngap = floor ((s[n-1] * (1.0 + GMT_CONV6_LIMIT)) / f->f_gap);	/* Allow 1 ppm noise since we use floor */
			gap = f->f_gap;	/* As given */
		}
		else {	/* Adjust so we start and end with a tick on each line */
			ngap = irint (s[n-1] / f->f_gap);
			gap = s[n-1] / ngap;	/* Adjust gap to fit line length */
			ngap++;
		}
		dist = f->f_off;	/* Start off at the offset distance [0] */
		GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Given gap: %g Adjusted gap: %g Number of front gaps: %d\n", f->f_gap, gap, ngap);
	}
	else {	/* Gave negative interval which means the # of ticks required */
		ngap = irint (fabs (f->f_gap));
		if (ngap == 0) {	/* Cannot happen but might as well leave the test in case of snafus */
			GMT_Report (GMT->parent, GMT_MSG_WARNING, "Number of front ticks reset from 0 to 1 (check your arguments)\n");
			ngap = 1;
		}
		if (ngap == 1)	/* Single centered tick */
			dist = 0.5 * s[n-1], gap = s[n-1];
		else		/* Equidistantly spaced tick starting at 1st point and ending at last */
			gap = s[n-1] / (ngap - 1);
		GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Given number of front gaps: %d Computed gap: %g\n", ngap, gap);
	}

	PSL_command (GMT->PSL, "V\n");	/* In case we change the graphic state regarding pens */

	len2 = 0.5 * f->f_len;
	len4 = 0.25 * f->f_len;
	len3 = 0.866025404 * f->f_len;
	if (f->f_sense == GMT_FRONT_CENTERED) len3 = len2;
	if (f->f_symbol) {	/* Temporarily use miter to get sharp points at slip vectors */
		tmp_join = PSL->internal.line_join;	PSL_setlinejoin (PSL, 0);
		tmp_limit = PSL->internal.miter_limit;	PSL_setmiterlimit (PSL, 0);
	}
	if (f->f_pen == 1) gmt_setpen (GMT, &f->pen);	/* Set alternate symbol pen */
	i = 0;
	while (i < n) {
		while ((s[i] - dist) > -GMT_CONV4_LIMIT) {	/* Time for tick */
			if (i > 0) {
				dx = x[i] - x[i-1];
				dy = y[i] - y[i-1];
			}
			else {
				dx = x[1] - x[0];
				dy = y[1] - y[0];
			}
			if (fabs (dist - s[i]) < GMT_CONV4_LIMIT) {
				x0 = x[i];
				y0 = y[i];
			}
			else {
				frac = (s[i] - dist) / (s[i] - s[i-1]);
				x0 = x[i] - dx * frac;
				y0 = y[i] - dy * frac;
			}
			angle = d_atan2 (dy, dx);
			skip = (GMT->current.map.is_world && fabs (dx) > gmt_half_map_width (GMT, y[i]));	/* Don't do ticks on jumps */
			if (skip) {
				dist += gap;	i++;
				continue;
			}

			switch (f->f_symbol) {
				case GMT_FRONT_TRIANGLE:	/* Triangle */
					switch (f->f_sense) {
						case GMT_FRONT_CENTERED:
							sincos (angle, &sina, &cosa);
							xx[0] = x0 + len2 * cosa;
							yy[0] = y0 + len2 * sina;
							xx[1] = x0 - len3 * sina;
							yy[1] = y0 + len3 * cosa;
							xx[2] = x0 - len2 * cosa;
							yy[2] = y0 - len2 * sina;
							xx[3] = x0 + len3 * sina;
							yy[3] = y0 - len3 * cosa;
							PSL_plotpolygon (PSL, xx, yy, 4);
							break;
						case GMT_FRONT_RIGHT:
							angle += M_PI;
							/* Intentionally fall through - after changing the angle */
						case GMT_FRONT_LEFT:
							sincos (angle, &sina, &cosa);
							xx[0] = x0 + len2 * cosa;
							yy[0] = y0 + len2 * sina;
							xx[1] = x0 - len3 * sina;
							yy[1] = y0 + len3 * cosa;
							xx[2] = x0 - len2 * cosa;
							yy[2] = y0 - len2 * sina;
							PSL_plotpolygon (PSL, xx, yy, 3);
							break;
					}
					break;

				case GMT_FRONT_ITRIANGLE:    /* Inverted triangle */
					a1 = angle + M_PI_3;	a2 = a1 + M_PI_3;
					sincos (a1, &sina1, &cosa1);
					sincos (a2, &sina2, &cosa2);
					switch (f->f_sense) {
						case GMT_FRONT_CENTERED:
							xx[0] = x0 + f->f_len * cosa1;   /* UR */
							yy[0] = y0 + f->f_len * sina1;
							xx[1] = x0 + f->f_len * cosa2;   /* UL */
							yy[1] = y0 + f->f_len * sina2;
							xx[2] = x0;   yy[2] = y0; 		 /* MC */
							xx[3] = x0 - f->f_len * cosa1;   /* LR */
							yy[3] = y0 - f->f_len * sina1;
							xx[4] = x0 - f->f_len * cosa2;   /* LL */
							yy[4] = y0 - f->f_len * sina2;
							xx[5] = x0;   yy[5] = y0; 		 /* MC */
							PSL_plotpolygon (PSL, xx, yy, 6);
							break;
						case GMT_FRONT_RIGHT:
							sina1 = -sina1;	cosa1 = -cosa1;	sina2 = -sina2;	cosa2 = -cosa2;
							/* Intentionally fall through - after changing the angle */
						case GMT_FRONT_LEFT:
							xx[0] = x0 + f->f_len * cosa1;   /* UR */
							yy[0] = y0 + f->f_len * sina1;
							xx[1] = x0 + f->f_len * cosa2;   /* UL */
							yy[1] = y0 + f->f_len * sina2;
							xx[2] = x0;   yy[2] = y0; 		 /* MC */
							PSL_plotpolygon (PSL, xx, yy, 3);
							break;
					}
					break;

				case GMT_FRONT_CIRCLE:	/* Circles */
					switch (f->f_sense) {
						case GMT_FRONT_CENTERED:
							PSL_plotsymbol (PSL, x0, y0, &(f->f_len), PSL_CIRCLE);
							break;
						case GMT_FRONT_RIGHT:
							angle += M_PI;
							/* Intentionally fall through - after changing the angle */
						case GMT_FRONT_LEFT:
							dir1 = R2D * angle;
							dir2 = dir1 + 180.0;
							if (dir1 > dir2) dir1 -= 360.0;
							dim[0] = len2, dim[1] = dir1, dim[2] = dir2;	dim[7] = 3;
							PSL_plotsymbol (PSL, x0, y0, dim, PSL_WEDGE);
							break;
					}
					break;

				case GMT_FRONT_BOX:	/* Squares */
					switch (f->f_sense) {
						case GMT_FRONT_CENTERED:	/* Full square centered on line */
							sincos (angle, &sina, &cosa);
							xx[0] = x0 + len2 * (cosa + sina);	/* LR */
							yy[0] = y0 + len2 * (sina - cosa);
							xx[1] = x0 + len2 * (cosa - sina);	/* UR */
							yy[1] = y0 + len2 * (sina + cosa);
							xx[2] = x0 + len2 * (-cosa - sina);	/* UL */
							yy[2] = y0 + len2 * (-sina + cosa);
							xx[3] = x0 + len2 * (-cosa + sina);	/* LL */
							yy[3] = y0 + len2 * (-sina - cosa);
							break;
						case GMT_FRONT_RIGHT:
							angle += M_PI;
							/* Intentionally fall through - after changing the angle */
						case GMT_FRONT_LEFT:
							/* Half square on the chosen side */
							sincos (angle, &sina, &cosa);
							xx[0] = x0 + len2 * (cosa);	/* LR */
							yy[0] = y0 + len2 * (sina);
							xx[1] = x0 + len2 * (cosa - sina);	/* UR */
							yy[1] = y0 + len2 * (sina + cosa);
							xx[2] = x0 + len2 * (-cosa - sina);	/* UL */
							yy[2] = y0 + len2 * (-sina + cosa);
							xx[3] = x0 + len2 * (-cosa);	/* LL */
							yy[3] = y0 + len2 * (-sina);
							break;
					}
					PSL_plotpolygon (PSL, xx, yy, 4);
					break;

				case GMT_FRONT_SLIP: /* draw strike-slip arrows */
					sincos (angle, &sina, &cosa);
					offx = GMT->current.setting.map_annot_offset[GMT_PRIMARY] * sina; /* get offsets from front line */
					offy = GMT->current.setting.map_annot_offset[GMT_PRIMARY] * cosa;
					/* sense == GMT_FRONT_LEFT == left-lateral, R_RIGHT = right lateral */
					/* arrow "above" line */
					sincos (angle + (f->f_sense * f->f_angle * D2R), &sa, &ca);
					xp = x0 + f->f_sense * offx;
					yp = y0 - f->f_sense * offy;
					xx[0] = xp - len2 * cosa;
					yy[0] = yp - len2 * sina;
					xx[1] = xp + len2 * cosa;
					yy[1] = yp + len2 * sina;
					xx[2] = xx[1] - len4 * ca;
					yy[2] = yy[1] - len4 * sa;
					PSL_plotline (PSL, xx, yy, 3, PSL_MOVE|PSL_STROKE);

					/* arrow "below" line */
					sincos (angle - (f->f_sense * (180.0 - f->f_angle) * D2R), &sa, &ca);
					xp = x0 - f->f_sense * offx;
					yp = y0 + f->f_sense * offy;
					xx[0] = xp + len2 * cosa;
					yy[0] = yp + len2 * sina;
					xx[1] = xp - len2 * cosa;
					yy[1] = yp - len2 * sina;
					xx[2] = xx[1] - len4 * ca;
					yy[2] = yy[1] - len4 * sa;
					PSL_plotline (PSL, xx, yy, 3, PSL_MOVE|PSL_STROKE);
					break;

				case GMT_FRONT_SLIPC: /* Draw curved strike-slip arrows a la USGS */
					PSL_command (PSL, "V ");	/* Place symbol under gsave/grestore since we will translate/rotate */
					PSL_setorigin (PSL, x0, y0, R2D * angle, PSL_FWD);
					offy = GMT->current.setting.map_annot_offset[GMT_PRIMARY]; /* get offset from front line */
					/* sense == GMT_FRONT_LEFT == left-lateral, R_RIGHT = right lateral */
					/* arrow "above" line */
					xx[0] = f->f_sense* len2;	xx[1] = -f->f_sense*len2;
					yy[0] = yy[1] = offy;
					PSL_plotline (PSL, xx, yy, 2, PSL_MOVE);
					PSL_plotarc (PSL, xx[1], offy + 0.4 * f->f_len, 0.4 * f->f_len, 270.0, 270.0+f->f_sense*f->f_angle, PSL_STROKE);

					/* arrow "below" line */
					xx[0] = -f->f_sense*len2;	xx[1] = f->f_sense*len2;
					yy[0] = yy[1] = -offy;
					PSL_plotline (PSL, xx, yy, 2, PSL_MOVE);
					PSL_plotarc (PSL, xx[1], -offy - 0.4 * f->f_len, 0.4 * f->f_len, 90.0, 90.0+f->f_sense*f->f_angle, PSL_STROKE);
					PSL_command (PSL, "U\n");
					break;

				case GMT_FRONT_FAULT:	/* Normal fault ticks */
					xx[0] = xx[1] = x0, yy[0] = yy[1] = y0;
					if (f->f_sense == GMT_FRONT_CENTERED) {
						angle -= M_PI_2;
						sincos (angle, &sina, &cosa);
						xx[0] += len2 * cosa;
						yy[0] += len2 * sina;
						xx[1] -= len2 * cosa;
						yy[1] -= len2 * sina;
					}
					else {
						angle += (f->f_sense * M_PI_2);
						sincos (angle, &sina, &cosa);
						xx[1] += len2 * cosa;
						yy[1] += len2 * sina;
					}
					PSL_plotline (PSL, xx, yy, 2, PSL_MOVE|PSL_STROKE);
					break;
			}
			dist += gap;
		}
		i++;
	}
	gmt_M_free (GMT, s);
	if (f->f_symbol) {	/* Switch line join style back */
		PSL_setlinejoin (PSL, tmp_join);
		PSL_setmiterlimit (PSL, tmp_limit);
	}
	PSL_command (GMT->PSL, "U\n");
}

void gmt_plane_perspective (struct GMT_CTRL *GMT, int plane, double level) {
	/* This routine write the PostScript code to change any following matter printed in the plane
	 * of the paper into a perspective view of that plane based on the GMT->current.proj.z_project
	 * parameters (azimuth and elevation).
	 * The plane is portrayed as a plane of constant X, Y, or Z.
	 * Input arguments:
	 * GMT	: The GMT struct
	 * PSL	: The PSL struct
	 * plane: The perspective plane if a constant X, Y, or Z (GMT_X = 0, GMT_Y = 1, GMT_Z = 2)
	 *        To indicate that the z-level is not in projected but "world" coordinates, add GMT_ZW = 3
	 *        To reset to normal printing, use -1.
	 * level: Level of X, Y, or Z in projected coordinates (inch).
	 */
	double a, b, c, d, e, f;
	struct PSL_CTRL *PSL= GMT->PSL;

	/* Only do this in 3D mode */
	if (!GMT->current.proj.three_D) return;

	/* Only do this at top module */
	if (GMT->hidden.func_level > GMT_TOP_MODULE) return;

	/* Nothing changed since last call, hence ignore */
	if (plane == GMT->current.proj.z_project.plane && gmt_M_eq(level,GMT->current.proj.z_project.level)) return;

	/* Store value of level (store plane at end) */
	GMT->current.proj.z_project.level = level;

	/* Concat contains the proper derivatives of these functions:
	x_out = - x * GMT->current.proj.z_project.cos_az + y * GMT->current.proj.z_project.sin_az + GMT->current.proj.z_project.x_off;
	y_out = - (x * GMT->current.proj.z_project.sin_az + y * GMT->current.proj.z_project.cos_az) *
		GMT->current.proj.z_project.sin_el + z * GMT->current.proj.z_project.cos_el + GMT->current.proj.z_project.y_off;
	*/

	a = b = c = d = e = f = 0.0;
	if (plane < 0)			/* Reset to original matrix */
		PSL_command (PSL, "PSL_GPP setmatrix\n");
	else {	/* New perspective plane: compute all derivatives and use full matrix */
		if (plane >= GMT_ZW) level = gmt_z_to_zz (GMT, level);	/* First convert world z coordinate to projected z coordinate */
		switch (plane % 3) {
			case GMT_X:	/* Constant x, Convert y,z to x',y' */
				a = GMT->current.proj.z_project.sin_az;
				b = -GMT->current.proj.z_project.cos_az * GMT->current.proj.z_project.sin_el;
				c = 0.0;
				d = GMT->current.proj.z_project.cos_el;
				e = GMT->current.proj.z_project.x_off - level * GMT->current.proj.z_project.cos_az;
				f = GMT->current.proj.z_project.y_off - level * GMT->current.proj.z_project.sin_az * GMT->current.proj.z_project.sin_el;
				break;
			case GMT_Y:	/* Constant y. Convert x,z to x',y' */
				a = -GMT->current.proj.z_project.cos_az;
				b = -GMT->current.proj.z_project.sin_az * GMT->current.proj.z_project.sin_el;
				c = 0.0;
				d = GMT->current.proj.z_project.cos_el;
				e = GMT->current.proj.z_project.x_off + level * GMT->current.proj.z_project.sin_az;
				f = GMT->current.proj.z_project.y_off - level * GMT->current.proj.z_project.cos_az * GMT->current.proj.z_project.sin_el;
				break;
			case GMT_Z:	/* Constant z. Convert x,y to x',y' */
				a = -GMT->current.proj.z_project.cos_az;
				b = -GMT->current.proj.z_project.sin_az * GMT->current.proj.z_project.sin_el;
				c = GMT->current.proj.z_project.sin_az;
				d = -GMT->current.proj.z_project.cos_az * GMT->current.proj.z_project.sin_el;
				e = GMT->current.proj.z_project.x_off;
				f = GMT->current.proj.z_project.y_off + level * GMT->current.proj.z_project.cos_el;
				break;
		}

		/* First restore the old matrix or save the old one when that was not done before */
		PSL_command (PSL, "%s [%.12g %.12g %.12g %.12g %.12g %.12g] concat\n",
			(GMT->current.proj.z_project.plane >= 0) ? "PSL_GPP setmatrix" : "/PSL_GPP matrix currentmatrix def",
			a, b, c, d, e * PSL->internal.x2ix, f * PSL->internal.y2iy);
	}

	/* Store value of plane */
	GMT->current.proj.z_project.plane = plane;
}

/* Creation of hidden PS0 filename used under modern GMT mode */

/*! . */
int gmt_set_psfilename (struct GMT_CTRL *GMT) {
	int k;
	/* Returns 0 or 1 if successful and -1 if failure */
	GMT->current.ps.figure = gmt_get_current_figure (GMT->parent);

	if (GMT->parent->gwf_dir == NULL) {	/* Use the established temp directory */
		GMT_Report (GMT->parent, GMT_MSG_ERROR, "GMT WorkFlow directory not set!\n");
		return GMT_NOTSET;
	}
	else
		snprintf (GMT->current.ps.filename, GMT_LEN256, "%s/gmt_%d.ps-", GMT->parent->gwf_dir, GMT->current.ps.figure);
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Use PS filename %s\n", GMT->current.ps.filename);
	k = 1 + access (GMT->current.ps.filename, W_OK);	/* 1 = File exists (must append) or 0 (must create) */
	GMT->current.ps.initialize = (k == 0);	/* False means it is an overlay and -R -J may come from history */
	return k;
}

/*! . */
int gmt_ps_append (struct GMT_CTRL *GMT, char *source, unsigned int mode, FILE *dest) {
	/* Append parts of the PS from source file to destination stream.  The part is decided by mode:
	 * mode & 1: Include the header. [-K]
	 * mode & 2: Include the trailer. [-O]
	 * The middle part is always included.
	 */
	FILE *fp = NULL;
	char buffer[GMT_BUFSIZ] = {""};
	bool go = true;
	if (mode == 0 || mode == 2)
		fprintf (dest, "/PSL_xorig 0 def /PSL_yorig 0 def\n");	/* Reset these since we did not use -K in making these pieces */

	if ((fp = fopen (source, "r")) == NULL) {
		GMT_Report (GMT->parent, GMT_MSG_ERROR, "Could not open PostScript file %s\n", source);
		return GMT_NOTSET;
	}
	go = (mode & 1);	/* Start output immediately if header was requested, else we wait until header has passed */
	while (fgets (buffer, GMT_BUFSIZ, fp)) {
		if (!strncmp (buffer, "%PSL_Begin_Trailer", 18U)) {	/* Beginning of the trailer */
			go = (mode & 2);	/* true if we want the trailer */
		}
		if (go) fprintf (dest, "%s", buffer);	/* Write this line to output stream */
		if (!strncmp (buffer, "%PSL_End_Header", 15U)) {	/* Now passed the header section */
			go = true;	/* The juicy middle is always delicious */
		}
	}
	fclose (fp);
	return GMT_NOERROR;
}

/* All functions involved in reading, writing, duplicating GMT_POSTSCRIPT structs and their PostScript content */

/*! . */
struct GMT_POSTSCRIPT * gmtlib_create_ps (struct GMT_CTRL *GMT, uint64_t length) {
	/* Makes an empty GMT_POSTSCRIPT struct - If length > 0 then we also allocate the string */
	struct GMT_POSTSCRIPT *P = gmt_get_postscript (GMT);
	struct GMT_POSTSCRIPT_HIDDEN *PH = gmt_get_P_hidden (P);
	PH->alloc_level = GMT->hidden.func_level;	/* Must be freed at this level. */
	PH->id = GMT->parent->unique_var_ID++;		/* Give unique identifier */
	if (length) {	/* Allocate a blank string */
		P->data = gmt_M_memory (GMT, NULL, length, char);
		PH->n_alloc = length;	/* But P->n_bytes = 0 since nothing was placed there */
		PH->alloc_mode = GMT_ALLOC_INTERNALLY;		/* Memory can be freed by GMT. */
	}
	return (P);
}

/*! . */
void gmtlib_free_ps_ptr (struct GMT_CTRL *GMT, struct GMT_POSTSCRIPT *P) {
	/* Free the memory allocated to hold a PostScript plot (which is pointed to by P->data) */
	unsigned k;
	struct GMT_POSTSCRIPT_HIDDEN *PH = gmt_get_P_hidden (P);
	if (P->data) {
		if (PH->alloc_mode == GMT_ALLOC_INTERNALLY)	/* Was allocated by GMT */
			gmt_M_free (GMT, P->data);
		/* Note: We never need to free the array allocated inside PSL since PSL always destroys it */
	}
	P->data = NULL;		/* Whatever we pointed to is now longer known to P */
	PH->n_alloc = P->n_bytes = 0;
	/* Use free() to free the headers since they were allocated with strdup */
	if (P->n_headers) {
		for (k = 0; k < P->n_headers; k++) gmt_M_str_free (P->header[k]);
		gmt_M_free (GMT, P->header);
	}
	gmt_M_free (GMT, P->hidden);
	P->mode = GMT_PS_EMPTY;
}

/*! . */
void gmtlib_free_ps (struct GMT_CTRL *GMT, struct GMT_POSTSCRIPT **P) {
	/* Free the memory allocated to hold a PostScript plot (which is pointed to by P->data) */
	gmtlib_free_ps_ptr (GMT, *P);
	gmt_M_free (GMT, *P);
	*P = NULL;
}

struct GMT_POSTSCRIPT * gmtlib_read_ps (struct GMT_CTRL *GMT, void *source, unsigned int source_type, unsigned int mode) {
	/* Opens and reads a PostScript file.
	 * Return the result as a GMT_POSTSCRIPT struct.
	 * source_type can be GMT_IS_[FILE|STREAM|FDESC]
	 * mode is not yet used.
	 */

	char ps_file[PATH_MAX] = {""}, buffer[GMT_LEN256] = {""};
	int c;
	bool close_file = false;
	size_t n_alloc = 0;
	struct GMT_POSTSCRIPT *P = NULL;
	struct GMT_POSTSCRIPT_HIDDEN *PH = NULL;
	FILE *fp = NULL;
	gmt_M_unused(mode);

	/* Determine input source */

	if (source_type == GMT_IS_FILE) {	/* source is a file name */
		struct stat buf;
		char path[PATH_MAX] = {""};
		strncpy (ps_file, source, PATH_MAX-1);
		if (!gmt_getdatapath (GMT, ps_file, path, R_OK)) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Cannot find PostScript file %s\n", ps_file);
			return (NULL);
		}
		if (stat (path, &buf)) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Cannot determine size of PostScript file %s\n", ps_file);
			return (NULL);
		}
		if ((fp = fopen (ps_file, "r")) == NULL) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Cannot open PostScript file %s\n", ps_file);
			return (NULL);
		}
		n_alloc = buf.st_size;	/* We know what to allocate here */
		close_file = true;	/* We only close files we have opened here */
	}
	else if (source_type == GMT_IS_STREAM) {	/* Open file pointer given, just copy */
		fp = (FILE *)source;
		if (fp == NULL) fp = GMT->session.std[GMT_IN];	/* Default input */
		if (fp == GMT->session.std[GMT_IN])
			strcpy (ps_file, "<stdin>");
		else
			strcpy (ps_file, "<input stream>");
	}
	else if (source_type == GMT_IS_FDESC) {		/* Open file descriptor given, just convert to file pointer */
		struct stat buf;
		int *fd = source;
		if (fstat (*fd, &buf)) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Cannot determine size of PostScript file give by file descriptor %d\n", *fd);
			return (NULL);
		}
		if ((fp = fdopen (*fd, "r")) == NULL) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Cannot convert PostScript file descriptor %d to stream in gmtlib_read_ps\n", *fd);
			return (NULL);
		}
		else
			close_file = true;	/* Since fdopen creates a FILE struct */
		n_alloc = buf.st_size;		/* We know what to allocate here */
		if (fp == GMT->session.std[GMT_IN])
			strcpy (ps_file, "<stdin>");
		else
			strcpy (ps_file, "<input file descriptor>");
	}
	else {
		GMT_Report (GMT->parent, GMT_MSG_ERROR, "Unrecognized source type %d in gmtlib_read_ps\n", source_type);
		return (NULL);
	}

	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "Reading PostScript from %s\n", ps_file);

	P = gmt_get_postscript (GMT);
	P->header = gmt_M_memory (GMT, NULL, 1, char *);
	snprintf (buffer, GMT_LEN256, "PostScript read from file: %s", ps_file);
	P->header[0] = strdup (buffer);
	P->n_headers = 1;
	if (n_alloc) P->data = gmt_M_memory (GMT, NULL, n_alloc, char);

	/* Start reading PostScript from fp */

	while ((c = fgetc (fp)) != EOF ) {
		if (P->n_bytes >= n_alloc) {
			n_alloc = (n_alloc == 0) ? GMT_INITIAL_MEM_ROW_ALLOC : n_alloc << 1;	/* Start at 2 Mb, then double */
			P->data = gmt_M_memory (GMT, P->data, n_alloc, char);
		}
		P->data[P->n_bytes++] = (char)c;
	}
	if (close_file) fclose (fp);
	if (P->n_bytes > n_alloc)
		P->data = gmt_M_memory (GMT, P->data, P->n_bytes, char);
	PH = gmt_get_P_hidden (P);
	PH->n_alloc = P->n_bytes;
	PH->alloc_mode = GMT_ALLOC_INTERNALLY;	/* So GMT can free the data array */
	/* Determine the mode by checking for typical starts and ends of PS */
	if (P->n_bytes > 4 && !strncmp (P->data, "%!PS", 4U))
		P->mode = GMT_PS_HEADER;	/* Found start of PS header */
	if (P->n_bytes > 10 && !strncmp (&P->data[P->n_bytes-10], "end\n%%EOF\n", 10U))
		P->mode += GMT_PS_TRAILER;	/* Found end of PS trailer */
	return (P);
}

int gmtlib_write_ps (struct GMT_CTRL *GMT, void *dest, unsigned int dest_type, unsigned int mode, struct GMT_POSTSCRIPT *P) {
	/* We write the PostScript file to fp [or stdout].
	 * dest_type can be GMT_IS_[FILE|STREAM|FDESC]
	 * mode is not used yet.
	 */

	bool close_file = false, append = false;
	char ps_file[PATH_MAX] = {""};
	static char *msg1[2] = {"Writing", "Appending"};
	FILE *fp = NULL;
	gmt_M_unused(mode);

	if (dest_type == GMT_IS_FILE && !dest) dest_type = GMT_IS_STREAM;	/* No filename given, default to stdout */

	if (dest_type == GMT_IS_FILE) {	/* dest is a file name */
		static char *msg2[2] = {"create", "append to"};
		strncpy (ps_file, dest, PATH_MAX-1);
		append = (ps_file[0] == '>');	/* Want to append to existing file */
		if ((fp = fopen (&ps_file[append], (append) ? "a" : "w")) == NULL) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Cannot %s PostScript file %s\n", msg2[append], &ps_file[append]);
			return (GMT_ERROR_ON_FOPEN);
		}
		close_file = true;	/* We only close files we have opened here */
	}
	else if (dest_type == GMT_IS_STREAM) {	/* Open file pointer given, just copy */
		fp = (FILE *)dest;
		if (fp == NULL) fp = GMT->session.std[GMT_OUT];	/* Default destination */
		if (fp == GMT->session.std[GMT_OUT])
			strcpy (ps_file, "<stdout>");
		else
			strcpy (ps_file, "<output stream>");
	}
	else if (dest_type == GMT_IS_FDESC) {		/* Open file descriptor given, just convert to file pointer */
		int *fd = dest;
		if (fd && (fp = fdopen (*fd, "a")) == NULL) {
			GMT_Report (GMT->parent, GMT_MSG_ERROR, "Cannot convert PostScript file descriptor %d to stream in gmtlib_write_ps\n", *fd);
			return (GMT_ERROR_ON_FDOPEN);
		}
		if (fd == NULL) fp = GMT->session.std[GMT_OUT];	/* Default destination */
		if (fp == GMT->session.std[GMT_OUT])
			strcpy (ps_file, "<stdout>");
		else
			strcpy (ps_file, "<output file descriptor>");
		close_file = true;	/* since fdopen allocates space */
	}
	else {
		GMT_Report (GMT->parent, GMT_MSG_ERROR, "Unrecognized source type %d in gmtlib_write_ps\n", dest_type);
		return (GMT_NOT_A_VALID_METHOD);
	}
	GMT_Report (GMT->parent, GMT_MSG_DEBUG, "%s PostScript to %s\n", msg1[append], &ps_file[append]);

	/* Start writing PostScript to fp */

	if (fwrite (P->data, 1U, P->n_bytes, fp) != P->n_bytes) {
		GMT_Report (GMT->parent, GMT_MSG_ERROR, "Failure while %s PostScript to %s\n", msg1[append], &ps_file[append]);
		if (close_file) fclose (fp);
		return (GMT_DATA_WRITE_ERROR);
	}
	if (close_file) fclose (fp);
	return (GMT_NOERROR);
}

void gmtplot_copy_ps (struct GMT_CTRL *GMT, struct GMT_POSTSCRIPT *P_copy, struct GMT_POSTSCRIPT *P_obj) {
	/* Just duplicate from P_obj into P_copy */
	struct GMT_POSTSCRIPT_HIDDEN *PH = gmt_get_P_hidden (P_copy);
	if (P_obj->n_bytes > PH->n_alloc) P_copy->data = gmt_M_memory (GMT, P_copy->data, P_obj->n_bytes, char);
	gmt_M_memcpy (P_copy->data, P_obj->data, P_obj->n_bytes, char);
	gmt_M_memcpy (P_copy->hidden, P_obj->hidden, 1, struct GMT_POSTSCRIPT_HIDDEN);
	P_copy->mode = P_obj->mode;
	PH->n_alloc = P_copy->n_bytes = P_obj->n_bytes;
	PH->alloc_mode = GMT_ALLOC_INTERNALLY;	/* So GMT can free the data array */
}

/*! . */
struct GMT_POSTSCRIPT * gmtlib_duplicate_ps (struct GMT_CTRL *GMT, struct GMT_POSTSCRIPT *P_from, unsigned int mode) {
	/* Duplicates a GMT_POSTSCRIPT structure.  Mode not used yet */
	struct GMT_POSTSCRIPT *P = gmtlib_create_ps (GMT, P_from->n_bytes);
	gmt_M_unused(mode);
	gmtplot_copy_ps (GMT, P, P_from);
	return (P);
}

struct GMT_POSTSCRIPT * gmt_get_postscript (struct GMT_CTRL *GMT) {
	struct GMT_POSTSCRIPT *P = NULL;
	P = gmt_M_memory (GMT, NULL, 1, struct GMT_POSTSCRIPT);
	P->hidden = gmt_M_memory (GMT, NULL, 1, struct GMT_POSTSCRIPT_HIDDEN);
	return (P);
}
