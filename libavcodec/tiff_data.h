/*
 * TIFF data tables
 * Copyright (c) 2011 Thomas Kuehnel
 *
 * This file is part of Libav.
 *
 * Libav is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * Libav is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Libav; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/**
 * @file
 * TIFF data tables
 * @author Thomas Kuehnel
 */

#ifndef AVCODEC_TIFF_DATA_H
#define AVCODEC_TIFF_DATA_H

#include "tiff.h"

#define TIFF_CONF_KEY_ID_OFFSET 1024
extern const Geokey ff_tiff_confkeys[3];

#define TIFF_GEOG_KEY_ID_OFFSET 2048
extern const Geokey ff_tiff_geogkeys[14];

#define TIFF_PROJ_KEY_ID_OFFSET 3072
extern const Geokey ff_tiff_projkeys[24];

#define TIFF_VERT_KEY_ID_OFFSET 4096
extern const Geokey ff_tiff_vertkeys[4];

#define TIFF_GEO_KEY_UNDEFINED    0
#define TIFF_GEO_KEY_USER_DEFINED 32767

#define TIFF_GTMODELTYPE_OFFSET 1
extern const char *const ff_tiff_gtmodeltypecodes[3];

#define TIFF_RASTERTYPE_OFFSET 1
extern const char *const ff_tiff_rastertypecodes[2];

#define TIFF_LINEARUNIT_OFFSET 9001
extern const char *const ff_tiff_linearunitcodes[15];

#define TIFF_ANGULARUNIT_OFFSET 9101
extern const char *const ff_tiff_angularunitcodes[8];

#define TIFF_GCSTYPE_OFFSET 4201
extern const char *const ff_tiff_gcstypecodes[133];

#define TIFF_GCSETYPE_OFFSET 4001
extern const char *const ff_tiff_gcsetypecodes[35];

#define TIFF_GEODETICDATUM_OFFSET 6201
extern const char *const ff_tiff_geodeticdatumcodes[120];

#define TIFF_GEODETICDATUME_OFFSET 6001
extern const char *const ff_tiff_geodeticdatumecodes[35];

#define TIFF_ELLIPSOID_OFFSET 7001
extern const char *const ff_tiff_ellipsoidcodes[35];

#define TIFF_PRIMEMERIDIAN_OFFSET 8901
extern const char *const ff_tiff_primemeridiancodes[11];

extern const Key_name ff_tiff_proj_cs_type_codes[978];

extern const Key_name ff_tiff_projection_codes[298];

#define COORD_TRANS_OFFSET 7001
extern const char *const ff_tiff_coord_trans_codes[27];
extern const Key_alias ff_tiff_coord_trans_aliases[9];

extern const Key_alias ff_tiff_key_aliases[4];

#endif
