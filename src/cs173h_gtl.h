/**
 * @file cs173h_gtl.h
 *
 * @section DESCRIPTION
 * 
 **/

#include "etree.h"

/** The configuration structure for the Vs30 map. */
typedef struct cs173h_vs30_map_config_t {
	/** Pointer to the e-tree file */
	etree_t *vs30_map;
	/** The type of map */
	char type[20];
	/** A description of the map */
	char description[50];
	/** The map's author */
	char author[30];
	/** The date the map was created */
	char date[10];
	/** The spacing in meters */
	double spacing;
	/** The map's schema */
	char schema[50];
	/** The projection string in Proj.4 format */
	char projection[128];
	/** The origin point */
	cs173h_point_t origin_point;
	/** The number of degrees the map was rotated around origin */
	double rotation;
	/** The X dimension of the map */
	double x_dimension;
	/** The Y dimension of the map */
	double y_dimension;
	/** The Z dimension of the map */
	double z_dimension;
	/** Number of e-tree ticks in the X direction */
	int x_ticks;
	/** Number of e-tree ticks in the Y direction */
	int y_ticks;
	/** Number of e-tree ticks in the Z direction */
	int z_ticks;
} cs173h_vs30_map_config_t;


/** Contains the Vs30 and surface values from the UCVM map. */
typedef struct cs173h_vs30_mpayload_t {
        /** Surface height in meters */
        float surf;
        /** Vs30 data from Wills and Wald */
        float vs30;
} cs173h_vs30_mpayload_t;

/** Location of the ucvm.e e-tree file. */
char cs173h_vs30_etree_file[128];

/** Holds the configuration parameters for the Vs30 map. */
cs173h_vs30_map_config_t *cs173h_vs30_map;

/** Proj.4 Vs30 map projection holder. */
projPJ cs173h_aeqd;

/** The cosine of the Vs30 map's rotation. */
double cs173h_cos_vs30_rotation_angle = 0;
/** The sine of the Vs30 map's rotation. */
double cs173h_sin_vs30_rotation_angle = 0;

// GTL related
/** Retrieves the vs30 value for a given point. */
int cs173h_get_vs30_based_gtl(cs173h_point_t *point, cs173h_properties_t *data);
/** Reads the specified Vs30 map from UCVM. */
int cs173h_read_vs30_map(char *filename, cs173h_vs30_map_config_t *map);
/** Gets the Vs30 value at a point */
double cs173h_get_vs30_value(double longitude, double latitude, cs173h_vs30_map_config_t *map);
