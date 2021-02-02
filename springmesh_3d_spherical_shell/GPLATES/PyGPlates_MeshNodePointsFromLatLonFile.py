
"""
    Copyright (C) 2014 The University of Sydney, Australia
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

from __future__ import print_function
from optparse import OptionParser
import sys

import pygplates


def parse_points(input_lat_lon_points_filename):
    """Parse ascii file containing a lat/lon point per line."""
    
    points = []
    with open(input_lat_lon_points_filename, 'r') as lat_lon_points_file:
            for line_number, line in enumerate(lat_lon_points_file):
                lat_lon = line.split()
                
                if len(lat_lon) != 2:
                    print('Line {0}: Ignoring point - line does not have exactly two white-space separated strings.'.format(line_number), file=sys.stderr)
                    continue
                try:
                    lat = float(lat_lon[0])
                    lon = float(lat_lon[1])
                except ValueError:
                    print('Line {0}: Ignoring point - cannot read lat/lon values.'.format(line_number), file=sys.stderr)
                    continue
                
                lat_lon_point = pygplates.LatLonPoint(lat, lon)
                point = pygplates.convert_lat_lon_point_to_point_on_sphere(lat_lon_point)
                points.append(point)
    
    return points


def create_mesh_node_feature_from_points(points):
    """Create a GPlates 'gpml:MeshNode' feature from a sequence of points (pygplates.PointOnSphere objects)."""
    
    # Create the new multipoint feature.
    # This type of feature 'MeshNode' will cause GPlates to automatically create a velocity layer upon loading the points.
    mesh_node_feature = pygplates.Feature(pygplates.FeatureType.create_gpml('MeshNode'))
    
    multipoint = pygplates.MultiPointOnSphere(points)
    mesh_node_feature.add(pygplates.PropertyName.create_gpml('meshPoints'),
    					  pygplates.GmlMultiPoint(multipoint))
    
    # Add time period property.
    mesh_node_feature.add(pygplates.PropertyName.create_gml('validTime'),
    					  pygplates.GmlTimePeriod(pygplates.GeoTimeInstant.create_distant_past(),
                            					  pygplates.GeoTimeInstant.create_distant_future()))
    
    # Add reconstruction plate id property (use plate id zero).
    mesh_node_feature.add(pygplates.PropertyName.create_gpml('reconstructionPlateId'),
    					  pygplates.GpmlConstantValue(pygplates.GpmlPlateId(0)))
    
    
    return mesh_node_feature


if __name__ == "__main__":

    __usage__ = "%prog [options] [-h --help] input_filename output_filename"
    __description__ = "Loads an ascii file (containing a latitude longitude point feature per line) and saves to a GPlates-compatible file."

    # Parse the command-line options.    
    parser = OptionParser(usage = __usage__,
                          description = __description__)

    # Parse command-line options.
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")
    
    input_lat_lon_points_filename = args[0]
    output_points_filename = args[1]
    
    # Parse the lat/lon points file and fill the points feature collection.
    points = parse_points(input_lat_lon_points_filename)
    mesh_node_feature = create_mesh_node_feature_from_points(points)
    
    # Create an empty points feature collection.
    points_feature_collection = pygplates.FeatureCollection()
    
    # A feature collection containing one multipoint feature.
    points_feature_collection.add(mesh_node_feature)

    file_registry = pygplates.FeatureCollectionFileFormatRegistry()
    
    # Write the points feature collection to disk.
    file_registry.write(points_feature_collection, output_points_filename)
