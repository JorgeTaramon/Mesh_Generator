#!/bin/bash
cd /Users/jorge/Downloads/Mesh_Generator-master/springmesh_3d_spherical_shell/GPLATES/
export PYTHONPATH=$PYTHONPATH://Users/jorge/Downloads/Mesh_Generator-master/springmesh_3d_spherical_shell/GPLATES/pygplates
python PyGPlates_MeshNodePointsFromLatLonFile.py surf_nodes.txt points.gpml