# Generation of unstructured meshes in 2-D, 3-D, and spherical geometries with embedded high-resolution sub-regions

### Authors:
Jorge M. Taramón (jorge.taramongomez.2014@live.rhul.ac.uk), Jason P. Morgan, Chao Shi, Jörg Hasenclever

### Description:
- Mesh generator for Cartesian 2-D/3-D and Earth-appropriate spherical coordinates
- Builds unstructured tetrahedral meshes with embedded high-resolution sub-regions
- Guide-mesh defines variable preferred element sizes in the mesh
- Tools for mesh quality improvement (removing tetrahedral slivers)
- Algorithm suitable for adaptive mesh refinement

### How to cite and background information:


### Running this code for the first time:

1. Download "SuiteSparse" from:
    http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.0.2.tar.gz
2. Install it running SuiteSparse_install.m
3. Download mutils-0.4-2 from:
    https://sourceforge.net/projects/milamin/files/mutils-0.4-2.zip
4. Install it running install.m (follow the instructions that will appear in the Command Window)
5. In order to set the path for mutils
    5.1. Write pctconfig in the Command Window
    5.2. Open SETUP_TEST/addpaths_mutils
    5.3. Create a new case for your computer:
	- Introduce the hostname you obtained in 5.1 as a new ‘case’, e.g., case 'fpdc462'
	- Introduce the path for the folder mutils in ‘path2mutils’
6. In order to set the OUTPUT folder in your computer:
    6.1. Open SETUP_TEST/data_storage_2d.m
    6.2. Create a new case for your computer:
        - Introduce the hostname you obtained in 5.1 as a new ‘case’, e.g., case 'fpdc462'
        - Introduce the path to save the output data in ‘path2data’
