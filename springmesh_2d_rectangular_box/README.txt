%==========================================================================
% Copyright (c) 2017, Jorge M. Taramon and Jason P. Morgan, RHUL
%==========================================================================

In order to run this code for the first time, follow the next steps: 
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

The springmesh_2d_rectangular_box folder contains:
- MESH_2D_SPRING.m (main script)
- mfiles_MESH. A folder that contains the routines for creating the mesh that are called by MESH_2D_SPRING.m
- SETUP _TEST. A folder that includes files to setup the output folder and the mesh parameters. 
	If you are going to load a mesh as a 1st guess or a guide-mesh they should be in this folder.
- TOOLS. Empty folder for adding new tools or test to the model.
