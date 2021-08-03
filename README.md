# XC_AC_PointSource_Pipeline
This is the pipeline that is used for the JCMT Transient Survey Calibration 2.0. 

To begin, create directories for each region with the names:

1. IC348
1. NGC1333
1. NGC2024
1. NGC2071
1. OPHCORE
1. OMC23
1. SERPM
1. SERPS
1. DR21C
1. DR21N
1. DR21S
1. M17
1. M17SWex
1. S255

Populate each directory with all the **850 micron** images you would like to process from the R3 directory on VOspace: (Example: DR21C_*_850_ER3.sdf)
The 450 images will be created automatically from the 850 micron names, so no need to place the 450 micron data in the directories, above.

Ensure an up-to-date Starlink installation is running on the computer.

*Note: If you have already run the pipeline on the full dataset, and you would just like the results for the newest observation, you can simply place
in these directories:*

1. The first observation of the region (ER3)
1. The newest observation of the region (ER3)

*The pipeline will automatically **not** re-reduce the first observation of the region, since that has already been done.*

The following is a chart describing the program: **Transient_Cal_Pipeline_v2_full.py**.

At the top of the program, you can select whether you would like to run the cross-correlation/alignment steps (which involve re-reducing the maps that have not already been reduced)
or if you just want to run the point source calibration (PScal) only (i.e. the weighted flux calibration and localised "source-skewering" calibration). The latter also
performs the pointing checks and generates the light curves.

![image info](./Pipeline.png)
