# XC_AC_PointSource_Pipeline
This is the pipeline that will be run to calibrate the transient 450 (and possibly 850) micron data - Dec 2020.

Plan:

1. Align with the cross correlation method (may need to iterate)
1. Flux Calibrate with the auto correlation method
1. Flux Calibrate with the Point Source Method (to avoid variables)

In order to work, the jcmt_transient_alignment directory will need to be available as it contains all the XC, AC, and (soon) Local point source analysis functions necessary for this pipeline.
