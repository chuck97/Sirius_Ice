#
#  Ani3D 3-X
#  Local Metric Recovery
#

add_library(lmr3D-${VERSION} STATIC
           "Func2MetricMAX.f"     "Func2GradMetricMAX.f"
           "Nodal2MetricVAR.f"    "Nodal2MetricZZ.f"       
	   "EdgeEst2MetricMAX.f"  "EdgeEst2GradMetricMAX.f" 
           "EdgeEst2MetricLS.f"   "CellEst2MetricZZ.f"     
           "Lp_norm.f"  "utils.f"  "error.f")

install(TARGETS lmr3D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




