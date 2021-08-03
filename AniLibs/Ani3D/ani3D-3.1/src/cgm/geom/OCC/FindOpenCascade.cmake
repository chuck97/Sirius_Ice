#
# ======================================================================================
# OpenCascade related stuff
# sets OPENCASCADE_INCLUDE_DIRS and OPENCASCADE_LIBRARIES
# ======================================================================================


find_path(OPENCASCADE_INC_DIR NAMES Standard_Version.hxx PATHS ${OPENCASCADE_DIR} /usr PATH_SUFFIXES include/opencascade include/oce include/occ include inc ros/inc)
#mark_as_advanced(OPENCASCADE_INC_DIR)
set(OPENCASCADE_INCLUDE_DIRS ${OPENCASCADE_INC_DIR})
set(OPENCASCADE_LIBRARIES)

set(OPENCASCADE_CORE_LIBS TKMesh TKTopAlgo TKGeomAlgo TKBRep TKGeomBase TKG3d TKG2d TKMath TKernel)
set(OPENCASCADE_GEOM_LIBS TKHLR TKOffset TKShHealing TKFillet TKFeat TKBool TKBO TKPrim)
set(OPENCASCADE_LCAF_LIBS TKBinL TKLCAF TKCDF TKCAF)
set(OPENCASCADE_BASE_LIBS TKXSBase)
set(OPENCASCADE_IGES_LIBS TKIGES)
set(OPENCASCADE_STEP_LIBS TKSTEP TKSTEP209 TKSTEPAttr TKSTEPBase)
set(OPENCASCADE_STL_LIBS TKSTL)


set(OPENCASCADE_LIB_DIR "" CACHE PATH "")
foreach(lib ${OPENCASCADE_STL_LIBS} ${OPENCASCADE_STEP_LIBS} ${OPENCASCADE_IGES_LIBS} ${OPENCASCADE_BASE_LIBS} ${OPENCASCADE_LCAF_LIBS} ${OPENCASCADE_GEOM_LIBS} ${OPENCASCADE_CORE_LIBS})
  find_library(OPENCASCADE_${lib}_LIBRARY NAMES ${lib} PATHS ${OPENCASCADE_DIR} /usr PATH_SUFFIXES ros/lib ros/${CMAKE_HOST_SYSTEM_NAME}/lib ros/${CMAKE_HOST_SYSTEM_NAME} lib)
  mark_as_advanced(OPENCASCADE_${lib}_LIBRARY)
  add_library(${lib} UNKNOWN IMPORTED)
  set_target_properties(${lib} PROPERTIES IMPORTED_LOCATION "${OPENCASCADE_${lib}_LIBRARY}")
  if(OPENCASCADE_${lib}_LIBRARY)
    set(OPENCASCADE_LIBRARIES ${OPENCASCADE_LIBRARIES} ${OPENCASCADE_${lib}_LIBRARY})
  endif()
endforeach(lib)

   #message(STATUS "OpenCascade lib dir = ${OPENCASCADE_LIB_DIR}")
   #message(STATUS "OpenCascade inc dir = ${OPENCASCADE_INC_DIR}")
   #message(STATUS "OpenCascade libs    = ${OPENCASCADE_LIBRARIES}")


