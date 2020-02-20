option (ENABLE_SDDSTOOLS "Compile conversion tools for SDDS to Opal and vice-versa")

if (ENABLE_SDDSTOOLS)
    find_package (SDDS REQUIRED)
    add_subdirectory (sdds2opal)
    add_subdirectory (opal2sdds)
endif ()

#add_subdirectory (SDDSReader)
option (ENABLE_MSLANG "Compile MSLang stand-alone compiler" OFF)

if (ENABLE_MSLANG)
    add_subdirectory (mslang)
endif ()

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
