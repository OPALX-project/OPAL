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

option (ENABLE_BANDRF "Compile BANDRF field conversion scripts" OFF)
if (ENABLE_BANDRF)
    add_subdirectory (BandRF)
endif  ()
