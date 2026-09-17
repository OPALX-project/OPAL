option (ENABLE_MSLANG "Compile MSLang stand-alone compiler" OFF)
if (ENABLE_MSLANG)
    add_subdirectory (mslang)
endif ()

option (ENABLE_BANDRF "Compile BANDRF field conversion scripts" OFF)
if (ENABLE_BANDRF)
    add_subdirectory (BandRF)
endif  ()
