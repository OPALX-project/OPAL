# gtest libraries are required for the tests
find_package(GTest)

## PythonExprTest is disabled
#find_package(PythonLibs)
#if (NOT PYTHONLIBS_FOUND)
#    message (STATUS "Python libraries not found.")
#endif()

if (NOT GTEST_FOUND)
    message (STATUS "gtest not found.")
    message (STATUS "No testing possible.")
    return()
endif()

include_directories (
    ${GTEST_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
#    ${PYTHON_INCLUDE_DIRS}
)

link_directories (
    ${GTEST_LIBRARY_DIRS}
#    ${PYTHON_LIBRARIES}
)

set (OPT_PILOT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/optimizer)
set (SDDSPARSER_SRC
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/array.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/associate.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/ast.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/column.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/data.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/description.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/file.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/include.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/parameter.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser/version.cpp
    ${OPT_PILOT_SOURCE_DIR}/Util/SDDSParser.cpp
)

set (ExpressionTest_SRC
    ${SDDSPARSER_SRC}
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
)

set (SumErrSqExpressionTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Expression/SumErrSq.cpp
    ${SDDSPARSER_SRC}
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
)

set (FromFileExpressionTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Expression/FromFile.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
)

#set (PythonExprTest_SRC
#    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
#    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
#)


set (HashNameGeneratorTest_SRC
)

set (CmdArgumentsTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Util/CmdArguments.cpp
)

set (MPIHelperTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Util/MPIHelper.cpp
)

set (ManagedIDsTest_SRC
)

set (IndividualTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
)

set (PopulationTest_SRC
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/expression.cpp
    ${OPT_PILOT_SOURCE_DIR}/Expression/Parser/evaluator.cpp
)

set (SDDSParserTest_SRC
    ${SDDSPARSER_SRC}
)

set (PILOT_TESTS
    MPIHelperTest
    SDDSParserTest
    ManagedIDsTest
    IndividualTest
    PopulationTest
    ExpressionTest
    SumErrSqExpressionTest
    FromFileExpressionTest
#    PythonExprTest
    CmdArgumentsTest
    HashNameGeneratorTest
)

set (TEST_LIBS
    ${GTEST_BOTH_LIBRARIES}
    pthread
    ${Boost_LIBRARIES}
    ${MPI_LIBRARIES}
)

#set (PYTHON_LIBS
#    ${PYTHON_LIBRARIES}
#    boost_python
#)

FOREACH(test ${PILOT_TESTS})
    SET( SRC ${${test}_SRC} )
    ADD_EXECUTABLE( ${test}.exe ${test}.cpp ${SRC} )
    TARGET_LINK_LIBRARIES( ${test}.exe ${TEST_LIBS}) # ${PYTHON_LIBS})
    ADD_TEST(${test} ${test}.exe)
ENDFOREACH(test ${PILOT_TESTS})

#link resources dir to test dir
if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/resources)
    EXECUTE_PROCESS(
        COMMAND ln -s ${CMAKE_CURRENT_SOURCE_DIR}/resources ${CMAKE_CURRENT_BINARY_DIR}
        OUTPUT_VARIABLE LN_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
        )
endif()

# vi: set et ts=4 sw=4 sts=4:
  
# Local Variables:
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:
