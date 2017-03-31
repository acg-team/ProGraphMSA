#
# Copyright (c) 2007-2011 ETH Zürich, Institute of Computational Science
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

FIND_PROGRAM(EXECUTABLE_JAVA
   NAMES java
   PATHS /bin /usr/bin /usr/local/bin
   DOC "path to java executable")

################################################################################
#                 _                        _   _                 
#   ___ __ _  ___| |__   ___    ___  _ __ | |_(_) ___  _ __  ___ 
#  / __/ _` |/ __| '_ \ / _ \  / _ \| '_ \| __| |/ _ \| '_ \/ __|
# | (_| (_| | (__| | | |  __/ | (_) | |_) | |_| | (_) | | | \__ \
#  \___\__,_|\___|_| |_|\___|  \___/| .__/ \__|_|\___/|_| |_|___/
#                                   |_|                          
################################################################################

SET(EIGEN_INCLUDE_DIR "/usr/include/eigen3" CACHE STRING "path to libEigen3 headers")
SET(TREKS_JAR_PATH "T-Reks.jar" CACHE STRING "path to T-Reks jar archive")
SET(WITH_DEBUG OFF CACHE BOOL "enable debug output")
SET(WITH_SSE ON CACHE BOOL "enable code using SSE2 instructions")
SET(WITH_LS_LOG ON CACHE BOOL "use fast approximation of logarithm")
SET(WITH_CODON ON CACHE BOOL "enable codon alignment")
SET(WITH_DNA ON CACHE BOOL "enable DNA alignment")

################################################################################
#                                       ___     _                       _       
#  ___  ___  _   _ _ __ ___ ___  ___   ( _ )   | |_ __ _ _ __ __ _  ___| |_ ___ 
# / __|/ _ \| | | | '__/ __/ _ \/ __|  / _ \/\ | __/ _` | '__/ _` |/ _ \ __/ __|
# \__ \ (_) | |_| | | | (_|  __/\__ \ | (_>  < | || (_| | | | (_| |  __/ |_\__ \
# |___/\___/ \__,_|_|  \___\___||___/  \___/\/  \__\__,_|_|  \__, |\___|\__|___/
#                                                            |___/              
################################################################################

INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

CONFIGURE_FILE(external_programs.cpp.in external_programs.cpp)

SET(SOURCES Alphabet.cpp CSProfile.cpp debug.cpp DistanceFactory.cpp
   DistanceFactoryAlign.cpp DistanceFactoryML.cpp Fasta.cpp FindRoot.cpp
   LeastSquares.cpp main.cpp ModelFactory.cpp ModelFactoryDarwin.cpp
   ModelFactoryEcm.cpp ModelFactoryWag.cpp newick.cpp PhyTree.cpp
   ProgressiveAlignment.cpp RepeatDetectionTReks.cpp Stockholm.cpp TempFile.cpp
   TreeNJ.cpp)

SET(DIST_SOURCES Alphabet.cpp dist.cpp DistanceFactoryAlign.cpp
   ModelFactoryDarwin.cpp ModelFactoryWag.cpp Fasta.cpp debug.cpp)

ADD_CUSTOM_COMMAND(
   OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/DateCompiled.cpp
   COMMAND /bin/sh
   ARGS    ${CMAKE_CURRENT_SOURCE_DIR}/DateCompiled.sh
   DEPENDS DateCompiled.sh ${SOURCES}
   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
   COMMENT "Updating build date")

ADD_EXECUTABLE(ProGraphMSA ${SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/DateCompiled.cpp ${CMAKE_CURRENT_BINARY_DIR}/external_programs.cpp)
ADD_EXECUTABLE(dist EXCLUDE_FROM_ALL ${DIST_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/DateCompiled.cpp)

################################################################################
#                                                              
#  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___  ___  _ __ 
# | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|/ _ \| '__|
# | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \ (_) | |   
# | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/\___/|_|   
# |_|            |_|                                           
#      _       __ _       _ _   _                 
#   __| | ___ / _(_)_ __ (_) |_(_) ___  _ __  ___ 
#  / _` |/ _ \ |_| | '_ \| | __| |/ _ \| '_ \/ __|
# | (_| |  __/  _| | | | | | |_| | (_) | | | \__ \
#  \__,_|\___|_| |_|_| |_|_|\__|_|\___/|_| |_|___/
#                                                 
################################################################################

IF(WITH_DEBUG)
   ADD_DEFINITIONS(-DDEBUG)
ENDIF(WITH_DEBUG)

IF(WITH_SSE)
   #TODO icc does only recognize -msse3
   ADD_DEFINITIONS(-msse2)
ENDIF(WITH_SSE)

IF(WITH_LS_LOG)
   ADD_DEFINITIONS(-DUSE_LS_LOG)
ENDIF(WITH_LS_LOG)

IF(WITH_CODON)
   ADD_DEFINITIONS(-DWITH_CODON)
ENDIF(WITH_CODON)

IF(WITH_DNA)
   ADD_DEFINITIONS(-DWITH_DNA)
ENDIF(WITH_DNA)

################################################################################
#  _           _        _ _       _   _             
# (_)_ __  ___| |_ __ _| | | __ _| |_(_) ___  _ __  
# | | '_ \/ __| __/ _` | | |/ _` | __| |/ _ \| '_ \ 
# | | | | \__ \ || (_| | | | (_| | |_| | (_) | | | |
# |_|_| |_|___/\__\__,_|_|_|\__,_|\__|_|\___/|_| |_|
#                                                   
################################################################################

INSTALL(TARGETS ProGraphMSA DESTINATION bin)

SET(CPACK_GENERATOR "STGZ;TGZ;TBZ2;DEB")
SET(CPACK_PACKAGE_CONTACT "adam.szalkowski@inf.ethz.ch")
INCLUDE(CPack)
