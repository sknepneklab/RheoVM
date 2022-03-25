# * ***************************************************************************
# *
# *  Copyright (C) 2017 University of Dundee
# *  All rights reserved. 
# *
# *  This file is part of RheoVM (Rheology of Vertex Model) program.
# *
# *  RheoVM is free software; you can redistribute it and/or modify
# *  it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation; either version 2 of the License, or
# *  (at your option) any later version.
# *
# *  RheoVM is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# *
# * ****************************************************************************/


macro(add_source_dir dir)
add_subdirectory(${dir})

# find all source files in this directory
file(GLOB RheoVM_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.cpp)
file(GLOB RheoVM_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.hpp)

# add them to the source group
source_group(${dir} FILES ${SAMoS_SRCS} ${SAMoS_HDRS})

# add them to the main samos list
list(APPEND _RheoVM_sources ${RheoVM_SRCS} ${RheoVM_HDRS})
endmacro(add_source_dir)
#################

#set(_RheoVM_sources "")


file(GLOB _RheoVM_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file(GLOB _RheoVM_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

SET(RheoVM_HDRS ${_RheoVM_HDRS} PARENT_SCOPE)
SET(RheoVM_SRCS ${_RheoVM_SRCS} PARENT_SCOPE)

#message(${RheoVM_SRCS})

#add_executable(RheoVM ${RheoVM_SRCS})



