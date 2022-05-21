# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Debug")
  file(REMOVE_RECURSE
  "CMakeFiles\\ModellingaSpacecraft_autogen.dir\\AutogenUsed.txt"
  "CMakeFiles\\ModellingaSpacecraft_autogen.dir\\ParseCache.txt"
  "ModellingaSpacecraft_autogen"
  )
endif()
