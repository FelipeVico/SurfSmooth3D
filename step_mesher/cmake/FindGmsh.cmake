set(_GMSH_HINT_INCLUDE_DIRS
  "$ENV{GMSH_DIR}/include"
  "/Library/Frameworks/Python.framework/Versions/3.12/include"
  "/opt/homebrew/include"
  "/opt/homebrew/opt/gmsh/include"
  "/usr/local/include"
)

set(_GMSH_HINT_LIBRARY_DIRS
  "$ENV{GMSH_DIR}/lib"
  "/Library/Frameworks/Python.framework/Versions/3.12/lib"
  "/opt/homebrew/lib"
  "/opt/homebrew/opt/gmsh/lib"
  "/usr/local/lib"
)

find_path(GMSH_INCLUDE_DIR
  NAMES gmsh.h
  HINTS ${_GMSH_HINT_INCLUDE_DIRS}
)

find_library(GMSH_LIBRARY
  NAMES gmsh libgmsh.4.12.dylib libgmsh.dylib
  HINTS ${_GMSH_HINT_LIBRARY_DIRS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gmsh
  REQUIRED_VARS GMSH_INCLUDE_DIR GMSH_LIBRARY
)

if(Gmsh_FOUND AND NOT TARGET Gmsh::Gmsh)
  get_filename_component(GMSH_LIBRARY_DIR "${GMSH_LIBRARY}" DIRECTORY)
  add_library(Gmsh::Gmsh UNKNOWN IMPORTED)
  set_target_properties(Gmsh::Gmsh PROPERTIES
    IMPORTED_LOCATION "${GMSH_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${GMSH_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(GMSH_INCLUDE_DIR GMSH_LIBRARY)
