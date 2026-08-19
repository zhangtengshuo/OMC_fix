#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
#***********************************************************************

set (sources
     main.F90
     exactemb_data.F90
     rdinp_exactemb.F90
     exactemb.F90
    )

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)
