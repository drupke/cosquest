-------------------------------------------------------------------------
OVERVIEW
-------------------------------------------------------------------------

Routines for COS QUEST paper.

-------------------------------------------------------------------------
REQUIREMENTS
-------------------------------------------------------------------------

IDL v8.0 or higher (tested with v8.3)

IDL libraries:
- IDL Astronomy User's Library, for various routines
  http://idlastro.gsfc.nasa.gov
- MPFIT, for non-linear least-squares fitting
  http://www.physics.wisc.edu/~craigm/idl/idl.html
- Coyote, for graphics
  http://www.idlcoyote.com/documents/programs.php#COYOTE_LIBRARY_DOWNLOAD
  [or from the subversion repository: https://code.google.com/p/idl-coyote/]
- PPXF
  http://www-astro.physics.ox.ac.uk/~mxc/software/#ppxf
- IDLUTILS, primarily for structure manipulation tasks (e.g.,
  STRUCT_ADDTAGS).
  http://www.sdss3.org/dr8/software/idlutils.php

Note that the IDL Astronomy User's Library ships with some Coyote
routines, and IDLUTILS ships with the IDL Astronomy User's Library and
MPFIT. However, it's not clear how well these libraries keep track of
each other, so it may be preferable to download each package
separately and delete the redundant routines that ship within other
packages.

-------------------------------------------------------------------------
QUESTIONS? BUGS? WANT TO MODIFY THE CODE?
-------------------------------------------------------------------------

Feel free to contact David Rupke at drupke@gmail.com with questions,
bug reports, etc.

Modifications are encouraged, but subject to the license.

-------------------------------------------------------------------------
LICENSE AND COPYRIGHT
-------------------------------------------------------------------------

Copyright (C) 2013--2016 Anthony D. To, David S. N. Rupke

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License or any
later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see http://www.gnu.org/licenses/.
