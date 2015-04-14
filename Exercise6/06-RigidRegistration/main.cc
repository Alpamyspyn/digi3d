//=============================================================================
//
//   Example code for the full-day course
//
//   M. Botsch, M. Pauly, C. Roessl, S. Bischoff, L. Kobbelt,
//   "Geometric Modeling Based on Triangle Meshes"
//   held at SIGGRAPH 2006, Boston, and Eurographics 2006, Vienna.
//
//   Copyright (C) 2006 by  Computer Graphics Laboratory, ETH Zurich,
//                      and Computer Graphics Group,      RWTH Aachen
//
//
//-----------------------------------------------------------------------------
//
//                                License
//
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License
//   as published by the Free Software Foundation; either version 2
//   of the License, or (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 51 Franklin Street, Fifth Floor,
//   Boston, MA  02110-1301, USA.
//
//=============================================================================
#include "RegistrationViewer.hh"
#include <fstream>

bool
readIndicesFile (const std::string &filename, std::vector<size_t> &indices)
{
  std::ifstream file (filename.c_str ());
  if (!file.is_open ())
  {
    printf ("Error reading indices file %s.\n", filename.c_str ());
    return (false);
  }
  size_t points_size;
  file >> points_size;
  for (size_t i = 0; i < points_size; ++i)
  {
    size_t v;
    file >> v;
    indices.push_back (v);
  }
  file.close ();
  return (true);
}

int main (int argc, char **argv)
{
  glutInit (&argc, argv);

  RegistrationViewer window ("Registration", 512, 512);

  std::vector<std::string> files;
  /// Template mesh
  files.push_back (std::string (argv[1]));
  /// Target mesh
  files.push_back (std::string (argv[2]));

  std::vector<size_t> indices;
  readIndicesFile (argv[3], indices);

  if (argc > 1)
  {
    window.open_meshes (files, indices);
    window.initialize ();
  }

  glutMainLoop ();
}
