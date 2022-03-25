/* ***************************************************************************
 *
 *  Copyright (C) 2017 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of RheoVM (Rheology of Vertex Model) program.
 *
 *  RheoVM is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  RheoVM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

/*!
 * \file dump.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief Dump class 
 */ 

#ifndef __DUMP_HPP__
#define __DUMP_HPP__

#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <exception>
#include <iomanip>
#include <memory>
#include <vector>
#include <utility>
#include <algorithm>
#include <cctype>
#include <sys/stat.h>
#include <regex>

#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolygon.h>
#include <vtkLine.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkZLibDataCompressor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>

#include "json.hpp"

#include "system.hpp"
#include "force_compute.hpp"

using std::back_inserter;
using std::cerr;
using std::copy;
using std::ifstream;
using std::istream_iterator;
using std::istringstream;
using std::map;
using std::move;
using std::ofstream;
using std::runtime_error;
using std::setprecision;
using std::setw;
using std::stod;
using std::stoi;
using std::string;
using std::stringstream;
using std::to_string;
using std::tolower;
using std::transform;
using std::unique_ptr;
using std::vector;

namespace RheoVM
{

  class Dump
  {
    public: 

      Dump(System& sys, ForceCompute& fc) : _sys{sys}, _force_compute{fc}, _sfc(0.95) { }

      void dump_cells(const string&, bool = false, bool = false);   
      void dump_junctions(const string&, bool = false);
      void dump_box(const string&, bool = false);
      void dump_mesh(const string&, bool = false);
      void dump_infiles(const string&);
      void dump_area(const string&);
      void dump_json(const string&);
      void dump_vertices(const list<int>&, const string&);
      void dump_stress(const string&, bool = false, const string& = "");
      void dump_energy(const string&);
      void dump_data(const string &, const list<string> &);

      void set_sfc(double sfc) { _sfc = sfc; }

    private: 

      System& _sys;
      ForceCompute& _force_compute;
      double _sfc;  // Scaling factor for junction output

  };

  void to_json(json&,  const HalfEdge<Property>&);
  void to_json(json&,  const Edge<Property>&);
  void to_json(json&,  const Vertex<Property>&);
  void to_json(json&,  const Face<Property>&);
  void to_json(json&,  const Box&);

  void from_json(const json&, HalfEdge<Property>&);
  void from_json(const json&, Edge<Property>&);
  void from_json(const json&, Vertex<Property>&);
  void from_json(const json&, Face<Property>&);

  vector<string> split(const std::string &, char);

  void export_Dump(py::module&);

}

#endif