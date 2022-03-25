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
 * \file measurement.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-May-2018
 * \brief Measurement class 
 */ 

#ifndef __MEASUREMENT_HPP__
#define __MEASUREMENT_HPP__

#include <vector>
#include <stdexcept>
#include <string>
#include <fstream>

#include "system.hpp"

using std::vector;
using std::runtime_error;
using std::string;
using std::ofstream;

namespace RheoVM
{

  struct Texture
  {
    Texture(double xx, double xy, double yy) : xx(xx), xy(xy), yy(yy) { }
    double xx;
    double xy;
    double yy;
  };

  struct Site 
  {
    Site(int id, Vec& r) : id(id), r(r) { }
    int id;
    Vec r;
    vector<int> neigh;
    vector<Vec> l;  // distance vectors to all neighbours
    vector< HEHandle<Property> > jc;  // junction references
    vector<Texture> m;
    vector<Texture> myom;
  };

  class Measurement 
  {
    public:

      Measurement(MyMesh& m, System& sys) : _mesh{m}, 
                                            _box{nullptr},
                                            _sys{sys}
                                            { }

      MyMesh& mesh()  { return _mesh; }
      vector<Site>& sites() { return _sites; }

      void set_box(const shared_ptr<Box>& box) { _box = box; }
      bool periodic() { return (_box != nullptr); }

      void populate();
      void populate_pairs(const std::vector<std::pair<int,int> >&);
      void compute_texture();
      void dump_texture(const string&);
      void dump_junctions(const std::vector<std::pair<int,int> >& , const string&);
      void dump_T1(const string&);


    private:

      MyMesh& _mesh;                             //tissue 
      System& _sys;                              // system
      shared_ptr<Box> _box;
      vector<Site> _sites;

  };

  void export_Texture(py::module&);
  void export_Site(py::module&);
  void export_Measurement(py::module&);

}

#endif
