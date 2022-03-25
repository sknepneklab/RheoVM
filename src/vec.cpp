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
 * \file vec.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Jan-2018
 * \brief 2D vectors 
 */ 


#include "vec.hpp"

namespace RheoVM
{
  void export_Vec(py::module& m)
  {
    py::class_<Vec>(m, "Vec")
          .def(py::init<double,double>())
          .def(py::init<double,double,const shared_ptr<Box>&>())
          .def_readwrite("x", &Vec::x)
          .def_readwrite("y", &Vec::y)
          .def("box", [](const Vec& v) { if (v.box) return v.box->get(); else return vector<double>({0,0,0,0}); })
          .def("__add__", [](const Vec& v1, const Vec& v2) { return Vec(v1 + v2); })
          .def("__sub__", [](const Vec& v1, const Vec& v2) { return Vec(v1 - v2); })
          .def("__mul__", [](const Vec& v, double c) { return c*v; })
          .def("__rmul__", [](const Vec& v, double c) { return c*v; })
          .def("__repr__", [](const Vec& v) { return "("+std::to_string(v.x)+","+std::to_string(v.y)+")";  })
          .def("unit", [](const Vec& v) { return Vec(v.x/v.len(), v.y/v.len(), v.box); });

  }
}