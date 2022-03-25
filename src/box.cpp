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
 * \file box.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Jan-2020
 * \brief Simulation Box
 */ 


#include "box.hpp"

namespace RheoVM
{
  void export_Box(py::module& m)
  {
    py::class_<Box, shared_ptr<Box>>(m, "Box")
          .def(py::init<double,double>())
          .def(py::init<double,double,double,double>())
          .def("a", [](const Box& b) { return vector<double>( {b.h._mxx, b.h._myx} ); } )
          .def("b", [](const Box& b) { return vector<double>( {b.h._mxy, b.h._myy} ); } );

  }
}