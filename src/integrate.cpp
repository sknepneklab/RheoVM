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
 * \file integrate.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-May-2019
 * \brief Integrate class 
 */ 

#include "integrate.hpp"

namespace RheoVM
{
  void export_Integrate(py::module& m)
  {
    py::class_<Integrate>(m, "Integrate")
      .def(py::init<System&, ForceCompute&, int>())
      .def("set_params", &Integrate::set_params)
      .def("set_type_params", &Integrate::set_type_params)
      .def("set_external_force", &Integrate::set_external_force)
      .def("set_radial_force", &Integrate::set_radial_force)
      .def("set_flag", &Integrate::set_flag)
      .def("enable", &Integrate::enable)
      .def("disable", &Integrate::disable)
      .def("set_dt", &Integrate::set_dt)
      .def("enable_constraint", &Integrate::enable_constraint)
      .def("converged", &Integrate::converged)
      .def("add", &Integrate::add_integrator);
  }
}