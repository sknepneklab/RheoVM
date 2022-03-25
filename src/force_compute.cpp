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
 * \file force_compute.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief ForceCompute class 
 */ 

#include "force_compute.hpp"

namespace RheoVM
{
  void export_ForceCompute(py::module& m)
  {
    py::class_<ForceCompute>(m, "Force")
        .def(py::init<System &>())
        .def("set_params", &ForceCompute::set_params)
        .def("set_vec_params", &ForceCompute::set_vec_params)
        .def("set_relative_params", &ForceCompute::set_relative_params)
        .def("set_flag", &ForceCompute::set_flag)
        .def("set_update", &ForceCompute::set_update)
        .def("set_compute_stress", &ForceCompute::set_compute_stress)
        .def("copy_type_params_to_cell", &ForceCompute::copy_type_param_to_cell)
        .def("add", &ForceCompute::add_force)
        .def("compute", &ForceCompute::compute_forces_and_stresses)
        .def("cell_force", [](ForceCompute &f, int i) -> Vec { return f.compute_face_force(i); });
  }
}