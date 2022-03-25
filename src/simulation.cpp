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
 * \file simulation.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Jan-2018
 * \brief Simulation class 
 */ 

#include "simulation.hpp"

namespace RheoVM
{
  void Simulation::run(int steps, bool topological_change, bool old_style)
  {
    double progress = 0.0;
    for (int i = sim_step; i < sim_step + steps; i++)
    {
      if (topological_change)
      {
        _topology.collapse();
        _topology.split();
      }
      _integ.apply();
      // reset topology change flag
      _sys.set_topology_change(false);

      if (this->print_freq > 0) 
      {
        if (old_style)
        {
          if (i % this->print_freq == 0)
            cout << "step : " << i << endl;
        }
        else 
        {
          progress_bar(progress, "\r");
          progress += 1.0 / steps;
        }
        
      }
    }
    if (this->print_freq > 0 && !old_style)
      progress_bar(progress, " ");
    
    sim_step += steps;
    if (this->print_freq > 0 && !old_style)
      std::cout << " --> Completed " << sim_step << " simlation steps " << std::endl;  
    
  }

  void Simulation::progress_bar(double progress, const string& end_of_line)
  {
    std::cout << "[";
    int pos = static_cast<int>(round(bar_width * progress));
    for (int i = 0; i < bar_width; ++i) 
    {
      if (i < pos) 
        std::cout << "=";
      else if (i == pos) 
        std::cout << ">";
      else 
        std::cout << " ";
    }
    std::cout << "] "  << static_cast<int>(round(progress * 100.0)) << "%" << end_of_line;
    std::cout.flush();
  }

  void export_Simulation(py::module& m)
  {
    py::class_<Simulation>(m, "Simulation")
          .def(py::init<System&,Integrate&,ForceCompute&,Topology&>())
          .def_readwrite("print_freq", &Simulation::print_freq)
          .def_readwrite("dump_cells_freq", &Simulation::dump_cells_freq)
          .def_readwrite("dump_junctions_freq", &Simulation::dump_junctions_freq)
          .def_readwrite("dump_mesh_freq", &Simulation::dump_mesh_freq)
          .def_readwrite("dump_area_freq", &Simulation::dump_area_freq)
          .def_readwrite("num_zeros", &Simulation::num_zeros)
          .def_readwrite("cell_base_name", &Simulation::cell_base_name)
          .def_readwrite("junction_base_name", &Simulation::junction_base_name)
          .def_readwrite("mesh_base_name", &Simulation::mesh_base_name)
          .def_readwrite("area_base_name", &Simulation::area_base_name)
          .def_readwrite("bar_width", &Simulation::bar_width)
          .def("run", &Simulation::run, py::arg("steps"), py::arg("topological_change") = true, py::arg("old_style") = true)
          .def("print_version", &Simulation::print_version);
  }  

}

PYBIND11_MODULE(RheoVM, m)
{
  RheoVM::export_Vec(m);
  RheoVM::export_Box(m);
  RheoVM::export_T1_stats(m);
  RheoVM::export_VertexProperty(m);
  RheoVM::export_HEProperty(m);
  RheoVM::export_EdgeProperty(m);
  RheoVM::export_FaceProperty(m);
  RheoVM::export_Vertex(m);
  RheoVM::export_Edge(m);
  RheoVM::export_HalfEdge(m);
  RheoVM::export_Face(m);
  RheoVM::export_Mesh(m);
  RheoVM::export_System(m);
  RheoVM::export_Texture(m);
  RheoVM::export_Site(m);
  RheoVM::export_Measurement(m);
  RheoVM::export_ForceCompute(m);
  RheoVM::export_Integrate(m);
  RheoVM::export_Topology(m);
  RheoVM::export_Dump(m);
  RheoVM::export_Simulation(m);
}

