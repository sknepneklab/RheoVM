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
 * \file simulation.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Jan-2018
 * \brief Simulation class 
 */ 

#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__

#define XSTR(s) STR(s)
#define STR(s) #s

#include <iostream>
#include <string>

#include "system.hpp"
#include "force_compute.hpp"
#include "integrate.hpp"
#include "measurement.hpp"
#include "topology.hpp"
#include "dump.hpp"
#include "version.hpp"

using std::cout;
using std::string;
using std::endl;
using std::to_string;

namespace RheoVM
{

  struct Simulation
  {
    Simulation(System& sys, Integrate& integ, ForceCompute& f, Topology& t) : _sys{sys}, 
                                                                              _integ{integ}, 
                                                                              _force_compute{f},
                                                                              _topology{t},
                                                                              print_freq{100},
                                                                              dump_cells_freq{100},
                                                                              dump_junctions_freq{0},
                                                                              dump_mesh_freq{0},
                                                                              dump_area_freq{0},
                                                                              num_zeros{8},
                                                                              cell_base_name{"cells_"},
                                                                              junction_base_name{"junctions_"},
                                                                              mesh_base_name{"mesh_"},
                                                                              area_base_name{"area_"},
                                                                              bar_width{40},
                                                                              sim_step{0}
                                                                              { }
    void run(int, bool = true, bool = true);
    const string print_version() { return  "branch : "+static_cast<string>(XSTR(GIT_BRANCH))+" commit : "+static_cast<string>(XSTR(GIT_COMMIT_HASH)); }
    void progress_bar(double, const string&);

    // variables and paramters
    System& _sys;
    Integrate& _integ;
    ForceCompute& _force_compute;
    Topology& _topology;

    //int rescale_freq;
    int print_freq;
    int dump_cells_freq;
    int dump_junctions_freq;
    int dump_mesh_freq;
    int dump_area_freq;
    int num_zeros;
    string cell_base_name;
    string junction_base_name;
    string mesh_base_name;
    string area_base_name;
    int bar_width;
    int sim_step;

  };

  void export_Simulation(py::module&);

}

#endif
