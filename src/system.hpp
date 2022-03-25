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
 * \file system.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 12-Jun-2017
 * \brief System class 
 */ 

#ifndef __SYSTEM_HPP__
#define __SYSTEM_HPP__

#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <exception>
#include <iomanip>
#include <memory>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/detail/xml_parser_writer_settings.hpp>

#include "json.hpp"

#include "rng.hpp"
#include "mesh.hpp"
#include "property.hpp"


using std::string;
using std::istringstream;
using std::istream_iterator;
using std::copy;
using std::back_inserter;
using std::stoi;
using std::stod;
using std::ifstream;
using std::map;
using std::cerr;
using std::runtime_error;
using std::setw;
using std::to_string;
using std::shared_ptr;
using std::make_shared;
using std::ofstream;

namespace pt = boost::property_tree;
using json = nlohmann::json;

namespace RheoVM
{

  typedef Mesh<Property> MyMesh;
  typedef map<string,int> type_data;
  typedef map<int, string> type_name;

  vector<string> split_string(const string&);
  string trim(const string&);

  enum SPLIT_TYPE { FORWARD, BACKWARD };

  typedef map<string, double> params_type;         // Used when we set a numerical value to a paramter, e.g., kappa = 1.0
  typedef map<string, Vec> vec_params_type;        // Used when we set a vectorial value to a paramter, e.g., n = Vec(1,0)
  typedef map<string, vector<double>> multi_params_type;   // Used when we need to at least two values to set a paramter, e.g., a paramters is drawn from a random distribution 

  bool operator<(const VertexHandle<Property>&, const VertexHandle<Property>&);

  class System
  {
    public:

      System(MyMesh& mesh) : _mesh{mesh}, 
                             _time_step{0},
                             _simulation_time{0.0},
                             _num_cell_types{0},
                             _num_vert_types{0},
                             _mesh_set{false},
                             _topology_changed{true}
                             { 
                               
                             }

      MyMesh& mesh()  { return _mesh; }

      void set_box(const shared_ptr<Box>& box) { _mesh.set_box(box); }
      const shared_ptr<Box> &box() const { return _mesh.box(); }
      bool periodic() { return (_mesh.box() != nullptr); }

      void read_input(const string&, bool = false);

      void set_simulation_time_step(int time_step) { _time_step = time_step; }
      int& time_step() { return _time_step; }
      double& simulation_time() { return _simulation_time; }
      
      
      type_data& cell_types() { return _cell_types; }
      type_data& vert_types() { return _vert_types; }
      const string get_cell_type_name(const int type_id) const { return _cell_types_map.at(type_id); }
      const string get_vert_type_name(const int type_id) const { return _vert_types_map.at(type_id); }
      void add_cell_type(const string& cell_type)
      {
        if (_cell_types.find(cell_type) == _cell_types.end())
        {
          _cell_types_map[_num_cell_types] = cell_type;
          _cell_types[cell_type] = _num_cell_types++;
        }
      }
      void add_vert_type(const string& vert_type)
      {
        if (_vert_types.find(vert_type) == _vert_types.end())
        {
          _vert_types_map[_num_vert_types] = vert_type;
          _vert_types[vert_type] = _num_vert_types++;
        }
      }

      
      void displace_vertices(const string &, const Vec&);

      void set_topology_change(bool flag) { _topology_changed = flag; }
      bool topology_change() { return _topology_changed;  }

    private:
      MyMesh &_mesh;
      type_data _cell_types;
      type_name _cell_types_map;
      type_data _vert_types;
      type_name _vert_types_map;
      int _time_step;
      double _simulation_time;
      int _num_cell_types;
      int _num_vert_types; 
      bool _mesh_set;
      bool _topology_changed;  // If true, mesh topology has changed
  };

  void export_T1_stats(py::module&);
  void export_VertexProperty(py::module&);
  void export_EdgeProperty(py::module&);
  void export_HEProperty(py::module&);
  void export_FaceProperty(py::module&);
  void export_Spoke(py::module &);
  void export_Vertex(py::module &);
  void export_Edge(py::module&);
  void export_HalfEdge(py::module&);
  void export_Face(py::module&);
  void export_Mesh(py::module&);
  void export_System(py::module&);

}

#endif
