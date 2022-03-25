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
 * \file topology.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com 
 * \date 21-May-2019
 * \brief Topology class 
*/


#include "system.hpp"
#include "force_compute.hpp"

#include "rng.hpp"

#include <exception>
#include <tuple>
#include <chrono>

using std::cerr;
using std::exception;
using std::make_tuple;
using std::tuple;
using namespace std::chrono;

namespace RheoVM
{
  class Topology
  {
    public:

      Topology(System& sys, ForceCompute& fc) : _sys{sys},
                                                _force_compute{fc},
                                                _min_edge_len{0.02},
                                                _new_edge_len{0.022},
                                                _rng{RNG(system_clock::now().time_since_epoch().count())}
                                                {  
                  
                                                }
      ~Topology() = default;

      void set_params(const params_type& params)  
      { 
        for (auto& p : params)
        {
          if (p.first == "min_edge_len")
            _min_edge_len = p.second;
          else if (p.first == "new_edge_len")
            _new_edge_len = p.second; 
          else if (p.first == "seed")
          {
            if (p.second >= 0)
              _rng = RNG(static_cast<int>(p.second));
            else
              throw runtime_error("Random Number Generator seed has to be positive.");
          }
        }
      };

      void set_type_params(const string& type, const params_type& params)  
      { 
        if (_sys.cell_types().find(type) == _sys.cell_types().end())
          throw runtime_error("Cell type " + type + " is not defined.");
      }

      void set_flag(const string& flag) 
      {  
      }

      void collapse();
      void split();
      

    private:

      System& _sys;
      ForceCompute&   _force_compute; 
      double _min_edge_len;
      double _new_edge_len;
      SPLIT_TYPE _split_type; 
      
      RNG _rng;

      // Auxiliary functions
      bool _unstable_vertex(VertexHandle<Property>&, Vec&, SPLIT_DIRECTION&);          // checks if a four-vertex is stable
      void _virtual_split(VertexHandle<Property>&, Vec&, double&, Vec&, SPLIT_TYPE&);
      void _shape_tensor(FaceHandle<Property>&, double&, double&, double&);
      void _split_direction(FaceHandle<Property>&, double&, HEHandle<Property>&, double&, HEHandle<Property>&);

  };

  void export_Topology(py::module&);

}

