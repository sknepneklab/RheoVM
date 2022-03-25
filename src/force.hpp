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
 * \file force.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 13-Jun-2017
 * \brief Force class 
 */ 

#ifndef __FORCE_HPP__
#define __FORCE_HPP__

#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <exception>
#include <chrono>

#include "system.hpp"
#include "rng.hpp"

double const SMALL_NUMBER = 1e-6;

using std::map;
using std::vector;
using std::string;
using std::runtime_error;
using std::cerr;
using std::exception;
using namespace std::chrono;

namespace RheoVM
{

  // Force on a vertex
  class Force 
  {
    public:
      // Orphaned model options
      // const param_type& gammaLen, const param_type& l0, const param_type& counter
      // _gammaLen(gammaLen), _l0(l0), _counter(counter)
                                                                                                                                
      Force(System& sys) : _sys{sys}, 
                           _use_cell_type{false} 
                           { }
      virtual ~Force() { }
        
      // computes force on vertex by a given edge
      virtual Vec compute(const VertexHandle<Property>&, const HEHandle<Property>&) = 0;  
      
      // compute edge tension for handling 4-vertices (and moves that computation out of integrator)
      // takes care of correct type of force law that way
      virtual double tension(const HEHandle<Property>&, double) = 0;

      // computes stress component due to give force type
      virtual void stress(FaceHandle<Property>&) = 0;

      // compute force on the face centre of a given type
      virtual Vec compute_face_force(FaceHandle<Property> &) = 0;

      // computed energy on of a given face
      virtual double energy(const FaceHandle<Property>&) = 0;
      
      // set all paramters for a given type
      virtual void set_params(const string&, const params_type&) = 0;

      // set all vector-valued paramters for a given type
      virtual void set_vec_params(const string&, const vec_params_type&) = 0;

      // Set relative parameters
      virtual void set_relative_paramters(const string&, const multi_params_type&) = 0;

      // set various compute flags
      virtual void set_flag(const string&) = 0;

      // Copy type parameters for cell paramters
      virtual void copy_type_param_to_cell() = 0;

    
    protected:

      System& _sys;    // Mesh
      bool _use_cell_type;  // if true, paramters are set per cell type (otherwise read them from face data)
      
  };

  
}

#endif
