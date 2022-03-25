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
 * \file integrator.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief Integrator class 
*/

#ifndef __INTEGRATOR_HPP__
#define __INTEGRATOR_HPP__

#include "rng.hpp"
#include "system.hpp"
#include "force.hpp"

#include "constrainer.hpp"
#include "constraint_none.hpp"
#include "constraint_fixed.hpp"

#include "force_compute.hpp"

#include <chrono>
#include <utility>
#include <map>

using namespace std::chrono;
using std::map;

namespace RheoVM
{

  class Integrator 
  {

    public:
    
      Integrator(System& sys, ForceCompute& fc) : _sys{sys}, 
                                                  _force_compute{fc}, 
                                                  _enabled{true},
                                                  _constraint_enabled{true},
                                                  _dt{0.01}
                                                   { }
      virtual ~Integrator() { }
      
      virtual void step() = 0;
      virtual void set_params(const params_type&) = 0;
      virtual void set_type_params(const string&, const params_type&) = 0;
      virtual void set_external_force(const string&, const Vec&) = 0;
      virtual void set_radial_force(const string &, double) = 0;
      virtual void set_flag(const string&) = 0;
      virtual bool converged() = 0;

      void set_dt(double dt) { _dt = dt; }
      void enable() { _enabled = true; }
      void disable() { _enabled = false; }
      bool is_enabled() { return _enabled; }
      void enable_constraint(bool enable) { _constraint_enabled = enable;  }

    protected:

      System& _sys;              // system
      ForceCompute&   _force_compute; 
      bool _enabled;
      bool _constraint_enabled;
      double _dt; // time step
  };

}

#endif

