/* ***************************************************************************
 *
 *  Copyright (C) 2020 University of Dundee
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
 * \file integrator_fore.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Sep-2020
 * \brief IntegratorFIRE class 
*/

#ifndef __INTEGRATOR_FIRE_HPP__
#define __INTEGRATOR_FIRE_HPP__

#include <cmath>

#include "constrainer.hpp"
#include "constraint_none.hpp"
#include "constraint_fixed.hpp"

#include "integrator.hpp"

using std::sqrt;
using std::fabs;
using std::min;
using std::cout;
using std::endl;

namespace RheoVM
{

  using ConstrainerType = unique_ptr<Constrainer>;

  /*! IntegratorFIRE class handles FIRE minimization. 
  *  \note No activity. Just minimization. 
  */
  class IntegratorFIRE : public Integrator
  {
  public:
    
    //! Constructor
    IntegratorFIRE(System& sys, ForceCompute& fc, int seed) : Integrator{sys, fc},
                                                              _converged{false},
                                                              _old_energy{1e15},
                                                              _last_neg{0},
                                                              _alpha{0.1},
                                                              _f_alpha{0.99},  
                                                              _F_tol{1e-4},
                                                              _E_tol{1e-6},        
                                                              _f_inc{1.1},
                                                              _f_dec{0.5},
                                                              _dt_max{1.0},
                                                              _N_min{5}

    {
      _constrainer = make_unique<Constrainer>();
      _constrainer->add<ConstraintNone>("none");
      _constrainer->add<ConstraintFixed>("fixed");
    }
    
    void step() override;
    void set_params(const params_type& params) override 
    { 
      for (auto& p : params)
      {
        if (p.first == "alpha")
          _alpha = p.second;
        if (p.first == "alpha_scale")
          _f_alpha = p.second;
        if (p.first == "f_tolerance")
          _F_tol = p.second;
        if (p.first == "e_tolerance")
          _E_tol = p.second;
        if (p.first == "dt_inc")
          _f_inc = p.second;
        if (p.first == "dt_dec")
          _f_dec = p.second;
        if (p.first == "dt_max")
          _dt_max = p.second;
        if (p.first == "min_neg_steps")
          _N_min = static_cast<int>(p.second);   
      }
    };
    void set_type_params(const string& type, const params_type& params) override { }
    void set_external_force(const string& vtype, const Vec& f) override { }
    void set_radial_force(const string& vtype, double f) override { }
    void set_flag(const string& flag) override {  }
    bool converged() override { return _converged; }
    
  private:

    double _alpha;                                   //!< alpha factor in the FIRE minimizer
    double _F_tol;                                   //!< Force tolerance for checking convergence 
    double _E_tol;                                   //!< Energy tolerance for checking convergence 
    int    _N_min;                                   //!< Number of steps energy change is negative before allowing m_alpha and _dt to adapt
    double _f_inc;                                   //!< Increase factor for _dt
    double _f_dec;                                   //!< Decrease factor for _dt
    double _alpha_init;                              //!< Initial value for m_alpha
    double _f_alpha;                                 //!< Decrease alpha by this much 
    int    _last_neg;                                //!< Last time (in terms of step number when force.velocity was negative
    double _old_energy;                              //!< Old value of energy
    bool   _converged;                               //!< Flag which tests if the method has converged
    double _dt_max;                                  //!< Maximum time step
    ConstrainerType _constrainer; // Apply various constraints

    
  };

  
}

#endif
