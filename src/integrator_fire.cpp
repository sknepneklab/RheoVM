/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

/*!
 * \file integrator_fire.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-Feb-2014
 * \brief Implementation of the FIRE integrator.
 */ 

#include "integrator_fire.hpp"

/*! This is the FIRE minimiser.  
**/
namespace RheoVM
{
  void IntegratorFIRE::step()
  {
    if (_converged)
      return; 
      
    double P = 0.0;
    double Fnorm = 0.0;
    double Vnorm = 0.0;
    double vx, vy, vz;
    double ax, ay, az;

    int N = _sys.mesh().vertices().size();
    double sqrt_ndof = sqrt(2*N);
    double dt_2 = 0.5*_dt;
    
    // Perform first half step for velocity
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
      {
        Vec f = vh->data().force;
        if (_constraint_enabled) 
          f = _constrainer->apply_vertex(vh, f);
        vh->data().vel += dt_2 * f;  //vh->data().vel += dt_2 * vh->data().force;
        vh->r += _dt*vh->data().vel;
      }

    // Compute forces    
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        _force_compute.compute(vh);

    // Perform second half step for velocity 
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
      {
        Vec f = vh->data().force;
        if (_constraint_enabled) 
          f = _constrainer->apply_vertex(vh, f);
        vh->data().vel += dt_2*f;  //vh->data().vel += dt_2*vh->data().force;
        // Compute FIRE related quantities
        P += dot(vh->data().vel, vh->data().force);
        Fnorm += vh->data().force.len2();
        Vnorm += vh->data().vel.len2();
      }

    Fnorm = sqrt(Fnorm);
    Vnorm = sqrt(Vnorm);  

    double E = 0.0;
    for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
      E += _force_compute.energy(fh);

    if (Fnorm/sqrt_ndof < _F_tol && fabs(E - _old_energy) < _E_tol)
      _converged = true;
    else
      _converged = false;
      
    
    double inv_Fnorm = 1.0/Fnorm;
    double fact_1 = 1.0 - _alpha;
    double fact_2 = _alpha * inv_Fnorm * Vnorm;

    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
      {
        Vec f = vh->data().force;
        if (_constraint_enabled) 
          f = _constrainer->apply_vertex(vh, f);
        vh->data().vel = fact_1*vh->data().vel + fact_2*f; //vh->data().vel = fact_1*vh->data().vel + fact_2*vh->data().force;
      }
    
    if (P > 0.0)
    {
      _last_neg++;
      if (_last_neg > _N_min)
      {
        _dt = min(_dt*_f_inc, _dt_max);
        _alpha *= _f_alpha;
      }
    }
    else
    {
      _dt *= _f_dec;
      _alpha = _alpha_init;
      _last_neg = 0;
      for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        vh->data().vel = Vec(0.0, 0.0);
    }
   
    _old_energy = E;
  }
}
