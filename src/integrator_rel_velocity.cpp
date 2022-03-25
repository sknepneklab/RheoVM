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
 * \file integrator_rel_velocity.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com 
 * \date 14-Apr-2021
 * \brief IntegratorRelVelocity class 
*/

#include "integrator_rel_velocity.hpp"

#ifdef HAS_EIGEN3

namespace RheoVM
{
  void IntegratorRelVelocity::step()
  {
    // Check is M needs an update
    if (_sys.topology_change())
      this->_update_matrix();

    // Compute force on each vertex
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        _force_compute.compute(vh);

    // If we are computing stress, add the term due to friction
    if (_force_compute.stress_compute())
      for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
      {
        if (!(fh->erased || fh->outer)) 
        {
          _force_compute.compute_stress(fh);
          HEHandle<Property> he = fh->he();
          HEHandle<Property> first = fh->he();
          Vec rc = _sys.mesh().get_face_centroid(fh);
          vector<double> stress(4,0.0);
          do
          {
            double z = he->from()->coordination;
            Vec ri = he->from()->r - rc;
            Vec fi = he->from()->data().f_type["friction"];
            stress[0] -= ri.x * fi.x / z;
            stress[1] -= 0.5*(ri.x*fi.y + ri.y*fi.x)/z;
            stress[2] -= 0.5*(ri.x*fi.y + ri.y*fi.x)/z;
            stress[3] -= ri.y*fi.y/z;
            he = he->next();
          } while (he != first);
          double fact = 1.0/_sys.mesh().area(fh);
          transform(stress.begin(),stress.end(),stress.begin(),[fact](double v) -> double {return fact*v; });
          transform(fh->data().stress.begin(),fh->data().stress.end(), stress.begin(), fh->data().stress.begin(), std::plus<double>());
          transform(fh->data().stress_v.begin(),fh->data().stress_v.end(), stress.begin(), fh->data().stress_v.begin(), std::plus<double>());
          // to here
        }
      }
    
    // This is actual integrator
    Eigen::VectorXd rhs_x(_M.rows()), rhs_y(_M.rows());
    Eigen::VectorXd fx(_M.rows()), fy(_M.rows());
    Eigen::VectorXd x(_M.rows()), y(_M.rows());
    Eigen::VectorXd x_sol(_M.rows()), y_sol(_M.rows());

    // Populate the right hand side of the system
    int i = 0;
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
    {
      if (!vh->erased)
      {
        Vec f(0, 0);
        
        if (_constraint_enabled && vh->data().constraint == "fixed")
          f = _constant_force[vh->data().vert_type];
        else
          f = vh->data().force;

        x[i] = vh->r.x;
        y[i] = vh->r.y;

        fx[i] = f.x;
        fy[i] = f.y;
        
        i++;
      }
    }

    rhs_x = _M*x + _dt * fx;
    rhs_y = _M*y + _dt * fy;

    x_sol = _solver.solve(rhs_x);
    y_sol = _solver.solve(rhs_y);

    // Update positions
    i = 0;
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
    {
      if (!vh->erased)
      {
        Vec r(x_sol[i], y_sol[i], vh->r.box);
        vh->data().vel = (1.0 / _dt) * (r - vh->r);
        vh->r = r;
        i++;
      }
    }
    
    // compute frictional force
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
    {
      vh->data().f_type["friction"] = -_gamma * vh->data().vel;  //Vec(0, 0); I think minus is OK here.
      HEHandle<Property> he = vh->he();
      HEHandle<Property> first = he;
      do 
      {
        VertexHandle<Property> vn = he->to();
        vh->data().f_type["friction"] += _zeta * (vn->data().vel - vh->data().vel);
        he = he->pair()->next();
      } while (he != first);
    }

  }

  // private methods
  void IntegratorRelVelocity::_update_matrix()
  {
    _vmap.clear();
    int N = 0;
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
      {
        _vmap[vh->id] = N;
        N++;
      }
    
    if (N != _M.rows())
    {
      _M.resize(N, N);
      _M.reserve(6 * N);
    }

    vector<Tripl> elems;

    if (_cell_centres)
    {
      // loop over all cells
      for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
        if (!fh->outer)
        {
          map<pair<int,int>, double> temp_elems;
          int Nc = fh->nsides;
          HEHandle<Property> he = fh->he();
          HEHandle<Property> first = fh->he();
          // lopp over neighbours of current cells
          do
          {
            FaceHandle<Property> fhn = he->pair()->face();
            int Ncn = fhn->nsides;
            HEHandle<Property> he_i = fh->he();
            HEHandle<Property> first_i = fh->he();
            // loop over vertices of the original cell : i-loop
            do 
            {
              VertexHandle<Property> vi = he_i->from();
              int i = vi->id;
              HEHandle<Property> he_j = fh->he();
              HEHandle<Property> first_j = fh->he();
              // loop over vertices of the original cell : first j loop
              do 
              {
                VertexHandle<Property> vj = he_j->from();
                int j = vj->id;
                if (temp_elems.find(make_pair(i,j)) == temp_elems.end())
                {
                  if (i == j)
                    temp_elems[make_pair(i,i)] = _gamma/_sys.mesh().coordination(vi) + _zeta/(Nc*Nc);
                  else 
                    temp_elems[make_pair(i,j)] = _zeta/(Nc*Nc);
                }
                else 
                  temp_elems[make_pair(i,j)] += _zeta/(Nc*Nc);
                he_j = he_j->next();
              } while (he_j != first_j); 
              // loop over vertices of the neigbouring cell : second j loop
              HEHandle<Property> hen_j = fhn->he();
              HEHandle<Property> firstn_j = fhn->he();
              do 
              {
                VertexHandle<Property> vnj = hen_j->from();
                int j = vnj->id;
                if (temp_elems.find(make_pair(i,j)) == temp_elems.end())
                {
                  if (i == j)
                    temp_elems[make_pair(i,i)] = _gamma/_sys.mesh().coordination(vi) - _zeta/(Nc*Ncn);
                  else 
                    temp_elems[make_pair(i,j)] = -_zeta/(Nc*Ncn);
                }
                else 
                  temp_elems[make_pair(i,j)] -= _zeta/(Nc*Ncn);
                hen_j = hen_j->next();
              } while (hen_j != firstn_j); 
              he_i = he_i->next();
            } while (he_i != first_i);
            he = he->next();
          } while (he != first);
          for (const auto& te : temp_elems)
            elems.push_back(Tripl(te.first.first, te.first.second, te.second));
        }
    }
    else
    {
      int i = 0;
      for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
        if (!vh->erased)
        {
          if (_constraint_enabled && vh->data().constraint == "fixed")
            elems.push_back(Tripl(i, i, 1.0));
          else
          {
            elems.push_back(Tripl(i,i, _gamma + _zeta*_sys.mesh().coordination(vh)));
            HEHandle<Property> he = vh->he();
            HEHandle<Property> first = vh->he();
            do
            {
              int j = _vmap[he->to()->id];
              elems.push_back(Tripl(i, j, -_zeta));
              he = he->pair()->next();
            } while (he != first);
          }
          i++;
        }
    }
    _M.setFromTriplets(elems.begin(), elems.end());
    _solver.analyzePattern(_M); 
    // Compute the numerical factorization 
    _solver.factorize(_M); 
  }

}

#endif