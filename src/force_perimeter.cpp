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
 * \file force_perimeter.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-May-2019
 * \brief ForceArea class 
*/ 

#include "force_perimeter.hpp"

namespace RheoVM 
{
  Vec ForcePerimeter::compute(const VertexHandle<Property>& vh, const HEHandle<Property>& he)
  {
    Vec l = he->to()->r - vh->r;                    // vector along the junction pointing away from the vertex
    FaceHandle<Property> fh   = he->face();         // cell to the right of the half edge
    FaceHandle<Property> fh_p = he->pair()->face(); // pair cell (opposite side of the same junction)
    double P1 = _sys.mesh().perim(fh);
    double P2 = _sys.mesh().perim(fh_p);
    double gamma_1, gamma_2;
    double lambda_1, lambda_2;
    
    if (_use_cell_type)
    {
      gamma_1 = (fh->outer)    ? 0.0 : _gamma[fh->data().face_type];
      gamma_2 = (fh_p->outer)  ? 0.0 : _gamma[fh_p->data().face_type];
      if (!_lambda_P0)
      {
        lambda_1 = (fh->outer)   ? 0.0 : _lambda[fh->data().face_type];
        lambda_2 = (fh_p->outer) ? 0.0 : _lambda[fh_p->data().face_type];
      }
    }
    else
    {
      gamma_1  = (fh->outer)   ? 0.0 : fh->data().gamma;
      gamma_2  = (fh_p->outer) ? 0.0 : fh_p->data().gamma;
      if (!_lambda_P0)
      {
        lambda_1 = (fh->outer)   ? 0.0 : fh->data().lambda;
        lambda_2 = (fh_p->outer) ? 0.0 : fh_p->data().lambda;
      }
    }

    if (_lambda_P0)
    {
      lambda_1 = gamma_1*fh->data().P0;
      lambda_2 = gamma_2*fh_p->data().P0;
    }

    if (_double_boundary_constants)
    {
      gamma_1  = (fh_p->outer) ? 2*gamma_1 : gamma_1;
      gamma_2  = (fh->outer)   ? 2*gamma_2 : gamma_2;
      lambda_1 = (fh_p->outer) ? 2*lambda_1 : lambda_1;
      lambda_2 = (fh->outer)   ? 2*lambda_2 : lambda_2;
    }

    double lambda = lambda_1 + lambda_2;
    double fedges = gamma_1*P1 + gamma_2*P2 - lambda;

    return fedges*l.unit();
  }

  double ForcePerimeter::tension(const HEHandle<Property>& he, double l0)
  {
    FaceHandle<Property> fh   = he->face();         // cell to the right of the half edge
    FaceHandle<Property> fh_p = he->pair()->face(); // pair cell (opposite side of the same junction)

    double gamma_1, gamma_2;
    double lambda_1, lambda_2;
    
    if (_use_cell_type)
    {
      gamma_1 = (fh->outer)    ? 0.0 : _gamma[fh->data().face_type];
      gamma_2 = (fh_p->outer)  ? 0.0 : _gamma[fh_p->data().face_type];
      if (!_lambda_P0)
      {
        lambda_1 = (fh->outer)   ? 0.0 : _lambda[fh->data().face_type];
        lambda_2 = (fh_p->outer) ? 0.0 : _lambda[fh_p->data().face_type];
      }
    }
    else
    {
      gamma_1  = (fh->outer)   ? 0.0 : fh->data().gamma;
      gamma_2  = (fh_p->outer) ? 0.0 : fh_p->data().gamma;
      if (!_lambda_P0)
      {
        lambda_1 = (fh->outer)   ? 0.0 : fh->data().lambda;
        lambda_2 = (fh_p->outer) ? 0.0 : fh_p->data().lambda;
      }
    }

    if (_lambda_P0)
    {
      lambda_1 = gamma_1*fh->data().P0;
      lambda_2 = gamma_2*fh_p->data().P0;
    }

    if (_double_boundary_constants)
    {
      gamma_1  = (fh_p->outer) ? 2*gamma_1 : gamma_1;
      gamma_2  = (fh->outer)   ? 2*gamma_2 : gamma_2;
      lambda_1 = (fh_p->outer) ? 2*lambda_1 : lambda_1;
      lambda_2 = (fh->outer)   ? 2*lambda_2 : lambda_2;
    }

    double lambda = lambda_1 + lambda_2;
  
    return gamma_1*_sys.mesh().perim(fh) + gamma_2*_sys.mesh().perim(fh_p) - lambda;
  }

  void ForcePerimeter::stress(FaceHandle<Property>& fh)
  {
    double P = _sys.mesh().perim(fh);
    double A = _sys.mesh().area(fh);
    double gamma, lambda;
    if (_use_cell_type)
    {
      gamma = (fh->outer)    ? 0.0 : _gamma[fh->data().face_type];
      if (!_lambda_P0)
        lambda = (fh->outer)   ? 0.0 : _lambda[fh->data().face_type];
    }
    else
    {
      gamma  = (fh->outer)   ? 0.0 : fh->data().gamma;
      if (!_lambda_P0)
        lambda = (fh->outer)   ? 0.0 : fh->data().lambda;
    }
    if (_lambda_P0)
      lambda = gamma*fh->data().P0;

    double fact = (gamma*P - lambda)/A;
    HEHandle<Property> he = fh->he();
    HEHandle<Property> first = fh->he();
    vector<double> stress(4, 0.0);
    do
    {
      VertexHandle<Property> vh = he->from();
      Vec r = vh->r;
      Vec rn = he->to()->r;
      Vec dr = rn - r;
      double l = dr.len();
      stress[0] -= dr.x*dr.x/l;
      stress[1] -= dr.x*dr.y/l;
      stress[2] -= dr.y*dr.x/l;
      stress[3] -= dr.y*dr.y/l;
      he = he->next();
    } while (he != first);
    transform(stress.begin(),stress.end(),stress.begin(),[fact](double v) -> double { return fact*v;});
    transform(fh->data().stress.begin(), fh->data().stress.end(), stress.begin(), fh->data().stress.begin(), std::plus<double>());
    transform(fh->data().stress_p.begin(), fh->data().stress_p.end(), stress.begin(), fh->data().stress_p.begin(), std::plus<double>());
  }

  Vec ForcePerimeter::compute_face_force(FaceHandle<Property>& fh)
  {
    if (fh->outer)
      return Vec(0, 0);

    Vec force(0, 0);
    double gamma, lambda;
    double P = _sys.mesh().perim(fh);

    if (_use_cell_type)
    {
      gamma = _gamma[fh->data().face_type];
      if (!_lambda_P0)
        lambda = _lambda[fh->data().face_type];
    }
    else
    {
      gamma  = fh->data().gamma;
      if (!_lambda_P0)
        lambda = fh->data().lambda;
    }

    if (_lambda_P0)
      lambda = gamma*fh->data().P0;

    double tc = gamma * P - lambda;

    HEHandle<Property> he = fh->he();
    HEHandle<Property> first = fh->he();
    do
    {
      VertexHandle<Property> vhp = he->prev()->from();
      VertexHandle<Property> vh = he->from();
      VertexHandle<Property> vhn = he->to();
      Vec l1 = vhp->r - vh->r;
      Vec l2 = vhn->r - vh->r;
      force += l1.unit() + l2.unit();
      he = he->next();
    } while (he != first);
    return tc * force;
  }

  double ForcePerimeter::energy(const FaceHandle<Property>& fh)
  {
    if (fh->outer || fh->erased)
      return 0.0;

    double P = _sys.mesh().perim(fh);
    
    double gamma, lambda = 0.0;
    if (_use_cell_type)
    {
      gamma = _gamma[fh->data().face_type];
      if (!_lambda_P0)
        lambda = _lambda[fh->data().face_type];
    }
    else
    {
      gamma = fh->data().gamma;
      if (!_lambda_P0)
        lambda = fh->data().lambda;
    }

    double P0;
    if (_lambda_P0)
      P0 = fh->data().P0;
    else
      P0 = lambda/gamma;
    double dP = P - P0;

    return 0.5*gamma*dP*dP;
    
  }

}
