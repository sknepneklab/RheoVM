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
 * \file force_area.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-May-2019
 * \brief ForceArea class 
*/ 

#include "force_area.hpp"

namespace RheoVM 
{
  Vec ForceArea::compute(const VertexHandle<Property>& vh, const HEHandle<Property>& he)
  {
    Vec l = he->to()->r - vh->r;                    // vector along the junction pointing away from the vertex
    FaceHandle<Property> fh   = he->face();         // cell to the right of the half edge
    FaceHandle<Property> fh_p = he->pair()->face(); // pair cell (opposite side of the same junction)
    double A1   = _sys.mesh().area(fh);
    double A2   = _sys.mesh().area(fh_p);
    double A0_1 = fh->data().A0;
    double A0_2 = fh_p->data().A0; 
    double kappa_1, kappa_2;
    
    if (_use_cell_type)
    {
      kappa_1 = (fh->outer)   ? 0.0 : _kappa[fh->data().face_type];
      kappa_2 = (fh_p->outer) ? 0.0 : _kappa[fh_p->data().face_type];
    }
    else
    {
      kappa_1 = (fh->outer)   ? 0.0 : fh->data().kappa;
      kappa_2 = (fh_p->outer) ? 0.0 : fh_p->data().kappa;
    }
    
    Vec farea_vec = 0.5*(kappa_1*(A1 - A0_1) - kappa_2*(A2 - A0_2))*l.ez_cross_v();

    return farea_vec;    
  }

  void ForceArea::stress(FaceHandle<Property>& fh)
  {
    double A   = _sys.mesh().area(fh);
    double A0  = fh->data().A0;
    double kappa;
    if (_use_cell_type)
      kappa = (fh->outer)   ? 0.0 : _kappa[fh->data().face_type];
    else
      kappa = (fh->outer)   ? 0.0 : fh->data().kappa;
    double pc = -kappa*(A-A0);
    fh->data().stress[0] += pc;
    fh->data().stress[3] += pc;
    fh->data().stress_a[0] += pc;
    fh->data().stress_a[3] += pc;
  }

  Vec ForceArea::compute_face_force(FaceHandle<Property>& fh)
  {
    if (fh->outer)
      return Vec(0, 0);

    Vec force(0, 0);
    double kappa;
    double A   = _sys.mesh().area(fh);
    double A0  = fh->data().A0;

    if (_use_cell_type)
      kappa = _kappa[fh->data().face_type];
    else
      kappa = fh->data().kappa;

    double PI = -0.5 * (A - A0);

    HEHandle<Property> he = fh->he();
    HEHandle<Property> first = fh->he();
    do
    {
      VertexHandle<Property> vhp = he->prev()->from();
      VertexHandle<Property> vh = he->from();
      VertexHandle<Property> vhn = he->to();
      Vec l1 = vhp->r - vh->r;
      Vec l2 = vhn->r - vh->r;
      Vec l = l1 - l2;
      force += l.ez_cross_v();
      he = he->next();
    } while (he != first);
    return PI * force;
  }

  double ForceArea::energy(const FaceHandle<Property>& fh)
  {
    double A   = _sys.mesh().area(fh);
    double A0 = fh->data().A0; 
    double kappa;
    
    if (fh->outer || fh->erased)
      return 0.0;

    if (_use_cell_type)
      kappa = _kappa[fh->data().face_type];
    else
      kappa = fh->data().kappa;
    
    double dA = A - A0;
    return 0.5*kappa*dA*dA;
  }
}
