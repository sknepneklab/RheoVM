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
 * \file force_self_propulsion.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 12-Jun-2019
 * \brief ForceSelfPropulsion class 
*/ 

#include "force_self_propulsion.hpp"

namespace RheoVM 
{
  Vec ForceSelfPropulsion::compute(const VertexHandle<Property>& vh, const HEHandle<Property>& he)
  {
    
    FaceHandle<Property> fh   = he->face();         // cell to the right of the half edge
    double v0;
    Vec n(0,0);
    
    if (_use_cell_type)
    {
      v0 = (fh->outer)   ? 0.0 : _v0[fh->data().face_type];
      n  = (fh->outer)   ? Vec(0,0) : _n[fh->data().face_type];
    }
    else
    {
      v0 = (fh->outer)   ? 0.0 : fh->data().v0;
      n  = (fh->outer)   ? Vec(0,0) : fh->data().n;
    }

    return (v0/vh->coordination)*n;
  }

  Vec ForceSelfPropulsion::compute_face_force(FaceHandle<Property> &fh) 
  {
    double v0;
    Vec n(0,0);
    
    if (_use_cell_type)
    {
      v0 = (fh->outer)   ? 0.0 : _v0[fh->data().face_type];
      n  = (fh->outer)   ? Vec(0,0) : _n[fh->data().face_type];
    }
    else
    {
      v0 = (fh->outer)   ? 0.0 : fh->data().v0;
      n  = (fh->outer)   ? Vec(0,0) : fh->data().n;
    }

    return v0*n;    
  }
}
