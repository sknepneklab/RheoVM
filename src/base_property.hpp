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
 * \file base_property.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 13-May-2019
 * \brief BaseProperty class 
 */ 

#ifndef __BASE_PROPERTY_HPP__
#define __BASE_PROPERTY_HPP__

struct BaseProperty 
{
  struct HEProperty
  {
    int he_type = -1;
  };
  struct VertexProperty
  {
    int vert_type = -1;
    VertexProperty& operator=(const VertexProperty& p)
    {
      if (this == &p) 
        return *this;
      this->vert_type = p.vert_type;
      return *this;
    }
  };
  struct EdgeProperty
  {
    int edge_type = -1;
  };
  struct FaceProperty
  {
    int face_type = -1;
    FaceProperty& operator=(const FaceProperty& p)
    {
      if (this == &p) 
        return *this;
      this->face_type = p.face_type;
      return *this;
    }
  };
};



#endif
