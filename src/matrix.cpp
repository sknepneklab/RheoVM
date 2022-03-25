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
 * \file matrix.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Jan-2020
 * \brief 2x2 matrix manipulations
 */ 

#include "matrix.hpp"

namespace RheoVM 
{
  Matrix Matrix::inv() 
  {
    double D = this->det();
    if (fabs(D) <= 1e-7)
      throw runtime_error("Encountered singular matrix.");
    return Matrix(_myy/D, -_mxy/D, -_myx/D, _mxx/D);
  }

  Matrix Matrix::operator*(const Matrix& m)
  {
    return Matrix(_mxx*m._mxx+_mxy*m._myx, 
                  _mxx*m._mxy+_mxy*m._myy,
                  _myx*m._mxx+_myy*m._myx,
                  _myx*m._mxy+_myy*m._myy
                 );
  }

}

