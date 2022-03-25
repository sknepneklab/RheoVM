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
 * \file matrix.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Jan-2020
 * \brief 2x2 matrix manipulations
 */ 

#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

#include <vector>
#include <stdexcept>
#include <cmath>

using std::vector;
using std::runtime_error;
using std::fabs;

namespace RheoVM 
{
  /*! Matrix class 
   * This class handles simple 2x2 matrices.
   */
  class Matrix 
  {
    public:
       
      Matrix() : _mxx{0}, _mxy{0}, _myx{0}, _myy{0} { }
      Matrix(double mxx, double mxy, double myx, double myy) : _mxx{mxx}, _mxy{mxy}, _myx{myx}, _myy{myy} { }
      Matrix(const vector<vector<double>>& m)
      {
        if (m.size() != 2)
          throw runtime_error("Matrix m has have 2 rows.");
        if (m[0].size() != 2 || m[1].size() != 2)
          throw runtime_error("Each row of matrix m has to have 2 columns.");
         _mxx = m[0][0];  _mxy = m[0][1];
         _myx = m[1][0];  _myy = m[1][1];
      }

      Matrix(const vector<double>& r1, const vector<double>& r2)
      {
        if (r1.size() != 2 || r2.size() != 2)
          throw runtime_error("Each row of matrix m has to have 2 columns.");
        _mxx = r1[0];  _mxy = r1[1];
        _myx = r2[0];  _myy = r2[1];
      }
      Matrix(const Matrix& m)
      {
        _mxx = m._mxx;   _mxy = m._mxy;
        _myx = m._myx;   _myy = m._myy;
      }

      double trace() { return (_mxx + _myy); }
      double trace() const { return (_mxx + _myy); }
      double det() { return (_mxx*_myy - _mxy*_myx); }
      double det() const { return (_mxx*_myy - _mxy*_myx); }
      Matrix inv();
      Matrix T()
      {
        return Matrix(_mxx, _myx, _mxy, _myy);
      }
      
      Matrix operator*(const Matrix&);

      double _mxx, _mxy, _myx, _myy;
     
  };

}
#endif