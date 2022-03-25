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
 * \file box.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 08-Jun-2017
 * \brief Simulation Box 
 */ 

#ifndef __BOX_HPP__
#define __BOX_HPP__

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "matrix.hpp"

using std::shared_ptr;

namespace py = pybind11;

namespace RheoVM
{

struct Box 
{
  Box(double lx, double ly) : h{lx,0.0,0.0,ly}, inv_h{h.inv()}, s_mat{1.0,0.0,0.0,1.0} { }
  Box(double ax, double ay, double bx, double by) : h{ax,bx,ay,by}, inv_h(h.inv()), s_mat{1.0,0.0,0.0,1.0} { }
  double area() { return h.det(); }
  vector<double> get() { return vector<double>({h._mxx,h._mxy,h._myx,h._myy}); }
  void stretch(double a, double b)
  {
    Matrix s(a,0.0,0.0,b);
    h = s*h;
    inv_h = h.inv();
    s_mat = s;
  }
  void shear(double gamma)
  {
    Matrix s(1.0,gamma,gamma,1.0);
    h = s*h; 
    inv_h = h.inv();
    s_mat = s;
  }
  void transform(double txx, double txy, double tyx, double tyy, bool undo = false)
  {
    Matrix s(txx, txy, tyx, tyy);
    if (undo)
      h = s_mat.inv()*h;
    h = s*h; 
    inv_h = h.inv();
    s_mat = s;
  }

  Matrix h;
  Matrix inv_h;
  Matrix s_mat;
  /*
  Box(double lx, double ly) : lx(lx), ly(ly) 
  { 
    xlo = -0.5*lx;   xhi = 0.5*lx;
    ylo = -0.5*ly;   yhi = 0.5*ly;
  }
  double xlo, ylo;
  double xhi, yhi;
  double lx, ly;
  */ 
};

void export_Box(py::module& m);

}

#endif