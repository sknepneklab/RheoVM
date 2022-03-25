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
 * \file vec.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 08-Jun-2017
 * \brief 2D vectors 
 */ 

#ifndef __VEC_HPP__
#define __VEC_HPP__

#include <cmath>
#include <memory>
#include <iostream>

#include <pybind11/pybind11.h>

#include "box.hpp"

using std::sqrt;
using std::sin;
using std::cos;
using std::shared_ptr;
using std::ostream;
using std::cout;
using std::endl;

namespace py = pybind11;

namespace RheoVM 
{
  /*! Vec class
  *  Handles vectors in 2d Euclidean space
  */
  class Vec
  {
  public:
    
    //! Default constructor
    Vec() : box{nullptr}, x{0.0}, y{0.0}, ix{0}, iy{0} { }
    //! Constructor for a vector object
    Vec(double x, double y) : box{nullptr}, x{x}, y{y}, ix{0}, iy{0} { }
    //! Constructor for a vector object
    Vec(double x, double y, const shared_ptr<Box>& box) : box{box}, ix{0}, iy{0}
    { 
      if (box)
      {
        double s_x = box->inv_h._mxx*x + box->inv_h._mxy*y;
        double s_y = box->inv_h._myx*x + box->inv_h._myy*y; 
        s_x -= rint(s_x);
        s_y -= rint(s_y);
        this->x = box->h._mxx*s_x + box->h._mxy*s_y;
        this->y = box->h._myx*s_x + box->h._myy*s_y;
      }
      else
      {
        this->x = x;
        this->y = y;
      }
      
    }
    //! Constructor for a vector object
    Vec(double x, double y, int ix, int iy, const shared_ptr<Box>& box) : x(x), y(y), ix(ix), iy(iy), box(box) { }
    //! Copy constructor  
    Vec(const Vec& v) : box(v.box) { x = v.x; y = v.y; ix = v.ix; iy = v.iy; }
    //! Assignment operator
    Vec& operator=(const Vec& rhs)
    {
      x = rhs.x;
      y = rhs.y;
      ix = rhs.ix;
      iy = rhs.iy;
      box = rhs.box;
      return *this;
    }
    
    //! Add two vectors
    Vec operator+(const Vec& v)
    {
      double xx = x + v.x, yy = y + v.y; 
      if (box || v.box)
      {
        const shared_ptr<Box>& sim_box = (box) ? box : v.box;
        Vec si = sim_box->inv_h*(*this);
        Vec sj = sim_box->inv_h*v;
        double sij_x = si.x + sj.x;
        double sij_y = si.y + sj.y;
        sij_x -= rint(sij_x);
        sij_y -= rint(sij_y);
        xx = sim_box->h._mxx*sij_x + sim_box->h._mxy*sij_y;
        yy = sim_box->h._myx*sij_x + sim_box->h._myy*sij_y;
      } 
      return Vec(xx, yy, box);
    }

    //! Add two vectors (constant version)
    Vec operator+(const Vec& v) const { return *(const_cast<Vec*>(this)) + v; }

    //! Subtract two vectors
    Vec operator-(const Vec& v)
    {
      double xx = x - v.x, yy = y - v.y;
      if (box || v.box)
      {
        const shared_ptr<Box>& sim_box = (box) ? box : v.box;
        Vec si = sim_box->inv_h*(*this);
        Vec sj = sim_box->inv_h*v;
        double sij_x = si.x - sj.x;
        double sij_y = si.y - sj.y;
        sij_x -= rint(sij_x);
        sij_y -= rint(sij_y);
        xx = sim_box->h._mxx*sij_x + sim_box->h._mxy*sij_y;
        yy = sim_box->h._myx*sij_x + sim_box->h._myy*sij_y;
      } 
      return Vec(xx, yy, box);
    }   

    //! Subtract two vectors (constant version)
    Vec operator-(const Vec& v) const { return *(const_cast<Vec*>(this)) - v; }
    
    //! Negate a vector
    Vec operator-()
    {
      return Vec(-x, -y, -ix, -iy, box);
    }
    
    //! Scale vector by a constant
    Vec operator*(const double c)
    {
      double xx = c*x, yy = c*y;
      if (box)
      {
        Vec s = box->inv_h*(*this);
        double s_x = c*s.x;
        double s_y = c*s.y;
        s_x -= rint(s_x);
        s_y -= rint(s_y);
        xx = box->h._mxx*s_x + box->h._mxy*s_y;
        yy = box->h._myx*s_x + box->h._myy*s_y;
      } 
      return Vec(xx, yy, box);
    }

    //! Scale vector by a constant (const version)
    Vec operator*(const double c) const { return (*(const_cast<Vec*>(this)))*c;  }
    
    //! Test equality
    bool operator==(const Vec& v)
    {
      return (x == v.x && y == v.y);
    }
    
    //! Add vector to current vector
    Vec& operator+=(const Vec& v)
    {
      double xx = x + v.x, yy = y + v.y;
      if (box || v.box)
      {
        const shared_ptr<Box>& sim_box = (box) ? box : v.box;
        Vec si = sim_box->inv_h*(*this);
        Vec sj = sim_box->inv_h*v;
        double sij_x = si.x + sj.x;
        double sij_y = si.y + sj.y;
        sij_x -= rint(sij_x);
        sij_y -= rint(sij_y);
        xx = sim_box->h._mxx*sij_x + sim_box->h._mxy*sij_y;
        yy = sim_box->h._myx*sij_x + sim_box->h._myy*sij_y;
      } 
      x = xx;  y = yy; 
      return *this;
    }
    
    //! Subtract vector from current vector
    Vec& operator-=(const Vec& v)
    {
      double xx = x - v.x, yy = y - v.y;
      if (box || v.box)
      {
        const shared_ptr<Box>& sim_box = (box) ? box : v.box;
        Vec si = sim_box->inv_h*(*this);
        Vec sj = sim_box->inv_h*v;
        double sij_x = si.x - sj.x;
        double sij_y = si.y - sj.y;
        sij_x -= rint(sij_x);
        sij_y -= rint(sij_y);
        xx = sim_box->h._mxx*sij_x + sim_box->h._mxy*sij_y;
        yy = sim_box->h._myx*sij_x + sim_box->h._myy*sij_y; 
      } 
      x = xx;  y = yy; 
      return *this;
    }
    
    //! Euclidean dot product with another vector
    double dot(const Vec& v)
    {
      return x*v.x + y*v.y;
    }
    
    //! Vector length 
    double len() const { return sqrt(x*x + y*y); }
    
    //! Vector length squared
    double len2() { return x*x + y*y; }
    
    //! Make the vector has unit length
    void normalise()
    {
      double len = this->len();
      if (len != double(0))
      {
        x /= len; y /= len;
      }
    }
    
    //! Return unit vector in the direction of this vector
    Vec unit()
    {
      double len = this->len();
      if (len != double(0))
        return Vec(x/len,y/len,box);
      return Vec(x,y,box);
    }

    Vec unit() const
    {
      double len = this->len();
      if (len != double(0))
        return Vec(x/len,y/len,box);
      return Vec(x,y,box);
    }
    
    //! Rotate vector by \f$ \phi \f$
    Vec rotate(const double phi)
    {
      double s = sin(phi);
      double c = cos(phi);
      
      double xx = c*x - s*y;
      double yy = s*x + c*y;
      
      if (box)
      {
        Vec s = box->inv_h*(*this);
        double s_x = s.x;
        double s_y = s.y;
        s_x -= rint(s_x);
        s_y -= rint(s_y);
        xx = box->h._mxx*s_x + box->h._mxy*s_y;
        yy = box->h._myx*s_x + box->h._myy*s_y; 
      } 
      
      return Vec(xx,yy,box);
    }

    //! Compute e_z x v (used for force)
    Vec ez_cross_v()   {  return Vec(-y,x,box);   }

    Vec ez_cross_v() const {  return Vec(-y,x,box);   }

    //! Fold back vector into the periodic box (used when changing simulation box to make sure all vertices are inside the box)
    void fold_back() 
    {
      if (box)
      {
        double s_x = box->inv_h._mxx*x + box->inv_h._mxy*y;
        double s_y = box->inv_h._myx*x + box->inv_h._myy*y; 
        s_x -= rint(s_x);
        s_y -= rint(s_y);
        this->x = box->h._mxx*s_x + box->h._mxy*s_y;
        this->y = box->h._myx*s_x + box->h._myy*s_y;
      }
    }

    // Scale vector 
    void scale(double a, double b)
    {
      if (box)
      {
        Vec s = box->inv_h*(*this);
        double s_x = a*s.x;
        double s_y = b*s.y; 
        s_x -= rint(s_x);
        s_y -= rint(s_y);
        this->x = box->h._mxx*s_x + box->h._mxy*s_y;
        this->y = box->h._myx*s_x + box->h._myy*s_y;
      }
      else
      {
        this->x = a*x;
        this->y = b*y;
      }
    }

    friend Vec operator*(const Matrix&, const Vec&);

    ///@{
    shared_ptr<Box>  box;
    double x, y;              //!< Position in the embedding 3d flat space
    int ix, iy;               //!< Image flags in periodic boundary conditions 
    //@}

  };

  
  //! Scale vector by a number
  inline Vec operator*(const double c, const Vec& v)
  {
    double xx = c*v.x, yy = c*v.y;
    if (v.box)
    {
      Vec s = v.box->inv_h*v;
      double s_x = c*s.x;
      double s_y = c*s.y;
      s_x -= rint(s_x);
      s_y -= rint(s_y);
      xx = v.box->h._mxx*s_x + v.box->h._mxy*s_y;
      yy = v.box->h._myx*s_x + v.box->h._myy*s_y;
    } 
    return Vec(xx,yy,v.box);
  }

  inline Vec operator*(const Matrix& m, const Vec& v)
  {
    return Vec(m._mxx*v.x+m._mxy*v.y,m._myx*v.x+m._myy*v.y, v.box);
  }
      
  inline double dot(const Vec& v1, const Vec& v2)
  {
    return (v1.x*v2.x + v1.y*v2.y);
  }

  inline ostream& operator<<(ostream& os, const Vec& v)
  {
    os << "(" << v.x << "," << v.y << ")";
    return os;
  }

  void export_Vec(py::module& m);

}

#endif
