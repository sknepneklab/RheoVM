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
 * \file topology.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com 
 * \date 21-May-2019
 * \brief Topology class 
*/

#include "topology.hpp"

namespace RheoVM
{
  void Topology::collapse()
  {
    for (EdgeHandle<Property> eh = _sys.mesh().edges().begin(); eh != _sys.mesh().edges().end(); eh++)
      if (!eh->erased && _sys.mesh().len(eh) < _min_edge_len)
      {
        
        HEHandle<Property>  he = eh->he();
        HEHandle<Property>  he_pair = he->pair();
        VertexHandle<Property> vh = he->from();
        VertexHandle<Property> vh_to = he->to();
	
        
        Vec r = vh_to->r - vh->r;  // vector from i towards j (vh towards vh_to)
        Vec r_faces =  _sys.mesh().get_face_centre(he->face()) - _sys.mesh().get_face_centre(he_pair->face());  // vector joining two centres of cells that share this edge
        
        // This is a hash based on face type that determines how many active cell there are.
        // NOTE: This relies on the assumption that active cells have different cell type.        
        int active_faces = he->face()->data().face_type + he_pair->face()->data().face_type  +
                           he->next()->pair()->face()->data().face_type + he->prev()->pair()->face()->data().face_type;
        T1_stats T1_stats(_sys.time_step(), std::atan2(r.y,r.x), active_faces, r_faces,vh->id, vh_to->id, he->face()->id, he_pair->face()->id, 0.0, 0.0, he->data().l0);        

        if (_sys.mesh().collapse(eh))
        {
          //store pre T1 data in a struct
          vh->data().pre_T1.push_back(T1_stats);
          _sys.set_topology_change(true);
        }
      }
  }

  void Topology::split()
  {
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
    {
      Vec F(0.0,0.0);
      SPLIT_DIRECTION split_direction;
      VertexHandle<Property> vh_new;
      EdgeHandle<Property> eh_new;
      
      if (!vh->erased && !vh->boundary)
        if (this->_unstable_vertex(vh, F, split_direction))
        {
          Vec l = _new_edge_len*F.unit();  
          _sys.mesh().split(vh, l, split_direction, vh_new, eh_new);	  
          HEHandle<Property> he = eh_new->he();
          HEHandle<Property> he_pair = eh_new->he()->pair();

          he->data().old_face_id = he->face()->id; 	  
          he_pair->data().old_face_id = he_pair->face()->id;

          _sys.set_topology_change(true);
        }
    }
  }


  bool Topology::_unstable_vertex(VertexHandle<Property>& vh, Vec& F, SPLIT_DIRECTION& split_dir)
  {
    
    if (vh->coordination < 4) return false;   // by definition, 3-vertex is always stable
    //cout << "Coordination of vertex to split: " << vh->coordination  << endl;
    _sys.mesh().order_star(vh);

    VertexHandle<Property> vh_new;
    EdgeHandle<Property>   eh_new;

    Vec F1(0.0,0.0); //Force on new vertex (except tension)
    Vec F2(0.0,0.0);

    Vec r12(0.0,0.0);
    Vec r34(0.0,0.0);
        
    double t_1_3 = 0.0, t_2_4 = 0.0;
    SPLIT_TYPE split_type_1_3;
    SPLIT_TYPE split_type_2_4;

    // Attempting a 1-3 split of 0 length in the system (topology 1)
    if (_sys.mesh().split(vh, Vec(0.0,0.0), SPLIT_1_3, vh_new, eh_new))
      this->_virtual_split(vh, F1, t_1_3, r12, split_type_1_3);
    
    if (_sys.mesh().split(vh, Vec(0.0,0.0), SPLIT_2_4, vh_new, eh_new))    
      this->_virtual_split(vh, F2, t_2_4, r34, split_type_2_4);
    
    
    // Here is our nice and central bit: 
    // The new edge will be split along the direction of either F1 or F2
    // Therefore, its tension t_1_3 (or t_2_4) will be along this new edge, and so of contractile (opposite) sign to the F1 or F2 *no matter what*
    // Then the scalar force along the new junction is always the one below
    double fa = F1.len() - 2.0*t_1_3;
    double fb = F2.len() - 2.0*t_2_4;

    // If either of them is positive, splitting the vertex is possible
    if (fa > 0 || fb > 0)
    {
      if (fa >= fb)               // in case fa is bigger, do a 13 spllit
      {
        F = F1 - 2.0*t_1_3*F1.unit(); 
        split_dir = SPLIT_1_3;            
        if (dot(F,r12) < 0.0)   // check if the force points in the way of pulling the vertex apart 
          return false;  // reject split that would lead to overlapping edges

        _split_type = split_type_1_3;
        return true;
      }
      else
      {
        F = F2 - 2.0*t_2_4*F2.unit();
        split_dir = SPLIT_2_4;          // and false for a 2 4 split
        if (dot(F,r34) < 0.0) 
          return false; // reject split that would lead to overlapping edges
        _split_type = split_type_2_4;
        return true;
      }
    }
    else // otherwise, keep the 4-vertex!
    {
      return false;
    }
  } 

  void Topology::_virtual_split(VertexHandle<Property>& vh, Vec& F, double& junction_tension, Vec& junction_vec, SPLIT_TYPE& split_type)
  {
          
    HEHandle<Property> he = vh->he(); // edge pointing from vh to vh->new 
    HEHandle<Property> first = he;
    EdgeHandle<Property> eh_new = he->edge();    
    VertexHandle<Property> vh_new=he->to();
  
    double l0 = 0; // l0 = 0 for forward split (this can also be set to minimum edge length)
    
    // Now we compute forces: Loop over the he's of the first vertex, F1 = f1 + f2, and t_1_3 is the edge tension on the infinitesimal junction
    do
    {
      if (he->to()->id != vh_new->id)
      {
        bool force_update = _force_compute.update();
        // if forward split : distribute myo_collapse to the two adjacent junctions
        _force_compute.set_update(false);
        F += _force_compute.compute(vh, he);
        _force_compute.set_update(force_update);
        junction_vec += he->to()->r - vh->r;	  
      }
      else // it must absolutely be the last junction to be reached otherwise myo isn't correctly computed
        junction_tension = _force_compute.tension(he, l0);	  
  
      he = he->pair()->next();
    } while (he != first);
    // Then loop over the other vertex, F1 = f1 + f2 - f3 -f4
    //he = vh_new->he(); 
    he = first->pair(); // I prefer to be sure
    first = he;

    // This is part of the force from other two vertices
    do
    {
      if (he->to()->id != vh->id)
      {
        bool force_update = _force_compute.update();
        _force_compute.set_update(false);
        F -= _force_compute.compute(vh_new, he);
        _force_compute.set_update(force_update);
      }
      he = he->pair()->next();
    } while (he != first);

   
    he = he->pair();
      
    _sys.mesh().collapse(eh_new);
    _sys.mesh().order_star(vh);
          
  }
  

  void Topology::_shape_tensor(FaceHandle<Property>& fh, double& sxx, double& sxy, double& syy)
  {
    sxx = 0.0; sxy = 0.0; syy = 0.0;
    Vec rc = _sys.mesh().get_face_centre(fh);
    HEHandle<Property> he = fh->he();
    HEHandle<Property> first = he;
    do 
    {
      Vec ri = he->from()->r;
      Vec dr = ri - rc;
      sxx += dr.x*dr.x;  sxy += dr.x*dr.y;   syy += dr.y*dr.y;
    } while ((he = he->next()) != first);
    sxx /= fh->nsides;
    sxy /= fh->nsides;
    syy /= fh->nsides;    
  }

  void Topology::_split_direction(FaceHandle<Property>& fh, double& t1, HEHandle<Property>& he1, double& t2, HEHandle<Property>& he2)
  {
    double sxx, sxy, syy;
    this->_shape_tensor(fh, sxx, sxy, syy);
    double D = sqrt(sxx*sxx + 4*sxy*sxy - 2*sxx*syy + syy*syy);
    double e1 = fabs(sxy) < 1e-10 ?  sxx : 0.5*(sxx + syy - D);
    double e2 = fabs(sxy) < 1e-10 ?  syy : 0.5*(sxx + syy + D);
    Vec v1, v2;
    if (fabs(sxy) < 1e-10)
    {
      v1 = Vec(0,1);
      v2 = Vec(1,0);
    }
    else
    {
      v1 = Vec(-(-sxx + syy + D)/(2*sxy),1);
      v2 = Vec(-(-sxx + syy - D)/(2*sxy),1);
    }
    
    double ca, sa;
    if (e1 < e2)
    {
      ca = v1.x;
      sa = v1.y;
    }
    else
    {
      ca = v2.x;
      sa = v2.y;
    }

    //Vec rc = _sys.mesh().get_face_centre(fh);
    Vec rc = _sys.mesh().get_face_centroid(fh);
    bool found = false;
    HEHandle<Property> he = fh->he();
    HEHandle<Property> first = he;
    do 
    {
      Vec r1 = he->from()->r;
      Vec r2 = he->to()->r;
      double denom = ca*(r2.y - r1.y) - sa*(r2.x - r1.x);
      double t = fabs(denom) >= 1e-10 ? (sa*(r1.x - rc.x) - ca*(r1.y - rc.y))/denom : -1.0; // Guard against div by zero
      if (fabs(t) < 1e-10) t = 0;
      else if (fabs(t-1.0) < 1e-10) t = 1.0;
      //cout << "denom = " << denom << "  t = " << t << endl;
      if ((t >= 0) && (t < 1.0))
      {
        if (!found)
        {
          t1 = t;
          he1 = he;
          found = true;
          //cout << "found first  one" << endl;
        }
        else
        {
          t2 = t;
          he2 = he;
          //cout << "found second one" << endl;
        }
      }
    } while ((he = he->next()) != first);
    //cout << "+++++++++++++" << endl;
    //cout << "t1 = " << t1 << " t2 = " << t2 << endl;
  }

  void export_Topology(py::module& m)
  {
    py::class_<Topology>(m, "Topology")
      .def(py::init<System&, ForceCompute&>())
      .def("set_params", &Topology::set_params)
      .def("set_type_params", &Topology::set_type_params)
      .def("set_flag", &Topology::set_flag);
  }
  
}
