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
 * \file property.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief Property class 
 */ 

#ifndef __PROPERTY_HPP__
#define __PROPERTY_HPP__

#include <map>
#include <vector>
#include <string> 

#include "vec.hpp"

#include "base_property.hpp"

using std::map;
using std::vector;
using std::string;

namespace RheoVM
{

  struct T1_stats
  {
    T1_stats(int ts, double a, int active_faces, Vec r_faces, int v1, int v2, int f1, int f2, double m1, double m2, double l00) : time_step(ts), 
                                                                                                                                  angle(a), 
                                                                                                                                  v1_id(v1), 
                                                                                                                                  v2_id(v2), 
                                                                                                                                  f1_id(f1), 
                                                                                                                                  f2_id(f2), 
                                                                                                                                  l00(l00),
                                                                                                                                  r_faces(r_faces),
                                                                                                                                  active_faces(active_faces) 
                                                                                                                                  {  }
    int time_step;     // time step at which the transition happened
    double angle;      // angle with the x-axis of the junction before/after the T1
    int v1_id, v2_id;  // ids of vertces that belong to the junction
    int f1_id, f2_id;  // ids of faces that belong to the junction
    double l00;         // l0 at collapse
    Vec r_faces; // vector joining the two face centre of the faces
    int active_faces; // number of faces which are active
  };


  struct Property : public BaseProperty 
  {  
    struct HEProperty : public BaseProperty::HEProperty
    {
      double tension = 0.0;  
      double l0;           // Native length of an edge
      double l0_collapse;   // l0 value before the collapse
      int old_face_id;      // the way to distinguish if a vertex split was actual T1 or a bounce back
      map<string, Vec> force_type;  // Force of a given type on the he->from() vertex due to along this half-edge
    };
    struct VertexProperty : public BaseProperty::VertexProperty
    {
      Vec vel;
      Vec force;
      string constraint;        // if "x" move only along x-axis, if "y" move only along y axis, if "radial", move along radius; otherwise ignore
      map<string, Vec> f_type;  // Force from a given interaction type (e.g., area term, perimeter term, etc)  
      vector<T1_stats> pre_T1;
      vector<T1_stats> post_T1;
      map<int,Vec>  he_force;
      map<int,int>  pre_T1_face_pair;
      map<int,int>  post_T1_face_pair;  
      string type_name;     // String with the actual name of the the vertex type 
      VertexProperty& operator=(const VertexProperty& p) 
      {
        if (this == &p)
          return *this;
        this->vert_type = p.vert_type;
        this->type_name = p.type_name;
        return *this;
      }
    };
    struct EdgeProperty : public BaseProperty::EdgeProperty
    {
      double tension = 0.0;  
      double l0;            // Native length of an edge
      double l0_collapse;   // l0 value before the collapse
    };
    struct FaceProperty : public BaseProperty::FaceProperty
    {
      double A0;
      double P0;
      double native_A0;
      double native_P0;
      double max_A0;       // Maximum A0 after which we divide 
      int original_face;    // one of the faces collapsed edge belonged to 
      double fa;          // self-propulsion (activity) for a given cell
      double phi;         // direction of the self-propulsion vector  
      double kappa;       // area modulus 
      double gamma;       // perimeter modulus 
      double lambda;      // line tension 
      double beta;        // activity     
      double alpha;       // activity in the alpha*l term
      double beta_a;      // area activity 
      double k;           // half-edge stiffness for harmonic force (total k = k[face_1] + k[face_2])
      double k_angle;     // stiffness constant for angle penalty
      double v0;          // Self-propulsion velocity on the cell centre
      int psi;            // Power of the term that couples beta to the cell's long axis direction
      Vec n;              // Self-propulsion direction
      vector<double> stress = {0.0, 0.0, 0.0, 0.0};    // Stress tensor, element 0 -> xx, element 1 -> xy, element 2 -> yx, element 3 -> yy
      vector<double> stress_a = {0.0, 0.0, 0.0, 0.0};    // Stress tensor, area component, element 0 -> xx, element 1 -> xy, element 2 -> yx, element 3 -> yy
      vector<double> stress_p = {0.0, 0.0, 0.0, 0.0};    // Stress tensor, perimeter component, element 0 -> xx, element 1 -> xy, element 2 -> yx, element 3 -> yy
      vector<double> stress_v = {0.0, 0.0, 0.0, 0.0};    // Stress tensor, velocity component, element 0 -> xx, element 1 -> xy, element 2 -> yx, element 3 -> yy
      double avg_tension;       // Average tension (sum of tensions in all junctions belonging to the cell divided by the number of edges)
      Vec rc;             // centre of the face
      vector<int> neighs; // indices of neighbouring faces 
      vector<Vec> junct_vec;   // vectors along each junction
      vector<Vec> neigh_vec;   // vectors connecting centres of neighbouring cells
      string type_name;     // String with the actual name of the the cell type
      FaceProperty& operator=(const FaceProperty& p) 
      {
        if (this == &p)
          return *this;
        this->face_type = p.face_type;
        this->A0 = p.A0;
        this->P0 = p.P0;
        this->native_A0 = p.native_A0;
        this->native_P0 = p.native_P0;
        this->max_A0 = p.max_A0;
        this->type_name = p.type_name;
        return *this;
      }
    };
  };

 
} // namespace RheoVM



#endif
