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
 * \file mesh.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 09-Jun-2017
 * \brief Mesh class 
 */ 

#ifndef __MESH_HPP__
#define __MESH_HPP__

#include <algorithm>
#include <vector>
#include <iterator>
#include <iostream>
#include <exception>
#include <memory>

#include <pybind11/stl.h>

#include "types.hpp"
#include "box.hpp"

using std::find_if;
using std::vector;
using std::next;
using std::prev;
using std::runtime_error;
using std::shared_ptr;

namespace py = pybind11;

namespace RheoVM
{

  template<typename Property>
  class Mesh 
  {

    public:

    Mesh() : _n_faces(0), _box(nullptr) { }

    //! Member functions 
    void add_vertex(const Vertex<Property>& v)  {  _vertices.push_back(v); }
    void add_edge(const Edge<Property>& e)  {  _edges.push_back(e); } 
    void add_halfedge(const HalfEdge<Property>& he)  {  _halfedges.push_back(he); }
    void add_face(const vector<int>&);

    list<HalfEdge<Property>>& halfedges() { return _halfedges; }
    list<Vertex<Property>>& vertices() { return _vertices; }
    list<Edge<Property>>& edges() { return _edges; }
    list<Face<Property>>& faces() { return _faces; }

    void wipe()
    {
      _halfedges.clear();
      _vertices.clear();
      _edges.clear();
      _faces.clear();
      _erased_edges.clear();
      _erased_halfedges.clear();
      _erased_vertices.clear();
      _n_faces = 0;
      _box = nullptr;
    }

    int num_vert() { return _vertices.size(); }
    int num_faces() { return _faces.size(); }

    Vertex<Property>& get_vertex(int);
    HalfEdge<Property>& get_halfedge(int);
    Edge<Property>& get_edge(int);
    Face<Property>& get_face(int);

    void move_vertex(int i, const Vec&);

    bool can_collapse(EdgeHandle<Property>&);
    bool collapse(EdgeHandle<Property>&);
    bool split(VertexHandle<Property>&, const Vec&, SPLIT_DIRECTION, VertexHandle<Property>&, EdgeHandle<Property>&);
    FaceHandle<Property> face_split(FaceHandle<Property>&, HEHandle<Property>&, HEHandle<Property>&, double, double);
    bool remove_edge(int, bool);

    double area(const FaceHandle<Property>&) const;
    double perim(const FaceHandle<Property>&) const;
    double len(const EdgeHandle<Property>&);
    int coordination(const VertexHandle<Property>&);
    int face_sides(const FaceHandle<Property>&);
    
    bool vertex_in_face(const VertexHandle<Property>&, const FaceHandle<Property>&);
    bool check_overlap(const HEHandle<Property>&);

    void tidyup();   // get the mesh in order (set boundary edges, outer faces, etc. )
    void order_star(VertexHandle<Property>&);  // orders start around a vertex
    Vec get_centre(); // compute geometric centre of the mesh
    Vec get_face_centre(const FaceHandle<Property>&);
    Vec get_face_direction(const FaceHandle<Property>&);
    Vec get_face_centroid(const FaceHandle<Property>&);
    void populate_face_neighbours();

    void set_box(const shared_ptr<Box>& box) { _box = box; }
    const shared_ptr<Box> &box() const { return _box; }

    void scale(double,double);
    void shear(double);
    void transform(double,double,double,double,bool=false);
    
    
    private:

      list<HalfEdge<Property>> _halfedges;
      list<Vertex<Property>>   _vertices;
      list<Edge<Property>>     _edges;
      list<Face<Property>>     _faces;

      list<VertexHandle<Property>> _erased_vertices;
      list<EdgeHandle<Property>> _erased_edges;
      list<HEHandle<Property>> _erased_halfedges;

      int _n_faces;  // number of faces;
      shared_ptr<Box> _box; // simulation box
  };

  template<typename Property>
  Vertex<Property>& Mesh<Property>::get_vertex(int i)
  {
    if ((i < 0) || (i > _vertices.size()))
      throw runtime_error("Vertex index out of bounds.");
    else
      return *(std::next(_vertices.begin(),i));
  }

  template<typename Property>
  HalfEdge<Property>& Mesh<Property>::get_halfedge(int i)
  {
    if ((i < 0) || (i > _halfedges.size()))
      throw runtime_error("Junction index out of bounds.");
    else
      return *(std::next(_halfedges.begin(),i));
  }

  template<typename Property>
  Edge<Property>& Mesh<Property>::get_edge(int i)
  {
    if ((i < 0) || (i > _edges.size()))
      throw runtime_error("Junction index out of bounds.");
    else
      return *(std::next(_edges.begin(), i));
  }

  template<typename Property>
  Face<Property>& Mesh<Property>::get_face(int i)
  {
    if ((i < 0) || (i > _faces.size()))
      throw runtime_error("Face index out of bounds.");
    else
      return *(find_if(_faces.begin(),_faces.end(),[i](const Face<Property>& f) -> bool { return (f.id == i); }));
  }

  template<typename Property>
  void Mesh<Property>::move_vertex(int i, const Vec& v)
  {
    if ((i < 0) || (i > _vertices.size()))
      throw runtime_error("Vertex index out of bounds.");
    VertexHandle<Property> vh = std::next(_vertices.begin(),i);
    vh->r += v;
  }

  template<typename Property>
  void Mesh<Property>::add_face(const vector<int>& vert_ids)
  {
    _faces.push_back(Face<Property>(_n_faces++));
    FaceHandle<Property> fh = prev(_faces.end());
    HEHandle<Property> he;
    HEHandle<Property> prev_he;
    HEHandle<Property> first_he;
    for (int i = 0; i < vert_ids.size(); i++)
    {
      int v1_id = vert_ids[i], v2_id = vert_ids[(i == (vert_ids.size()-1)) ? 0 : i + 1];
      VertexHandle<Property> vh_from = find_if(_vertices.begin(), _vertices.end(), [v1_id](const Vertex<Property>& v)->bool { return (v.id == v1_id); });
      VertexHandle<Property> vh_to   = find_if(_vertices.begin(), _vertices.end(), [v2_id](const Vertex<Property>& v)->bool { return (v.id == v2_id); });
      _halfedges.push_back(HalfEdge<Property>());
      _halfedges.back().set_idx(_halfedges.size()-1);
      he = prev(_halfedges.end());
      if (i == 0) first_he = he;
      he->from() = vh_from;
      he->to()   = vh_to;
      vh_from->he()   = he;
      EdgeHandle<Property> eh = find_if(_edges.begin(), _edges.end(), [v1_id,v2_id](const Edge<Property>& e)->bool { return (v1_id == e.i && v2_id == e.j) || (v1_id == e.j && v2_id == e.i); });

      if (eh == _edges.end())
      {
        _edges.push_back(Edge<Property>(v1_id,v2_id));
        _edges.back().set_idx(_edges.size() - 1);
        eh = prev(_edges.end());
        eh->he() = he;
      }
      else
      {
        eh->he()->pair() = he;
        he->pair() = eh->he();
        eh->boundary = false;
      }
      he->edge() = eh;
      he->face() = fh;
      if (i > 0) 
      {
        he->prev() = prev_he;
        prev_he->next() = he;
      }
      prev_he = he;
    }
    first_he->prev() = he;
    he->next() = first_he;
    fh->he() = first_he;  // Make sure that the first half edge is the face half-edge. This makes life easier when reading in "per edge" data.
    fh->nsides = this->face_sides(fh);
  }

  template<typename Property>
  bool Mesh<Property>::can_collapse(EdgeHandle<Property>& eh)
  {
    if (eh->boundary) 
    {
      //cout << "edge is a boundary "<<endl;
      return false;   // cannot collapse boundary edge
    }
    HEHandle<Property>  he_1 = eh->he();
    HEHandle<Property>  he_2 = he_1->pair();
    VertexHandle<Property> vh_from = he_1->from();
    VertexHandle<Property> vh_to   = he_1->to();
    FaceHandle<Property>   fh_1    = he_1->face();
    FaceHandle<Property>   fh_2    = he_2->face();
    if (vh_to->boundary || vh_from->boundary)  return false;
    
    if ((this->coordination(vh_to) > 3) || (this->coordination(vh_from) > 3))    return false;
       
    if ((vh_from->boundary) || (vh_to->boundary)) return false;  // cannot collapse edge that has boundary vertex
    if (fh_1->nsides < 4 || fh_2->nsides < 4) return false;      // cannot collapse an edge if either of faces is a triangle
    return true;
  }

  template<typename Property>
  bool Mesh<Property>::collapse(EdgeHandle<Property>& eh)
  {
    if (!this->can_collapse(eh)) return false;

    HEHandle<Property>  he_1 = eh->he();
    HEHandle<Property>  he_2 = he_1->pair();
    VertexHandle<Property> vh_from = he_1->from();
    VertexHandle<Property> vh_to   = he_1->to();
    FaceHandle<Property>   fh_1    = he_1->face();
    FaceHandle<Property>   fh_2    = he_2->face();

    Vec r = 0.5*(vh_to->r-vh_from->r);
    vh_from->r += r;
    
    //if (vh_from->he()->to()->id == eh->i || vh_from->he()->to()->id == eh->j)  // shift he of vh_from if its original was erased
    if (vh_from->he() == he_1)
      vh_from->he() = he_1->pair()->next();
    else if (vh_from->he() == he_2)
      vh_from->he() = he_2->pair()->next();

    he_1->next()->prev() = he_1->prev();
    he_2->next()->prev() = he_2->prev();

    he_1->prev()->next() = he_1->next();
    he_2->prev()->next() = he_2->next();

    he_1->next()->from() = vh_from;
    he_2->prev()->to()   = vh_from;

    he_1->next()->pair()->to()   = vh_from;
    he_2->prev()->pair()->from() = vh_from;
    
    if (fh_1->he() == he_1) fh_1->he() = he_1->next();
    if (fh_2->he() == he_2) fh_2->he() = he_2->next();

    fh_1->nsides = this->face_sides(fh_1);
    fh_2->nsides = this->face_sides(fh_2);

    vh_from->coordination = this->coordination(vh_from);

    vh_to->erased = true;
    eh->erased = true;
    he_1->erased = true;
    he_2->erased = true;

    _erased_vertices.push_back(vh_to);
    _erased_edges.push_back(eh);
    _erased_halfedges.push_back(he_1);
    _erased_halfedges.push_back(he_2);

    return true;
  }

  // last two arguments are handles of the new vertex and new edge
  template<typename Property>
  bool Mesh<Property>::split(VertexHandle<Property>& vh, const Vec& v, SPLIT_DIRECTION direction, VertexHandle<Property>& vh_new, EdgeHandle<Property>& eh_new)
  {
    if (vh->erased) 
    {
      //cout<<"erased"<<endl;
      return false;   // cannot split erased vertex
    }
    if (this->coordination(vh) != 4) 
    {
      //cout<<"coordination: "<<this->coordination(vh)<<" " <<vh->r<<endl;
      return false;  // can only split vertices with coordination greater or equal to 4
    }
    if (vh->boundary) 
    {
      //cout<<"boundary"<<endl;
      return false;  // cannot split boundary vertex
    }

    Vec r1 = vh->r + 0.5*v;
    Vec r2 = vh->r - 0.5*v;

    vh->r = r1;

    if (_erased_vertices.size() != 0)
    {
      vh_new = *(_erased_vertices.begin());
      _erased_vertices.pop_front();
    }
    else 
    {
      this->add_vertex(Vertex<Property>(_vertices.size(),r2));
      vh_new = prev(_vertices.end());
    }
    vh_new->erased = false;
    vh_new->r = r2;

    /* 
       This part of code might be unnecessary, but we keep it here for the time being to ensure that the
       edge he_1 always points to the vertex with lowest index.
       Sometimes it may happen that the edge collapse removes the outgoing half-edge of vertex vh. 
       In this case, half-edges are relabelled, which may lead to the four-vertex being split the 
       wrong way. This piece of code should prevent this from happening (hopefully).  
    */
    this->order_star(vh);
    /* ************************************************************************************************** */

    HEHandle<Property> he_1 = vh->he();
    HEHandle<Property> he_2 = he_1->pair()->next();
    HEHandle<Property> he_3 = he_2->pair()->next();
    HEHandle<Property> he_4 = he_3->pair()->next();

    HEHandle<Property> he_n1, he_n2;
    if (_erased_halfedges.size() != 0)
    {
      he_n1 = *(_erased_halfedges.begin());
      _erased_halfedges.pop_front();
    }
    else
    {
      _halfedges.push_back(HalfEdge<Property>());
      _halfedges.back().set_idx(_halfedges.size()-1);
      he_n1 = prev(_halfedges.end());
    }
    he_n1->erased = false;
    if (_erased_halfedges.size() != 0)
    {
      he_n2 = *(_erased_halfedges.begin());
      _erased_halfedges.pop_front();
    }
    else
    {
      _halfedges.push_back(HalfEdge<Property>());
      _halfedges.back().set_idx(_halfedges.size()-1);
      he_n2 = prev(_halfedges.end());
    }
    he_n2->erased = false;

    if (_erased_edges.size() != 0)
    {
      eh_new = *(_erased_edges.begin());
      _erased_edges.pop_front();
    }
    else
    {
      _edges.push_back(Edge<Property>(vh->id, vh_new->id));
      eh_new = prev(_edges.end());
    }
    eh_new->erased = false;
    eh_new->boundary = false;
    eh_new->i = vh->id;
    eh_new->j = vh_new->id;
    //std::cout << "Making new vertex " << vh_new->id << endl;
    //std::cout << "and new edge from " << vh->id << " to " << vh_new->id << endl;
    
    if (direction == SPLIT_1_3)
    {
      eh_new->he() = he_n1;

      he_1->pair()->next() = he_n1;

      he_2->from() = vh_new;
      he_2->pair()->to() = vh_new;
      he_2->prev() = he_n1;

      he_3->from() = vh_new;
      he_3->pair()->to() = vh_new;
      he_3->pair()->next() = he_n2;

      he_4->prev() = he_n2;

      he_n1->from() = vh;
      he_n1->to() = vh_new;
      he_n1->prev() = he_1->pair();
      he_n1->next() = he_2;
      he_n1->face() = he_2->face();
      he_n1->pair() = he_n2; 
      he_n1->edge() = eh_new;

      he_n2->from() = vh_new;
      he_n2->to() = vh;
      he_n2->prev() = he_3->pair();
      he_n2->next() = he_4;
      he_n2->face() = he_4->face();
      he_n2->pair() = he_n1;
      he_n2->edge() = eh_new;
    }
    else if (direction == SPLIT_2_4)
    { 
      eh_new->he() = he_n1;

      he_1->prev() = he_n2;

      he_2->pair()->next() = he_n1;

      he_3->from() = vh_new;
      he_3->pair()->to() = vh_new;
      he_3->prev() = he_n1;

      he_4->from() = vh_new;
      he_4->pair()->to() = vh_new;
      he_4->pair()->next() = he_n2;

      he_n1->from() = vh;
      he_n1->to() = vh_new;
      he_n1->prev() = he_2->pair();
      he_n1->next() = he_3;
      he_n1->face() = he_3->face();
      he_n1->pair() = he_n2;
      he_n2->edge() = eh_new;

      he_n2->from() = vh_new;
      he_n2->to() = vh;
      he_n2->prev() = he_4->pair();
      he_n2->next() = he_1;
      he_n2->face() = he_1->face();
      he_n2->pair() = he_n1;
      he_n2->edge() = eh_new;
    }
    else
    {
      throw runtime_error("Unknown vertex split direction.");
    }
    
    
    vh->he() = he_n1;
    vh_new->he() = he_n2;

    vh->coordination = this->coordination(vh);
    vh_new->coordination = this->coordination(vh_new);

    he_n1->face()->nsides = this->face_sides(he_n1->face());
    he_n2->face()->nsides = this->face_sides(he_n2->face());

    if (he_n1->face()->outer)
    {
      he_n1->edge()->boundary = true;
      he_n1->from()->boundary = true;
      he_n1->to()->boundary = true;
      vh_new->boundary = true;
      //std::cout << "new vertex " << vh_new->id << " is boundary" << endl;
    }
    else if (he_n2->face()->outer)
    {
      he_n2->edge()->boundary = true;
      he_n2->from()->boundary = true;
      he_n2->to()->boundary = true;
      vh_new->boundary = true;
      //std::cout << "new vertex " << vh_new->id << " is boundary" << endl;
    }

    return true;
  }

  /*! Implementation of the code that splits a face (this is used for cell division)
   *  \param fh Handle of the face that will be split
   *  \param he1 Handle of the first half-edge that will be split
   *  \param he2 Handle of the second half-edge that will be split
   *  \param t1 position of the split on he1 (0 < t1 < 1)
   *  \param t2 position of the split on he2 (0 < t2 < 1)
   */
  template<typename Property>
  FaceHandle<Property> Mesh<Property>::face_split(FaceHandle<Property>& fh, HEHandle<Property>& he1, HEHandle<Property>& he2, double t1, double t2)
  {
    if (fh->nsides < 4)
      return _faces.end();
    if ((t1 < 0.0) || (t1 > 1.0))
      throw runtime_error("Face split position t1 has to be between 0 and 1.");
    if ((t2 < 0.0) || (t2 > 1.0))
      throw runtime_error("Face split position t2 has to be between 0 and 1.");
    if ((he1->face()->id != fh->id) || (he2->face()->id != fh->id))
      throw runtime_error("One or both half-edges do not belong to the face.");
    int fnew_id = _faces.size();
    int vn1_id = _vertices.size();
    int vn2_id = vn1_id + 1;

    // Steps 1-6
    // get four vertices involved
    VertexHandle<Property> v1 = he1->from();
    VertexHandle<Property> v2 = he1->to();
    VertexHandle<Property> v3 = he2->from();
    VertexHandle<Property> v4 = he2->to();

    // calculate position of two new vertices
    Vec r1 = v1->r + t1*(v2->r - v1->r);
    Vec r2 = v3->r + t2*(v4->r - v3->r);

    // create and append those two new vertices
    _vertices.push_back(Vertex<Property>(vn1_id,Vec(r1.x,r1.y)));
    VertexHandle<Property> vn1 = prev(_vertices.end());
    _vertices.push_back(Vertex<Property>(vn2_id,Vec(r2.x,r2.y)));
    VertexHandle<Property> vn2 = prev(_vertices.end());

    // Copy vertex properties
    vn1->data() = v1->data();
    vn2->data() = v3->data();

    // create and append three new edge
    _edges.push_back(Edge<Property>(vn1_id,vn2_id));
    EdgeHandle<Property> enew = prev(_edges.end()); 
    _edges.push_back(Edge<Property>(v2->id,vn1_id));
    EdgeHandle<Property> en1 = prev(_edges.end());
    _edges.push_back(Edge<Property>(v4->id,vn2_id));
    EdgeHandle<Property> en2 = prev(_edges.end());

    // create new half-edges
    _halfedges.push_back(HalfEdge<Property>());
    _halfedges.back().set_idx(_halfedges.size()-1);
    HEHandle<Property> hen = prev(_halfedges.end());
    _halfedges.push_back(HalfEdge<Property>());
    _halfedges.back().set_idx(_halfedges.size()-1);
    HEHandle<Property> henp = prev(_halfedges.end());
    _halfedges.push_back(HalfEdge<Property>());
    _halfedges.back().set_idx(_halfedges.size()-1);
    HEHandle<Property> hen1 = prev(_halfedges.end());
    _halfedges.push_back(HalfEdge<Property>());
    _halfedges.back().set_idx(_halfedges.size()-1);
    HEHandle<Property> hen2 = prev(_halfedges.end());
    _halfedges.push_back(HalfEdge<Property>());
    _halfedges.back().set_idx(_halfedges.size()-1);
    HEHandle<Property> henp1 = prev(_halfedges.end());
    _halfedges.push_back(HalfEdge<Property>());
    _halfedges.back().set_idx(_halfedges.size()-1);
    HEHandle<Property> henp2 = prev(_halfedges.end());

    // Set boundary
    if (v1->boundary && v2->boundary)
      vn1->boundary = true;
    
    if (v3->boundary && v4->boundary)
      vn2->boundary = true;

    // Add new face
    // If the last face is not last face, simply add it
    FaceHandle<Property> flast = prev(_faces.end());
    FaceHandle<Property> fnew;
    if (!flast->outer)
    {
      _faces.push_back(Face<Property>(fnew_id));
      fnew = prev(_faces.end());
    }
    else
    {
      _faces.insert(flast, Face<Property>(flast->id));
      flast->id++;
      fnew = prev(_faces.end(),2);
    }
    
    // Step 7
    hen->face() = fh;
    henp->face() = fnew;
    hen1->face() = fh;
    hen2->face() = fnew;
    henp1->face() = he1->pair()->face();
    henp2->face() = he2->pair()->face();

    // Step 8
    enew->he() = hen;
    en1->he() = hen1;
    en2->he() = hen2; 

    // Step 9
    hen->pair() = henp;
    henp->pair() = hen;

    // Step 10
    henp1->pair() = he1;
    hen1->pair() = he1->pair();
    henp2->pair() = he2;
    hen2->pair() = he2->pair();

    // Step 11
    hen->from() = vn2;
    hen->to() = vn1;
    henp->from() = vn1;
    henp->to() = vn2; 

    // Step 12
    he1->to() = vn1;
    hen1->from() = vn1;
    hen1->to() = v2;
    hen1->pair()->to() = vn1;
    henp1->from() = vn1; 
    henp1->to() = v1;
    he2->to() = vn2;
    hen2->from() = vn2;
    hen2->to() = v4; 
    hen2->pair()->to() = vn2; 
    henp2->from() = vn2;
    henp2->to() = v3; 

    // Step 13
    hen1->prev() = hen;
    hen1->next() = he1->next();
    he1->next() = henp;
    hen->next() = hen1;
    henp->prev() = he1;
    henp1->prev() = hen1->pair();
    henp1->next() = hen1->pair()->next();
    hen1->pair()->next() = henp1;
    henp1->next()->prev() = henp1;
    hen1->next()->prev() = hen1;

    // Step 14
    hen1->pair()->pair() = hen1;
    he1->pair() = henp1; 

    // Step 15
    hen2->prev() = henp;
    hen2->next() = he2->next();
    hen->prev() = he2;
    henp->next() = hen2;
    he2->next() = hen; 
    henp2->prev() = hen2->pair();
    henp2->next() = hen2->pair()->next();
    hen2->pair()->next() = henp2;
    henp2->next()->prev() = henp2;
    hen2->next()->prev() = hen2;

    // Step 16
    hen2->pair()->pair() = hen2;
    he2->pair() = henp2;

    // Step 17
    fh->he() = hen;
    fnew->he() = henp;

    // Step 18
    vn1->he() = henp;
    vn2->he() = hen;

    // Step 19
    hen->edge() = enew;
    henp->edge() = enew;
    hen1->edge() = en1;
    hen1->pair()->edge() = en1;
    henp1->edge() = he1->edge();
    hen2->edge() = en2;
    hen2->pair()->edge() = en2;
    henp2->edge() = he2->edge();

    // Step 20
    HEHandle<Property> he = fnew->he();
    HEHandle<Property> first = he;
    do
    {
      he->face() = fnew;
    } while ((he = he->next()) != first);

    // Step 21
    fh->nsides = this->face_sides(fh);
    fnew->nsides = this->face_sides(fnew);

    // Update coordination of new vertices
    vn1->coordination = this->coordination(vn1);
    vn2->coordination = this->coordination(vn2);


    // Split successful
    return fnew;
   
  }
  
  template<typename Property>
  bool Mesh<Property>::remove_edge(int eid, bool make_hole)
  {
    if ((eid < 0) || (eid > _edges.size()))
      throw runtime_error("Junction index out of bounds."); 
    
    EdgeHandle<Property> eh = std::next(_edges.begin(),eid);
    HEHandle<Property> he = eh->he();
    HEHandle<Property> hep = he->pair();
    VertexHandle<Property> vh_from = he->from();
    VertexHandle<Property> vh_to = he->to();

    // Abort cut is this would leave the mesh with dangling vertices
    if ((vh_from->coordination < 3) || (vh_to->coordination < 3))
      return false;

    // or id the edge has already been erased
    if (eh->erased)
      return false;

    FaceHandle<Property> f = he->face();
    FaceHandle<Property> fp = hep->face();

    // Make sure that HEHandle of the surviving face is not erased
    if (f->he() == he)
      f->he() = he->prev();

    // Make sure that vertices do not have erased edge as their half edges
    if (vh_from->he() == he)
      vh_from->he() = he->prev()->pair();

    if (vh_to->he() == hep)
      vh_to->he() = hep->prev()->pair(); 
  
    // update coordination
    vh_from->coordination--;
    vh_to->coordination--;

    // Make sure half-edges of the old face point to the surviving face
    HEHandle<Property> first = hep;
    HEHandle<Property> he_run = hep;
    do 
    {
      he_run->face() = f;
    } while ((he_run = he_run->next()) != first);

    // Fix the connectivity
    he->prev()->next() = hep->next();
    hep->next()->prev() = he->prev();
    he->next()->prev() = hep->prev();
    hep->prev()->next() = he->next();

    // Mark edges and half-edges erased
    he->erased = true;
    hep->erased = true;
    eh->erased = true;

    // Update number of sides
    f->nsides = f->nsides + fp->nsides - 2;

    // Mark the surviving face as outer
    if (make_hole)
      f->outer = true;

    // Erase one of the faces
    fp->erased = true;

    return true;
  }
  
  template<typename Property>
  double Mesh<Property>::area(const FaceHandle<Property>& fh) const
  {
    if (fh->outer) return 0.0;
    HEHandle<Property> he = fh->he();
    HEHandle<Property> first = he;
    VertexHandle<Property> vh0 = he->from();
    double A = 0.0;
    do 
    {
      VertexHandle<Property> vh_from = he->from();
      VertexHandle<Property> vh_to = he->to();
      Vec r1 = vh_from->r - vh0->r; // this takes care of the boundary conditions 
      Vec r2 = vh_to->r - vh0->r;
      A += r1.x*r2.y - r2.x*r1.y;
    } while ((he = he->next()) != first);
    if (A < 0.0) 
      fh->outer = true; 
    else 
      fh->outer = false;
    return 0.5*fabs(A);
  }

  template<typename Property>
  double Mesh<Property>::perim(const FaceHandle<Property>& fh) const
  {
    if (fh->outer) return 0.0;
    HEHandle<Property> he = fh->he();
    HEHandle<Property> first = he;
    double P = 0.0;
    int i = 0;
    do 
    {
      VertexHandle<Property> vh_from = he->from();
      VertexHandle<Property> vh_to = he->to();
      P += (vh_to->r - vh_from->r).len();
      if (i++ > 50) throw runtime_error("Something wrong with the mesh.");
    } while ((he = he->next()) != first);
    return P;
  }

  template<typename Property>
  double Mesh<Property>::len(const EdgeHandle<Property>& eh)
  {
    HEHandle<Property> he = eh->he();
    VertexHandle<Property> v_from = he->from();
    VertexHandle<Property> v_to = he->to();
    return (v_to->r - v_from->r).len();
  }

  template<typename Property>
  int Mesh<Property>::coordination(const VertexHandle<Property>& vh)
  {
    if (vh->erased) return -1;
    HEHandle<Property> he = vh->he();
    HEHandle<Property> first = he;
    int i = 0;
    do 
    {
      i++;
      if (i > 50) throw runtime_error("Something wrong with the mesh.");
      he = he->pair()->next();
    } while (he != first);
    return i;
  }

  template<typename Property>
  int Mesh<Property>::face_sides(const FaceHandle<Property>& fh)
  {
    HEHandle<Property> he = fh->he();
    HEHandle<Property> first = he;
    int i = 0;
    do 
    {
      i++;
      if (i > 100000) throw runtime_error("Something wrong with the mesh.");
      he = he->next();
    } while (he != first);
    return i;
  }

  /* Check if a given vertex is inside a given face. */
  template<typename Property>
  bool Mesh<Property>::vertex_in_face(const VertexHandle<Property>& vh, const FaceHandle<Property>& fh)
  {

    HEHandle<Property> he = fh->he();
    HEHandle<Property> first = he;
    bool inside = false;
    Vec r = vh->r;
    do
    {
      Vec rv = he->from()->r;
      Vec rn = he->to()->r;
      // First make sure that the vertex is not on one of the edges 
      double dx = rn.x - rv.x, dy = rn.y - rv.y;
      if ( (r.x - rv.x)*dy == (r.y - rv.y)*dx)
        return false;
      else if ( ((rv.y > r.y) != (rn.y > r.y)) &&  (r.x < (rn.x - rv.x) * (r.y-rv.y) / (rn.y - rv.y) + rv.x) )
        inside = !inside;
      he = he->next();
    } while (he != first);
    if (fh->outer)   // One have to be carefull as for the outer face, 'inside' is actually 'outside'.
      return (!inside);
    else
      return inside;
  }

  /* For a given half-edge make check if we have a situation where the vertices
     crossed each other. This means that a vertex is inside the face of of one 
     of the faces that belong to other vertex.
  */
  template<typename Property>
  bool Mesh<Property>::check_overlap(const HEHandle<Property>& he)
  {
    VertexHandle<Property> vfrom = he->from();
    VertexHandle<Property> vto = he->to();
    FaceHandle<Property> fh = he->prev()->pair()->face();
    return this->vertex_in_face(vto,fh);
  }

  template<typename Property>
  void Mesh<Property>::tidyup()
  {      
    for (VertexHandle<Property> vh = _vertices.begin(); vh != _vertices.end(); vh++)
      vh->coordination = this->coordination(vh);

    for (EdgeHandle<Property> eh = _edges.begin(); eh != _edges.end(); eh++)
    {
      HEHandle<Property> he = eh->he();
      VertexHandle<Property> vh_from = he->from();
      VertexHandle<Property> vh_to   = he->to();
      eh->boundary = false;
      if (vh_from->boundary && vh_to->boundary)
        eh->boundary = true;
    }
    // sign of area determines if a face is outer
    for (FaceHandle<Property> fh = _faces.begin(); fh != _faces.end(); fh++) 
      this->area(fh);
  }

  // This is an auxiliary function that orders edges around a vertex in such a way that the half-edge
  // associated with the vertex is the one pointing to the vertex with the lowest id
  // This is useful as it ensures that edges the first half-edge in the vertex star is always the same
  // which makes handling things like force dependent vertex split easier to implement as those
  // require comparisons between precisely defined pairs of vertices. 
  template<typename Property>
  void Mesh<Property>::order_star(VertexHandle<Property>& vh)
  {
    HEHandle<Property> he = vh->he();
    int min_id = he->to()->id;
    HEHandle<Property> first = he;
    HEHandle<Property> he_tmp = he;
    int i = 0;
    do
    {
      if (he_tmp->to()->id < min_id)
      {
        min_id = he_tmp->to()->id;
        he = he_tmp;
      }
      he_tmp = he_tmp->pair()->next();
      if (i++ > 50) throw runtime_error("Something's wrong with the mesh.");
    } while (he_tmp != first);
    vh->he() = he;
  }

  // Compute geometric center of the mesh by tracing positions of boundary vertices
  template<typename Property>
  Vec Mesh<Property>::get_centre()
  {
    Vec cm(0.0,0.0);
    for (VertexHandle<Property> vh = _vertices.begin(); vh != _vertices.end(); vh++)
      cm += vh->r;
    return static_cast<double>(1.0/_vertices.size())*cm;
  }

  // Compute centre of a face 
  template<typename Property>
  Vec Mesh<Property>::get_face_centre(const FaceHandle<Property>& fh)
  {
    HEHandle<Property> he = fh->he();
    VertexHandle<Property> vh0 = he->from();
    HEHandle<Property> first = he;
    Vec rc(0.0,0.0);
    do
    {    
      Vec dr = (he->from()->r - vh0->r);
      rc += Vec(dr.x, dr.y);
      he = he->next();      
    } while (he != first);
    Vec Rc = (1.0/fh->nsides)*rc + vh0->r;
    return Vec(Rc.x, Rc.y, this->_box);
  }

  // Compute direction (unit-length) of the eigenvector corresponding to the larger eigenvalue
  // of the tensor of intertia
  template<typename Property>
  Vec Mesh<Property>::get_face_direction(const FaceHandle<Property>& fh)
  {
    double A = 0.0, B = 0.0, C = 0.0;
    Vec rc = this->get_face_centre(fh);
    HEHandle<Property> he = fh->he();
    HEHandle<Property> first = he;
    do
    {    
      Vec dr = he->from()->r - rc;
      A += dr.x*dr.x;  B += dr.x*dr.y;  C += dr.y*dr.y;
      he = he->next();      
    } while (he != first);
    A /= fh->nsides;  B /= fh->nsides;  C /= fh->nsides;
    double l1 = 0.5*(A + C + sqrt((A-C)*(A-C) + 4*B*B));
    double len = sqrt((l1-A)*(l1-A) + B*B);
    return Vec(B/len, (l1-A)/len);
  }

  // Compute position of the face centroid
  template<typename Property>
  Vec Mesh<Property>::get_face_centroid(const FaceHandle<Property>& fh)
  {
    if (fh->outer)
      return Vec(0.0, 0.0);
    HEHandle<Property> he = fh->he();
    HEHandle<Property> first = fh->he();
    VertexHandle<Property> vh0 = he->from();
    Vec rc(0.0, 0.0);
    do
    {    
      Vec ri = he->from()->r - vh0->r;
      Vec rj = he->to()->r - vh0->r;
      double fact = ri.x*rj.y - ri.y*rj.x;
      rc.x += (ri.x + rj.x)*fact;
      rc.y += (ri.y + rj.y)*fact;
      he = he->next();      
    } while (he != first);
    return (1.0/(6*this->area(fh)))*rc + vh0->r;
  }

  // Populate face neighbours info for output
  template<typename Property>
  void Mesh<Property>::populate_face_neighbours()
  {
    for (FaceHandle<Property> fh = _faces.begin(); fh != _faces.end(); fh++) 
    {
      if (!fh->erased)
      {
        Vec rc = this->get_face_centre(fh);
        fh->data().rc = rc;
        fh->data().neighs.clear();
        fh->data().junct_vec.clear();
        fh->data().neigh_vec.clear();
        HEHandle<Property> first = fh->he();
        HEHandle<Property> he = fh->he();
        do
        {
          fh->data().neighs.push_back(he->pair()->face()->id);
          fh->data().junct_vec.push_back(he->to()->r - he->from()->r);
          fh->data().neigh_vec.push_back(this->get_face_centre(he->pair()->face()) - rc);
          he = he->next();
        } while (he != first);
      }
    }
  }

  // Scale the entire system by a factor in each direction
  template<typename Property>
  void Mesh<Property>::scale(double a, double b)
  {
    this->_box->stretch(a, b);
    for (VertexHandle<Property> vh = _vertices.begin(); vh != _vertices.end(); vh++)
      if (!vh->erased)
        vh->r.scale(a, b);
  }
  
  // Shear the entire system
  template<typename Property>
  void Mesh<Property>::shear(double gamma)
  {
    Matrix s(1.0,gamma,gamma,1.0);
    this->_box->shear(gamma);
    for (VertexHandle<Property> vh = _vertices.begin(); vh != _vertices.end(); vh++)
      if (!vh->erased)
        vh->r = s*vh->r;
  }

  // Transform the entire mesh according to an arbitrary transformation matrix
  template<typename Property>
  void Mesh<Property>::transform(double txx, double txy, double tyx, double tyy, bool undo)
  {
    Matrix s(txx, txy, tyx, tyy);
    if (this->_box)
    {
      Matrix old_s = this->_box->s_mat;
      Matrix inv_old_s = old_s.inv();
      this->_box->transform(txx, txy, tyx, tyy, undo);
      for (VertexHandle<Property> vh = _vertices.begin(); vh != _vertices.end(); vh++)
        if (!vh->erased)
        {
          if (undo)
          {
            vh->r = s*inv_old_s*vh->r;
          }
          else
          {
            vh->r = s*vh->r;  
          }
        }  
    }
    else 
    {
      for (VertexHandle<Property> vh = _vertices.begin(); vh != _vertices.end(); vh++)
        if (!vh->erased)
          vh->r = s*vh->r;
    }
    
  }


};

#endif
