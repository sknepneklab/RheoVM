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
 * \file types.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 09-Jun-2017
 * \brief Data types for the Mesh 
 */ 

#ifndef __TYPES_HPP__
#define __TYPES_HPP__

//#define VERSION

#include <list>

#include <pybind11/pybind11.h>

#include "vec.hpp"


namespace py = pybind11;
using std::list;



namespace RheoVM
{
  
  template<typename Property> class HalfEdge;
  template<typename Property> class Vertex;
  template<typename Property> class Edge;
  template<typename Property> class Face;

  template<typename Property> using HEHandle     = typename list<HalfEdge<Property>>::iterator;
  template<typename Property> using VertexHandle = typename list<Vertex<Property>>::iterator;
  template<typename Property> using EdgeHandle   = typename list<Edge<Property>>::iterator;
  template<typename Property> using FaceHandle   = typename list<Face<Property>>::iterator;

  template<typename Property> using HECHandle     = typename list<HalfEdge<Property>>::const_iterator;
  template<typename Property> using VertexCHandle = typename list<Vertex<Property>>::const_iterator;
  template<typename Property> using EdgeCHandle   = typename list<Edge<Property>>::const_iterator;
  template<typename Property> using FaceCHandle   = typename list<Face<Property>>::const_iterator;
  
  //! HalfEdge class
  template<typename Property>
  class HalfEdge
  {

    public:

      // Constructors
      HalfEdge() : _idx{-1}, 
                   _property{typename Property::HEProperty()}, 
                   erased{false}  
                   { 

                   }
      //HalfEdge(const Property& p) : _property(p), erased(false)  { }

      // Member functions
      int idx() { return _idx; }  
      void set_idx(int idx) { _idx = idx; }    

      typename Property::HEProperty& data() { return _property; }
      typename Property::HEProperty  data() const { return _property; }

      VertexHandle<Property>& from() {  return _from; }
      VertexHandle<Property>& to()   {  return _to; }

      VertexCHandle<Property> from() const  {  return _from; }
      VertexCHandle<Property> to()   const  {  return _to; }
 
      EdgeHandle<Property>&  edge()  {  return _edge;  }
      EdgeCHandle<Property>  edge() const  {  return _edge;  }

      FaceHandle<Property>&  face()  {  return _face; }
      FaceCHandle<Property>  face()  const {  return _face; }

      HEHandle<Property>&    pair()  {  return _pair; }
      HEHandle<Property>&    next()  {  return _next; }
      HEHandle<Property>&    prev()  {  return _prev; }

      HECHandle<Property>     pair() const {  return _pair; }
      HECHandle<Property>     next() const {  return _next; }
      HECHandle<Property>     prev() const {  return _prev; }

      Vec direction() { return _to->r - _from->r; }

      // Public members
      bool erased;
      
    private:

      int _idx;                               // Half-edge index (used for debugging)
      typename Property::HEProperty          _property;

      VertexHandle<Property>  _from;         // vertex it starts from (that vertex will have this halfedge as its he)
      VertexHandle<Property>  _to;           // vertex it points to (its pair will have this vertex as its he)

      EdgeHandle<Property>   _edge;          // edge this he is part of

      FaceHandle<Property>    _face;         // face to the left of it, when looking in the direction this he points to

      HEHandle<Property>   _pair;            // its pair half edge (undefined for boundary edges)
      HEHandle<Property>   _next;            // next he in the same face
      HEHandle<Property>   _prev;            // previous he in the same face
      
  };

  
  //!< Vertex class
  template<typename Property>
  class Vertex
  {
    public:

      // Constructors    
      Vertex()  : id(0), r(0.0,0.0), _property(typename Property::VertexProperty()), erased(false), boundary(false) { }
      Vertex(int id, const Vec& r) : id(id), r(r), _property(typename Property::VertexProperty()), erased(false), boundary(false) { }
      Vertex(int id, const Vec& r, bool bnd) : id(id), r(r), _property(typename Property::VertexProperty()), erased(false), boundary(bnd) { }
      //Vertex(int id, const Vec& r, const Property& p) : id(id), r(r), _property(p), erased(false), boundary(false) { }
      
      // Member functions
      typename Property::VertexProperty& data() { return _property; }

      typename Property::VertexProperty data() const { return _property; }

      HEHandle<Property>& he()  { return _he; }
      
      HECHandle<Property> he() const { return _he; }

      // Public members 
      Vec   r;              // position
      int   id;             // unique id
      bool  erased;         // marks vertices that are not connected to the rest of the mesh, but are still in memory 
      bool  boundary;       // if true, vertex is on boundary 
      int   coordination;   // number of neighbours this vertex has

    private:

      typename Property::VertexProperty  _property;
      HEHandle<Property>   _he;          // outgoing half edge
         
  };

  
  //!< Edge class
  template<typename Property>
  class Edge
  {
    public:

      // Constructors    
      Edge()  : _idx{0}, i{0}, j{0}, _property{typename Property::EdgeProperty()}, boundary{true}, erased{false} { }
      //Edge(int i, int j, const Property& p) : i(i), j(j), _property(p), boundary(true), erased(false) { }
      Edge(int i, int j) : _idx{0}, i{i}, j{j}, _property{typename Property::EdgeProperty()}, boundary{true}, erased{false} { }
      
      // Member functions
      int idx() const { return _idx;  }
      void set_idx(int idx) { _idx = idx;  }
      typename Property::EdgeProperty& data() { return _property; }
      typename Property::EdgeProperty  data() const { return _property; }

      HEHandle<Property>& he()  { return _he; }
      
      HECHandle<Property> he() const { return _he; }

      // Public members
      int i, j;               // indices of two vertices
      bool boundary;          // if true, edge is a boundary edge
      bool erased;            // markes all erased edges that are still in memory

    private:

      int _idx;  // Unique edge idx
      typename Property::EdgeProperty _property;
      HEHandle<Property>   _he;     // one of the two half edges
      
  };

  //!< Face class
  template<typename Property>
  class Face
  {
    public:

      // Constructors    
      Face()  : id(0), _property{typename Property::FaceProperty()}, outer{false}, erased{false} { }
      Face(int id)  : id(id), _property{typename Property::FaceProperty()}, outer{false}, erased{false} { }
      //Face(const Property& p) : id(0), _property(p), outer(false) { }
      
      // Member functions
      typename Property::FaceProperty& data() { return _property; }
      typename Property::FaceProperty  data() const { return _property; }

      HEHandle<Property>& he()  { return _he; }
      
      HECHandle<Property> he() const { return _he; }

      // public members 
      int id;        // face id
      bool outer;    // if true, face is a ghost outer face
      int nsides;     // number of sides face has
      bool erased;    // if true, face is marked as erased

    private:

      typename Property::FaceProperty      _property;
      HEHandle<Property>   _he;           // one of its half edges
      
  };

  enum SPLIT_DIRECTION { SPLIT_1_3, SPLIT_2_4 };
  
}

#endif
