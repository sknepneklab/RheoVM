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
 * \file measurement.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-May-2018
 * \brief Mmeasurement class 
 */ 

#include "measurement.hpp"

using std::pair;
using std::min;
using std::max;


namespace RheoVM
{
  void Measurement::populate()
  {
    _sites.clear();
    for (FaceHandle<Property> fh = _mesh.faces().begin(); fh != _mesh.faces().end(); fh++)
    {
      if (!fh->outer)
      {
        int id = fh->id;
        Vec r = _mesh.get_face_centre(fh);
        _sites.push_back(Site(id,r));
      }
    }

    for(FaceHandle<Property> fh = _mesh.faces().begin(); fh != _mesh.faces().end(); fh++)
    {
      if (!fh->outer)
      {
        int id = fh->id;
        HEHandle<Property> he = fh->he();
        HEHandle<Property> first = he;
        do
        {
          FaceHandle<Property> f = he->pair()->face();
          if (!f->outer)
            _sites[id].neigh.push_back(f->id);
          he = he->next();
        } while (he != first);
      }
    }
  }

  void Measurement::populate_pairs(const vector< pair<int,int> > & faces_pair)
  { //same function but for selected faces, also store junction(s)
    _sites.clear();
    for (FaceHandle<Property> fh = _mesh.faces().begin(); fh != _mesh.faces().end(); fh++)
    {
      if (!fh->outer)
      {
        int id = fh->id;
        Vec r = _mesh.get_face_centre(fh);
        _sites.push_back(Site(id,r));
      }
    }

    for (vector< pair<int,int> >::const_iterator it=faces_pair.begin();it != faces_pair.end();it++)
    {
      Face<Property> f1 = _mesh.get_face(it->first);
      if (!f1.outer)
      {
        int id1 = it->first;
        int id2 = it->second;
        HEHandle<Property> he = f1.he();
        HEHandle<Property> first = he;
        do
        {
          FaceHandle<Property> f2 = he->pair()->face();
          if ( f2->id == id2  && !f2->outer  ){
            Vec l = he->to()->r - he->from()->r;
            _sites[id1].neigh.push_back(f2->id);
            _sites[id1].jc.push_back(he);
          }
          he = he->next();
        } while (he != first);
      }
    }
  }


void Measurement::compute_texture()
{
  if (_sites.size() == 0) 
    throw runtime_error("Measurement data structure has not been populated. Please use populate() function to do it.");
  for (Site& s : _sites)
  {
    Vec r1 = s.r;
    double x2 = 0.0, xy = 0.0, y2 = 0.0;
    for (int i : s.neigh)
    {
      Vec r2 = _sites[i].r;
      Vec l = r2 - r1;
      x2 = l.x*l.x;  xy = l.x*l.y;  y2 = l.y*l.y;
      s.l.push_back(l);
      s.m.push_back(Texture(x2,xy,y2));
    }
  }
}


void Measurement::dump_texture(const string& filename)
{
  ofstream out(filename.c_str(),ofstream::app);
  //out << "# id  xc yc  <xx>  <xy>  <yy>" << std::endl;
  for (Site& s : _sites)
  {
    double XX = 0.0, XY = 0.0, YY = 0.0;
    for (int i = 0; i < s.m.size(); i++)
    {
      XX += s.m[i].xx;  XY += s.m[i].xy;  YY += s.m[i].yy; 
    }
    int Ntot = s.m.size();
    XX /= Ntot;  XY /= Ntot;  YY /= Ntot;

    if (Ntot > 0)
    {
      out << _sys.time_step() <<  "  " << s.id << "  " << s.r.x << "  " << s.r.y << "  " << XX << "  " << XY << "  " << YY << "  " << Ntot;
      if(s.jc.size() > 0){
        for(auto he : s.jc)          
        {
          Vec l = he->to()->r - he->from()->r;
          
          out << "  " << l.x << "  " << l.y << "  ";
        }
      }      
      out << std::endl;
    }
  }
  out.close();
}


void Measurement::dump_junctions(const vector< pair<int,int> > & faces_pair, const  string &  filename)
{
  ofstream out(filename, ofstream::app);

  for (vector< pair<int,int> >::const_iterator it=faces_pair.begin();it != faces_pair.end();it++)
  {
    int id = it->first;
    Face<Property> fh = _mesh.get_face(id);
    
    if (!fh.outer )
    {          
      int id2 = it->second;
      HEHandle<Property> he = fh.he();
      HEHandle<Property> first = he;
      do
      { // if the two faces are neighbors find common junction and dump it 
        FaceHandle<Property> f = he->pair()->face();
        if (f->id == id2 && !f->outer )
        {
          Vec l = he->to()->r - he->from()->r;
          out << _sys.time_step() << "  " << l.x << "  " << l.y << "  " << endl;
        }
        he = he->next();
      } while (he != first);
    }
  }
}

void Measurement::dump_T1(const string& filename) 
{ // this function must be called at the end of the simulation otherwise 4-vertices will be output several time 
  ofstream out(filename.c_str(),std::ios_base::app);

  if( out.tellp() == 0) //check if file empty
    out << "#vert_id|vec_pos|time_step_pre|num_active_faces|f1|f2|vec_faces_pre|time_step_post|f3|f4|vec_faces_post|kind" << endl; 
  
  
  for (VertexHandle<Property> vh = _mesh.vertices().begin(); vh != _mesh.vertices().end(); vh++)
  {
    for (std::vector<T1_stats>::iterator it_pre = vh->data().pre_T1.begin() ,  it_post = vh->data().post_T1.begin() ; it_pre != vh->data().pre_T1.end(); it_pre++)
    {		
      out << vh->id << "|" << vh->r << "|";
      out << it_pre->time_step << "|" ;
      out << it_pre->active_faces << "|" ; 	      // number of active faces involved:
      out << it_pre->f1_id << "|" << it_pre->f2_id << "|" << it_pre->r_faces << "|" ;
      if(it_post != vh->data().post_T1.end())
      {
        pair<int,int> faces_pre  = std::make_pair(min(it_pre->f1_id,it_pre->f2_id),max(it_pre->f1_id,it_pre->f2_id));
        pair<int,int> faces_post = std::make_pair(min(it_post->f1_id,it_post->f2_id),max(it_post->f1_id,it_post->f2_id));
        int kind = (faces_pre == faces_post ? 0 : 1) ;		  
        out << it_post->time_step << "|";
        out << it_post->f1_id << "|" << it_post->f2_id << "|" << it_post->r_faces << "|";
        out << kind;
        it_post++;
      }
      else
      {
        out  << "-1|-1|-1|(-1,-1)|-1";		  
      }

      out << endl;
    }
    vh->data().pre_T1.erase(vh->data().pre_T1.begin(),vh->data().pre_T1.begin()+vh->data().post_T1.size());
    vh->data().post_T1.clear();	  
  }
}

void export_Texture(py::module& m)
{
  py::class_<Texture>(m,"Texture")
        .def(py::init<double,double,double>())
        .def_readonly("xx", &Texture::xx)
        .def_readonly("xy", &Texture::xy)
        .def_readonly("yy", &Texture::yy);
}

void export_Site(py::module& m)
{
  py::class_<Site>(m,"Site")
        .def(py::init<int,Vec&>())
        .def_readonly("id", &Site::id)
        .def_readonly("r", &Site::r)
        .def_readonly("neigh", &Site::neigh)
        .def_readonly("l", &Site::l)
        .def_readonly("m", &Site::m);        
}

void export_Measurement(py::module& m)
{
  py::class_<Measurement>(m, "Measurement")
    //          .def(py::init<MyMesh&>())
    .def(py::init<MyMesh&,System&>())
        .def("sites", &Measurement::sites)
        .def("populate", &Measurement::populate)
        .def("populate_pairs", &Measurement::populate_pairs)
        .def("dump_junctions", &Measurement::dump_junctions)
        .def("texture", &Measurement::compute_texture)
        .def("dump_texture", &Measurement::dump_texture)
        .def("dump_T1", &Measurement::dump_T1);
}


}
