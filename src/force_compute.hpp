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
 * \file force_compute.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief ForceCompute class 
 */ 

#ifndef __FORCE_COMPUTE_HPP__
#define __FORCE_COMPUTE_HPP__

#include <exception>
#include <algorithm>
#include <map>
#include <string>

#include "system.hpp"
#include "class_factory.hpp"
#include "force.hpp"
#include "force_area.hpp"
#include "force_perimeter.hpp"
#include "force_self_propulsion.hpp"

 
using std::runtime_error;
using std::transform;
using std::map;
using std::string;


namespace RheoVM
{

  class ForceCompute : public ClassFactory<Force>
  {
    public:

      ForceCompute(System& sys) : _sys{sys},
                                  _update_flag{false},
                                  _compute_stress{false}
                                  { }

      ~ForceCompute() = default; 

      ForceCompute(const ForceCompute&) = delete;

      void compute_forces_and_stresses()
      {
        for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
          if (!vh->erased)
            this->compute(vh);
        if (this->_compute_stress)
          for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
              if (!(fh->erased || fh->outer)) 
                this->compute_stress(fh);
      }
      
      void compute(VertexHandle<Property> &vh)
      {
        vh->data().force = Vec(0.0,0.0);
        if (_update_flag)
          for (auto f_type : vh->data().f_type)
            vh->data().f_type[f_type.first] = Vec(0.0,0.0);
        HEHandle<Property> he = vh->he();
        HEHandle<Property> first = he;
        int i = 0;
        do 
        {
          vh->data().force += this->compute(vh, he);
          he = he->pair()->next();
          if (i++ > 50) throw runtime_error("Something's wrong with the mesh."); // This is some consistency check of the mesh being a regular mesh
        } while (he != first);
      }

      Vec compute(VertexHandle<Property> &vh, const HEHandle<Property> &he)
      {
        Vec ftot(0,0);
        for (auto& f : this->factory_map)
        {
          Vec force = f.second->compute(vh, he);
          if (_update_flag)
          {
            vh->data().f_type[f.first] += force;
            he->data().force_type[f.first] = force;
          }
          ftot += force;
        }
        if (_update_flag)
          he->data().tension = dot(ftot, he->direction().unit());
        return ftot;
      }

      void compute_stress(FaceHandle<Property> &fh)
      {
        fh->data().stress.assign(4,0.0);
        fh->data().stress_a.assign(4,0.0);
        fh->data().stress_p.assign(4,0.0);
        fh->data().stress_v.assign(4,0.0);
        for (auto& f : this->factory_map)
          f.second->stress(fh);
      }

      Vec compute_face_force(int i)
      {
        if (i < 0 || i >= _sys.mesh().faces().size())
          throw runtime_error("Wrong face id.");

        FaceHandle<Property> fh = find_if(_sys.mesh().faces().begin(), _sys.mesh().faces().end(), [i](const Face<Property> &f) -> bool { return (f.id == i); });
        Vec ftot(0, 0);
        for (auto& f : this->factory_map)
          ftot += f.second->compute_face_force(fh);
        return ftot;
      }

      double tension(HEHandle<Property>& he,  double l0)
      {
        double T = 0.0;
        for (auto& f : this->factory_map)
          T += f.second->tension(he, l0);
        return T;
      }

      double energy(const FaceHandle<Property>& fh)
      {
        double E = 0.0;
        for (auto& f : this->factory_map)
          E += f.second->energy(fh);
        return E;
      }

      map<string, double> energy_comp(const FaceHandle<Property>& fh)
      {
        map<string, double> E;
        for (auto& f : this->factory_map)
          E[f.first] = f.second->energy(fh);
        return E;
      }

      void set_params(const string& fname, const string& type, const params_type& params)
      {
        if (this->factory_map.find(fname) != this->factory_map.end())
          this->factory_map[fname]->set_params(type, params);
        else
          throw runtime_error("set_params: Force type " + fname + " is not used in this simulation.");
      }

      void set_vec_params(const string& fname, const string& type, const vec_params_type& params)
      {
        if (this->factory_map.find(fname) != this->factory_map.end())
          this->factory_map[fname]->set_vec_params(type, params);
        else
          throw runtime_error("set_vec_params: Force type " + fname + " is not used in this simulation.");
      }

      void set_relative_params(const string& fname, const string& type, const multi_params_type& params)
      {
        if (this->factory_map.find(fname) != this->factory_map.end())
          this->factory_map[fname]->set_relative_paramters(type, params);
        else
          throw runtime_error("set_relative_paramters: Force type " + fname + " is not used in this simulation.");
      }

      void set_flag(const string& fname, const string& flag)
      {
        if (this->factory_map.find(fname) != this->factory_map.end())
          this->factory_map[fname]->set_flag(flag);
        else
          throw runtime_error("set_flag: Force type " + fname + " is not used in this simulation.");
      }

      void copy_type_param_to_cell() 
      {
        cout << "Warning! Overwriting cell paramters with parameters set by type. This action is not reversible." << endl;
        for (auto& f : this->factory_map)
          f.second->copy_type_param_to_cell();
      }

      void set_update(bool flag) { _update_flag = flag; }
      void set_compute_stress(bool flag) { _compute_stress = flag; }

      bool update() { return _update_flag; }
      bool stress_compute() { return _compute_stress; }

      void add_force(const string& fname)
      {
        string name = fname; 
        transform(name.begin(), name.end(), name.begin(), ::tolower);
        if (name == "area")
          this->add<ForceArea,System&>(name, _sys);
        else if (name == "perimeter")
          this->add<ForcePerimeter,System&>(name, _sys);
        else if (name == "self-propulsion")
          this->add<ForceSelfPropulsion,System&>(name, _sys);
        else 
          throw runtime_error("Unknown force type : " + name + ".");
      }

    private: 

      System& _sys;
      bool _update_flag;
      bool _compute_stress;

  };

  void export_ForceCompute(py::module&);

}
#endif
