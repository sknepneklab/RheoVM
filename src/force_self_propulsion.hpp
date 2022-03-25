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
 * \file force_self_propulsion.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 12-Jun-2019
 * \brief ForceSelfPropulsion class 
 */ 

#ifndef __FORCE_SELF_PROPULSION_HPP__
#define __FORCE_SELF_PROPULSION_HPP__

#include "force.hpp"


namespace RheoVM
{

  // Force on a vertex
  class ForceSelfPropulsion : public Force
  {
    public:
      // Orphaned model options
      // const param_type& gammaLen, const param_type& l0, const param_type& counter
      // _gammaLen(gammaLen), _l0(l0), _counter(counter)
                                                                                                                                
      ForceSelfPropulsion(System& sys) : Force(sys) 
                               { 
                                  _v0.resize(_sys.cell_types().size(), 0.0);
                                  _n.resize(_sys.cell_types().size(), Vec(0,0));
                               }
      virtual ~ForceSelfPropulsion() { }
        
      // computes force on vertex by a given edge
      Vec compute(const VertexHandle<Property>&, const HEHandle<Property>&) override;  
      
      // compute edge tension for handling 4-vertices (and moves that computation out of integrator)
      // takes care of correct type of force law that way
      double tension(const HEHandle<Property>& he, double l0) override
      {
        return 0.0;
      }

      // Stress calculation not implemented yet
      void stress(FaceHandle<Property>& fh) override 
      {
        //throw runtime_error("Stress calculation has not been implemented for self-propulsion force.");
      }

      // We do not compute face force for this force type
      Vec compute_face_force(FaceHandle<Property> &) override;

      // Energy of self propulsion is not defined, so we return 0
      double energy(const FaceHandle<Property>& fh) override
      {
        return 0.0;   
      }
      
      // set all paramters for a given type
      void set_params(const string& cell_type, const params_type& params) override
      {
        for (auto p : params)
          if (p.first != "v0")
            throw runtime_error("Unknown parameter " + p.first + ".");
            
        if (params.find("v0") == params.end())
          throw runtime_error("SelfPropulsion force requires parameter v0.");

        try 
        {
          if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
            throw runtime_error("Fore self-propulsion: Cell type " + cell_type + " is not defined.");
          int ct = _sys.cell_types()[cell_type];
          if (ct + 1 > _v0.size())
            _v0.push_back(params.at("v0"));
          else
            _v0[ct] = params.at("v0");
        } 
        catch(const exception& e)
        {
          cerr << "Problem with setting self-propulsion force paramters. Exception: " << e.what() << '\n';
          throw;
        }
      }

      // set all vector-valued paramters for a given type
      void set_vec_params(const string& cell_type, const vec_params_type& params) override 
      { 
        for (auto p : params)
          if (p.first != "n")
            throw runtime_error("Unknown parameter " + p.first + ".");

        if (params.find("n") == params.end())
          throw runtime_error("SelfPropulsion force requires vector parameter n.");

        try 
        {
          int ct = _sys.cell_types()[cell_type];
          if (ct + 1 > _n.size())
            _n.push_back(params.at("n").unit());
          else
            _n[ct] = params.at("n").unit();
        } 
        catch(const exception& e)
        {
          cerr << "Problem with setting self-propulsion force paramters. Exception: " << e.what() << '\n';
          throw e;
        }

      };

      void set_relative_paramters(const string& cell_type, const multi_params_type& params) override { };

      void set_flag(const string& flag) override 
      { 
        if (flag == "use_cell_type") 
        {
          _use_cell_type = true;
          cout << "Warning! Setting use_cell_type flag in self-propulsion force will override all parameters read from the input JSON file (if the read_params flag in the read_input function was set to True)." << endl;
        }
        else
          throw runtime_error("Unknown flag : " + flag + ".");
      }

      void copy_type_param_to_cell() override
      {
        for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
          fh->data().v0 = _v0[fh->data().face_type];
      }
    
    private:

      vector<double> _v0; 
      vector<Vec> _n; 
      
      
  };

  
}

#endif
