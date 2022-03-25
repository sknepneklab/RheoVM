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
 * \file force_area.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-May-2019
 * \brief ForceArea class 
 */ 

#ifndef __FORCE_AREA_HPP__
#define __FORCE_AREA_HPP__

#include "force.hpp"


namespace RheoVM
{

  // Force on a vertex
  class ForceArea : public Force
  {
    public:
      // Orphaned model options
      // const param_type& gammaLen, const param_type& l0, const param_type& counter
      // _gammaLen(gammaLen), _l0(l0), _counter(counter)
                                                                                                                                
      ForceArea(System& sys) : Force(sys) 
                               { 
                                  _kappa.resize(_sys.cell_types().size(), 0.0);
                               }
      virtual ~ForceArea() { }
        
      // computes force on vertex by a given edge
      Vec compute(const VertexHandle<Property>&, const HEHandle<Property>&) override;  
      
      // compute edge tension for handling 4-vertices (and moves that computation out of integrator)
      // takes care of correct type of force law that way
      double tension(const HEHandle<Property>& he, double l0) override
      {
        return 0.0;
      }

      void stress(FaceHandle<Property>&) override;

      // We do not compute face force for this force type
      Vec compute_face_force(FaceHandle<Property> &) override;

      // Energy calculation 
      double energy(const FaceHandle<Property>&) override;
      
      // set all paramters for a given type
      void set_params(const string& cell_type, const params_type& params) override
      {
        for (auto p : params)
          if (p.first != "kappa")
            throw runtime_error("Unknown parameter "+p.first+".");
            
        if (params.find("kappa") == params.end())
          throw runtime_error("Area force requires parameter kappa.");

        try 
        {
          if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
            throw runtime_error("Fore area: Cell type " + cell_type + " is not defined.");
          int ct = _sys.cell_types()[cell_type];
          if (ct + 1 > _kappa.size())
            _kappa.push_back(params.at("kappa"));
          else
            _kappa[ct] = params.at("kappa");
        } 
        catch(const exception& e)
        {
          cerr << "Problem with setting area force paramters. Exception: " << e.what() << '\n';
          throw;
        }
      }

      // set all vector-valued paramters for a given type
      void set_vec_params(const string& cell_type, const vec_params_type& params) override { };

      void set_relative_paramters(const string& cell_type, const multi_params_type& params) override
      {
        if (_use_cell_type)
          throw runtime_error("Using the \"use_cell_type\" flag is not compatible with setting relative values of per-cell paramters.");
        for (auto p : params)
          if (p.first != "kappa" && p.first != "seed" && p.first != "A0")
            throw runtime_error("Unknown parameter " + p.first + ".");
        
        bool has_kappa = params.find("kappa") != params.end();
        bool has_A0 = params.find("A0") != params.end();

        try 
        {
          unsigned int seed = system_clock::now().time_since_epoch().count();
          if (params.find("seed") != params.end())
            seed = static_cast<unsigned int>(params.at("seed")[0]);
          RNG rng(seed);
          bool setall = false; 
          if (cell_type == "all")
            setall = true;
          else
          {
            if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
              throw runtime_error("Fore area: Cell type " + cell_type + " is not defined.");
          }
          
          for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
          {
            if (has_kappa)
            {
              double mean_kappa = params.at("kappa")[0];
              double std_kappa = 0.0; 
              if (params.at("kappa").size() > 1) std_kappa = params.at("kappa")[1];
              double scale_kappa = (fabs(std_kappa) < 1e-6) ? mean_kappa : rng.normal(mean_kappa, std_kappa);
              if (setall)
                fh->data().kappa *= scale_kappa;
              else
                if (fh->data().type_name == cell_type)
                  fh->data().kappa *= scale_kappa;
            }
            if (has_A0)
            {
              double mean_A0 = params.at("A0")[0];
              double std_A0 = 0.0; 
              if (params.at("A0").size() > 1) std_A0 = params.at("A0")[1];
              double scale_A0 = (fabs(std_A0) < 1e-6) ? mean_A0 : rng.normal(mean_A0, std_A0);
              if (setall)
                fh->data().A0 *= scale_A0;
              else
                if (fh->data().type_name == cell_type)
                  fh->data().A0 *= scale_A0;
            }
          }
        } 
        catch(const exception& e)
        {
          cerr << "Problem with setting relative paramters for area force. Exception: " << e.what() << '\n';
          throw;
        }  

      }

      void set_flag(const string& flag) override 
      { 
        if (flag == "use_cell_type") 
        {
          _use_cell_type = true;
          cout << "Warning! Setting use_cell_type flag in area force will override all parameters read from the input JSON file (if the read_params flag in the read_input function was set to True)." << endl;
        }
        else
          throw runtime_error("Unknown flag : " + flag + ".");
      }

      void copy_type_param_to_cell() override
      {
        for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
          fh->data().kappa  = _kappa[fh->data().face_type];
      }

    
    private:

      vector<double> _kappa; 
      
      
  };

  
}

#endif
