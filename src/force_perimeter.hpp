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
 * \file force_perimeter.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-May-2019
 * \brief ForcePerimeter class 
 */ 

#ifndef __FORCE_PERIMETER_HPP__
#define __FORCE_PERIMETER_HPP__

#include "force.hpp"

#include <algorithm>

using std::transform;

namespace RheoVM
{

  // Force on a vertex
  class ForcePerimeter : public Force
  {
    public:
      // Orphaned model options
      // const param_type& gammaLen, const param_type& l0, const param_type& counter
      // _gammaLen(gammaLen), _l0(l0), _counter(counter)
                                                                                                                                
      ForcePerimeter(System& sys) : Force{sys},
                                    _double_boundary_constants{false},
                                    _lambda_P0{false}
                                    { 
                                      _gamma.resize(_sys.cell_types().size(), 0.0);
                                      _lambda.resize(_sys.cell_types().size(), 0.0);
                                    }
      virtual ~ForcePerimeter() { }
        
      // computes force on vertex by a given edge
      Vec compute(const VertexHandle<Property>&, const HEHandle<Property>&) override;  
      
      // compute edge tension for handling 4-vertices (and moves that computation out of integrator)
      // takes care of correct type of force law that way
      double tension(const HEHandle<Property>&, double) override;

      void stress(FaceHandle<Property>&) override;

      // We do not compute face force for this force type
      Vec compute_face_force(FaceHandle<Property> &) override;
      
      // Energy calculation 
      double energy(const FaceHandle<Property>&) override;

      // set all paramters for a given type
      void set_params(const string& cell_type, const params_type& params) override
      {
        for (auto p : params)
          if (p.first != "gamma" && p.first != "lambda")
            throw runtime_error("Unknown parameter "+p.first+".");

        if (params.find("gamma") == params.end())
          throw runtime_error("Perimeter force requires parameter gamma.");
        if (params.find("lambda") == params.end())
          throw runtime_error("Perimeter force requires parameter lambda.");

        try
        {
          if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
            throw runtime_error("Fore perimeter: Cell type " + cell_type + " is not defined.");
          int ct = _sys.cell_types()[cell_type];
          if (ct + 1 > _gamma.size())
            _gamma.push_back(params.at("gamma"));
          else
            _gamma[ct] = params.at("gamma");
          if (ct + 1 > _lambda.size())
            _lambda.push_back(params.at("lambda"));
          else
            _lambda[ct] = params.at("lambda");
        }
        catch(const exception& e)
        {
          cerr << "Problem with setting perimeter force paramters. Exception: " << e.what() << '\n';
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
          if (p.first != "gamma" && p.first != "lambda" && p.first != "seed")
            throw runtime_error("Unknown parameter " + p.first + ".");
        
        bool has_gamma = params.find("gamma") != params.end();
        bool has_lambda = params.find("lambda") != params.end();

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
            if (has_gamma)
            {
              double mean_gamma = params.at("gamma")[0];
              double std_gamma = 0.0;
              if (params.at("gamma").size() > 1) std_gamma = params.at("gamma")[1];
              double scale_gamma = (fabs(std_gamma) < 1e-6) ? mean_gamma : rng.normal(mean_gamma, std_gamma);
              if (setall)
                fh->data().gamma *= scale_gamma;
              else
                if (fh->data().type_name == cell_type)
                  fh->data().gamma *= scale_gamma;
            }
            if (has_lambda)
            {
              double mean_lambda = params.at("lambda")[0];
              double std_lambda = 0.0;
              if (params.at("lambda").size() > 1) std_lambda = params.at("lambda")[1];
              double scale_lambda = (fabs(std_lambda) < 1e-6) ? mean_lambda : rng.normal(mean_lambda, std_lambda);
              if (setall)
                fh->data().lambda *= scale_lambda;
              else
                if (fh->data().type_name == cell_type)
                  fh->data().lambda *= scale_lambda;
            }
          }
        } 
        catch(const exception& e)
        {
          cerr << "Problem with setting relative paramters for perimeter force. Exception: " << e.what() << '\n';
          throw;
        }  

      }

      void set_flag(const string& flag) override 
      { 
        if (flag == "use_cell_type")  
        {
          _use_cell_type = true;
          cout << "Warning! Setting use_cell_type flag in perimeter force will override all parameters read from the input JSON file (if the read_params flag in the read_input function was set to True)." << endl;
        }
        else if (flag == "double_boundary_constants") 
          _double_boundary_constants = true;
        else if (flag == "use_P0")
          _lambda_P0 = true;
        else
          throw runtime_error("Unknown flag : " + flag + ".");
      }

      void copy_type_param_to_cell() override
      {
        for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
        {
          fh->data().gamma = _gamma[fh->data().face_type];
          fh->data().lambda  = _lambda[fh->data().face_type];
        }
      }

    
    private:

      vector<double> _gamma; 
      vector<double> _lambda; 
      bool _double_boundary_constants;   // If true, gamma and lambda for cells bordering a hole are doubled
      bool _lambda_P0;                   // If true, lambda will be computed as gamma*P0 where P0 is read from the input configuration 
      
  };

  
}

#endif
