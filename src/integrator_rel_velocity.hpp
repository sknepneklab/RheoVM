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
 * \file integrator_brownian.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 14-Apr-2021
 * \brief IntegratorRelVelocityclass 
*/

#ifndef __INTEGRATOR_REL_VELOCITY_HPP__
#define __INTEGRATOR_REL_VELOCITY_HPP__
//#define HAS_EIGEN3
#ifdef HAS_EIGEN3

#include "rng.hpp"

#include "constrainer.hpp"
#include "constraint_none.hpp"
#include "constraint_fixed.hpp"

#include "integrator.hpp"

#include <chrono>
#include <utility>
#include <map>
#include <memory>
#include <vector>

#include <Eigen/Sparse>


using namespace std::chrono;
using std::make_unique;
using std::map;
using std::vector;
using std::pair;
using std::make_pair;


namespace RheoVM
{
  using ConstrainerType = unique_ptr<Constrainer>;

  using Tripl = Eigen::Triplet<double>;
  using SpMat = Eigen::SparseMatrix<double>;
  using SolType = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;

  class IntegratorRelVelocity : public Integrator 
  {

    public:

      IntegratorRelVelocity(System& sys, ForceCompute& fc, int seed) : Integrator{sys, fc},
                                                                       _gamma{0.0},
                                                                       _zeta{1.0},
                                                                       _cell_centres{false}
      { 
        map<string,int>& vert_types = _sys.vert_types();
        map<string,int>& cell_types = _sys.cell_types();
        for (int i = 0; i < vert_types.size(); i++)  _constant_force.push_back(Vec(0.0,0.0));
        _constrainer = make_unique<Constrainer>();
        _constrainer->add<ConstraintNone>("none");
        _constrainer->add<ConstraintFixed>("fixed");
        int N = _sys.mesh().vertices().size();
        _M.resize(N, N);
        _M.reserve(4 * N);
        this->_update_matrix();
        for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
          vh->data().f_type["friction"] = Vec(0, 0);
      }

      void step() override;
      void set_params(const params_type& params) override 
      { 
        for (auto& p : params)
        {
          if (p.first == "gamma")
            _gamma = p.second;
          if (p.first == "zeta")
            _zeta = p.second;
        }
      };
      void set_type_params(const string& type, const params_type& params) override { }
      void set_external_force(const string& vtype, const Vec& f) override 
      { 
        _constant_force[_sys.vert_types()[vtype]] = f;
      }
      void set_radial_force(const string& vtype, double f) override { }
      void set_flag(const string& flag) override 
      {  
        if (flag == "cell_centres")
          _cell_centres = true; 
      }

      bool converged() override { return true; }

      
    private:

      
      ConstrainerType _constrainer; // Apply various constraints
      double _gamma;                // external friction 
      double _zeta;                 // internal friction 
      bool _cell_centres;           // If true, construct fiction matix using cell centres
      vector<Vec> _constant_force;
      SpMat _M;   // Connectivity matrix
      SolType _solver; // Sparse solver
      map<int, int> _vmap;

      void _update_matrix();
  };

}
#else 
 #warning "Eigen3 library needs to be installed to compile the Relative Velocity integrator."
#endif

#endif

