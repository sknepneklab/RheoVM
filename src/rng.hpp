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
 * \file rng.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 24-Oct-2013
 * \brief Class RNG provides wrappers for the GSL random number generate
 */ 

#ifndef __RNG_H__
#define __RNG_H__

#include <random>

namespace RheoVM
{

  /*! Class handles random numbers in the system */
  class RNG
  {
  public:
    
    //! Constructor (initialize random number generator)
    RNG(unsigned int seed) : _generator(seed), _uniform_distribution(0.0,1.0), _normal_distribution(0.0,1.0) { }
    
    //! Destructor
    ~RNG() { }
    
    //! Return random number between 0 and 1
    double drnd();
    
    //! Return random integer between 0 and N
    int lrnd(int);
    
    //! Return a Gaussian distributed number with a given standard deviation
    double gauss_rng();

    //! Return a Gaussian random number from a normal distribution with a given mean and standard deviation
    double normal(double, double);

  private:
    
    std::mt19937_64 _generator;  //!< Mersenne Twister engine 
    std::uniform_real_distribution<double> _uniform_distribution;  // Uniform random numbers
    std::normal_distribution<double> _normal_distribution; // Gaussian distribution zero mean, unit variance
    
  };

}

#endif
