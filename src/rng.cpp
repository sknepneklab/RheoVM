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
 * \file rng.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 24-Oct-2013
 * \brief Wrappers for the GSL random number generate
 */ 

#include "rng.hpp"

namespace RheoVM
{
  //! Get a random number between 0 and 1 drawn from an uniform distribution
  //! \return random number between 0 and 1
  double RNG::drnd()
  {
    return _uniform_distribution(_generator);
  }

  //! Return a random number from a Gaussian distribution with a given standard deviation 
  //! \param sigma standard deviation 
  double RNG::gauss_rng()
  {
    return _normal_distribution(_generator);
  }

  //! Return a random number from a Gaussian distribution with a given mean and standard deviation 
  double RNG::normal(double avg, double std)
  {
    return avg + std*this->gauss_rng();
  }

  //! Get an integer random number between 0 and N drawn from an uniform distribution
  //! \param N upper bound for the interval
  //! \return integer random number between 0 and N
  int RNG::lrnd(int N)
  {
    return static_cast<int>(N*drnd());
  }
}


