#ifndef __CONSTRAINT_HPP__  
#define __CONSTRAINT_HPP__ 

#include "system.hpp"

namespace RheoVM
{

  class Constraint 
  {
    public: 

      Constraint() { }
      virtual ~Constraint() { }
      virtual Vec apply(const VertexHandle<Property>&, const Vec&) = 0;
      virtual Vec apply(const Vec&) = 0;

  };

}
#endif