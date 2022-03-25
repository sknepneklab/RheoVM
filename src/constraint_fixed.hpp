#ifndef __CONSTRAINT_FIXED_HPP__  
#define __CONSTRAINT_FIXED_HPP__ 

#include "constraint.hpp"

namespace RheoVM
{

  class ConstraintFixed : public Constraint
  {
    public: 

      ConstraintFixed() { }
      Vec apply(const VertexHandle<Property>& vh, const Vec& f) override
      {
        return Vec(0.0,0.0);
      }
      Vec apply(const Vec& f) override
      {
        return Vec(0.0,0.0);
      }
  };

}
#endif