#ifndef __CONSTRINER_HPP__
#define __CONSTRINER_HPP__

#include "class_factory.hpp"

#include "class_factory.hpp"
#include "constraint.hpp"

namespace RheoVM
{

  class Constrainer : public ClassFactory<Constraint>
  {
    public:

      Constrainer()  = default; 
      ~Constrainer() = default; 

      Constrainer(const Constrainer&) = delete;

      Vec apply_vertex(const VertexHandle<Property> &vh, const Vec& fc)
      {
        if (this->factory_map.find(vh->data().constraint) != this->factory_map.end())
          return this->factory_map[vh->data().constraint]->apply(vh, fc);
        else 
          return Vec(fc);
      }

      Vec apply_vector(const VertexHandle<Property> &vh, const Vec& fc)
      {
        if (this->factory_map.find(vh->data().constraint) != this->factory_map.end())
          return this->factory_map[vh->data().constraint]->apply(fc);
        else 
          return Vec(fc);
      }
  };

}
#endif
