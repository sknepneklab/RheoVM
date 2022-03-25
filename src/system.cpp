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
 * \file system.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 12-Jun-2017
 * \brief System class 
 */ 

#include "system.hpp"

namespace RheoVM
{
  vector<string> split_string(const string& input) 
  { 
    istringstream buffer(input);
    vector<std::string> ret;
    copy(istream_iterator<string>(buffer), istream_iterator<string>(), back_inserter(ret));
    return ret;
  }

  string trim(const string& str)
  {
      size_t first = str.find_first_not_of(' ');
      if (string::npos == first) return "";
      size_t last = str.find_last_not_of(' ');
      return str.substr(first, (last - first + 1));
  }

  // Read input from a JSON file
  void System::read_input(const string& json_file, bool read_params)
  {
    if (_mesh_set)
    {
      cout << "Warning! Mesh has already been set. Overwriting it." << endl;
      _mesh.wipe();
    }
    ifstream inp(json_file.c_str());
    json j;
    inp >> j;
    inp.close();
    // Check if simulation box exists
    if (j["mesh"].find("box") != j["mesh"].end())
      if (j["mesh"]["box"]["periodic"])
      {
        if (j["mesh"]["box"].find("a") != j["mesh"]["box"].end() && j["mesh"]["box"].find("b") != j["mesh"]["box"].end())
        {
          vector<double> a = j["mesh"]["box"]["a"];
          vector<double> b = j["mesh"]["box"]["b"];
          this->set_box(make_shared<Box>(a[0], a[1], b[0], b[1]));
        }
        else if (j["mesh"]["box"].find("lx") != j["mesh"]["box"].end() && j["mesh"]["box"].find("ly") != j["mesh"]["box"].end())
        {
          double lx = j["mesh"]["box"]["lx"];
          double ly = j["mesh"]["box"]["ly"];
          this->set_box(make_shared<Box>(lx, ly));
        }
        else 
          throw runtime_error("There is a problem with the box field in the input JSON file.");
        cout << "Setting periodic simulation box." << endl;
      }

    // Check for time step
    if (j["mesh"].find("time_step") != j["mesh"].end())
      _time_step = j["mesh"]["time_step"];

    // Populate vertices
    for (int i = 0; i < j["mesh"]["vertices"].size(); i++)
    {
      int id = j["mesh"]["vertices"][i]["id"];
      double x = j["mesh"]["vertices"][i]["r"][0];
      double y = j["mesh"]["vertices"][i]["r"][1];
      bool boundary = j["mesh"]["vertices"][i]["boundary"];
      _mesh.add_vertex(Vertex<Property>(id, Vec(x, y, _mesh.box()), boundary));
      Vertex<Property> &v = _mesh.vertices().back();
      this->add_vert_type(j["mesh"]["vertices"][i]["type"]);
      v.data().vert_type = _vert_types[j["mesh"]["vertices"][i]["type"]];
      v.data().type_name = get_vert_type_name(v.data().vert_type);
      v.data().constraint = j["mesh"]["vertices"][i]["constraint"];
    }
    cout << "Finished reading vertices." << endl;
    // Populate faces
    for (int i = 0; i < j["mesh"]["faces"].size(); i++)
    {
      _mesh.add_face(j["mesh"]["faces"][i]["vertices"]);
      Face<Property>& f = _mesh.faces().back();
      this->add_cell_type(j["mesh"]["faces"][i]["type"]);
      f.data().face_type = _cell_types[j["mesh"]["faces"][i]["type"]];
      f.data().type_name = get_cell_type_name(f.data().face_type);
      f.outer = j["mesh"]["faces"][i]["outer"];
      if (j["mesh"]["faces"][i].find("A0") != j["mesh"]["faces"][i].end())
        f.data().A0 = j["mesh"]["faces"][i]["A0"];
      if (j["mesh"]["faces"][i].find("P0") != j["mesh"]["faces"][i].end())
        f.data().P0 = j["mesh"]["faces"][i]["P0"];
      HEHandle<Property> he = f.he();
      HEHandle<Property> first = f.he();
      do
      {
        he->data().old_face_id = f.id;
        he = he->next();
      } while (he != first);
    }
    cout << "Finished reading faces." << endl;
    _mesh.tidyup();
    cout << "Finished mesh setup." << endl;
    cout << "Mesh has " << _mesh.vertices().size() << " vertices " << _mesh.edges().size() << " edges and " << _mesh.faces().size() << " faces." << endl;
    // We now need to loop over all edges and set the native length
    if (j["mesh"].find("l0") != j["mesh"].end())
    {
      for (EdgeHandle<Property> eh = _mesh.edges().begin(); eh != _mesh.edges().end(); eh++)
        eh->data().l0 = j["mesh"]["l0"];
      cout << "Finished setting l0." << endl;
    }
    // If read_params flag is set, we loop over faces (i.e., cells) again and read in cell-specific paramters
    if (read_params)
      for (int i = 0; i < j["mesh"]["faces"].size(); i++)
      {
        Face<Property>& f = _mesh.get_face(i);
        if (j["mesh"]["faces"][i].find("kappa") != j["mesh"]["faces"][i].end())
          f.data().kappa = j["mesh"]["faces"][i]["kappa"];
        if (j["mesh"]["faces"][i].find("gamma") != j["mesh"]["faces"][i].end())
          f.data().gamma = j["mesh"]["faces"][i]["gamma"];
        if (j["mesh"]["faces"][i].find("lambda") != j["mesh"]["faces"][i].end())
          f.data().lambda = j["mesh"]["faces"][i]["lambda"];
        if (j["mesh"]["faces"][i].find("beta") != j["mesh"]["faces"][i].end())
          f.data().beta = j["mesh"]["faces"][i]["beta"];
        if (j["mesh"]["faces"][i].find("beta_a") != j["mesh"]["faces"][i].end())
          f.data().beta_a = j["mesh"]["faces"][i]["beta_a"];
        if (j["mesh"]["faces"][i].find("alpha") != j["mesh"]["faces"][i].end())
          f.data().alpha = j["mesh"]["faces"][i]["alpha"];
        if (j["mesh"]["faces"][i].find("k") != j["mesh"]["faces"][i].end())
          f.data().k = j["mesh"]["faces"][i]["k"];
        if (j["mesh"]["faces"][i].find("l0") != j["mesh"]["faces"][i].end())
        {
          HEHandle<Property> he = f.he();
          for (int k = 0; k < j["mesh"]["faces"][i]["l0"].size(); k++)
          {
            he->edge()->data().l0 = j["mesh"]["faces"][i]["l0"][k];
            he = he->next();
          }
        }
        if (j["mesh"]["faces"][i].find("maxA0") != j["mesh"]["faces"][i].end())
          f.data().max_A0 = j["mesh"]["faces"][i]["maxA0"];
        if (j["mesh"]["faces"][i].find("nativeA0") != j["mesh"]["faces"][i].end())
          f.data().native_A0 = j["mesh"]["faces"][i]["nativeA0"];
        else
          f.data().native_A0 = f.data().A0;
      }
      _mesh_set = true;
  }

    
  // Displace vertices of a given type by a given amount
  void System::displace_vertices(const string & vtype, const Vec& dr)
  {
    if (_vert_types.find(vtype) == _vert_types.end())
      throw runtime_error("Unknown vertex type " + vtype + " in displace vertices.");
    int vert_type = _vert_types[vtype];
    for (VertexHandle<Property> vh = _mesh.vertices().begin(); vh != _mesh.vertices().end(); vh++)
      if (vh->data().vert_type == vert_type)
        _mesh.move_vertex(vh->id, dr);
  }

  // Used to be able to make a map of MyoStore
  bool operator<(const VertexHandle<Property>& lv, const VertexHandle<Property>& rv) 
  {
    return (lv->id < rv->id);
  }

  // Python exports

  void export_T1_stats(py::module& m)
  {
    py::class_<T1_stats>(m,"T1Stats")
      .def_readonly("time_step", &T1_stats::time_step)
      .def_readonly("angle", &T1_stats::angle)
      .def_readonly("v1", &T1_stats::v1_id)
      .def_readonly("v2", &T1_stats::v2_id)
      .def_readonly("f1", &T1_stats::f1_id)
      .def_readonly("f2", &T1_stats::f2_id);
  }


  void export_VertexProperty(py::module& m)
  {
    py::class_<Property::VertexProperty>(m, "VertexProperty")
      .def_readonly("type", &Property::VertexProperty::vert_type)
      .def_readonly("vel", &Property::VertexProperty::vel)
      .def_readonly("force", &Property::VertexProperty::force)
      .def_readonly("pre_T1", &Property::VertexProperty::pre_T1)
      .def_readonly("post_T1", &Property::VertexProperty::post_T1)
      .def_readonly("pre_T1_face_pair", &Property::VertexProperty::pre_T1_face_pair)
      .def_readonly("post_T1_face_pair", &Property::VertexProperty::post_T1_face_pair);
  }

  void export_EdgeProperty(py::module& m)
  {
    py::class_<Property::EdgeProperty>(m, "EdgeProperty")
      .def_readwrite("tension", &Property::EdgeProperty::tension)
      .def_readwrite("l0", &Property::EdgeProperty::l0);
  }

  void export_HEProperty(py::module& m)
  {
    py::class_<Property::HEProperty>(m, "HEProperty")
      .def_readonly("tension", &Property::HEProperty::tension)
      .def_readonly("l0", &Property::HEProperty::l0)
      .def_readonly("l0_collapse", &Property::HEProperty::l0_collapse)
      .def_readonly("force_type", &Property::HEProperty::force_type);
  }

  
  void export_FaceProperty(py::module& m)
  {
    py::class_<Property::FaceProperty>(m, "CellProperty")
      .def_readonly("type", &Property::FaceProperty::face_type)
      .def_readonly("type_name", &Property::FaceProperty::type_name)
      .def_readwrite("A0", &Property::FaceProperty::A0)
      .def_readwrite("P0", &Property::FaceProperty::P0)
      .def_readwrite("native_A0", &Property::FaceProperty::native_A0)
      .def_readwrite("native_P0", &Property::FaceProperty::native_P0)
      .def_readwrite("max_A0", &Property::FaceProperty::max_A0)
      .def_readwrite("fa", &Property::FaceProperty::fa)
      .def_readwrite("phi", &Property::FaceProperty::phi)
      .def_readwrite("cell_type", &Property::FaceProperty::face_type)
      .def_readwrite("kappa", &Property::FaceProperty::kappa)
      .def_readwrite("gamma", &Property::FaceProperty::gamma)
      .def_readwrite("lam", &Property::FaceProperty::lambda)
      .def_readwrite("beta", &Property::FaceProperty::beta)
      .def_readwrite("v0", &Property::FaceProperty::v0)
      .def_readwrite("n", &Property::FaceProperty::n);
  }

  
  void export_Vertex(py::module& m)
  {
    py::class_<Vertex<Property>>(m, "Vertex")
      .def(py::init<>())
      .def_readwrite("r", &Vertex<Property>::r)
      .def_readonly("id", &Vertex<Property>::id)
      .def_readonly("erased", &Vertex<Property>::erased)
      .def_readwrite("boundary", &Vertex<Property>::boundary)
      .def_readonly("coordination", &Vertex<Property>::coordination)
      .def("he", [](Vertex<Property>& v) { return *(v.he()); })
      .def("force", [](Vertex<Property>& v) { return v.data().force; })
      .def("he_force", [](Vertex<Property>& v,int idx) { return v.data().he_force[idx]; })     
      .def("force_type", [](Vertex<Property>& v,const std::string& type) { return v.data().f_type[type]; })
      .def("vel", [](Vertex<Property>& v) { return v.data().vel; })
      .def("collapsed", [](Vertex<Property>& v) { return v.data().pre_T1.size(); })          
      .def("property", (Property::VertexProperty& (Vertex<Property>::*)()) &Vertex<Property>::data, py::return_value_policy::reference);
  }

  void export_Edge(py::module& m)
  {
    py::class_<Edge<Property>>(m, "Edge")
      .def(py::init<>())
      .def_readonly("i", &Edge<Property>::i)
      .def_readonly("j", &Edge<Property>::j)
      .def_readonly("erased", &Edge<Property>::erased)
      .def_readonly("boundary", &Edge<Property>::boundary)
      .def_property_readonly("id", &Edge<Property>::idx)
      .def("property", (Property::EdgeProperty& (Edge<Property>::*)()) &Edge<Property>::data, py::return_value_policy::reference);
  }

  void export_HalfEdge(py::module& m)
  {
    py::class_<HalfEdge<Property>>(m, "Junction")
      .def(py::init<>())
      .def("vfrom", [](HalfEdge<Property>& he) { return *(he.from()); })
      .def("vto", [](HalfEdge<Property>& he) { return *(he.to()); })
      .def("edge", [](HalfEdge<Property>& he) { return *(he.edge()); })
      .def("id", &HalfEdge<Property>::idx)
      .def("next", [](HalfEdge<Property>& he) { return *(he.next()); })
      .def("prev", [](HalfEdge<Property>& he) { return *(he.prev()); })
      .def("pair", [](HalfEdge<Property>& he) { return *(he.pair()); })
      .def("face", [](HalfEdge<Property>& he) { return *(he.face()); })
      .def("property", (Property::HEProperty& (HalfEdge<Property>::*)()) &HalfEdge<Property>::data, py::return_value_policy::reference);
  }

  void export_Face(py::module& m)
  {
    py::class_<Face<Property>>(m, "Cell")
      .def(py::init<>())
      .def_readonly("id", &Face<Property>::id)
      .def_readonly("neighbours", &Face<Property>::nsides)
      .def_readonly("outer", &Face<Property>::outer)
      .def("type", [](Face<Property>& f) { return f.data().face_type; })
      .def("A0", [](Face<Property>& f) { return f.data().A0; })
      .def("P0", [](Face<Property>& f) { return f.data().P0; })
      .def("maxA0", [](Face<Property>& f) { return f.data().max_A0; })
      // .def("he", (HEHandle<Property>& (Face<Property>::*)()) &Face<Property>::he, py::return_value_policy::reference)
      .def("he", [](Face<Property>& f) { return *(f.he()); })
      .def("property", (Property::FaceProperty& (Face<Property>::*)()) &Face<Property>::data, py::return_value_policy::reference);
  }

  void export_Mesh(py::module& m)
  {
    py::class_<Mesh<Property>>(m, "Tissue")
      .def(py::init<>())
      .def("num_vert", &Mesh<Property>::num_vert)
      .def("num_cells", &Mesh<Property>::num_faces)
      .def("tidyup", &Mesh<Property>::tidyup)
      .def("set_cell_type", [](Mesh<Property>& m, int i, int type) { m.get_face(i).data().face_type = type; })
      .def("set_cell_A0", [](Mesh<Property>& m, int i, double A0) { m.get_face(i).data().A0 = A0; })
      .def("set_cell_P0", [](Mesh<Property>& m, int i, double P0) { m.get_face(i).data().P0 = P0; })
      .def("get_vertex", &Mesh<Property>::get_vertex , py::return_value_policy::reference)
      .def("get_junction", &Mesh<Property>::get_halfedge, py::return_value_policy::reference)
      .def("get_cell", &Mesh<Property>::get_face, py::return_value_policy::reference)
      .def("vertices", &Mesh<Property>::vertices, py::return_value_policy::reference)
      .def("junctions", &Mesh<Property>::edges, py::return_value_policy::reference)
      .def("halfedges", &Mesh<Property>::halfedges, py::return_value_policy::reference)
      .def("cells", &Mesh<Property>::faces, py::return_value_policy::reference)
      .def("move_vertex", &Mesh<Property>::move_vertex)
      .def("scale", &Mesh<Property>::scale)
      .def("shear", &Mesh<Property>::shear)
      .def("transform", &Mesh<Property>::transform, py::arg("txx"), py::arg("txy"), py::arg("tyx"), py::arg("tyy"), py::arg("undo") = false)
      .def("get_cell_centre", [](Mesh<Property>& m, int i) -> Vec { return m.get_face_centre(find_if(m.faces().begin(),m.faces().end(),[i](const Face<Property>& f) -> bool { return (f.id == i); })); })
      .def("get_cell_centroid", [](Mesh<Property>& m, int i) -> Vec { return m.get_face_centroid(find_if(m.faces().begin(),m.faces().end(),[i](const Face<Property>& f) -> bool { return (f.id == i); })); })
      .def("get_centre", &Mesh<Property>::get_centre)
      .def("cell_area", [](Mesh<Property> &m, int i) -> double { return m.area(find_if(m.faces().begin(), m.faces().end(), [i](const Face<Property> &f) -> bool { return (f.id == i); })); })
      .def("cell_perim", [](Mesh<Property> &m, int i) -> double { return m.perim(find_if(m.faces().begin(), m.faces().end(), [i](const Face<Property> &f) -> bool { return (f.id == i); })); });
  }

  void export_System(py::module& m)
  {
    py::class_<System>(m, "System")
      .def(py::init<MyMesh&>())
      .def("read_input", &System::read_input, py::arg("input_file"), py::arg("read_params") = false)
      .def("mesh", &System::mesh)
      .def("cell_types", &System::cell_types)
      .def("time_step", &System::time_step)
      .def("simulation_time", &System::simulation_time)
      .def("set_simulation_time_step", &System::set_simulation_time_step)
      .def("get_cell_type_name", &System::get_cell_type_name)
      .def("get_vert_type_name", &System::get_vert_type_name)
      .def("displace_vertices", &System::displace_vertices);
  }

  
}

