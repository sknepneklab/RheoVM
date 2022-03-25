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
 * \file dump.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief Dump class 
 */ 

#include "dump.hpp"

namespace RheoVM
{
  void Dump::dump_cells(const string& vtk_file, bool binary_output, bool draw_periodic)
  {
    vtkSmartPointer<vtkPolyData> polydata =  vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkDoubleArray> forces =  vtkSmartPointer<vtkDoubleArray>::New();
    std::map<std::string, vtkSmartPointer<vtkDoubleArray> > force_comp;
    for (auto& fc : _force_compute.factory())
      force_comp[fc.first] = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkPolygon> face =  vtkSmartPointer<vtkPolygon>::New();
    vtkSmartPointer<vtkIntArray> ids =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> face_ids =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> true_ids =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> vert_type =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkDoubleArray> areas =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> perims =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkIntArray> cell_types =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkDoubleArray> stress =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> time =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> num_neigh =  vtkSmartPointer<vtkDoubleArray>::New();
    
    ids->SetName("Id");
    ids->SetNumberOfComponents(1);
    true_ids->SetName("TrueId");
    true_ids->SetNumberOfComponents(1);
    areas->SetName("Area");
    areas->SetNumberOfComponents(1);
    perims->SetName("Perimeter");
    perims->SetNumberOfComponents(1);
    face_ids->SetName("FaceId");
    face_ids->SetNumberOfComponents(1);
    vert_type->SetName("VertType");
    vert_type->SetNumberOfComponents(1);
    forces->SetName("Force");
    forces->SetNumberOfComponents(3);
    stress->SetName("Stress");
    stress->SetNumberOfComponents(9);
    time->SetName("time");
    time->SetNumberOfComponents(1);
    time->InsertNextValue(_sys.simulation_time());
    num_neigh->SetName("NumNeigh");
    num_neigh->SetNumberOfComponents(1);
    
    for (auto& fc : _force_compute.factory())
    {
      force_comp[fc.first]->SetName(("force_"+fc.first).c_str());
      force_comp[fc.first]->SetNumberOfComponents(3);
    }
    
    cell_types->SetName("CellTypes");
    cell_types->SetNumberOfComponents(1);
    
    int id = 0;
    map<int,int> id_map;
    for (auto vh : _sys.mesh().vertices())
    {
      if (!vh.erased)
      {
        points->InsertNextPoint(vh.r.x, vh.r.y, 0.0);
        ids->InsertNextValue(id);
        true_ids->InsertNextValue(vh.id);
        vert_type->InsertNextValue(vh.data().vert_type);
        double f[3] = {vh.data().force.x, vh.data().force.y, 0.0};
        forces->InsertNextTuple(f);
        for (auto& fc : _force_compute.factory())
        {
          Vec fcomp = vh.data().f_type[fc.first];
          double fcomponent[3] = {fcomp.x, fcomp.y, 0.0};
          force_comp[fc.first]->InsertNextTuple(fcomponent);
        }
        id_map[vh.id] = id;
        id++;
      }
    }
    
    for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
    {
      if (!fh->erased)
      {
        double alpha=0;
        bool omit = false;
        bool crosses_pbc = false;
        vector<int> face_verts;
        if (_sys.periodic())
        {
          HEHandle<Property> he = fh->he();
          HEHandle<Property> first = he;
          Vec r0 = _sys.mesh().get_face_centre(fh);
          Vec s0 = r0.box->inv_h*r0;
          do
          {
            Vec s = r0.box->inv_h*he->from()->r;
            double sx = s.x - s0.x, sy = s.y - s0.y;
            if (fabs(rint(sx)) != 0.0 || fabs(rint(sy)) != 0.0) 
            {
              if (!draw_periodic)
              {
                omit = true;
                break;
              }
              else
              {
                crosses_pbc = true;
                Vec r = r0.box->h*Vec(s.x - rint(sx), s.y - rint(sy));
                points->InsertNextPoint(r.x, r.y, 0.0);
                ids->InsertNextValue(id);
                true_ids->InsertNextValue(he->from()->id);
                vert_type->InsertNextValue(he->from()->data().vert_type);
                double f[3] = {0.0, 0.0, 0.0};
                forces->InsertNextTuple(f);
                for (auto& fc : _force_compute.factory())
                {
                  double fcomponent[3] = {0.0, 0.0, 0.0};
                  force_comp[fc.first]->InsertNextTuple(fcomponent);
                }
                face_verts.push_back(id);
                id++;
              }
            }
            else
              face_verts.push_back(he->from()->id);
            he = he->next();
          } while (he != first);
        }
        if (!(fh->outer || omit))
        {
          face->GetPointIds()->SetNumberOfIds(_sys.mesh().face_sides(fh));
          HEHandle<Property> he = fh->he();
          HEHandle<Property> first = he;
          
          int i = 0;
          do
          {
            double le = (he->from()->r-he->to()->r).len();
            if (!crosses_pbc)
              face->GetPointIds()->SetId(i++,id_map[he->from()->id]);
            else 
            {
              face->GetPointIds()->SetId(i,face_verts[i]);
              i++;
            }
            VertexHandle<Property> vh = he->from();
            
            he = he->next();
          } while (he != first);
          faces->InsertNextCell(face);
          areas->InsertNextValue(_sys.mesh().area(fh));
          perims->InsertNextValue(_sys.mesh().perim(fh));
          face_ids->InsertNextValue(fh->id);
          cell_types->InsertNextValue(fh->data().face_type);
          if (_force_compute.stress_compute())
            stress->InsertNextTuple9(fh->data().stress[0], fh->data().stress[1], 0.0, fh->data().stress[2], fh->data().stress[3], 0.0, 0.0, 0.0, 0.0);
          num_neigh->InsertNextValue(fh->nsides);
        }
      }
    }

    polydata->SetPoints(points);
    polydata->GetPointData()->AddArray(ids);
    polydata->GetPointData()->AddArray(true_ids);
    polydata->GetPointData()->AddArray(vert_type);
    polydata->GetPointData()->AddArray(forces);
    
    for (auto& fc : _force_compute.factory())
      polydata->GetPointData()->AddArray(force_comp[fc.first]);
    
    polydata->SetPolys(faces);
    polydata->GetCellData()->AddArray(areas);
    polydata->GetCellData()->AddArray(perims);
    polydata->GetCellData()->AddArray(face_ids);
    polydata->GetCellData()->AddArray(cell_types);
    polydata->GetCellData()->AddArray(num_neigh);
    if (_force_compute.stress_compute())
      polydata->GetCellData()->AddArray(stress);
    polydata->GetFieldData()->AddArray(time);
    

    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =  vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(vtk_file.c_str());
    #if VTK_MAJOR_VERSION <= 5
      writer->SetInput(polydata);
    #else
      writer->SetInputData(polydata);
    #endif
    if (binary_output)
    {
      writer->SetDataModeToBinary();
      vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();   
      compressor->SetCompressionLevel(9); // max compression level
      writer->SetCompressor(compressor);
    }
    else
      writer->SetDataModeToAscii();
    writer->Write();

  }

  void Dump::dump_junctions(const string& vtk_file, bool binary_output)
  {
    vtkSmartPointer<vtkPolyData> polydata =  vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines =  vtkSmartPointer<vtkCellArray>::New();
    
    vtkSmartPointer<vtkLine> edge =  vtkSmartPointer<vtkLine>::New();
    vtkSmartPointer<vtkIntArray> ids =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> he_id =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> edge_id =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkDoubleArray> lens =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> tension =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> l0 =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkIntArray> true_ids =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkDoubleArray> time =  vtkSmartPointer<vtkDoubleArray>::New();
    
    ids->SetName("Id");
    ids->SetNumberOfComponents(1);
    he_id->SetName("HE_Id");
    he_id->SetNumberOfComponents(1);
    edge_id->SetName("Edge_Id");
    edge_id->SetNumberOfComponents(1);
    lens->SetName("Length");
    lens->SetNumberOfComponents(1);
    tension->SetName("Tension");
    tension->SetNumberOfComponents(1);
    l0->SetName("l0");
    l0->SetNumberOfComponents(1);
    true_ids->SetName("TrueId");
    true_ids->SetNumberOfComponents(1);
    time->SetName("time");
    time->SetNumberOfComponents(1);
    time->InsertNextValue(_sys.simulation_time());
    
    int id = 0;
    map<int,int> id_map;
    map<int, bool> omit;
    
    // build list of points belonging to faces, with shift (for display purpose)
    for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
    {
      if (!fh->erased)
      {
        HEHandle<Property>   he = fh->he();
        HEHandle<Property>   first = he;
        Vec fc = _sys.mesh().get_face_centre(fh); // face centre coordinates

        do
        {
          EdgeHandle<Property> eh = he->edge();
          if (!eh->erased)
          {		
            VertexHandle<Property> vh1 = he->from();
            VertexHandle<Property> vh2 = he->to();
            Vec r1_shift = fc + _sfc*(vh1->r - fc); 
            Vec r2_shift = fc + _sfc*(vh2->r - fc);
            if (fh->outer)
            {
              r1_shift = vh1->r; 
              r2_shift = vh2->r;
            }
            omit[he->idx()] = false;
            if (_sys.periodic())
            {
              Vec s1 = r1_shift.box->inv_h*r1_shift;
              Vec s2 = r2_shift.box->inv_h*r2_shift;
              double sx = s1.x - s2.x, sy = s1.y - s2.y;
              if (fabs(rint(sx)) != 0.0 || fabs(rint(sy)) != 0.0)  
                omit[he->idx()] = true;
            }
            if (!(omit[he->idx()]))
            {
              points->InsertNextPoint(r1_shift.x, r1_shift.y, 0.0);
              ids->InsertNextValue(id);
              true_ids->InsertNextValue(vh1->id);
              id_map[vh1->id] = id;
              id++;
          
              points->InsertNextPoint(r2_shift.x, r2_shift.y, 0.0);
              ids->InsertNextValue(id);
              true_ids->InsertNextValue(vh2->id);
              id_map[vh2->id] = id;
              id++;
            }	       
          }
          he = he->next();
        } while (he !=first);
      }
    }

	  polydata->SetPoints(points);
	  polydata->GetPointData()->AddArray(ids);
	  polydata->GetPointData()->AddArray(true_ids);

	  //create connectivity 
	  id = 0;
	  for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
    {
      if (!fh->erased)
      {
        HEHandle<Property> he = fh->he();
        HEHandle<Property> first = he;
        do
        {
          EdgeHandle<Property> eh = he->edge();
          if (!(eh->erased || omit[he->idx()]))
          {		      
            edge->GetPointIds()->SetId(0,id);
            id++;
            edge->GetPointIds()->SetId(1, id);
            id++;
            he_id->InsertNextValue(he->idx());
            edge_id->InsertNextValue(eh->idx());
            lines->InsertNextCell(edge);
            lens->InsertNextValue(_sys.mesh().len(eh));
            double he_tension = he->data().tension;
            tension->InsertNextValue(he_tension);
            l0->InsertNextValue(eh->data().l0);		  		      
          }
          he = he->next();
        } while (he != first);
      }
    }
    
    polydata->SetLines(lines);
    polydata->GetCellData()->AddArray(lens);
    polydata->GetCellData()->AddArray(he_id);
    polydata->GetCellData()->AddArray(edge_id);
    polydata->GetCellData()->AddArray(tension);
    polydata->GetCellData()->AddArray(l0);
    polydata->GetFieldData()->AddArray(time);


    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =  vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(vtk_file.c_str());
    #if VTK_MAJOR_VERSION <= 5
      writer->SetInput(polydata);
    #else
      writer->SetInputData(polydata);
    #endif
    if (binary_output)
    {
      writer->SetDataModeToBinary();
      vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();   
      compressor->SetCompressionLevel(9); // max compression level
      writer->SetCompressor(compressor);
    }
    else
      writer->SetDataModeToAscii();
    writer->Write();

  }

  void Dump::dump_box(const string& vtk_file, bool binary_output)
  {
    if (!_sys.box())
      throw runtime_error("Simulation box not defined.");
    
    vtkSmartPointer<vtkPolyData> polydata =  vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines =  vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkLine> edge =  vtkSmartPointer<vtkLine>::New();

    Vec r1 = _sys.box()->h*Vec(1,0,nullptr);
    Vec r2 = _sys.box()->h*Vec(0,1,nullptr);
    double xc = 0.5*(r1.x + r2.x), yc = 0.5*(r1.y + r2.y);


    points->InsertNextPoint(-xc, -yc, 0.0);
    points->InsertNextPoint(-xc+r1.x, -yc+r1.y, 0.0);
    points->InsertNextPoint(-xc+r1.x+r2.x, -yc+r1.y+r2.y, 0.0);
    points->InsertNextPoint(-xc+r2.x, -yc+r2.y, 0.0);

    polydata->SetPoints(points);

    edge->GetPointIds()->SetId(0,0);
    edge->GetPointIds()->SetId(1,1);
    lines->InsertNextCell(edge);
    edge->GetPointIds()->SetId(0,1);
    edge->GetPointIds()->SetId(1,2);
    lines->InsertNextCell(edge);
    edge->GetPointIds()->SetId(0,2);
    edge->GetPointIds()->SetId(1,3);
    lines->InsertNextCell(edge);
    edge->GetPointIds()->SetId(0,3);
    edge->GetPointIds()->SetId(1,0);
    lines->InsertNextCell(edge);

    polydata->SetLines(lines);
    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =  vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(vtk_file.c_str());
    #if VTK_MAJOR_VERSION <= 5
      writer->SetInput(polydata);
    #else
      writer->SetInputData(polydata);
    #endif
    if (binary_output)
    {
      writer->SetDataModeToBinary();
      vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();   
      compressor->SetCompressionLevel(9); // max compression level
      writer->SetCompressor(compressor);
    }
    else
      writer->SetDataModeToAscii();
    writer->Write();
  }

  
  void Dump::dump_mesh(const string& mesh_file, bool copy_params)
  {
    if (copy_params)
      _force_compute.copy_type_param_to_cell();
    vector<string> strs = split(mesh_file, '.');
    string ext = strs[strs.size()-1];
    transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c){ return tolower(c); });
    pt::ptree out;
    pt::ptree mesh;
    pt::ptree vertices;
    pt::ptree faces;
    if (_sys.periodic())
    {
      pt::ptree box;
      pt::ptree lx, ly;
      pt::ptree A, B;
      pt::ptree periodic;
      pt::ptree mxx, mxy, myx, myy;
      mxx.put("", _sys.mesh().box()->h._mxx);
      myx.put("", _sys.mesh().box()->h._myx);
      A.push_back(std::make_pair("", mxx));
      A.push_back(std::make_pair("", myx));
      mxy.put("", _sys.mesh().box()->h._mxy);
      myy.put("", _sys.mesh().box()->h._myy);
      B.push_back(std::make_pair("", mxy));
      B.push_back(std::make_pair("", myy));
      lx.put("", _sys.mesh().box()->h._mxx);
      ly.put("", _sys.mesh().box()->h._myy);
      if (fabs(_sys.mesh().box()->h._mxy) >= 1e-6 || fabs(_sys.mesh().box()->h._myx) >= 1e-6)
      {
        box.add_child("a", A);
        box.add_child("b", B);
      }
      else
      {
        box.add_child("lx", lx);
        box.add_child("ly", ly);
      }
      periodic.put("", "true");
      box.add_child("periodic", periodic);
      mesh.add_child("box", box);
    }
    for (auto vh : _sys.mesh().vertices())
    {
      pt::ptree vertex;
      pt::ptree id;
      pt::ptree boundary;
      pt::ptree x;
      pt::ptree y;
      pt::ptree type;
      pt::ptree r;
      pt::ptree constraint;
      pt::ptree erased;
      id.put("", vh.id);
      boundary.put("", vh.boundary);
      type.put("", vh.data().type_name);
      x.put("", vh.r.x);
      y.put("", vh.r.y);
      r.push_back(std::make_pair("", x));
      r.push_back(std::make_pair("", y));
      constraint.put("", vh.data().constraint);
      erased.put("", vh.erased);
      vertex.add_child("id", id);
      vertex.add_child("boundary", boundary);
      vertex.add_child("type", type);
      vertex.add_child("r", r);
      vertex.add_child("constraint", constraint);
      vertex.add_child("erased", erased);
      vertices.push_back(std::make_pair("", vertex));
    }
    mesh.add_child("vertices", vertices);
    _sys.mesh().populate_face_neighbours();
    int face_id = 0;
    for (auto fh : _sys.mesh().faces())
    {
      if (!fh.erased)
      {
        vector<int> verts;
        HECHandle<Property> he = fh.he();
        HECHandle<Property> first = he;
        vector<double> l0,  tension;
        do
        {
          verts.push_back(he->from()->id);
          l0.push_back(he->edge()->data().l0);
          Vec er = he->to()->r - he->from()->r;
          Vec fr = he->pair()->face()->data().rc - fh.data().rc;
          tension.push_back(he->data().tension);
        } while ((he = he->next()) != first);
        pt::ptree face;
        pt::ptree id;
        pt::ptree original_id;
        pt::ptree outer;
        pt::ptree nsides;
        pt::ptree type;
        pt::ptree A0;
        pt::ptree P0;
        pt::ptree fverts;
        pt::ptree fneighs;
        pt::ptree rc;
        pt::ptree rc_x, rc_y;
        pt::ptree kappa, gamma, lambda;
        pt::ptree beta, beta_a, alpha;
        pt::ptree k;
        pt::ptree fl0;
        pt::ptree ftension;

        id.put("", face_id++);
        original_id.put("", fh.id);
        outer.put("", fh.outer);
        nsides.put("", fh.nsides);
        type.put("", fh.data().type_name);
        A0.put("", fh.data().A0);
        P0.put("", fh.data().P0);
        for (int vid : verts)
        {
          pt::ptree pvid;
          pvid.put("", vid);
          fverts.push_back(std::make_pair("", pvid));
        }
        for (int nid : fh.data().neighs)
        {
          pt::ptree pnid;
          pnid.put("", nid);
          fneighs.push_back(std::make_pair("", pnid));
        }
        rc_x.put("", fh.data().rc.x);
        rc_y.put("", fh.data().rc.y);
        rc.push_back(std::make_pair("", rc_x));
        rc.push_back(std::make_pair("", rc_y));
        kappa.put("", fh.data().kappa);
        gamma.put("", fh.data().gamma);
        lambda.put("", fh.data().lambda);
        beta.put("", fh.data().beta);
        beta_a.put("", fh.data().beta_a);
        alpha.put("", fh.data().alpha);
        k.put("", fh.data().k);

        face.add_child("id", id);
        face.add_child("original_id", original_id);
        face.add_child("outer", outer);
        face.add_child("nsides", nsides);
        face.add_child("type", type);
        face.add_child("A0", A0);
        face.add_child("P0", P0);
        face.add_child("vertices", fverts);
        face.add_child("neighbours", fneighs);
        face.add_child("rc", rc);
        face.add_child("kappa", kappa);
        face.add_child("gamma", gamma);
        face.add_child("lambda", lambda);
        face.add_child("beta", beta);
        face.add_child("beta_a", beta_a);
        face.add_child("alpha", alpha);
        face.add_child("k", k);
        for (double ll0 : l0)
        {
          pt::ptree pl0;
          pl0.put("", ll0);
          fl0.push_back(std::make_pair("", pl0));
        }
        face.add_child("l0", fl0);
        for (double t : tension)
        {
          pt::ptree pt;
          pt.put("", t);
          ftension.push_back(std::make_pair("", pt));
        }
        face.add_child("tension", ftension);
        faces.push_back(std::make_pair("", face));
      }
    }
    mesh.add_child("faces", faces);

    // save some system properties in the json file
    pt::ptree system;
    pt::ptree time;
    pt::ptree time_step;
    time.put("", _sys.simulation_time());
    system.add_child("time", time);
    time_step.put("", _sys.time_step());
    system.add_child("time_step", time_step);
    out.add_child("system",system);

    out.add_child("mesh", mesh);

    if (ext == "json")
    {
      std::ostringstream oss;
      pt::write_json(oss, out);
      std::regex reg("\\\"([+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?))\\\"");
      //std::regex reg("\\\"(\\-{0,1}[0-9]+(\\.[0-9]+){0,1})\\\"");
      std::regex reg_bool("\\\"(true|false)\\\"");
      std::string result = std::regex_replace(oss.str(), reg, "$1");
      result = std::regex_replace(result, reg_bool, "$1");

      std::ofstream file;
      file.open(mesh_file);
      file << result;
      file.close();
    }
    else if (ext == "xml")
      pt::write_xml(mesh_file, out, std::locale(), pt::xml_writer_make_settings<string>(' ', 4));
    else
      pt::write_json(mesh_file,out);

  }

  void Dump::dump_infiles(const string& save_file)
  {
    ofstream myverts;
    myverts.open ("vertices_"+save_file+"_save.dat");

    for (auto vh : _sys.mesh().vertices())
    {
      myverts << vh.id << " " << vh.r.x <<  " " << vh.r.y << endl;
    }
    myverts.close();
    
    ofstream mycellverts;
    ofstream myboundverts;    
    ofstream mytypes;    
    mycellverts.open ("cells_"+save_file+"_save.dat");
    myboundverts.open ("boundary_"+save_file+"_save.dat");
    mytypes.open ("types_"+save_file+"_save.dat");
    
    for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
    {
      mytypes << fh->data().face_type << endl;
      HEHandle<Property> he = fh->he();
      HEHandle<Property> first = he;
      do
      {
        mycellverts << he->from()->id << " ";
        if(fh == --_sys.mesh().faces().end())
            myboundverts << he->from()->id << endl;
        he = he->next();
      } while (he != first);
      mycellverts << endl << " ";
    }
    mycellverts.close();
    myboundverts.close();
    mytypes.close();

  }
  
  
  void Dump::dump_area(const string& save_file)
  {
    ofstream myarea;
    myarea.open(save_file);
    for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
    {
      // not the outside face, though!
      // since there does not appear to be a clean label for this,
        if (!fh->outer) 
        {
          myarea << _sys.mesh().area(fh) << "     " << _sys.mesh().perim(fh) << endl;
        }
    }
    
    myarea.close();

  }
  
  
  // JSON output 
  void Dump::dump_json(const string& json_file)
  {
    ofstream jout(json_file.c_str());
    json j;
    _sys.mesh().populate_face_neighbours();
    if (_sys.periodic())
      j["mesh"]["box"] = *(_sys.mesh().box());
    j["mesh"]["time_step"] = _sys.time_step();
    j["mesh"]["vertices"] = _sys.mesh().vertices();
    //j["mesh"]["edges"] = _sys.mesh().edges();
    j["mesh"]["faces"] = _sys.mesh().faces();
    //j["mesh"]["half_edges"] = _sys.mesh().halfedges();
    jout << setw(2) << j << endl;
    jout.close();
  }

  // Dump vertices
  void Dump::dump_vertices(const list<int>& vert_list, const string& save_file)
  {
    ofstream vert_file;
    vert_file.open(save_file);
    vert_file << "# vert_id  x  y" << endl;
    for (auto v: _sys.mesh().vertices())
      if (find(vert_list.begin(), vert_list.end(), v.id) != vert_list.end())
        vert_file << v.id << "   " << v.r.x << "   " << v.r.y << endl;
    vert_file.close();
  }

  void Dump::dump_stress(const string& stress_file, bool average, const string& cell_type)
  {
    if (!_force_compute.stress_compute())
      throw runtime_error("Compute stress flag needs to be set to true in order to compute and print cell stresses.");
    ofstream stressf;
    bool write_header = false;
    struct stat buffer;   
    if (stat(stress_file.c_str(), &buffer) != 0)
      write_header = true;
    stressf.open(stress_file, std::ofstream::app);
    if (write_header)
    {
      if (average)
        stressf << "# time_step s_xx  s_xy  s_yx  s_yy  sa_xx  sa_xy  sa_yx  sa_yy  sp_xx  sp_xy  sp_yx  sp_yy  sv_xx  sv_xy  sv_yx  sv_yy " << endl;
      else
        stressf << "# face_id  xc  yc  s_xx  s_xy  s_yx  s_yy" << endl;
    }
    if (cell_type != "")
      if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
          throw runtime_error("Dump stress: Cell type " + cell_type + " is not defined.");
    int tp = cell_type != "" ? _sys.cell_types()[cell_type] : -1;
    double area = 0.0;
    double Sxx = 0.0, Sxy = 0.0, Syx = 0.0, Syy = 0.0;
    double Saxx = 0.0, Saxy = 0.0, Sayx = 0.0, Sayy = 0.0;
    double Spxx = 0.0, Spxy = 0.0, Spyx = 0.0, Spyy = 0.0;
    double Svxx = 0.0, Svxy = 0.0, Svyx = 0.0, Svyy = 0.0;
    bool include_cell = true;
    for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
    {
      if (cell_type != "" && fh->data().face_type != tp)
        include_cell = false;
      else
        include_cell = true;
      if (!(fh->erased || fh->outer))
      {
        Vec rc = _sys.mesh().get_face_centre(fh);
        if (fh->data().stress.size() == 0)
        {
          if (average && include_cell)
          {
            double Ac = _sys.mesh().area(fh);
            area += Ac;
          }
          else if (!average)
          {
            stressf << fh->id << "   " << rc.x << "   " << rc.y << "  " << 0.0 << "  " << 0.0 << "  " << 0.0 << "  " << 0.0 << endl;
          }
        }
        else
        {
          if (average && include_cell)
          {
            double Ac = _sys.mesh().area(fh);
            area += Ac;
            Sxx += Ac*fh->data().stress[0];
            Sxy += Ac*fh->data().stress[1];
            Syx += Ac*fh->data().stress[2];
            Syy += Ac*fh->data().stress[3];
            Saxx += Ac*fh->data().stress_a[0];
            Saxy += Ac*fh->data().stress_a[1];
            Sayx += Ac*fh->data().stress_a[2];
            Sayy += Ac*fh->data().stress_a[3]; 
            Spxx += Ac*fh->data().stress_p[0];
            Spxy += Ac*fh->data().stress_p[1];
            Spyx += Ac*fh->data().stress_p[2];
            Spyy += Ac*fh->data().stress_p[3];
            Svxx += Ac*fh->data().stress_v[0];
            Svxy += Ac*fh->data().stress_v[1];
            Svyx += Ac*fh->data().stress_v[2];
            Svyy += Ac*fh->data().stress_v[3];  
          }
          else if (!average)
          {
            stressf << fh->id << "   " << rc.x << "   " << rc.y << "  " << fh->data().stress[0] << "  " << fh->data().stress[1] << "  " << fh->data().stress[2] << "  " << fh->data().stress[3] << endl;
          }
        }
      }
    }
    if (average)
      stressf << _sys.time_step() << "  " << Sxx/area << "  " << Sxy/area << "  " << Syx/area << "  " << Syy/area  
                                  << "  " << Saxx/area << "  " << Saxy/area << "  " << Sayx/area << "  " << Sayy/area
                                  << "  " << Spxx/area << "  " << Spxy/area << "  " << Spyx/area << "  " << Spyy/area
                                  << "  " << Svxx/area << "  " << Svxy/area << "  " << Svyx/area << "  " << Svyy/area
                                  << endl;
    stressf.close();
  }

  void Dump::dump_energy(const string& energy_file)
  {
    ofstream engf;
    struct stat buffer;
    bool write_header = false;   
    if (stat(energy_file.c_str(), &buffer) != 0)
      write_header = true;
    engf.open(energy_file, std::ofstream::app);
    if (write_header)
    {
      engf << "# time_step energy ";
      FaceHandle<Property> fh = _sys.mesh().faces().begin();
      map<string,double> eng = _force_compute.energy_comp(fh);
      for (const auto& kv : eng)
        engf << kv.first << " ";
      engf << endl;
    }   
    double totE = 0.0;
    map<string, double> tot_eng;
    for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
    {
      if (!(fh->outer || fh->erased))
      {
        totE += _force_compute.energy(fh);
        map<string, double> eng = _force_compute.energy_comp(fh);
        for (const auto& kv : eng)
          if (tot_eng.find(kv.first) == tot_eng.end())
            tot_eng[kv.first] = kv.second;
          else
            tot_eng[kv.first] += kv.second;
      }
    }
    engf << std::setw(7) << _sys.time_step() << "  " << std::setprecision(16) << std::scientific << totE << "  ";
    for (const auto& kv : tot_eng)
      engf << std::setprecision(16) << std::scientific << kv.second << " ";
    engf << endl;
    engf.close();
  }

  // Dump data
  void Dump::dump_data(const string& data_file_name, const list<string>& to_print)
  {
    ofstream data_file;
    struct stat buffer;
    bool write_header = false;   
    if (stat(data_file_name.c_str(), &buffer) != 0)
      write_header = true;
    data_file.open(data_file_name, std::ofstream::app);
    if (write_header)
    {
      data_file << "# ";
      for (auto quant : to_print)
        data_file << quant << " ";
      data_file << endl;
    }
    for (auto v: _sys.mesh().vertices())
    {
      for (auto p : to_print)
      {
        if (p == "id")
          data_file << v.id << " ";
        if (p == "x")
          data_file << setprecision(10) << v.r.x << " ";
        if (p == "y")
          data_file << setprecision(10) << v.r.y << " ";
        if (p == "boundary")
          data_file << v.boundary << " ";  
        if (p == "type")
          data_file << _sys.get_vert_type_name(v.data().vert_type) << " ";
        if (p == "vx")
          data_file << setprecision(10) << v.data().vel.x << " ";
        if (p == "vy")
          data_file << setprecision(10) << v.data().vel.y << " ";
        if (p == "fx")
          data_file << setprecision(10) << v.data().force.x << " ";
        if (p == "fy")
          data_file << setprecision(10) << v.data().force.y << " ";
        if (p == "constraint")
          data_file << v.data().constraint << " ";
        if (p == "frict_x")
          data_file << v.data().f_type["friction"].x << " ";
        if (p == "frict_y")
          data_file << v.data().f_type["friction"].y << " ";
      }
      data_file << endl;
    }

    data_file.close();
  }

  // JSON conversions
  // Conversation to and from JSON

  // Half-edge
  void to_json(json& j, const HalfEdge<Property>& he) 
  {
    j = json{{"from", he.from()->id}, {"to", he.to()->id}};
  }

  void from_json(const json& j, HalfEdge<Property>& he)
  {
    
  } 

  // Edge
  void to_json(json& j, const Edge<Property>& e) 
  {
    j = json{{"i", e.i}, {"j", e.j}, {"boundary", e.boundary}, {"l0", e.data().l0} };
  }

  void from_json(const json& j, Edge<Property>& e)
  {
    e.i = j.at("i").get<int>();
    e.j = j.at("j").get<int>();
    e.boundary = j.at("boundary").get<bool>();
  }

  // Vertex
  void to_json(json& j, const Vertex<Property>& v) 
  {
    vector<int> neigh;
    HECHandle<Property> he = v.he();
    HECHandle<Property> first = v.he();
    do
    {
      neigh.push_back(he->to()->id);
      he = he->pair()->next();
    } while (he != first);
    string vert_type;
    j = json{
              {"id", v.id}, 
              {"r", {v.r.x,v.r.y}}, 
              {"type", v.data().type_name},
              {"erased", v.erased}, 
              {"boundary", v.boundary}, 
              {"constraint", v.data().constraint},
              {"coordination", v.coordination},
              {"force", {v.data().force.x, v.data().force.y}},
              {"neighbours", neigh}
            };
  }

  void from_json(const json& j, Vertex<Property>& v) 
  {
    v.id   = j.at("id").get<int>();
    v.r.x  = j.at("r").get<vector<double>>()[0];
    v.r.y  = j.at("r").get<vector<double>>()[1];
    v.data().vert_type   = j.at("type").get<int>();
    v.erased = j.at("erased").get<bool>();
    v.boundary = j.at("boundary").get<bool>();
    v.data().constraint = j.at("constraint").get<string>();
    v.coordination = j.at("coordination").get<int>();
    v.data().force.x = j.at("force").get<vector<double>>()[0];
    v.data().force.y = j.at("force").get<vector<double>>()[1];
  }

  // Face
  void to_json(json& j, const Face<Property>& f) 
  {
    vector<int> verts;
    HECHandle<Property> he = f.he();
    HECHandle<Property> first = he;
    vector<double> lx, ly, ndx, ndy, l0,  tension;
    do
    {
      verts.push_back(he->from()->id);
      l0.push_back(he->edge()->data().l0);
      Vec er = he->to()->r - he->from()->r;
      Vec fr = he->pair()->face()->data().rc - f.data().rc;
      lx.push_back(er.x);
      ly.push_back(er.y);
      ndx.push_back(fr.x);
      ndy.push_back(fr.y);
      tension.push_back(he->data().tension);
    } while ((he = he->next()) != first);
    j = json{
        {"id", f.id},
        {"outer", f.outer},
        {"nsides", f.nsides},
        {"type", f.data().type_name},
        {"A0", f.data().A0},
        {"P0", f.data().P0},
        {"vertices", verts},
        {"neighbours", f.data().neighs},
        {"rc", {f.data().rc.x, f.data().rc.y}},
        /*
        {"lx", lx},
        {"ly", ly},
        {"ndx", ndx},
        {"ndy", ndy},
        */
        {"kappa", f.data().kappa},
        {"gamma", f.data().gamma},
        {"lambda", f.data().lambda},
        {"beta", f.data().beta},
        {"beta_a", f.data().beta_a},
        {"alpha", f.data().alpha},
        {"k", f.data().k},
        {"l0", l0},
        {"tension", tension}};
  }

  void from_json(const json& j, Face<Property>& f) 
  {
    f.id = j.at("id").get<int>();
    f.outer = j.at("outer").get<bool>();
    f.nsides = j.at("nsides").get<int>();
    f.data().face_type = j.at("type").get<int>();
    f.data().kappa = j.at("kappa").get<double>();
    f.data().gamma = j.at("gamma").get<double>();
    f.data().lambda = j.at("lambda").get<double>();
    f.data().beta = j.at("beta").get<double>();
  }

  // Box
  void to_json(json& j, const Box& b)
  {
    vector<double> A = {b.h._mxx, b.h._myx};
    vector<double> B = {b.h._mxy, b.h._myy};
    j = json{
          {"lx", b.h._mxx},
          {"ly", b.h._myy},
          {"a", A},
          {"b", B},
          {"periodic", true}
        };
  }

  vector<string> split(const string &s, char delim)
  {
    stringstream ss{s};
    string item;
    vector<string> elems;
    while (std::getline(ss, item, delim)) 
      elems.push_back(move(item)); 
    return elems;
  }

  void export_Dump(py::module& m)
  {
    py::class_<Dump>(m, "Dump")
        .def(py::init<System &, ForceCompute &>())
        .def("dump_cells", &Dump::dump_cells, py::arg("vtk_file"), py::arg("binary_output") = false, py::arg("draw_periodic") = false)
        .def("dump_junctions", &Dump::dump_junctions, py::arg("vtk_file"), py::arg("binary_output") = false)
        .def("dump_box",&Dump::dump_box, py::arg("vtk_file"), py::arg("binary_output") = false)
        .def("dump_mesh", &Dump::dump_mesh, py::arg("mesh_file"), py::arg("copy_params") = false)
        .def("dump_area", &Dump::dump_area)
        .def("dump_json", &Dump::dump_json)
        .def("dump_vertices", &Dump::dump_vertices)
        .def("dump_stress", &Dump::dump_stress, py::arg("stress_file"), py::arg("average") = false, py::arg("cell_type") = "")
        .def("dump_energy", &Dump::dump_energy)
        .def("dump_data", &Dump::dump_data)
        .def("set_sfc", &Dump::set_sfc);
  }
}

