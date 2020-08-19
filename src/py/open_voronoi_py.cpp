/*
 *  Copyright (c) 2010-2012 Anders Wallin (anders.e.e.wallin "at" gmail.com).
 *
 *  This file is part of Openvoronoi
 *  (see https://github.com/aewallin/openvoronoi).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

//#include <boost/python.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>


//#include "voronoidiagram_py.hpp"
#include "common/point.hpp"
#include "common/numeric.hpp"
#include "graph.hpp"
//#include "medial_axis_walk_py.hpp"
//#include "offset_py.hpp"
//#include "offset_sorter_py.hpp"

#include "utility/vd2svg.hpp"
#include "version.hpp"


#include "voronoidiagram_py.hpp"

#include "offset.hpp"
#include "offset_py.hpp"
#include "offset_sorter_py.hpp"

// filters:
#include "polygon_interior_filter.hpp"
#include "island_filter.hpp"
#include "medial_axis_filter.hpp"
#include "medial_axis_pocket_py.hpp"
#include "medial_axis_walk_py.hpp"


//#include "medial_axis_pocket_py.hpp"

/*
 *  pybind11 wrapping of voronoi diagram and related classes.
 */

//using namespace ovd;
//using namespace pyovd;
/*!
 * \namespace ovd::pyovd
 * \brief Python wrappers for OpenVoronoi
 */

namespace py = pybind11;

namespace ovd {
namespace pyovd {

PYBIND11_MODULE(openvoronoi, m) {
    m.def("version", version);
    m.def("build_type", build_type);
    m.def("vd2svg", vd2svg);

    py::class_<HEGraph>(m, "Graph")
        .def(py::init<>())
    ;

    py::class_<VoronoiDiagram_py>(m, "VoronoiDiagram")
        .def(py::init<double, unsigned int>())
        .def("addVertexSite",  &VoronoiDiagram_py::insert_point_site1) // (point)
        //.def("addVertexSite",  &VoronoiDiagram_py::insert_point_site ) // (point, step)
        .def("addLineSite",  &VoronoiDiagram_py::insert_line_site2 ) // takes two arguments
        .def("addLineSite",  &VoronoiDiagram_py::insert_line_site3 ) // takes three arguments (idx1, idx2, step)
        .def("addArcSite",  &VoronoiDiagram_py::insert_arc_site ) // arc-site (idx1,idx2, center, cw?, step)
        .def("addArcSite",  &VoronoiDiagram_py::insert_arc_site4 ) // arc-site (idx1,idx2, center, cw?, step)
        .def("getGenerators",  &VoronoiDiagram_py::getGenerators)
        .def("getEdgesGenerators",  &VoronoiDiagram_py::getEdgesGenerators)
        .def("getVoronoiVertices",  &VoronoiDiagram_py::getVoronoiVertices)
        .def("getFaceVertices",  &VoronoiDiagram_py::get_face_vertices)
        .def("getFarVoronoiVertices",  &VoronoiDiagram_py::getFarVoronoiVertices)
        .def("getFarRadius",  &VoronoiDiagram_py::get_far_radius)
        .def("getVoronoiEdges",  &VoronoiDiagram_py::getVoronoiEdges)
        .def("getVoronoiEdgesOffset",  &VoronoiDiagram_py::getVoronoiEdgesOffset)
        .def("numPointSites", &VoronoiDiagram_py::num_point_sites)
        .def("numLineSites", &VoronoiDiagram::num_line_sites)
        .def("numArcSites", &VoronoiDiagram_py::num_arc_sites)
        .def("numVertices", &VoronoiDiagram_py::num_vertices)
        .def("numFaces", &VoronoiDiagram_py::num_faces)
        .def("numSplitVertices", &VoronoiDiagram_py::num_split_vertices)
        .def("__str__", &VoronoiDiagram_py::print)
        .def("reset_vertex_count", &VoronoiDiagram_py::reset_vertex_count)
        .def("setEdgePoints", &VoronoiDiagram_py::set_edge_points)
        .def("setEdgeOffset", &VoronoiDiagram_py::set_null_edge_offset)
        .def("debug_on", &VoronoiDiagram_py::debug_on)
        .def("set_silent", &VoronoiDiagram_py::set_silent)
        .def("check", &VoronoiDiagram_py::check)
        .def_static("reset_vertex_count_", &VoronoiDiagram_py::reset_vertex_count)
        .def("getStat", &VoronoiDiagram_py::getStat)
        .def("filterReset", &VoronoiDiagram_py::filter_reset)
        .def("filter_graph", &VoronoiDiagram_py::filter) // "filter" is a built-in function in Python!
        .def("getFaceStats", &VoronoiDiagram_py::getFaceStats)
        .def("getGraph", &VoronoiDiagram_py::get_graph_reference, py::return_value_policy::reference)
    ;

    py::enum_<VertexStatus>(m, "VertexStatus")
        .value("OUT", OUT)
        .value("IN", IN)
        .value("UNDECIDED", UNDECIDED)
        .value("NEW", NEW)
    ;
    py::enum_<VertexType>(m, "VertexType")
        .value("OUTER", OUTER)
        .value("NORMAL", NORMAL)
        .value("POINTSITE", POINTSITE)
        .value("ENDPOINT", ENDPOINT)
        .value("SEPPOINT", SEPPOINT)
        .value("APEX", APEX)
        .value("SPLIT", SPLIT)
    ;
    py::enum_<VoronoiFaceStatus>(m, "VoronoiFaceStatus")
        .value("INCIDENT", INCIDENT)
        .value("NONINCIDENT", NONINCIDENT)
    ;
    py::enum_<EdgeType>(m, "EdgeType")
        .value("LINE", LINE)
        .value("LINELINE", LINELINE)
        .value("PARA_LINELINE", PARA_LINELINE)
        .value("OUTEDGE", OUTEDGE)
        .value("PARABOLA", PARABOLA)
        .value("ELLIPSE", ELLIPSE)
        .value("HYPERBOLA", HYPERBOLA)
        .value("SEPARATOR", SEPARATOR)
        .value("LINESITE", LINESITE)
        .value("ARCSITE", ARCSITE)
        .value("NULLEDGE", NULLEDGE)
    ;

    py::class_<Point>(m, "Point")
        .def(py::init<>())
        .def(py::init<double, double>())
        .def(py::init<Point>())
        .def(double() * py::self)
        .def(py::self *  double())
        .def(py::self -= py::self)
        .def(py::self -  py::self)
        .def(py::self += py::self)
        .def(py::self +  py::self)
        .def("norm", &Point::norm)
        .def("normalize", &Point::normalize)
        .def("dot", &Point::dot)
        .def("cross", &Point::cross)
        .def("is_right", &Point::is_right)
        .def("xy_perp", &Point::xy_perp)
        .def("__str__", &Point::str)
        .def_readwrite("x", &Point::x)
        .def_readwrite("y", &Point::y)
        .def(py::pickle(
            [](Point const& p) {
                return py::make_tuple( p.x, p.y );
            },[](py::tuple t){
                if (t.size() != 2)
                    return Point();
                return Point(t[0].cast<double>(), t[1].cast<double>());
            }))
    ;

    // Offsetting
    py::class_<Offset_py>(m, "Offset")
        .def(py::init<HEGraph&>())
        .def("str", &Offset_py::print )
        .def("offset", &Offset_py::offset )
        .def("offset_loop_list", &Offset_py::offset_loop_list )
    ;

    py::class_<OffsetVertex>(m, "OffsetVertex")
        .def_readwrite("p", &OffsetVertex::p)
        .def_readwrite("r", &OffsetVertex::r)
        .def_readwrite("c", &OffsetVertex::c)
        .def_readwrite("cw", &OffsetVertex::cw)
        .def_readwrite("f", &OffsetVertex::f)
    ;

    py::class_<OffsetLoop>(m, "OffsetLoop")
        .def(py::init<>())
        .def_readwrite("offset_distance", &OffsetLoop::offset_distance)
        .def("__iter__", [](const OffsetLoop &self) {
            return py::make_iterator(self.vertices.begin(), self.vertices.end());
        })
        .def("__len__", [](const OffsetLoop &self) { return self.vertices.size(); })
        .def("__getitem__", [](const OffsetLoop &self, py::slice slice) {
            ssize_t start, stop, step, slicelength;
            if (!slice.compute(self.vertices.size(), &start, &stop, &step, &slicelength))
                throw py::error_already_set();
            int istart = static_cast<int>(start);
            int istop =  static_cast<int>(stop);
            int istep =  static_cast<int>(step);
            return std::make_tuple(istart, istop, istep);
        })
    ;

    py::class_<OffsetSorter_py>(m, "OffsetSorter")
        .def(py::init<HEGraph&>())
        .def("add_loop", &OffsetSorter_py::add_loop )
        .def("sort_loops", &OffsetSorter_py::sort_loops )
        .def("get_loops", &OffsetSorter_py::offset_list_py )
    ;

// Filters
    py::class_<Filter>(m, "Filter_base"); // pure virtual base class!

    py::class_<polygon_interior_filter, Filter>(m, "PolygonInterior")
        .def(py::init<>())
        .def(py::init<bool>())
    ;
    py::class_<island_filter,  Filter>(m, "IslandFilter")
        .def(py::init<>())
    ;

    py::class_<medial_axis_filter, Filter>(m, "MedialAxis")
        .def(py::init<>())
        .def(py::init<double>())
    ;

    py::class_<MedialAxisWalk_py>(m, "MedialAxisWalk")
        .def(py::init<HEGraph&>())
        .def(py::init<HEGraph&, int>())
        .def("walk", &MedialAxisWalk_py::walk_py)
    ;

    py::class_<medial_axis_pocket_py>(m, "MedialAxisPocket")
        .def(py::init<HEGraph&>())
        .def("run", &medial_axis_pocket_py::run)
        //.def("run2", &medial_axis_pocket_py::run2)
        //.def("get_mic_list", &medial_axis_pocket_py::py_get_mic_list)
        .def("get_mic_components", &medial_axis_pocket_py::py_get_mic_components)
        .def("set_width", &medial_axis_pocket_py::set_width)
        .def("debug", &medial_axis_pocket_py::set_debug)
    ;
}

} // pyovd namespace
} // ovd namespace
