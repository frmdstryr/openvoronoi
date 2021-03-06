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
#pragma once

#include "offset_sorter.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;


namespace ovd
{


/// Python wrapper for OffsetSorter
class OffsetSorter_py : public OffsetSorter {
public:
    /// create walk
    OffsetSorter_py(HEGraph& gi): OffsetSorter(gi) { }

    /// return list of offsets as python list
    py::list offset_list_py() {
        py::list py_offsets;
        BOOST_FOREACH( MGVertex v, boost::vertices(g) ) { // loop through each loop
            py::list py_loop;
            bool first = true;
            int vdeg =boost::out_degree( v, g );
            /*
            BOOST_FOREACH( Edge e, boost::out_edges( v, g ) ) {
                vdeg++;
            }*/


            BOOST_FOREACH( OffsetVertex lpt, g[v].vertices ) { //loop through each line/arc
                py::list py_lpt;
                double offset_distance = g[v].offset_distance;
                if (first) {
                    first = false;
                    py_lpt.append( lpt.p ); // 0
                    py_lpt.append( -1 ); // 1
                    py_lpt.append( offset_distance ); // 2
                    py_lpt.append( vdeg ); // 3
                } else {
                    py_lpt.append( lpt.p ); // 0, position
                    py_lpt.append( lpt.r ); // 1, radius
                    py_lpt.append( lpt.c ); // 2, center
                    py_lpt.append( lpt.cw ); // 3, cw or ccw
                    py_lpt.append( lpt.f ); // 4, face
                    py_lpt.append( offset_distance ); // 5
                    py_lpt.append( vdeg ); // 6
                }
                py_loop.append( py_lpt );
            }
            py_offsets.append( py_loop );
        }
        return py_offsets;
    }
private:
    OffsetSorter_py(); // don't use.
};


} // end ovd namespace
// end file offset_sorter.hpp
