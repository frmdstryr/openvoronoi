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

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "offset.hpp"

namespace ovd {
namespace pyovd {

/// \brief python wrapper for Offset
class Offset_py : public Offset {
public:
    /// offset of graph \a gi
    Offset_py(HEGraph& gi): Offset(gi) { }

    /// return list of offsets at given offset distance \a t
    py::list offset_py(double t) {
        offset(t);
        py::list py_offsets;
        BOOST_FOREACH( OffsetLoop loop, offset_list ) { // loop through each loop
            py::list py_loop;
            bool first = true;
            BOOST_FOREACH( OffsetVertex lpt, loop.vertices ) { //loop through each line/arc
                py::list py_lpt;
                double offset_distance = loop.offset_distance;
                if (first) {
                    first = false;
                    py_lpt.append( lpt.p );
                    py_lpt.append( -1 );
                    py_lpt.append( offset_distance ); // 2
                } else {
                    py_lpt.append( lpt.p ); // 0, position
                    py_lpt.append( lpt.r ); // 1, radius
                    py_lpt.append( lpt.c ); // 2, center
                    py_lpt.append( lpt.cw ); // 3, cw or ccw
                    py_lpt.append( lpt.f ); // 4, face
                    py_lpt.append( offset_distance ); // 5
                }
                py_loop.append( py_lpt );
            }
            py_offsets.append( py_loop );
        }
        return py_offsets;
    }
    /// return a python-list of OffsetLoop objects
    py::list offset_loop_list(double t) {
        offset(t);
        py::list out;
        BOOST_FOREACH( OffsetLoop loop, offset_list ) {
            out.append( loop );
        }
        return out;
    }

private:
    Offset_py(); // don't use.
};

} // pyovd
} // end ovd namespace
// end offset_py.hpp
