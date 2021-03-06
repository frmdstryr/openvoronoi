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

#include "medial_axis_pocket.hpp"

namespace ovd {
namespace pyovd {

/// python wrapper for medial_axis_pocket
class medial_axis_pocket_py : public medial_axis_pocket {
public:
    /// create medial_axis_pocket object
    medial_axis_pocket_py(HEGraph& gi): medial_axis_pocket(gi) {}
    /// return list of MICs \todo replaced by get_mic_components()
    py::list py_get_mic_list() {
        py::list out;
        BOOST_FOREACH(MIC mic, mic_list) {
            py::list m;
            m.append( mic.c2 );          //0
            m.append( mic.r2 );          //1
            m.append( mic.t1 );         //2
            m.append( mic.t2 );         //3
            m.append( mic.t3 );         //4
            m.append( mic.t4 );         //5
            m.append( mic.c1 );          //6
            m.append( mic.r1 );          //7
            m.append( mic.new_branch ); //8
            m.append( mic.c_prev );     //9
            m.append( mic.r_prev );     //10
            out.append(m);
        }
        return out;
    }

    /// return list of MICs. one list per connected component in the graph
    py::list py_get_mic_components() {
        py::list out;
        BOOST_FOREACH(MICList miclist, ma_components) {
            py::list mic_out;
            BOOST_FOREACH(MIC mic, miclist) {
                py::list m;
                m.append( mic.c2 );          //0
                m.append( mic.r2 );          //1
                m.append( mic.t1 );         //2
                m.append( mic.t2 );         //3
                m.append( mic.t3 );         //4
                m.append( mic.t4 );         //5
                m.append( mic.c1 );          //6
                m.append( mic.r1 );          //7
                m.append( mic.new_branch ); //8
                m.append( mic.c_prev );     //9
                m.append( mic.r_prev );     //10
                mic_out.append(m);
            }
            out.append(mic_out);
            //return out;

        }
        return out;
    }
};

} // pyovd
} // end ovd namespace
// end file
