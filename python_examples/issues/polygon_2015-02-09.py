import openvoronoi as ovd
import ovdvtk
import time
import vtk
import datetime
import math
import random
import os
import sys
import pickle
import gzip
import ovdgenerators as gens


def drawLine(myscreen, pt1, pt2, lineColor):
    myscreen.addActor(ovdvtk.Line(p1=(pt1.x, pt1.y, 0), p2=(pt2.x, pt2.y, 0), color=lineColor))


def drawArc(myscreen, pt1, pt2, r, arcColor):
    myscreen.addActor(ovdvtk.Line(p1=(pt1.x, pt1.y, 0), p2=(pt2.x, pt2.y, 0), color=arcColor))


def drawOffsets(myscreen, ofs):
    # draw loops
    nloop = 0
    lineColor = ovdvtk.green
    arcColor = ovdvtk.grass
    for lop in ofs:
        n = 0
        N = len(lop)
        first_point = []
        previous = []
        for p in lop:
            # p[0] is the Point
            # p[1] is -1 for lines, and r for arcs
            if n == 0:  # don't draw anything on the first iteration
                previous = p[0]
                # first_point = p[0]
            else:
                r = p[1]
                p = p[0]
                if r == -1:
                    drawLine(myscreen, previous, p, lineColor)
                else:
                    drawArc(myscreen, previous, p, r, arcColor)
                # myscreen.addActor( ovdvtk.Line(p1=(previous.x,previous.y,0),p2=(p.x,p.y,0),color=loopColor) )
                previous = p
            n = n + 1
        print("rendered loop %s with %s points" % (nloop, len(lop)))
        nloop = nloop + 1


poly_points = [(-0.2567719874411157, -0.4983049800651602),
               (0.12205285479992212, -0.640371712930281),
               (-0.25972854724944455, -0.5143879072702902),
               (-0.34168692840153536, -0.6418861147966213),
               (-0.5288215108461576, 0.18480346369654843),
               (-0.35263585687204546, -0.50735692278175),
               (-0.4821854389417177, 0.46463421861462373)]

if __name__ == "__main__":
    # w=2500
    # h=1500

    # w=1920
    # h=1080
    w = 1024
    h = 1024
    myscreen = ovdvtk.VTKScreen(width=w, height=h)
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version())

    scale = 1
    myscreen.render()
    random.seed(42)
    far = 1
    camPos = far
    zmult = 3
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos / float(1000), zmult * camPos)
    myscreen.camera.SetClippingRange(-(zmult + 1) * camPos, (zmult + 1) * camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)

    vd = ovd.VoronoiDiagram(far, 120)
    print(ovd.version())

    # for vtk visualization
    vod = ovdvtk.VD(myscreen, vd, float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()

    vod.textScale = 0.02
    vod.vertexRadius = 0.0031
    vod.drawVertices = 0
    vod.drawVertexIndex = 1
    vod.drawGenerators = 1
    vod.offsetEdges = 0
    vd.setEdgeOffset(0.05)

    """
    p1=ovd.Point(-0.1,-0.2)
    p2=ovd.Point(0.2,0.1)
    p3=ovd.Point(0.4,0.2)
    p4=ovd.Point(0.6,0.6)
    p5=ovd.Point(-0.6,0.3)

    pts = [p1,p2,p3,p4,p5]
    """
    pts = []
    for p in poly_points:
        pts.append(ovd.Point(p[0], p[1]))

    # t_after = time.time()
    # print ".done in {0:.3f} s.".format( t_after-t_before )
    times = []
    id_list = []
    m = 0
    t_before = time.time()
    for p in pts:
        pt_id = vd.addVertexSite(p)
        id_list.append(pt_id)
        print("%s added vertex %s at %s" % (m, pt_id, p))
        m = m + 1

    t_after = time.time()
    times.append(t_after - t_before)
    # exit()

    # print "   ",2*Nmax," point-sites sites took {0:.3f}".format(times[0])," seconds, {0:.2f}".format( 1e6*float( times[0] )/(float(2*Nmax)*float(math.log10(2*Nmax))) ) ,"us/n*log(n)"
    print("all point sites inserted. ")
    print("VD check: %s" % vd.check())

    print("now adding line-segments.")
    t_before = time.time()
    for n in [0]:  # range(len(id_list)):
        if n == len(id_list) - 1:
            vd.addLineSite(id_list[n], id_list[n + 1])
            print("%s added segment %s to %s" % (n, n, n + 1))
        else:
            vd.addLineSite(id_list[n], id_list[0])
            print("%s added final segment %s to %s" % (n, n, 0))

    # vd.addLineSite( id_list[1], id_list[2])
    # vd.addLineSite( id_list[2], id_list[3])
    # vd.addLineSite( id_list[3], id_list[4])
    # vd.addLineSite( id_list[4], id_list[0])
    vd.check()

    t_after = time.time()
    line_time = t_after - t_before
    if line_time < 1e-3:
        line_time = 1
    times.append(line_time)

    # of = ovd.Offset( vd.getGraph() ) # pass the created graph to the Offset class
    # of.str()
    # ofs = of.offset(0.123)
    # print ofs
    # drawOffsets(myscreen, ofs)

    pi = ovd.PolygonInterior(True)
    vd.filter_graph(pi)

    of = ovd.Offset(vd.getGraph())  # pass the created graph to the Offset class

    ofs = of.offset(0.123)
    # print ofs
    ovdvtk.drawOffsets(myscreen, ofs)

    # of.offset(0.125)

    vod.setVDText2(times)
    vod.setAll()
    print("PYTHON All DONE.")
    myscreen.render()
    myscreen.iren.Start()
