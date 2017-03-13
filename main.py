
import pydcel
from pydcel.dcel import vertex, hedge, face, DCEL
from sets import Set

# 0 => horizontal partition
# 1 => grid partition
horizontal_or_grid_part = int(raw_input())

# vertex count
vertex_count = int(raw_input())

# vertex list, in CCW order
vertex_list = []
for i in xrange(vertex_count):
    s = raw_input().split()
    x = int(s[0])
    y = int(s[1])

    vertex_list.append((x, y))

# hole count, each hole described later
hole_count = int(raw_input())

# description for each of the holes
holes = []
for i in xrange(hole_count):
    # hole vertex list, in CW order
    hole_vertex_count = int(raw_input())
    hole_vertex_list = []
    for j in xrange(hole_vertex_count):
        s = raw_input().split()
        x = int(s[0])
        y = int(s[1])

        hole_vertex_list.append((x, y))

    holes.append(hole_vertex_list)

# build DCEL
d = DCEL()

hedges = {} # he_id: (v_origin, v_end), f, nextedge, prevedge

vertex_counter = 0
vertex_start = vertex_counter
vertex_end = vertex_counter + vertex_count
j = 0
for v_idx in xrange(vertex_start, vertex_end):
    if v_idx == vertex_start:
        hedges[j] = (v_idx, v_idx+1), 0, j+1, vertex_end-1
    elif v_idx < vertex_end-1:
        hedges[j] = (v_idx, v_idx+1), 0, j+1, j-1
    else:
        hedges[j] = (v_idx, vertex_start), 0, vertex_start, j-1

    j += 1

vertex_counter += vertex_count

d.createFace() # only one face
infinite_face = d.createInfFace() # outside face

for (x, y) in vertex_list:
    d.createVertex(x, y, 0)

for e in xrange(len(hedges)):
    d.createHedge()

inf_edge = None
for this_edge, value in hedges.iteritems():
    v, face, nextedge, prevedge = value
    v_origin, v_end = v
    
    e = d.hedgeList[this_edge]

    # def setTopology(self, newOrigin, newTwin, newIncindentFace, newNext, newPrevious):
    e_twin = d.createHedge()
    e_twin.setTopology(d.vertexList[v_end], e, infinite_face, None, None) # oops, forgetting to set something here...
    inf_edge = e_twin

    d.faceList[face].setTopology(e)
    
    e.setTopology(d.vertexList[v_origin], e_twin, d.faceList[face], d.hedgeList[nextedge], d.hedgeList[prevedge])
    e.origin.setTopology(e)

# now fix prev/next refs for all edges incident to inf face
infinite_face.innerComponent = inf_edge
current_edge = last_correct_edge = inf_edge

while inf_edge.previous == None:
    current_edge = last_correct_edge
    while current_edge.twin.incidentFace != infinite_face:
        current_edge = current_edge.twin.previous
    current_edge = current_edge.twin

    last_correct_edge.next = current_edge
    current_edge.previous = last_correct_edge
    last_correct_edge = current_edge

for hole in holes:
    hedges = {}
    
    vertex_start = vertex_counter
    vertex_end = vertex_counter + len(hole)
    edge_count = len(hole)
    edge_list = []
    for e in xrange(edge_count):
        edge_list.append(d.createHedge())

    j = 0
    for v_idx in xrange(vertex_start, vertex_end):
        print j, len(edge_list)
        if v_idx == vertex_start:
            hedges[j] = (v_idx, v_idx+1), 0, edge_list[j+1], edge_list[-1]
        elif v_idx < vertex_end-1:
            hedges[j] = (v_idx, v_idx+1), 0, edge_list[j+1], edge_list[j-1]
        else:
            hedges[j] = (v_idx, vertex_start), 0, edge_list[0], edge_list[j-1]

        j += 1
    
    vertex_counter += len(hole)

    for (x, y) in hole:
        d.createVertex(x, y, 0)

    inf_edge = None
    for this_edge, value in hedges.iteritems():
        v, face, nextedge, prevedge = value
        v_origin, v_end = v

        e = edge_list[this_edge]

        # def setTopology(self, newOrigin, newTwin, newIncindentFace, newNext, newPrevious):
        e_twin = d.createHedge()
        e_twin.setTopology(d.vertexList[v_end], e, infinite_face, None, None) # oops, forgetting to set something here...
        inf_edge = e_twin

        d.faceList[face].setTopology(e)
        
        e.setTopology(d.vertexList[v_origin], e_twin, d.faceList[face], nextedge, prevedge)
        e.origin.setTopology(e)

    # now fix prev/next refs for all edges incident to inf face
    infinite_face.innerComponent = inf_edge
    current_edge = last_correct_edge = inf_edge

    while inf_edge.previous == None:
        current_edge = last_correct_edge
        while current_edge.twin.incidentFace != infinite_face:
            print (current_edge, current_edge.twin.previous)
            current_edge = current_edge.twin.previous

        current_edge = current_edge.twin

        last_correct_edge.next = current_edge
        current_edge.previous = last_correct_edge
        last_correct_edge = current_edge


gui = pydcel.dcelVis(d)
gui.mainloop()





