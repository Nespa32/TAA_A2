
import sys
import pydcel
from pydcel.dcel import vertex, hedge, face, DCEL
from sets import Set
from time import sleep

# enums
# hedge directions, -x, +x, -y, +y
XMinus, XPlus, YMinus, YPlus = range(4)

def getOppositeDirection(dir):
    if dir == XMinus:
        return XPlus
    elif dir == XPlus:
        return XMinus
    elif dir == YMinus:
        return YPlus

    assert(dir == YPlus)
    return YMinus

def isHorizontalEdge(e):
    return e.origin.x == e.next.origin.x

def isInternalEdge(e, d):
    return e.incidentFace != d.infiniteFace
    
def getHedgeDirection(e):
    v1 = e.origin
    v2 = e.next.origin
    if v2.x > v1.x:
        return XPlus
    elif v2.x < v1.x:
        return XMinus
    elif v2.y > v1.y:
        return YPlus

    assert(v2.y < v1.y)
    return YMinus

def canHedgeVertexExpandTo(e, dir):
    return dir != getHedgeDirection(e) and getOppositeDirection(dir) != getHedgeDirection(e.previous)

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
            current_edge = current_edge.twin.previous

        current_edge = current_edge.twin

        last_correct_edge.next = current_edge
        current_edge.previous = last_correct_edge
        last_correct_edge = current_edge

# at this point, the dcel is built
# now let's compute the horizontal or grid partition
# horizontal partition must be computed in both cases

def fixupNearbyHedges(e):
    # fixup previous/next edge
    e.next.setTopology(e.next.origin, e.next.twin, e.next.incidentFace, e.next.next, e)
    e.previous.setTopology(e.previous.origin, e.previous.twin, e.previous.incidentFace, e, e.previous.previous)

# d = dcel, hes = source hedge, het = target hedge
# connects source hedge origin to target hedge, adjusts dcel
# returns dest vertex
def connectHedgeTo(d, hes, het, dir):
    hetDir = getHedgeDirection(het)
    
    # check if the direction makes sense
    if ((dir == XMinus and hetDir != YMinus) or
        (dir == XPlus and hetDir != YPlus) or
        (dir == YMinus and hetDir != XPlus) or
        (dir == YPlus and hetDir != XMinus)):
        return None

    # first vertex
    v1 = hes.origin
    # compute second vertex location
    destX, destY = hes.origin.x, hes.origin.y
    if dir == XMinus or dir == XPlus:
        destX = het.origin.x
    else:
        destY = het.origin.y

    # can't form a zero-len hedge
    assert(hes.origin.x != destX or hes.origin.y != destY)

    # get or create second vertex
    v2 = None
    edgeTarget, twinTarget = None, None
    if het.origin.x == destX and het.origin.y == destY:
        v2 = het.origin
        edgeTarget, twinTarget = het.previous, het
    elif het.next.origin.x == destX and het.next.origin.y == destY:
        v2 = het.next.origin
        edgeTarget, twinTarget = het, het.next
    else:
        # need to create a vertex
        # we also need to split the edge into 2 edges, because of this new vertex
        v2 = d.createVertex(destX, destY, 0)

        # 2 new edges
        e = d.createHedge()
        e_twin = d.createHedge()
        
        # need to take into account that sweep lines progress through YMinus/XMinus
        # meaning that the XMinus/YMinus (or 'lower') part of any edge needs to remain intact,
        #  and the XPlus/YPlus (or 'upper') part is the one to mutate
        mutateUpperPart = (hetDir == XPlus or hetDir == YPlus)
        if mutateUpperPart:
            e.setTopology(v2, e_twin, het.incidentFace, het.next, het)
            e_twin.setTopology(het.next.origin, e, het.twin.incidentFace, het.twin, het.twin.previous)
            e_twin.next.origin = v2
            edgeTarget, twinTarget = het, e
        else:
            e.setTopology(het.origin, e_twin, het.incidentFace, het, het.previous)
            e_twin.setTopology(v2, e, het.incidentFace, het.twin.next, het.twin)
            het.origin = v2
            edgeTarget, twinTarget = e, het

        e.origin.setTopology(e)
        e_twin.origin.setTopology(e_twin)

        assert(e.origin != e.next.origin and e.origin != e.previous.origin)

        fixupNearbyHedges(e)
        fixupNearbyHedges(e_twin)

    # now that we have v1 and v2, we connect them
    e = d.createHedge()
    e_twin = d.createHedge()

    e.setTopology(v2, e_twin, hes.incidentFace, hes, edgeTarget)
    e_twin.setTopology(v1, e, hes.incidentFace, twinTarget, hes.previous)

    fixupNearbyHedges(e)
    fixupNearbyHedges(e_twin)
    return v2

# returns he2 dist from he1 for the direction axis (i.e can be negative if he2 is 'behind' he1 for dir)
def getHedgeDirDist(he1, he2, dir):
    if dir == XMinus:
        return he1.origin.x - he2.origin.x
    elif dir == XPlus:
        return -(he1.origin.x - he2.origin.x)
    elif dir == YMinus:
        return he1.origin.y - he2.origin.y
    else:
        assert(dir == YPlus)
        return -(he1.origin.y - he2.origin.y)
    
def processSweepLine(he, activeHedges, dir, stopAtFirst=True):
    assert(dir == XMinus or dir == YMinus)

    if dir == XMinus:
        dirMinus, dirPlus = YMinus, YPlus
        hitMinus, hitPlus = XPlus, XMinus
    else:
        dirMinus, dirPlus = XMinus, XPlus
        hitMinus, hitPlus = YMinus, YPlus
    
    # filters hedges and returns only the ones that are in the requested direction
    def getHedgesInDir(he, hedges, dir, hitDir):
        if not canHedgeVertexExpandTo(he, dir): return []
        l = [ahe for ahe in hedges if getHedgeDirDist(he, ahe, dir) > 0]
        l = [ahe for ahe in l if getHedgeDirection(ahe) == hitDir]
        return sorted(l, key = lambda ahe: getHedgeDirDist(he, ahe, dir))

    closestMinus = getHedgesInDir(he, activeHedges, dirMinus, hitMinus)
    closestPlus = getHedgesInDir(he, activeHedges, dirPlus, hitPlus)

    def getHedgesOfVertex(v, nextHe):
        l = []
        while nextHe not in l:
            l.append(nextHe)
            if nextHe.origin == v:
                nextHe = nextHe.previous
            else:
                nextHe = nextHe.twin

        return l


    # after connecting, need to keep going in the same direction until
    #  we hit a hole/outside
    if len(closestMinus) > 0:
        closestHe = closestMinus[0]
        v = connectHedgeTo(d, he, closestHe, dirMinus)
        if not stopAtFirst and v:
            adjHedges = [ahe for ahe in getHedgesOfVertex(v, closestHe) if ahe.origin == v and getHedgeDirection(ahe) == getHedgeDirection(he)]
            if len(adjHedges) > 0:
                processSweepLine(adjHedges[0], activeHedges, dir, stopAtFirst)

    if len(closestPlus) > 0:
        closestHe = closestPlus[0]
        v = connectHedgeTo(d, he, closestHe, dirPlus)
        if not stopAtFirst and v:
            adjHedges = [ahe for ahe in getHedgesOfVertex(v, closestHe) if ahe.origin == v and getHedgeDirection(ahe) == getHedgeDirection(he)]
            if len(adjHedges) > 0:
                processSweepLine(adjHedges[0], activeHedges, dir, stopAtFirst)

# gather relevant hedges
hedges = [he for he in d.hedgeList if isInternalEdge(he, d)] # only internal edges
# sort for sweep
hedges = sorted(hedges, key = lambda he: (he.origin.y, he.origin.x), reverse=True)

activeHedges = [] # edges active for sweep line
# for each internal, we have the following possible events
while len(hedges) > 0:
    he = hedges.pop(0)

    # start hedge, handle if:
    # hedge moving in sweep direction
    # previous hedge moving opposite to sweep direction
    dir = getHedgeDirection(he)
    prevDir = getHedgeDirection(he.previous)
    assert(dir != getOppositeDirection(prevDir))
    if dir == YMinus:
        activeHedges.append(he)
    elif prevDir == YPlus:
        activeHedges.append(he.previous)

    # close hedge, handle if:
    # hedge moving opposite to sweep direction
    if dir == YPlus:
        activeHedges.remove(he)
    elif prevDir == YMinus:
        activeHedges.remove(he.previous)

    # line extend to closest edges
    # need sorted active hedges, but best would be a tree-like structure
    processSweepLine(he, activeHedges, YMinus)

# need to build grid partition
if horizontal_or_grid_part == 1:
    # gather relevant hedges
    hedges = [he for he in d.hedgeList if isInternalEdge(he, d)] # only internal edges
    # sort for sweep
    hedges = sorted(hedges, key = lambda he: (he.origin.x, he.origin.y), reverse=True)

    activeHedges = [] # edges active for sweep line
    # for each internal, we have the following possible events
    while len(hedges) > 0:
        he = hedges.pop(0)
        # start hedge, handle if:
        # hedge moving in sweep direction
        # previous hedge moving opposite to sweep direction
        dir = getHedgeDirection(he)
        prevDir = getHedgeDirection(he.previous)
        assert(dir != getOppositeDirection(prevDir))
        if dir == XMinus:
            activeHedges.append(he)
        elif prevDir == XPlus:
            activeHedges.append(he.previous)

        # close hedge, handle if:
        # hedge moving opposite to sweep direction
        if dir == XPlus:
            activeHedges.remove(he)
        elif prevDir == XMinus:
            activeHedges.remove(he.previous)

        # line extend to closest edges
        # need sorted active hedges, but best would be a tree-like structure
        processSweepLine(he, activeHedges, XMinus, stopAtFirst=False)

sys.stdout.flush()

gui = pydcel.dcelVis(d)
gui.mainloop()





