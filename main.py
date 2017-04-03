
import sys
import dep.pydcel.pydcel as pydcel
from dep.pydcel.pydcel.dcel import vertex, hedge, face, DCEL
from sets import Set
from time import sleep
import argparse

# enums
# hedge directions, -x, +x, -y, +y
# ordered clockwise
YPlus, XPlus, YMinus, XMinus = range(4)

# get next clockwise direction
def getNextQuadrant(dir):
    return (dir + 1) % 4

# get the previous clockwise direction
def getPreviousQuadrant(dir):
    return (dir + 4 - 1) % 4

def getOppositeDirection(dir):
    return (dir + 2) % 4

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

def getHedgesOfVertex(v, nextHe):
    l = []
    while nextHe not in l:
        if nextHe.origin == v:
            l.append(nextHe)
            nextHe = nextHe.previous
        else:
            nextHe = nextHe.twin

    return l

def canHedgeVertexExpandTo(e, dir):
    return len([he for he in getHedgesOfVertex(e.origin, e) if getHedgeDirection(he) == dir]) == 0

def printMessage(s):
    print(s)
    sys.stdout.flush()

parser = argparse.ArgumentParser()
parser.add_argument('--display_init', action='store_true',
                    help='display initial DCEL and exit')
parser.add_argument('--display_interval', action='store', metavar='INTERVAL',
                    type=float, default=0.0,
                    help='set display interval in seconds. Non-zero values shows partition step progress')
parser.add_argument('--skip_horizontal', action='store_true')
parser.add_argument('--compute_visibility', action='store', metavar=('X', 'Y'), nargs=2, type=int,
                    help='compute visibility after partition for vertex at (X, Y)')

args = parser.parse_args()

# 0 => horizontal partition
# 1 => grid partition
horizontal_or_grid_part = int(raw_input())

# vertex count
vertex_count = int(raw_input())

# vertex list, in CCW order
vertex_list = []
for i in xrange(vertex_count):
    s = raw_input().split()
    x = float(s[0])
    y = float(s[1])

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
        x = float(s[0])
        y = float(s[1])

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

gui = pydcel.dcelVis(d)

# show dcel and exit
if args.display_init:
    gui.mainloop()
    sys.exit()

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
def connectHedgeTo(d, hes, het, dir, display):
    hetDir = getHedgeDirection(het)
    
    # check if the direction makes sense
    if ((dir == XMinus and hetDir != YMinus) or
        (dir == XPlus and hetDir != YPlus) or
        (dir == YMinus and hetDir != XPlus) or
        (dir == YPlus and hetDir != XMinus)):
        return None

    if getHedgeDirection(hes) != getOppositeDirection(hetDir):
        if len([xhe for xhe in getHedgesOfVertex(hes.origin, hes) if getHedgeDirection(xhe) == getOppositeDirection(hetDir)]) > 0:
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

    # only connect in the direction of the sweep line if there's no existing vertex on the other side
    # otherwise, we disturb the active hedge list (and even when the vertex is created, this is taken into account)
    needsNewVertex = (dir == XMinus or dir == YMinus)
    
    # get or create second vertex
    v2 = None
    edgeTarget, twinTarget = None, None
    if het.origin.x == destX and het.origin.y == destY:
        v2 = het.origin
        edgeTarget, twinTarget = het.previous, het
        if needsNewVertex:
            return None

    elif het.next.origin.x == destX and het.next.origin.y == destY:
        v2 = het.next.origin
        edgeTarget, twinTarget = het, het.next
        if needsNewVertex:
            return None

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
            e_twin.setTopology(v2, e, het.twin.incidentFace, het.twin.next, het.twin)
            het.origin = v2
            edgeTarget, twinTarget = e, het

        e.origin.setTopology(e)
        e_twin.origin.setTopology(e_twin)

        assert(e.origin != e.next.origin and e.origin != e.previous.origin)

        fixupNearbyHedges(e)
        fixupNearbyHedges(e_twin)
        
        if display and args.display_interval > 0:
            gui.wm_title("Inserting vertex")
            gui.draw_dcel()
            # gui.explain_hedge(e)
            gui.update()
            sleep(args.display_interval)

    # don't create edges inside holes
    if not isInternalEdge(hes, d) and not isInternalEdge(edgeTarget, d):
        return v2

    # now that we have v1 and v2, we connect them
    e = d.createHedge()
    e_twin = d.createHedge()

    e.setTopology(v2, e_twin, hes.incidentFace, hes, edgeTarget)
    e_twin.setTopology(v1, e, hes.incidentFace, twinTarget, hes.previous)

    fixupNearbyHedges(e)
    fixupNearbyHedges(e_twin)

    if display and args.display_interval > 0:
        gui.wm_title("Adding intersection hedge")
        gui.draw_dcel()
        gui.explain_hedge(e)
        gui.update()
        sleep(args.display_interval)

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
    
def processSweepLine(he, activeHedges, dir, stopAtFirst=True, display=True):
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

    # after connecting, need to keep going in the same direction until
    #  we hit a hole/outside
    if len(closestMinus) > 0:
        closestHe = closestMinus[0]
        v = connectHedgeTo(d, he, closestHe, dirMinus, display)
        if not stopAtFirst and v:
            adjHedges = [ahe for ahe in getHedgesOfVertex(v, closestHe) if getHedgeDirection(ahe) == getHedgeDirection(he)]
            if len(adjHedges) > 0:
                processSweepLine(adjHedges[0], activeHedges, dir, stopAtFirst)

    if len(closestPlus) > 0:
        closestHe = closestPlus[0]
        v = connectHedgeTo(d, he, closestHe, dirPlus, display)
        if not stopAtFirst and v:
            adjHedges = [ahe for ahe in getHedgesOfVertex(v, closestHe) if getHedgeDirection(ahe) == getHedgeDirection(he)]
            if len(adjHedges) > 0:
                processSweepLine(adjHedges[0], activeHedges, dir, stopAtFirst)

# gather relevant hedges
hedges = [he for he in d.hedgeList]
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
    processSweepLine(he, activeHedges, YMinus, display=not args.skip_horizontal)

# need to build grid partition
if horizontal_or_grid_part == 1:
    # gather relevant hedges
    hedges = [he for he in d.hedgeList]
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

# last but not least, fixup the dcel's face info, run some sanity checks
hedges = [he for he in d.hedgeList if isInternalEdge(he, d)]
# consume one polygon for face0 (which exists by default after dcel creation)
if len(hedges) > 0:
    he = hedges[0]
    while he in hedges:
        he.incidentFace.setTopology(he)
        hedges.remove(he)
        he = he.next

# now consume the rest of the polygons
while len(hedges) > 0:
    face = d.createFace()
    he = hedges[0]
    while he in hedges:
        he.incidentFace = face
        face.setTopology(he)
        hedges.remove(he)
        he = he.next

# https://stackoverflow.com/questions/328107/how-can-you-determine-a-point-is-between-two-other-points-on-a-line-segment
epsilon = 0.0001
def isBetween(a, b, c):
    crossproduct = (c.y - a.y) * (b.x - a.x) - (c.x - a.x) * (b.y - a.y)
    if abs(crossproduct) > epsilon : return False   # (or != 0 if using integers)

    dotproduct = (c.x - a.x) * (b.x - a.x) + (c.y - a.y)*(b.y - a.y)
    if dotproduct < 0 : return False

    squaredlengthba = (b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y)
    if dotproduct > squaredlengthba: return False

    return True

# https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines-in-python
def calcVertexLineHedgeIntersect(v1, v2, he):
    xdiff = (v1.x - v2.x, he.origin.x - he.next.origin.x)
    ydiff = (v1.y - v2.y, he.origin.y - he.next.origin.y)

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
        return None # no intersection

    # @todo: refactor this
    d = (det([v1.x, v1.y], [v2.x, v2.y]), det([he.origin.x, he.origin.y], [he.next.origin.x, he.next.origin.y]))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div

    from dep.pydcel.pydcel.dcel import vertex
    v = vertex(x, y, 0, 0)
    if isBetween(v1, v2, v) and isBetween(he.origin, he.next.origin, v):
        return (x, y)

    return None

def vertexLineIntersectsHedge(v1, v2, he):
    return calcVertexLineHedgeIntersect(v1, v2, he) is not None

def vertexLineIntersectsAtVertex(v1, v2, he):
    res = calcVertexLineHedgeIntersect(v1, v2, he)
    if res is None:
        return False

    (x, y) = res
    return ((he.origin.x == x and he.origin.y == y) or
            (he.next.origin.x == x and he.next.origin.y == y))

# @todo: rename
def getVertexRayDirection(v1, v2):
    dirs = []
    if v2.x > v1.x:
        dirs.append(XPlus)
    elif v2.x < v1.x:
        dirs.append(XMinus)

    if v2.y > v1.y:
        dirs.append(YPlus)
    elif v2.y < v1.y:
        dirs.append(YMinus)

    # handle unidirectional ray
    if len(dirs) == 1:
        if dirs[0] == XMinus:
            dirs.append(YPlus)
        elif dirs[0] == XPlus:
            dirs.append(YMinus)
        elif dirs[0] == YMinus:
            dirs.append(XMinus)
        elif dirs[0] == YPlus:
            dirs.append(XPlus)

    assert(len(dirs) > 0)
    return dirs

def isVertexVisible(v1, v2, d):
    if v1 == v2:
        return True

    # compute start quadrant
    dirs = []
    if v2.x > v1.x:
        dirs.append(XPlus)
    elif v2.x < v1.x:
        dirs.append(XMinus)
    if v2.y > v1.y:
        dirs.append(YPlus)
    elif v2.y < v1.y:
        dirs.append(YMinus)

    obj = {}
    for he in getHedgesOfVertex(v1, v1.incidentEdge):
        obj[getHedgeDirection(he)] = he

    startHe = None
    if len(dirs) == 1: # ray is moving along an axis
        # which means there are 2 possible directions to move
        dir = dirs[0]
        if dir in obj and isInternalEdge(obj[dir], d):
            startHe = obj[dir]
        elif getNextQuadrant(dir) in obj and isInternalEdge(obj[getNextQuadrant(dir)], d):
            startHe = obj[getNextQuadrant(dir)]
    else:
        # ray not moving along an axis, can only move to corresponding quadrant
        dir = dirs[0]
        if getNextQuadrant(dir) in dirs:
            dir = getNextQuadrant(dir)

        # can't go in the opposite direction
        invalidDir = getPreviousQuadrant(dir)
        while dir not in obj and dir != invalidDir:
            dir = getNextQuadrant(dir)

        if dir in obj and isInternalEdge(obj[dir], d):
            startHe = obj[dir]
        else: # if the quadrant is faced to infinite, then the ray is outside the polygon
            return False

    if startHe is None:
        return False

    dirs = getVertexRayDirection(v1, v2)
    def isEntryHedge(he, dirs, checkNext=True):
        return ((getHedgeDirection(he) in dirs and getHedgeDirection(he.next) in dirs) or
                (checkNext and getHedgeDirection(he.next) in dirs and getHedgeDirection(he.next.next) in dirs))

    def isExitHedge(he, dirs):
        return ((getOppositeDirection(getHedgeDirection(he)) in dirs and getOppositeDirection(getHedgeDirection(he.next)) in dirs) or
                (getOppositeDirection(getHedgeDirection(he.next)) in dirs and getOppositeDirection(getHedgeDirection(he.next.next)) in dirs))

    he = startHe

    while True:
        # at this point, he is an 'entry' hedge to a polygon
        # we need to loop through to find the 'exit' hedge
        nextHe = he.next
        while nextHe != he:
            if vertexLineIntersectsHedge(v1, v2, nextHe) and isExitHedge(nextHe, dirs):
                break

            nextHe = nextHe.next

        assert(nextHe != he) # there *must* be a different exit hedge

        if nextHe.origin == v2 or nextHe.next.origin == v2:
            return True

        he = nextHe.twin
        if not isInternalEdge(he, d):
            return False

    assert(False)
    return False

def isFaceVisible(v1, f, d):
    return all([isVertexVisible(v1, v2, d) for v2 in f.loopOuterVertices()])

# compute visibility if required
if args.compute_visibility is not None:
    x = args.compute_visibility[0]
    y = args.compute_visibility[1]

    v = [vertex for vertex in d.vertexList if vertex.x == x and vertex.y == y]
    if len(v) == 0:
        printMessage("Can't compute visibility - vertex (%d, %d) does not exist" % (x, y))
        exit(1)

    assert(len(v) == 1) # only one vertex at (x, y) can exist
    v = v[0]

    def checkHedge(he):
        return he.incidentFace in d.faceList or he.twin.incidentFace in d.faceList
    def checkVertex(ve):
        return len([he for he in getHedgesOfVertex(ve, ve.incidentEdge) if checkHedge(he)]) > 0

    d.faceList = [face for face in d.faceList if isFaceVisible(v, face, d)]
    d.vertexList = [ve for ve in d.vertexList if checkVertex(ve)]
    d.hedgeList = [he for he in d.hedgeList if checkHedge(he)]

gui.wm_title("Final DCEL")
gui.draw_dcel()
gui.mainloop()





