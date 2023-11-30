#!/usr/bin/env sage -python
from sage.all import *
from sage.geometry.polyhedron.plot import cyclic_sort_vertices_2d
import itertools


def orientation(a,b,c):
    return sign(a[0]*b[1]-a[0]*c[1]-a[1]*b[0]+a[1]*c[0]+b[0]*c[1]-b[1]*c[0])

def orientation2(i,j,k,coords):
    return orientation(coords[i],coords[j],coords[k])

def compute_convex_hull(pointset):
    hull=Polyhedron(vertices=pointset)
    cycle=cyclic_sort_vertices_2d(hull.vertices())
    convex_hull=[pointset.index(list(vect.vector())) for vect in cycle]
    return convex_hull

def compute_empty_triangles(points_x,coords): #compute empty triangulations, where points are indexed by increasing x-coordinate, polynomial pre-processing
    n=len(coords)
    empty_triangles_top=set()
    empty_triangles_bot=set()
    for triangle in itertools.combinations(range(n),3):
        p1,p2,p3=triangle
        for i in range(n):
            if orientation2(points_x[p1],points_x[p2],points_x[i],coords)==orientation2(points_x[p2],points_x[p3],points_x[i],coords)==orientation2(points_x[p3],points_x[p1],points_x[i],coords): #if i inside p1p2p3 the trangle is not empty we skip it
                break
        else: #if the loop was never broken the triangle is empty
            if orientation2(points_x[p1],points_x[p3],points_x[p2],coords) > 0: #p2 above p1p3
                empty_triangles_top.add(frozenset(triangle))
            else: #p2 below p1p3
                empty_triangles_bot.add(frozenset(triangle))

    return (empty_triangles_top,empty_triangles_bot)

def decompose_hull(points_x,coords): #returns the upper-most chain U and bottom-most chain B, where points are indexed by increasing x-coordinate
    convex_hull=compute_convex_hull(coords)
    convex_hull_x=[]
    n=len(coords)
    B=[]
    U=[]
    for i in range(n):
        if points_x[i] in convex_hull:
            convex_hull_x.append(i)
            orientation_i = orientation2(points_x[0],points_x[n-1],points_x[i],coords)
            if orientation_i==0:
                B.append(i)
                U.append(i)
            elif orientation_i>0:
                U.append(i)
            else:
                B.append(i)
    return tuple(U),tuple(B),convex_hull_x

def compute_successors(C,m,empty_triangles,convex_hull,x): #e_m = [C[m-1],C[m]]
    l=[]
    for i in range(m,len(C)): #e_i=[pi-1 pi]

        #advances of type 1
        for j in range(C[i-1]+1,C[i]):
            if {C[i-1],j,C[i]} in empty_triangles[0]:
                monom=1
                if j in convex_hull:
                    monom*=x[j]**2
                if C[i-1] in convex_hull:
                    monom*=x[C[i-1]]
                if C[i] in convex_hull:
                    monom*=x[C[i]]

                C2=C[0:i]+(j,)+C[i:]
                l.append((C2,i,monom)) #copy could be avoided at this step and be treated later

        #advances of type 2
        if i-2>=0:
            if {C[i-2],C[i-1],C[i]} in empty_triangles[1]:
                monom=1
                if C[i-2] in convex_hull:
                    monom*=x[C[i-2]]
                if C[i] in convex_hull:
                    monom*=x[C[i]]
                C2=C[0:i-1]+C[i:]
                l.append((C2,i-1,monom))
    return l

def reciproc_list(l):
    n=len(l)
    l2=[0]*n
    for i in range(n):
        l2[l[i]]=i;
    return l2

def init_variables(n_points):
    rang=range(n_points)
    x = list(var('x_%i' % i) for i in rang)
    return x

def convert_indices(pol,points_x,x): #convert indices in increasing x-coordinate to indices of the points given in the original coordinate array
    substitution=dict((x[i],x[points_x[i]]) for i in range(len(x)))
    return pol.subs(substitution)

def coords_to_triangulations_pol(coords,proxies=None):
    #initialisation
    n_points=len(coords)
    points_x = [t[0] for t in sorted(enumerate(coords), key=lambda t: t[1][0])]; # sorting by x-coordinate
    reciproc = reciproc_list(points_x)
    U,B,convex_hull=decompose_hull(points_x,coords) # U chain up, B chain bottom, indexed by x
    empty_triangles=compute_empty_triangles(points_x,coords)
    monom=1
    x=init_variables(n_points)
    if proxies is not None:
        convex_hull=[reciproc[x] for x in proxies]

    #at this point convex_hull only contains the points of the convex hull of interest, using the numerotation of increasing x-coordinate
    for p in B:
        if p in convex_hull:
            monom*=x[p]
            if 0<p<n_points-1: #interior points have arity 2
                monom*=x[p]

    d_level1=dict()
    d_level1[(B,1)]=monom
    U_polynom=0
    while(d_level1):
        d_next=dict()
        for chaine, pol in d_level1.items():
            C,m=chaine
            successors=compute_successors(C,m,empty_triangles,convex_hull,x)
            for successor in successors:
                C2,l,monom_s=successor
                if C2==U:
                    U_polynom=U_polynom+pol*monom_s
                else:
                    if (C2,l) in d_next:
                        d_next[(C2,l)]+=pol*monom_s
                    else:
                        d_next[(C2,l)]=pol*monom_s
        d_level1=d_next
    U_polynom_rename = convert_indices(U_polynom,points_x,x)
    if proxies is not None:
        return U_polynom_rename
    else:
        return U_polynom_rename,[points_x[i] for i in convex_hull]


def compute_pol_triangulations(coords):
    #initialisation
    n_points=len(coords)
    points_x = [t[0] for t in sorted(enumerate(coords), key=lambda t: t[1][0])]; # sorting by x-coordinate
    U,B,convex_hull=decompose_hull(points_x,coords) # U chain up, B chain bottom, indexed by x
    empty_triangles=compute_empty_triangles(points_x,coords)
    monom=1
    x=init_variables(n_points)
    for p in B:
        monom*=x[p]
        if 0<p<n_points-1: #les points intérieurs sont d'arité 2
            monom*=x[p]

    d_level1=dict()
    d_level1[(B,1)]=monom
    U_polynom=0
    while(d_level1):
        d_next=dict()
        for chaine, pol in d_level1.items():
            C,m=chaine
            successors=compute_successors(C,m,empty_triangles,convex_hull,x)
            for successor in successors:
                C2,l,monom_s=successor
                if C2==U:
                    U_polynom=U_polynom+pol*monom_s
                else:
                    if (C2,l) in d_next:
                        d_next[(C2,l)]+=pol*monom_s
                    else:
                        d_next[(C2,l)]=pol*monom_s
        d_level1=d_next
    U_polynom_rename = convert_indices(U_polynom,points_x,x)
    return U_polynom_rename,[points_x[i] for i in convex_hull]

def transform_pol(pol,kept_indices,coords): #remove unwanted variables
    pol2=pol
    n_points=len(coords)
    x=init_variables(n_points)
    subs=dict((x[i],x[i]) if i in kept_indices else (x[i],1) for i in range(n_points))
    pol2=pol2.subs(subs)
    return pol2.expand()
