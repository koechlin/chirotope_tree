#!/usr/bin/env sage -python
from sage.all import *
from sage.geometry.polyhedron.plot import cyclic_sort_vertices_2d
from sweep_polynom_triangulations import coords_to_triangulations_pol,init_variables,orientation2
from time import perf_counter

def sum_lz(k):
    R=PolynomialRing(ZZ,'z,l,r')#['z,l,r']
    z,l,r=R.gens()
    return sum([l**i*z**(k-i) for i in range(k+1)])

@cached_function
def Rab_z(a,b):
    R=PolynomialRing(ZZ,'z,l,r')#['z,l,r']
    z,l,r=R.gens()
    if a==0 or b==0:
        return 0;
    elif b==1:
        return z**(a-1);
    else:
        return Rab_z(a-1,b)+Rab_z(a,b-1);

@cached_function
def Rab_lz(a,b):
    R=PolynomialRing(ZZ,'z,l,r')#['z,l,r']
    z,l,r=R.gens()
    if a==0 or b<=1:
        return 0
    elif a==1:
        return 1;
    elif b==2:
        return sum_lz(a-1)
    else:
        return Rab_lz(a-1,b)+Rab_lz(a,b-1)

def empty_cache(): #for benchmarks
    Rab_z.clear_cache()
    Rab_lz.clear_cache()

def unary_merge(P,Q): #P is merge into Q, P(z) and Q(z,l), returns ~Q(z)
    R=PolynomialRing(ZZ,'z,l,r')
    z,l,r=R.gens()
    d1=P.degree(z)
    d2=Q.degree(l)
    return sum([Q.coefficient(l**b)*sum([P.coefficient(z**a)*Rab_z(a,b) for a in range(2,d1+1)]) for b in range(2,d2+1)])

def binary_merge_right(P,Q): #P *right* child is merge into Q, P(z) and Q(z,l,r), returns unary node ~Q(z,l)
    R=PolynomialRing(ZZ,'z,l,r')
    z,l,r=R.gens()
    d1=P.degree(z)
    d2=Q.degree(r)
    return sum([Q.coefficient(r**b)*sum([P.coefficient(z**a)*Rab_lz(a,b) for a in range(2,d1+1)]) for b in range(2,d2+1)])

def binary_merge_left(P,Q): #P *left* child is merge into Q, P(z) and Q(z,l,r), ~Q(z,r) returns unary node ~Q(z,l)
    R=PolynomialRing(ZZ,'z,l,r')
    z,l,r=R.gens()
    return binary_merge_right(P,Q.substitute({l:r,r:l}))

def poltree_to_pol(tree): #compute the triangulation polynomials of the chirotope tree, given the triangulation polynomials of the nodes
    R=PolynomialRing(ZZ,'z,l,r')#['z,l,r']
    z,l,r=R.gens()
    pnode,ltree,rtree=tree
    pnode=R(pnode)
    if ltree is None and rtree is None:
        return pnode.substitute({l:1,r:1}); #should be a polynom in z normally
    elif ltree is not None and rtree is not None : #binary tree
        l_p=poltree_to_pol(ltree) #polynom in z
        r_p=poltree_to_pol(rtree) #polynom in z
        if l_p.degree(z) <= r_p.degree(z):
            return unary_merge(r_p,binary_merge_left(l_p,pnode))
        else:
            return unary_merge(l_p,binary_merge_right(r_p,pnode))
    else: #unary tree, unique node is in left child is not None, but in case we check
        if rtree is not None:
            l_p=poltree_to_pol(rtree)
        else:
            l_p=poltree_to_pol(ltree)
        return unary_merge(l_p,pnode)

def coordtree_to_poltree(tree): #compute the triangulation polynomials of the nodes of the tree according to the coordinates indicated in the nodes
    z,l,r=var('z l r')
    node,ltree,rtree=tree
    #z=is the root, l left, and r right
    coords,proxies,index=node
    n=len(coords)
    x=init_variables(n)
    i_root=proxies[0]
    if ltree is None and rtree is None: #leaf
        p=coords_to_triangulations_pol(coords, [i_root]).subs({x[i_root]:z})
    elif ltree is not None and rtree is not None : #binary tree
        i_left=proxies[1]
        i_right=proxies[2]
        p=coords_to_triangulations_pol(coords, proxies).subs({x[i_root]:z,x[i_left]:l,x[i_right]:r})
        ltree=coordtree_to_poltree(ltree)
        rtree=coordtree_to_poltree(rtree)
    else: #unary tree, unique node is in left child is not None, but in case we check
        i_left=proxies[1]
        if rtree is not None:
            ltree=coordtree_to_poltree(rtree)
            rtree=None
        else:
            ltree=coordtree_to_poltree(ltree)
        p=coords_to_triangulations_pol(coords,[i_root,i_left]).substitute({x[i_root]:z,x[i_left]:l})
    return (p,ltree,rtree)

def compute_triangulations(coordtree):
    tic=perf_counter()
    poltree=coordtree_to_poltree(coordtree)
    t1=perf_counter()-tic
    tic=perf_counter()
    pol=poltree_to_pol(poltree)
    t2=perf_counter()-tic
    nb_triangulations=pol.substitute(z=1)
    return nb_triangulations,t1,t2
