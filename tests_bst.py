#!/usr/bin/env sage -python
from sage.all import *
from random_gen_tree import fill_tree_otypes,bstrandomtree
from tree_to_polynom import compute_triangulations, empty_cache
from parsedb import *

size_tree=146
npoints=9 #number of points in nodes, files are provided in ./db for npoints between 4 and 9
totalnbpoints=npoints*size_tree-2*(size_tree-1)

print(size_tree, "nodes,", npoints, "points per node,", totalnbpoints, "points in total")


nb_iterations=100

number_chirotopes_triangle=[0,0,0,1,1,1,6,49,1178,55235,4876476] #number of chirotopes of size n with  triangular convex hull for n<=10
z,l,r=var('z l r')

points=parsedb(npoints, "otypes3filter",2*npoints) # loading of the database of chirotopes of size npoints, with triangular convex hull, each read of the database returns an array of 2*npoints corresponding of the xy coordinates of the npoints elements of the consideres chirotope

for i in range(nb_iterations):
    tree=bstrandomtree(size_tree) #create an empty chirotope tree
    nbchirotopes=number_chirotopes_triangle[npoints]
    fill_tree_otypes(tree,points,nbchirotopes) #fill the tree with random chirotopes in the db points
    empty_cache()
    nb_triangulations,t1,t2=compute_triangulations(tree) # returns the number of triangulation of the chirotope tree tree, and the time t1 spent on Alvarez and Seidel precomputation, and the time t2 spent on the actual computation from the tree
    root1=real_nth_root(nb_triangulations,totalnbpoints)
    print(root1.n(digits=15),t1,t2)
