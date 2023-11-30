# chirotope_tree
## Description
This repository gathers the Sagemath sourcecode that was developed for experimenting on counting triangulations of chirotope trees
## List of files
The repository contains:
- parsedb.py
      which is used to parse the db of chirotopes from Aichholzer's webpage http://www.ist.tugraz.at/aichholzer/research/rp/triangulations/ordertypes/
- ./db
      which contains a filtered version of Aichholzer's database for order types of size from 4 to 9, with triangular convex hull. More precisely, the files of name otypes3filter* contain the coordinates of such chirotopes, and the files of name otypes3index* contain the index of such chirotopes in the full database of Aichholzer
- random_gen_tree.py which focuses on generating random chirotope trees
- sweep_polynom_triangulations which implements Alvarez and Seidel's algorithm to compute the triangulation polynomial of a given chirotope
- tree_to_polynom.py which implements the computation of the number of triangulation of a chirotope tree

And finally two test scripts:
- tests_bst.py which generate random chirotope trees and compute their number of triangulations
- example.py which computes the number of triangulations on an example

## Running the examples
In order to run the examples, Sagemath must be installed on your computer. The two examples can be then run with sage using the commands
```sage -python tests_bst.py```
and
```sage -python example.py```

## Description of a chirotope tree

The main function for the example is `compute_triangulation()` from  `tree_to_polynom.py` which takes a chirotope tree as input and returns the number of triangulations of the chirotope tree, as well as the time t1 spent on running Alvarez and Seidel's algorithm on nodes for precomputation, and the time t2 spent on computing the number of triangulation for the full tree.

A *chirotope tree* is implemented as unary-binary rooted tree of the form:
- either `None`, and it is considered as empty
- or `(node,left_tree, right_tree)` where `left_tree`, `right_tree` are chirotope trees, and `node` is of the form `node = [coords, proxies, id]` where:
  - `coords` lists the coordinates of the points
  - `proxies` lists the position in coords of the proxies, in the order `[proxy towards root, proxy towards left child, proxy towards right child]`
  - `id` is just a name for recognizing the node, and is actually not useful for the computation. In `example.py`, it contains the position of the chirotope in the filtered database of chirotopes of length 9 with triangular convex hull. The actual position in the full database of chirotopes of size 9 can be retrieved from the file `otypes3index09.b32`.

Leaves are trees with two empty children, and unary nodes are implemented as binary nodes with an empty right child.

The function `is_well_built()` from `random_gen_tree.py` takes a chirotope tree as inputs and checks that it is well built.
