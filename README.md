# chirotope_tree

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



A coord-tree is a rooted tree of the form:
    — either None, and it is considered as empty
    — or (node,left_tree, right_tree) where
        left_tree, right_tree are coord-trees
        and node is of the form node = [coords, proxies, id] where
            -coords lists the coordinates of the points
            -proxies lists the position in coords of the proxies, in the order [proxy towards root, proxy towards left child, proxy towards right child]
            - id is actually not useful for the computation, it is here the position of the chirotope in the filtered database of chirotopes of length 9 with triangular convex hull. The actual position in the full database of chirotopes of size 9 can be retrieved from the file otypes3index09.b32.
