# chirotope_tree

A coord-tree is a rooted tree of the form:
    — either None, and it is considered as empty
    — or (node,left_tree, right_tree) where
        left_tree, right_tree are coord-trees
        and node is of the form node = [coords, proxies, id] where
            -coords lists the coordinates of the points
            -proxies lists the position in coords of the proxies, in the order [proxy towards root, proxy towards left child, proxy towards right child]
            - id is actually not useful for the computation, it is here the position of the chirotope in the filtered database of chirotopes of length 9 with triangular convex hull. The actual position in the full database of chirotopes of size 9 can be retrieved from the file otypes3index09.b32.
