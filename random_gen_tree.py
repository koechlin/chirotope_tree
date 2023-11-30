#!/usr/bin/env sage -python
from random import randint, choice, choices,sample
import sys
from sweep_polynom_triangulations import orientation, compute_convex_hull
from parsedb import *


def count(n): #count number of empty trees of size <= n
    T = [0,1]
    debut = 2
    for i in range(debut, n+1):
        T.append( ( 3*(i-2)*T[i-2] + (2*i-1)*T[i-1] )//(1+i) )
    return T

def bstrandomtree(n):
    if n==0 :
        print("Error: no tree of size 0")
        sys.exit(1)
    if n == 1:
        return ([],None,None)
    if n==2:
        return([],([],None,None),None)
    childs = randint(1,4) #draw number between 1 and ...
    if childs==1:
        return ([],bstrandomtree(n-1),None)
    else:
        k=randint(1,n-2) #draw number between 1 and n-2
        return ([], bstrandomtree(k),bstrandomtree(n-1-k))

def urandomtree(n):
    T=count(n)
    return urandomtree_rec(n,T)

def urandomtree_rec(n,T): #generate an empty *uniform* unary-binary tree of size n
    if n==0 or len(T)<=n:
        print("Error: no tree of size 0 or sizes not computed")
        sys.exit(1)
    if n == 1:
        return ([],None,None)
    dice = randint(1,T[n]) #draw number between 1 and T[n]
    if dice <= T[n-1]: #the node is unary
        return ([],urandomtree_rec(n-1,T),None)
    dice -= T[n-1]
    k = 0
    while dice > 0:
        k += 1
        dice -=  T[k]*T[n-1-k]
    return ([], urandomtree_rec(k,T),urandomtree_rec(n-1-k,T))


def print_tree(tree):
    if tree is None:
        return
    (node,ltree,rtree)=tree
    print(node)
    print_tree(ltree)
    print_tree(rtree)
    return

def fill_tree_otypes(tree,chirotopedb,len_chirotope): #fill tree with random uniform chirotopes from chirotopedb
    if tree is None:
        return
    (node,ltree,rtree)=tree
    index=randint(0,len_chirotope-1)
    coords=chirotopedb.getelem(index,1)
    convex_hull=compute_convex_hull(coords)
    # print(convex_hull)
    if ltree is None and rtree is None: #leaf
        proxies=sample(convex_hull,1)
    elif ltree is None or rtree is None: #unary node
        proxies=sample(convex_hull,2)
    else:
        proxies=sample(convex_hull,3)
        if orientation(coords[proxies[0]],coords[proxies[1]],coords[proxies[2]])==-1:
            proxies[1],proxies[2]=proxies[2],proxies[1]
    node.append(coords)
    node.append(proxies)
    node.append(index)
    fill_tree_otypes(ltree,chirotopedb,len_chirotope)
    fill_tree_otypes(rtree,chirotopedb,len_chirotope)

def is_well_built(tree):
    if tree is None:
        return True
    (node,ltree,rtree)=tree
    coords,proxies,index=node
    i_root=proxies[0]
    convex_hull=compute_convex_hull(coords)
    if not i_root in convex_hull:
        print("Proxy 0 (", i_root, ") in node", node, "should be extreme")
        return False
    if ltree is None and rtree is None: #leaf
        if len(proxies)!=1:
            print("The leaf", node, "should have one proxy, but has", len(proxies), "instead");
            return False;
        return True
    elif ltree is not None and rtree is not None : #binary tree
        if len(proxies)!=3:
            print("The binary node ", node, "should have three proxies, but has", len(proxies), "instead");
            return False;
        i_left=proxies[1]
        i_right=proxies[2]
        if orientation(coords[proxies[0]],coords[proxies[1]],coords[proxies[2]])==-1:
            print("The proxies of the node ", node, "are not oriented counterclockwise.");
            return False;

        if set(convex_hull) != set([i_root,i_left,i_right]):
            print("The proxies of the node ", node, "do not form the convex hull of the chirotope");
            return False;
        return is_well_built(ltree) and is_well_built(rtree)
    else: #unary tree, unique node is in left child is not None, but in case we check
        if rtree is not None:
            print("The unary node ", node, " is supposed to have empty right child.");
            return False;
        if len(proxies)!=2:
            print("The unary node ", node, "should have two proxies, but has", len(proxies), "instead");
            return False;
        i_left=proxies[1]
        if not i_left in convex_hull:
            print("Proxy 1 (", i_root, ") in node", node, "should be extreme")
            return False
        return is_well_built(ltree)



if __name__ == '__main__':

    tree=([[[65174, 18626], [55817, 16017], [59425, 22530], [61385, 21294], [63421, 20586], [61721, 22788], [63741, 23489], [64307, 26250], [63839, 30589]], [0, 8, 1], 16857], ([[[4801, 37057], [65372, 63586], [61686, 42059], [55910, 39095], [21406, 36112], [46619, 33927], [34098, 25491], [45327, 19373], [60033, 1950]], [1, 0, 8], 7404], ([[[18213, 19968], [13347, 18118], [13605, 18664], [15562, 19247], [16528, 19932], [14469, 19949], [15086, 20360], [15607, 21018], [14409, 22519]], [1, 0], 28842], ([[[48485, 13508], [46673, 21344], [47074, 21010], [48066, 16565], [48501, 17526], [48733, 17993], [49086, 20115], [49554, 19035], [50091, 21148]], [1, 0, 8], 50641], ([[[61982, 766], [11138, 718], [15438, 12458], [33998, 11918], [27538, 15438], [40478, 28438], [33798, 37358], [34298, 44058], [28438, 64818]], [8, 1, 0], 31166], ([[[65529, 24669], [677, 3688], [4838, 23277], [49944, 27704], [38012, 30649], [17435, 36881], [27828, 35931], [22261, 46590], [8163, 61847]], [8, 1], 29845], ([[[33940, 48786], [32020, 46365], [31917, 46582], [32476, 47466], [32931, 48067], [31727, 47654], [32752, 48449], [32189, 48621], [30805, 48589]], [1, 0, 8], 21498], ([[[58721, 38074], [62967, 38504], [63091, 37321], [62192, 37397], [61454, 37441], [60598, 37584], [62016, 35611], [61853, 35211], [63611, 33579]], [0, 8, 1], 15658], ([[[15627, 14124], [3166, 47754], [37377, 49303], [33963, 42139], [31111, 35209], [42271, 49657], [38593, 44317], [39493, 40883], [55317, 51411]], [8], 38115], None, None), ([[[17991, 22065], [21083, 20842], [20467, 19958], [20086, 20158], [19674, 20522], [19082, 20811], [19189, 19658], [19034, 19049], [19536, 17339]], [8], 11634], None, None)), ([[[5486, 63922], [61594, 56939], [34735, 33594], [30451, 36170], [25005, 42085], [21869, 35372], [6511, 56379], [8590, 35890], [9612, 330]], [8, 1, 0], 3075], ([[[28272, 29586], [27517, 26293], [27310, 26884], [27534, 27979], [27004, 27871], [27493, 29244], [26730, 28927], [26857, 29522], [25026, 30345]], [1], 24229], None, None), ([[[14467, 49774], [16404, 41661], [13917, 40500], [14146, 44470], [13176, 45184], [12449, 42646], [11750, 42980], [11389, 43780], [8072, 37672]], [0], 28968], None, None))), None), ([[[54111, 8150], [10370, 28244], [18170, 36596], [27218, 29444], [31634, 33716], [39338, 35108], [45710, 31016], [48950, 34412], [43370, 63920]], [1, 8], 23558], ([[[65154, 301], [1673, 10178], [5378, 14543], [18308, 12413], [35213, 10223], [32198, 24038], [25478, 39863], [44108, 30878], [23753, 62048]], [8], 41164], None, None), None)), ([[[52933, 19923], [50252, 11423], [49527, 11486], [49406, 12154], [49599, 12700], [49657, 14700], [48728, 16227], [47970, 16110], [41850, 11756]], [8, 0], 2460], ([[[64317, 34064], [41868, 8172], [37305, 14854], [32027, 23109], [38813, 25501], [39671, 30181], [43233, 33704], [38618, 41062], [1815, 57364]], [1], 556], None, None), None)), None), ([[[59126, 1208], [6410, 39490], [16095, 37091], [46567, 12300], [44604, 22236], [46347, 28529], [49329, 51027], [55993, 26765], [51415, 64327]], [1, 0, 8], 45091], ([[[32653, 6121], [30610, 12794], [34545, 13458], [34142, 11581], [35331, 12081], [35473, 11283], [36722, 12011], [36230, 10636], [39328, 14444]], [0], 19086], None, None), ([[[12524, 56174], [53012, 32037], [48171, 29999], [45776, 28613], [25750, 42677], [34571, 31143], [21430, 32153], [13782, 32657], [13822, 9361]], [0, 1], 33727], ([[[44470, 48975], [63520, 18730], [43775, 25012], [36484, 17929], [36119, 19752], [32281, 20751], [28741, 21871], [24561, 22306], [2016, 16561]], [8, 1], 52693], ([[[55361, 795], [11456, 9014], [15464, 18170], [16550, 37472], [31592, 23492], [41510, 17486], [25436, 41018], [28682, 42626], [16160, 65336]], [8], 45828], None, None), None), None))), ([[[21953, 6171], [512, 56816], [24002, 48536], [24722, 36638], [25874, 26648], [31004, 16568], [50480, 35576], [42002, 21176], [65024, 34298]], [1, 0, 8], 11457], ([[[57027, 50268], [60269, 30656], [16937, 19069], [23142, 23982], [17556, 21427], [17400, 21343], [7718, 16486], [15126, 21568], [5266, 15267]], [0, 8, 1], 43988], ([[[11277, 13976], [12877, 17289], [13301, 16965], [12728, 15595], [12860, 15320], [12365, 14862], [15369, 15432], [13058, 14114], [16863, 14373]], [0], 40076], None, None), ([[[20519, 3059], [1555, 44598], [17879, 30571], [20066, 9025], [19920, 11918], [19874, 14185], [19655, 18517], [22359, 3805], [48499, 12957]], [8], 20365], None, None)), ([[[56913, 30204], [13187, 1619], [14693, 6587], [23195, 13043], [28925, 24995], [33515, 26471], [38819, 33089], [28535, 35081], [27491, 63917]], [1, 8], 26323], ([[[6898, 36992], [10504, 39042], [10209, 38156], [8797, 37367], [9837, 37534], [10670, 37542], [10079, 37221], [8447, 36778], [10847, 36398]], [8, 1, 0], 49777], ([[[64218, 43130], [6469, 1471], [16839, 15836], [30201, 23673], [22823, 23129], [29402, 32224], [4820, 38990], [40758, 47745], [1607, 64065]], [8, 0], 45013], ([[[7412, 52319], [63467, 54383], [36659, 34781], [20723, 44327], [20747, 42281], [12527, 46145], [21299, 34853], [21215, 34199], [2069, 8543]], [0, 8], 21253], ([[[51061, 46970], [43251, 49469], [44479, 52621], [45009, 52693], [45122, 53051], [47783, 50539], [46892, 53803], [48149, 55187], [46639, 59652]], [0], 8832], None, None), None), None), ([[[26178, 44349], [55663, 43320], [42908, 39123], [28215, 43643], [29950, 35133], [28082, 39475], [26689, 39331], [25120, 39275], [9872, 21187]], [8, 1, 0], 23129], ([[[64532, 18706], [39904, 13550], [32580, 21346], [37792, 21306], [49028, 20970], [37692, 25822], [42612, 26730], [44360, 28322], [3716, 51986]], [0, 1], 6477], ([[[32289, 8507], [2179, 51527], [15283, 50695], [27724, 27087], [30935, 47263], [38605, 31871], [42349, 27659], [53269, 42492], [63357, 47874]], [8], 24020], None, None), None), ([[[29833, 50423], [25264, 48700], [25870, 50344], [26719, 50435], [26134, 51322], [28840, 50792], [27494, 51344], [28088, 51130], [25913, 52471]], [1, 8], 34541], ([[[44449, 62590], [56160, 2945], [54060, 5547], [23521, 33168], [44271, 62402], [41195, 59762], [22560, 43962], [21980, 46792], [9375, 40367]], [8], 46510], None, None), None))), None)))

    print(is_well_built(tree))
