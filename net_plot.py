# Create circular plot of network based on sensitivity.
# Draw edges to show set membership and damage.

import argparse
from graphics import *
from math import *
from csv import DictReader
from collections import defaultdict
import edge_sorter

# calculate position of node centers in circular ordering
def calc_centers():
    dth = 2*pi/N
    th = 3*pi/2
    centers = dict()
    disol = 4*n_radius
    isolx = margin
    isoly = int(margin+2*g_radius+base)
    
    for q_n,n in sorted_nodes:
        if n not in isol:
            cent = Point(cos(th)*g_radius+g_center,sin(th)*g_radius+g_center)
        else:
            cent = Point(isolx,isoly)
            isolx += disol
            
        centers[n] = cent
        th += dth

    return centers

# draw the nodes
def draw_nodes():
    for q_n,n in sorted_nodes:
        cir = Circle(centers[n],get_node_radius(n))
        #cir.setFill(get_node_fill(q_n))
        cir.setFill(get_node_fill(n))
        cir.setOutline(get_node_outline(n))
        cir.setWidth(2)
        cir.draw(win)

# determine node fill color based on set membership
def get_node_fill(n):
    if nodes[n]['GOUTSTAR'] == 'true' and args.star != 0:
        return 'red'
    elif nodes[n]['DAMAGED'] == 'true':
        return dcolor
    elif nodes[n]['GOUT'] == 'true' and args.star == 0:
        return 'black'
    else:
        return 'white'

# determine node outline color
def get_node_outline(n):
##    if nodes[n]['GOUTSTAR'] == 'true' and args.star != 0:
##        return 'red'
##    elif nodes[n]['DAMAGED'] == 'true':
##        return dcolor
##    elif nodes[n]['GOUT'] == 'true' and args.star == 0:
##        return 'black'
##    else:
##        return 'black'
    return 'black'

# determine node radius
def get_node_radius(n):
    if args.star == 0:
        if nodes[n]['GOUT'] == 'true':
            return floor(n_radius*2)
        else:
            return n_radius

    else:
        if nodes[n]['GOUTSTAR'] == 'true':
            return floor(n_radius*2)
        else:
            return n_radius

# draw isolated nodes
def draw_isol():
    for n in isol:
        cir = Circle(centers[n],get_node_radius(n)+2)
        cir.setOutline('black')
        cir.setWidth(2)
        cir.draw(win)

# draw a single edge from s to d
def draw_edge(s,d,l):
    p_s = centers[s]
    p_d = centers[d]

    lin = Line(p_s,p_d)
    lin.setWidth(get_edge_width(l))
    lin.setOutline(get_edge_color(l))

    lin.draw(win)
    
# draw all edges
def draw_edges():
    #print(edges)

    e = order_edges()

    for s,d,l in e:
        draw_edge(s,d,l)


# order the edges so that "important" edges are drawn last (on top)
def order_edges():
    e_dict = defaultdict(list)

    for s,d,l in edges:
        if l!= 's':
            e_color = get_edge_color(l)
            if e_color == dcolor and l == 'gd':
                e_color = 'gorange'
                
            e_dict[e_color].append((s,d,l))

    e = []
    if args.star == 0:
        e.extend(e_dict[gray])
        e.extend(e_dict[dcolor])
        e.extend(e_dict['black'])
        e.extend(e_dict['gorange'])

    else:
        e.extend(e_dict[gray])
        e.extend(e_dict[dcolor])
        e.extend(e_dict['gorange'])
        e.extend(e_dict['red'])

    return e


# determine edge thickness based on nodes
def get_edge_width(l):
    thick = 3
    thin = 1
    
    if args.star == 0:
        if 'go' in l:
            return thick
        else:
            return thin

    else:
        if 'gs' in l:
            return thick
        else:
            return thin

# determine edge color based on damage/membership of connecting nodes
def get_edge_color(l):
    if 'gs' in l:
        if args.star == 0:
            return dcolor
        else:
            return 'red'
    elif 'd' in l:
        return dcolor
    elif 'go' in l and args.star == 0:
        return 'black'
    else:
        return gray

##    if args.star == 0:
##        if 'd' in l or 'gs' in l:
##            return 'red'
##        else:
##            if 'g' in l:
##                return 'black'
##            else:
##                return color_rgb(gray,gray,gray)
##
##    else:
##        if 'gs' in l:
##            return 'red'
##        elif 'd' in l:
##            return 'orange'
##        else:
##            return color_rgb(gray,gray,gray)



            
# read file
parser = argparse.ArgumentParser(description='Graph drawing options')

# match file numbering
parser.add_argument('--mode', type=int, default=0, help="Edge addition rule")
parser.add_argument('--m', type=int, default=2, help="Edge competition subset")
parser.add_argument('--tag', type=int, default=1, help="Data set")

# 0 for GOUT-focused style, 1 for GOUT*-focused style
parser.add_argument('--star', type=int, default=0, help="GOUT* format")

args = parser.parse_args()

# base folder location of file
folder = "/Users/bralex1/eclipse-workspace/EvoNets/"
file_root = "data_dir_" + str(args.mode) + "_m" + str(args.m) + "_l0_dmg0_tag"
file_root += str(args.tag)

margin = 20 # spacing from edge of figure
g_radius = 300 # total graph radius
n_radius = 4 # individual node radius
base = 40 # node spacing

gray_rgb = 178
gray = color_rgb(gray_rgb,gray_rgb,gray_rgb)

dcolor = color_rgb(0, 221, 255)

# total image size
w = 2*margin+2*g_radius
h = w + base

win = GraphWin('',w,h)

edges = []

# read file and parse node and edge data
with open(folder + file_root + "_struct.csv", 'r') as f:
    next(f)
    for line in f:
        line_list = line.split(',')
        #print(line_list)
 #       edges.append((int(line_list[0]), int(line_list[1]), line_list[2].strip()))
        edges.append((int(line_list[0]), int(line_list[1].strip())))
        
node_data = list(DictReader(open(folder + file_root + "_node.csv", 'r')))


node_data,edges = edge_sorter.sort_edges(node_data,edges)


N = len(node_data)
all_n = list(range(N))
nodes = defaultdict(dict)
sorted_nodes = []


for n in all_n:
    nodes[n] = node_data[n]
    sorted_nodes.append((float(node_data[n]['Q']),n))

sorted_nodes.sort()

isol = all_n.copy()

for s,d,l in edges:
    if s in isol:
        isol.remove(s)
    if d in isol:
        isol.remove(d)


#print(edges)
#print(node_data)

##for e in edges:
##    if get_edge_color(e[2]) == dcolor:
##        print(e)
##
##for n in range(N):
##    if get_node_fill(n) == dcolor:
##        print(n)


# draw graph
g_center = margin+g_radius
centers = calc_centers()
draw_edges()
#draw_isol()
draw_nodes()

# save file
win.postscript(file="figs/"+file_root+".ps",colormode="color")

#win.getMouse()
win.close()
