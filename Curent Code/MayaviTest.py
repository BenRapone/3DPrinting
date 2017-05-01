from numpy import *
import sys
# import os
# os.environ['ETS_TOOLKIT'] = 'qt4'
# orig_stdout = sys.stdout


NumLayers = 100;
Nodes = {};
NodeIndices = {};
NodePath = [];
x = [];
y= [];
z = [];
s = [];

i=0;
# for j in range(1,NumLayers+1):
#   with open('Layer'+str(j)+'.node') as f:
#     NumNodes = int(open('Layer'+str(j)+'.node').readlines()[0].split()[0]);
#     for line in f.readlines()[1:NumNodes+1]:
#       i=i+1;
#       node = [float(line.split()[1]),float(line.split()[2]),float(2*(j-1))];
#       Nodes[i]=node;
#       NodeIndices[node[0],node[1],node[2]]=i;

# NumNodes = len(Nodes);

# pz = 0.0;

with open('VertexPath.txt') as f:
  for line in f.readlines():
    if float(line.split()[2]) == 40:
        if float(line.split()[3]) > 0:
    #       NodePath.append([NodeIndices[float(line.split()[0]),float(line.split()[1]),float(line.split()[2])],7]);
          x.append(float(line.split()[0]));
          y.append(float(line.split()[1]));
          z.append(float(line.split()[2]));
          s.append(7);

        else:
    #       NodePath.append([NodeIndices[float(line.split()[0]),float(line.split()[1]),float(line.split()[2])],1]);
          x.append(float(line.split()[0]));
          y.append(float(line.split()[1]));
          z.append(float(line.split()[2]));
          s.append(1);

    # pz = float(line.split()[2]);

import pylab as pl
from matplotlib import collections  as mc

lines = [];
c = [];

for i in range(0,len(x)-1):
    lines.append([(x[i],y[i]),(x[i+1],y[i+1])])
    if s[i] == 1:
        c.append((1, 0, 0, 1))
    else:
        c.append((0, 1, 0, 1))

# lines = [[(0, 1,2), (1, 1,2)], [(2, 3,2), (3, 3,2)], [(1, 2,2), (1, 3,2)]]
# c = array([(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)])
c = array(c)

lc = mc.LineCollection(lines, colors=c, linewidths=2)
fig, ax = pl.subplots()
ax.add_collection(lc)
ax.autoscale()
ax.margins(0.1)

pl.show()



#####################################
#########  MayaviTest #############
####################################


# from mayavi import mlab
# # View it.
#
# l = mlab.plot3d(x, y, z, s)
# mlab.show()
