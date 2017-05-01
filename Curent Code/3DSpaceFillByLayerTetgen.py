# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 12:28:49 2015

@author: bastinux
"""
###############################################################################################################
###               Notes: Program requires 'Layer'+str(i)+'.  Node, edge and ele (face) files;
###############################################################################################################

import time
from numpy import *
import sys
import scipy.sparse as sp
from scipy.sparse import lil_matrix
from scipy import linalg
import scipy.spatial.distance
import cplex
from cplex.exceptions import CplexError
import sys
orig_stdout = sys.stdout
tti=time.clock();
import triangle
# import triangle.plot as plot
# import matplotlib.pyplot as plt
from Boundary_Info import InteriorEdgeList_Grab

#########################################################################################################
########         Number of Layers must be Manually assigned as is the gap between layers        #########
#########################################################################################################

NumLayers = 3;
LayerGap = 2.0;


###################################################################################
###    Initialize/Define Layer information and vars: Dictionaries, lists etc.  ####
###################################################################################

Collection = {};

for i in range(1,NumLayers+1):

#### Original Mesh Elements
  Collection["NumNodesLayer{0}".format(i)] = int(open('Layer'+str(i)+'.node').readlines()[0].split()[0]);   #### Extract number of nodes from each layer
  Collection["NumEdgesLayer{0}".format(i)] = int(open('Layer'+str(i)+'.edge').readlines()[0].split()[0]);   #### Extract number of edges from each layer
  Collection["NumFacesLayer{0}".format(i)] = int(open('Layer'+str(i)+'.ele').readlines()[0].split()[0]);    #### Extract number of faces/triangles from each layer
  Collection["NodesLayer{0}".format(i)] = {};             #### Initialize Dictionary for the Nodes of each layer
  Collection["NodeListLayer{0}".format(i)] = [];          #### Initialize List of Nodes for each layer
  Collection["NodeIndicesLayer{0}".format(i)] = {};       #### Initialize Dictionary of Node Indices for each layer
  Collection["EdgesLayer{0}".format(i)] = {};             #### Initialize Dictionary of Edges for each layer
  Collection["EdgeIndicesLayer{0}".format(i)] = {};       #### Initialize Dictionary of Edge Indices for each layer
  Collection["EdgeListLayer{0}".format(i)] = [];          #### Initialize List of Edges for each layer
  Collection["FacesLayer{0}".format(i)] = {};             #### Initialize Dictionary of Faces for each layer

#### Boundary Elements
  Collection["BNumNodesLayer{0}".format(i)] = int(open('Layer'+str(i)+'.node').readlines()[0].split()[0]);  #### Initialy set the number of Nodes on the boundary of each layer as the number of nodes for each layer
  Collection["BNumEdgesLayer{0}".format(i)] = int(open('Layer'+str(i)+'.edge').readlines()[0].split()[0]);  #### Initialy set the number of Edges on the boundary of each layer as the number of edges for each layer
  Collection["BNodesLayer{0}".format(i)] = {};            #### Initialize Dictionary for the Nodes on the Boundary of each layer
  Collection["BEndingNodesLayer{0}".format(i)] = [];      #### Initialize List for the potential nodes on the boundary that may be ended on at each layer
  Collection["BNodeOutsLayer{0}".format(i)] = {};         #### Initialize Dictionary that stores the indices of terminating nodes connected to edges emanating from a node for each node on the boundary
  Collection["BNodeInsLayer{0}".format(i)] = {};          #### Initialize Dictionary that stores the indices of starting nodes connected to edges terminating at a node for each node on the boundary
  Collection["BNodeListLayer{0}".format(i)] = [];         #### Initialize List of Nodes on the Boundary of each layer
  Collection["BEdgeListLayer{0}".format(i)] = [];         #### Initialize List of Edges on the Boundary of each layer
  Collection["BEdgesLayer{0}".format(i)] = {};            #### Initialize Dictionary for the Edges on the Boundary of each layer
  Collection["BNodeIndicesLayer{0}".format(i)] = {};      #### Initialize Dictionary for the Indices of each node on the boundary of each layer
  Collection["BEdgeIndicesLayer{0}".format(i)] = {};      #### Initialize Dictionary for the Indices of each edge on the boundary of each layer

#### Interior Elements
  Collection["INumNodesLayer{0}".format(i)] = int(open('Layer'+str(i)+'.node').readlines()[0].split()[0]);  #### Initialy set number of interior nodes to the number of nodes on each layer
  Collection["INumEdgesLayer{0}".format(i)] = int(open('Layer'+str(i)+'.edge').readlines()[0].split()[0]);  #### Initial set the number of interior edges to the number of edges on each layer
  Collection["INodesLayer{0}".format(i)] = {};            #### Initialize Dictionary for the interior nodes on each layer
  Collection["INodeOutsLayer{0}".format(i)] = {};         #### Initialize Dictionary that stores the indices of terminating nodes connected to edges emanating from a node for each node in the interior
  Collection["INodeInsLayer{0}".format(i)] = {};          #### Initialize Dictionary that stores the indices of starting nodes connected to edges terminating at a node for each node in the interior
  Collection["INodeListLayer{0}".format(i)] = [];         #### Initialize List of Nodes in the interior of each layer
  Collection["IEdgesLayer{0}".format(i)] = {};            #### Initialize Dictionary of Edges in the interior of each layer
  Collection["INodeIndicesLayer{0}".format(i)] = {};      #### Initialize Dictionary of Indicies for the nodes in the interior of each layer
  Collection["IEdgeIndicesLayer{0}".format(i)] = {};      #### Initialize Dictionary of Indicies for the edges in the interior of each layer
  Collection["IEdgeListLayer{0}".format(i)] = [];         #### Initialize List of Edges in the interior of each layer


######################################################################################
####################      Define Global Variables and Restrictions     ###############
######################################################################################


timelim = 1500.0;                                         #### Set Cplex Time limit
InteriorGaplim = 0.05;                                            #### Set the gap limit to optimal solution for Cplex
BoundaryGap = 0.01;
Solim = 1;                                                #### Set the Number of Solutions necessary to obtain before Cplex terminates
NodeWeight = 0;                                           #### Set the weight placed on each node in the LP
PotentialRadius = 2*LayerGap;                             #### Set the radius for which to search for possible starting nodes on the next layer from the ending node on the current layer
BlueEdgePenalty = -8;
LayerMovePenalty = -0.25;

# Thickness Constraint and Overlap Incentive#

VtxRad = LayerGap+1;                                       #### Set the maximal accepatable distance a node can be away from the path without being on it
OverlapIncentive = 10;                                     #### Set the weight on each edge that overlaps an edge used in the previous layer
AnticipateIncentive = 3;                                   #### Set the weight on the edges in the current layer that are overlapped by existing edges in the next layer
CompConectsIncentive = 10;                                 #### Set the weight for utilizing edges that are added to connect disjoint components (these edges are not printed) very important to have a high weight if there are relatively small components


############# Simple Defs ################

def substract_lists(a, b):                                 #### Input is two vectors of the same length and output a vector of the same length whose entries are the difference in the corresponding entries of the input vectors (subtracts 1st - 2nd)
    if len(a) != len(b):
      print "Lists not length compatible";
      quit();
    else:
      val = zeros(len(a));
      for i in range(0,len(a)):
        val[i] = abs(a[i] - b[i]);
    return val

def CosAngleBtwConEdges(u1,u2,u3,LU,LV):                   #### Input is 3 vertices representing two connected edges with u2 the vertex in common, and LU (the length of edge u1,u2) and LV (the length of edge u2,u3)
  if len(u1) != len(u2) or len(u3) != len(u2):             #### Output is cosine of the angel between the two edges
    print "Lists not length compatible";
    quit();
  else:
    U1 = zeros(len(u1));
    V1 = zeros(len(u1));
    U2 = zeros(len(u1));
    V2 = zeros(len(u1));
    for i in range(0,len(u1)):
      U1[i] = u1[i] - u2[i];
      V1[i] = u3[i] - u2[i];
      U2[i] = u2[i] - u1[i];
      V2[i] = u2[i] - u3[i];
    U1V1= dot(U1,V1)/(LU*LV);
    U2V2= dot(U2,V2)/(LU*LV);
    return max(U1V1,U2V2)


##############################################################################
###      Original Mesh: Import Node Information into dictionay format      ###
##############################################################################
###    Stores vertices by vertex number as [float_1, float_2, float_3]     ###
##############################################################################

for j in range(1,NumLayers+1):
  i=0;
  with open('Layer'+str(j)+'.node') as f:                 #### open file titled Layerj.node for j = 1, 2, 3,...
    for line in f.readlines()[1:Collection['NumNodesLayer'+str(j)]+1]:    #### read in the line
      i=i+1;
#       node = [float(line.split()[1]),float(line.split()[2]),float(line.split()[3])];  #### When meshes are given and z coord defined
      node = [float(line.split()[1]),float(line.split()[2]),float(LayerGap*(j-1))];  #### split the line and for Demo Purposes when z coord is undefined use the layer*gap
      Collection['NodesLayer'+str(j)][i] = [node[0],node[1],node[2]];                #### store node coordinate information by index
      Collection['NodeIndicesLayer'+str(j)][node[0],node[1],node[2]] = [i];          #### store node index information by coordinates
      Collection['NodeListLayer'+str(j)].append([node[0],node[1],node[2]]);          #### store nodes into searchable list


###############################################################################################################
################################  Begin Boundary Information Grab #############################################
###############################################################################################################

for j in range(1,NumLayers+1):
  i = 0;
  l = 0;
  InteriorEdgeList=InteriorEdgeList_Grab(j);
  BNodeList = [];
  with open('Layer'+str(j)+'.edge') as f:    #####  Go through the list of edges and if they are not in the interior list then they are on the boundary and so extract and build Boundary Node info
    for line in f.readlines()[1::]:
      edge = sorted([int(line.split()[1]),int(line.split()[2])])
      if edge not in InteriorEdgeList:
        for k in range(0,2):
          if edge[k] not in BNodeList:
            BNodeList.append(edge[k]);
            i = i+1;
            node = Collection['NodesLayer'+str(j)][edge[k]];
            Collection['BNodesLayer'+str(j)][i] = [node[0],node[1],node[2]];
            Collection['BNodeOutsLayer'+str(j)][i]= [];
            Collection['BNodeInsLayer'+str(j)][i]= [];
            Collection['BNodeIndicesLayer'+str(j)][node[0],node[1],node[2]] = [i];
            Collection['BNodeListLayer'+str(j)].append([node[0],node[1],node[2]]);
            Collection['BEndingNodesLayer'+str(j)].append(i);
            i = i+1;
            node = Collection['NodesLayer'+str(j)][edge[k]];
            node = [node[0],node[1],-j];
            Collection['BNodesLayer'+str(j)][i] = [node[0],node[1],node[2]];
            Collection['BNodeOutsLayer'+str(j)][i]= [];
            Collection['BNodeInsLayer'+str(j)][i]= [];
            Collection['BNodeIndicesLayer'+str(j)][node[0],node[1],node[2]] = [i];
            Collection['BNodeListLayer'+str(j)].append([node[0],node[1],node[2]]);
            Collection['BEndingNodesLayer'+str(j)].append(i);

          if k == 1:
            l=l+1;
            Len = abs(scipy.spatial.distance.pdist([array(Collection['NodesLayer'+str(j)][edge[0]]),array(Collection['NodesLayer'+str(j)][edge[1]])]));

            nodea=[Collection['NodesLayer'+str(j)][edge[0]][0],Collection['NodesLayer'+str(j)][edge[0]][1],Collection['NodesLayer'+str(j)][edge[0]][2]];
            nodeb=[Collection['NodesLayer'+str(j)][edge[1]][0],Collection['NodesLayer'+str(j)][edge[1]][1],Collection['NodesLayer'+str(j)][edge[0]][2]];
            Collection['BEdgesLayer'+str(j)][l]=[Collection['BNodeIndicesLayer'+str(j)][nodea[0],nodea[1],nodea[2]][0],Collection['BNodeIndicesLayer'+str(j)][nodeb[0],nodeb[1],nodeb[2]][0],Len[0]];

            l = l+1;
            nodea=[Collection['NodesLayer'+str(j)][edge[0]][0],Collection['NodesLayer'+str(j)][edge[0]][1],-j];
            nodeb=[Collection['NodesLayer'+str(j)][edge[1]][0],Collection['NodesLayer'+str(j)][edge[1]][1],-j];
            Collection['BEdgesLayer'+str(j)][l]=[Collection['BNodeIndicesLayer'+str(j)][nodea[0],nodea[1],nodea[2]][0],Collection['BNodeIndicesLayer'+str(j)][nodeb[0],nodeb[1],nodeb[2]][0],Len[0]];

  Collection["BNumNodesLayer"+str(j)] = int(i);  ####  Set the number of boundary nodes to be the number of boundary nodes found
  Collection["BNumEdgesLayer"+str(j)] = int(l);  ####  Set the number of boundary edges to be the number of boundary edges found
  Collection["BONumEdgesLayer"+str(j)] = int(l);  ####  Set the number of boundary edges original to the mesh

############################ Plotting and convex fill
##### utilize package "triangle" to take boundary nodes and quickly mesh them without adding more vertices. We store them in a list csegments

  Collection['BConvexEdgeIndicesListLayer'+str(j)] = [];
  cvertices = {};
  cvertices['vertices']=[];
  for line in Collection['BNodesLayer'+str(j)]:
    line = Collection['BNodesLayer'+str(j)][line];
    if [float(line[0]),float(line[1])] not in cvertices['vertices']:
      cvertices['vertices'].append([float(line[0]),float(line[1])]);
  cvertices['vertices'] = array(cvertices['vertices']);


  csegments=[];
  t = triangle.triangulate(cvertices, 'c');
  for tri in t['triangles']:
    if sorted([tri[0],tri[1]]) not in csegments:
      csegments.append(sorted([tri[0],tri[1]]));
    if sorted([tri[2],tri[1]]) not in csegments:
      csegments.append(sorted([tri[1],tri[2]]));
    if sorted([tri[0],tri[2]]) not in csegments:
      csegments.append(sorted([tri[0],tri[2]]));

  # triangle.plot.plot(plt.axes(), vertices=cvertices['vertices'], segments=csegments)    ##### for plotting purposes and not part of the main code
  # plt.show()
  # quit()


############################# Calculating number of new edges ##########################################

  Numnewsegments = int(4*(len(csegments)-Collection["BNumEdgesLayer"+str(j)]/2.0)+Collection["BNumNodesLayer"+str(j)]);   #### we calculate the number of new edges we will add to the boundary edges in addition to edges in the opposite direction of the current edges

################## Boundary Edge Information extraction and Convex Fill Incorporation

  i=0;

  for line in range(1,Collection["BNumEdgesLayer"+str(j)]+1):  #### We now add the edges going in the opposite direction to the edges we already have to our dictionaries and lists
    edge = Collection['BEdgesLayer'+str(j)][line];
    i=i+1;

    Collection['BEdgeIndicesLayer'+str(j)][edge[0],edge[1]] = [i];
    Collection['BEdgeIndicesLayer'+str(j)][edge[1],edge[0]] = [int(Collection['BNumEdgesLayer'+str(j)]+Numnewsegments/2.0+i)]; ##### we use Numnewsegments/2.0 since half the number of new segments includes edges going in both directions

    Collection['BEdgesLayer'+str(j)][int(Collection['BNumEdgesLayer'+str(j)]+Numnewsegments/2.0+i)]=[edge[1],edge[0],edge[2]];

    Collection['BEdgeListLayer'+str(j)].append([edge[0],edge[1]]);
    Collection['BEdgeListLayer'+str(j)].append([edge[1],edge[0]]);

    Nout = Collection['BNodeOutsLayer'+str(j)][edge[0]];
    Nout.append(edge[1]);
    Collection['BNodeOutsLayer'+str(j)][edge[0]]= Nout;

    Nout = Collection['BNodeOutsLayer'+str(j)][edge[1]];
    Nout.append(edge[0]);
    Collection['BNodeOutsLayer'+str(j)][edge[1]]= Nout;

    Nin = Collection['BNodeInsLayer'+str(j)][edge[1]];
    Nin.append(edge[0]);
    Collection['BNodeInsLayer'+str(j)][edge[1]]= Nin;

    Nin = Collection['BNodeInsLayer'+str(j)][edge[0]];
    Nin.append(edge[1]);
    Collection['BNodeInsLayer'+str(j)][edge[0]]= Nin;

  Collection['BNoPrintEdgeIndLayer'+str(j)] = [];
  count =0;

  for cedge in csegments:   ##### we now add the new edges from the convex fill to our information dictionaries and lists
    cedge = [int(cedge[0]),int(cedge[1])];

    a = Collection['BNodeIndicesLayer'+str(j)][cvertices['vertices'][cedge[0]][0],cvertices['vertices'][cedge[0]][1],float(LayerGap*(j-1))][0];
    b = Collection['BNodeIndicesLayer'+str(j)][cvertices['vertices'][cedge[1]][0],cvertices['vertices'][cedge[1]][1],float(LayerGap*(j-1))][0];

    c = Collection['BNodeIndicesLayer'+str(j)][cvertices['vertices'][cedge[0]][0],cvertices['vertices'][cedge[0]][1],-j][0];
    d = Collection['BNodeIndicesLayer'+str(j)][cvertices['vertices'][cedge[1]][0],cvertices['vertices'][cedge[1]][1],-j][0];


#     cedge = [Collection['BNodeIndicesLayer'+str(j)][cvertices['vertices'][cedge[0]][0],cvertices['vertices'][cedge[0]][1],float(LayerGap*(j-1))][0],Collection['BNodeIndicesLayer'+str(j)][cvertices['vertices'][cedge[1]][0],cvertices['vertices'][cedge[1]][1],float(LayerGap*(j-1))][0]];
    cedge = sorted([a,b])

    if cedge not in Collection['BEdgeListLayer'+str(j)]:
      count = count +1;
      i = i+1;
      Collection['BEdgeListLayer'+str(j)].append([cedge[0],cedge[1]]);
      Collection['BEdgeIndicesLayer'+str(j)][cedge[0],cedge[1]] = [i];

      Len = abs(scipy.spatial.distance.pdist([array(Collection['NodesLayer'+str(j)][cedge[0]]),array(Collection['NodesLayer'+str(j)][cedge[1]])]));

      Collection['BEdgesLayer'+str(j)][i]=[cedge[0],cedge[1], BlueEdgePenalty*Len]; #### Here we add a penalty for choosing non-vertical blue edges
      Collection['BEdgeListLayer'+str(j)].append([cedge[1],cedge[0]]);

      #### adding new edge in opposite direction
      Collection['BEdgeIndicesLayer'+str(j)][cedge[1],cedge[0]] = [int(Collection['BNumEdgesLayer'+str(j)]+i+Numnewsegments/2.0)];
      Collection['BEdgesLayer'+str(j)][int(Collection['BNumEdgesLayer'+str(j)]+i+Numnewsegments/2.0)]=[cedge[1],cedge[0],BlueEdgePenalty*Len];

      Collection['BNoPrintEdgeIndLayer'+str(j)].append(i);
      Collection['BNoPrintEdgeIndLayer'+str(j)].append(int(Collection['BNumEdgesLayer'+str(j)]+i+Numnewsegments/2.0));

      Nout = Collection['BNodeOutsLayer'+str(j)][cedge[0]];
      Nout.append(cedge[1]);
      Collection['BNodeOutsLayer'+str(j)][cedge[0]]= Nout;

      Nout = Collection['BNodeOutsLayer'+str(j)][cedge[1]];
      Nout.append(cedge[0]);
      Collection['BNodeOutsLayer'+str(j)][cedge[1]]= Nout;

      Nin = Collection['BNodeInsLayer'+str(j)][cedge[1]];
      Nin.append(cedge[0]);
      Collection['BNodeInsLayer'+str(j)][cedge[1]]= Nin;

      Nin = Collection['BNodeInsLayer'+str(j)][cedge[0]];
      Nin.append(cedge[1]);
      Collection['BNodeInsLayer'+str(j)][cedge[0]]= Nin;


      ######## adding duplicate edge information
      count = count +1
      i = i+1;
      Collection['BEdgeListLayer'+str(j)].append([c,d]);
      Collection['BEdgeIndicesLayer'+str(j)][c,d] = [i];
      Collection['BEdgesLayer'+str(j)][i]=[c,d,BlueEdgePenalty*Len];
      Collection['BEdgeListLayer'+str(j)].append([c,d]);


      #### adding new edge in opposite direction
      Collection['BEdgeIndicesLayer'+str(j)][d,c] = [int(Collection['BNumEdgesLayer'+str(j)]+i+Numnewsegments/2.0)];
      Collection['BEdgesLayer'+str(j)][int(Collection['BNumEdgesLayer'+str(j)]+i+Numnewsegments/2.0)]=[d,c,BlueEdgePenalty*Len];

      Collection['BNoPrintEdgeIndLayer'+str(j)].append(i);
      Collection['BNoPrintEdgeIndLayer'+str(j)].append(int(Collection['BNumEdgesLayer'+str(j)]+i+Numnewsegments/2.0));


      Nout = Collection['BNodeOutsLayer'+str(j)][c];
      Nout.append(d);
      Collection['BNodeOutsLayer'+str(j)][c]= Nout;

      Nout = Collection['BNodeOutsLayer'+str(j)][d];
      Nout.append(c);
      Collection['BNodeOutsLayer'+str(j)][d]= Nout;

      Nin = Collection['BNodeInsLayer'+str(j)][d];
      Nin.append(c);
      Collection['BNodeInsLayer'+str(j)][d]= Nin;

      Nin = Collection['BNodeInsLayer'+str(j)][c];
      Nin.append(d);
      Collection['BNodeInsLayer'+str(j)][c]= Nin;


################################## Add in edges between duplicate nodes, all duplicate nodes are listed 1 apart

  for num in range(0,Collection["BNumNodesLayer"+str(j)]/2):
    cedge = [int(2*num+1),int(2*num+2)];

    if cedge not in Collection['BEdgeListLayer'+str(j)]:
      count = count +1;
      i = i+1;
      Collection['BEdgeListLayer'+str(j)].append([cedge[0],cedge[1]]);
      Collection['BEdgeIndicesLayer'+str(j)][cedge[0],cedge[1]] = [i];
      Collection['BEdgesLayer'+str(j)][i]=[cedge[0],cedge[1],LayerMovePenalty];
      Collection['BEdgeListLayer'+str(j)].append([cedge[1],cedge[0]]);

      #### adding new edge in opposite direction
      Collection['BEdgeIndicesLayer'+str(j)][cedge[1],cedge[0]] = [int(Collection['BNumEdgesLayer'+str(j)]+i+Numnewsegments/2.0)];
      Collection['BEdgesLayer'+str(j)][int(Collection['BNumEdgesLayer'+str(j)]+i+Numnewsegments/2.0)]=[cedge[1],cedge[0],LayerMovePenalty];

      Collection['BNoPrintEdgeIndLayer'+str(j)].append(i);
      Collection['BNoPrintEdgeIndLayer'+str(j)].append(int(Collection['BNumEdgesLayer'+str(j)]+i+Numnewsegments/2.0));

      Nout = Collection['BNodeOutsLayer'+str(j)][cedge[0]];
      Nout.append(cedge[1]);
      Collection['BNodeOutsLayer'+str(j)][cedge[0]]= Nout;

      Nout = Collection['BNodeOutsLayer'+str(j)][cedge[1]];
      Nout.append(cedge[0]);
      Collection['BNodeOutsLayer'+str(j)][cedge[1]]= Nout;

      Nin = Collection['BNodeInsLayer'+str(j)][cedge[1]];
      Nin.append(cedge[0]);
      Collection['BNodeInsLayer'+str(j)][cedge[1]]= Nin;

      Nin = Collection['BNodeInsLayer'+str(j)][cedge[0]];
      Nin.append(cedge[1]);
      Collection['BNodeInsLayer'+str(j)][cedge[0]]= Nin;


  Collection['BNumEdgesLayer'+str(j)] = len(Collection['BEdgesLayer'+str(j)]);   ###### This number of boundary edges includes all the edges (convex fill, forward, backwards and duplicates)


###############################################################################################
########### End Boundary Information Grab ############################
###################################################################################



#############################################################################################
###   Interior Node and Edge Information into dictionay format   ###
#############################################################################################
#### We repeat the same process for the interior nodes and edges except that we no longer have to add duplicate nodes or edges

  INodemap = {};
  i=0;
  count = 0;
  for node in Collection['NodesLayer'+str(j)]:
    count = count+1;
    if Collection['NodesLayer'+str(j)][node] not in Collection['BNodeListLayer'+str(j)]:
      node=[Collection['NodesLayer'+str(j)][node][0],Collection['NodesLayer'+str(j)][node][1],float(LayerGap*(j-1))];
      i = i+1;
      INodemap[count]=i;
      Collection['INodesLayer'+str(j)][i] = [node[0],node[1],node[2]];
      Collection['INodeOutsLayer'+str(j)][i]= [];
      Collection['INodeInsLayer'+str(j)][i]= [];
      Collection['INodeIndicesLayer'+str(j)][node[0],node[1],node[2]] = [i];
      Collection['INodeListLayer'+str(j)].append([node[0],node[1],node[2]]);

  Collection["INumNodesLayer"+str(j)] = int(i);

######## Interior edges utilizing reassignment of node indice information

  i = 0;
  with open('Layer'+str(j)+'.edge') as f:
    for line in f.readlines()[1:Collection['NumEdgesLayer'+str(j)]+1]:
      edge1= sort([int(line.split()[1]), int(line.split()[2])]);
      if [Collection['NodesLayer'+str(j)][edge1[0]][0],Collection['NodesLayer'+str(j)][edge1[0]][1],float(LayerGap*(j-1))] in Collection['INodeListLayer'+str(j)] and [Collection['NodesLayer'+str(j)][edge1[1]][0],Collection['NodesLayer'+str(j)][edge1[1]][1],float(LayerGap*(j-1))] in Collection['INodeListLayer'+str(j)]:

        i=i+1;
        nodea=[Collection['NodesLayer'+str(j)][edge1[0]][0],Collection['NodesLayer'+str(j)][edge1[0]][1],float(LayerGap*(j-1))];
        nodeb=[Collection['NodesLayer'+str(j)][edge1[1]][0],Collection['NodesLayer'+str(j)][edge1[1]][1],float(LayerGap*(j-1))];
        Collection['IEdgesLayer'+str(j)][i]=[Collection['INodeIndicesLayer'+str(j)][nodea[0],nodea[1],nodea[2]][0],Collection['INodeIndicesLayer'+str(j)][nodeb[0],nodeb[1],nodeb[2]][0],Len[0]];


  Collection["INumEdgesLayer"+str(j)] = int(i);
  Collection["IONumEdgesLayer"+str(j)] = int(i);


############################ Plotting and convex fill

  Collection['IConvexEdgeIndicesListLayer'+str(j)] = [];
  cvertices = {};
  cvertices['vertices']=[];
  for line in Collection['INodesLayer'+str(j)]:
    line = Collection['INodesLayer'+str(j)][line];
    if [float(line[0]),float(line[1])] not in cvertices['vertices']:
      cvertices['vertices'].append([float(line[0]),float(line[1])]);
  cvertices['vertices'] = array(cvertices['vertices']);

  cvertices['segments']=[];
  for line in Collection['IEdgesLayer'+str(j)]:
    line = Collection['IEdgesLayer'+str(j)][line];
    cvertices['segments'].append(sorted([int(line[0])-1,int(line[1])-1]));
  csegments=cvertices['segments'];
  cvertices['segments'] = array(cvertices['segments']);

  with open('Layer'+str(j)+'.ele') as f:
    cvertices['triangles']=[];
    for line in f.readlines()[1::]:
      if int(line.split()[1]) in INodemap and int(line.split()[2]) in INodemap and int(line.split()[3]) in INodemap:
        line = [INodemap[int(line.split()[1])],INodemap[int(line.split()[2])],INodemap[int(line.split()[3])]];
        cvertices['triangles'].append(sorted([int(line[0])-1,int(line[1])-1,int(line[2])-1]));
    cvertices['triangles'] = array(cvertices['triangles']);

  # csegments=[]

  t = triangle.triangulate(cvertices, 'c');

  for tri in t['triangles']:
    Len = abs(scipy.spatial.distance.pdist([array(Collection['INodesLayer'+str(j)][tri[0]+1]),array(Collection['INodesLayer'+str(j)][tri[1]+1])])[0]);
    # print Len, LayerGap*sqrt(2);
    if sorted([tri[0],tri[1]]) not in csegments and Len > LayerGap*sqrt(2)+0.01:
      csegments.append(sorted([tri[0],tri[1]]));
      # print sorted([tri[0],tri[1]])

    Len = abs(scipy.spatial.distance.pdist([array(Collection['INodesLayer'+str(j)][tri[2]+1]),array(Collection['INodesLayer'+str(j)][tri[1]+1])])[0]);
    # print Len, LayerGap*sqrt(2);
    if sorted([tri[2],tri[1]]) not in csegments and Len > LayerGap*sqrt(2)+0.01:
      csegments.append(sorted([tri[1],tri[2]]));
      sorted([tri[0],tri[1]])

    Len = abs(scipy.spatial.distance.pdist([array(Collection['INodesLayer'+str(j)][tri[0]+1]),array(Collection['INodesLayer'+str(j)][tri[2]+1])])[0]);
    # print Len, LayerGap*sqrt(2);
    if sorted([tri[0],tri[2]]) not in csegments and Len > LayerGap*sqrt(2)+0.01:
      csegments.append(sorted([tri[0],tri[2]]));
      sorted([tri[0],tri[1]])


  # triangle.plot.plot(plt.axes(), vertices=cvertices['vertices'], segments=csegments)   ##### for plotting purposes and not connected with main code
  # plt.show()
  # triangle.plot.plot(plt.axes(), vertices=cvertices['vertices'], triangles=cvertices['triangles'])
  # plt.show()
  # triangle.plot.plot(plt.axes(), vertices=t['vertices'], triangles=t['triangles'])
  # plt.show()
  # quit()

############################# Calculating number of new edges ##########################################

  Numnewsegments = len(csegments)-Collection["INumEdgesLayer"+str(j)];   ##### this number does not take into account the edges that will go in the opposite direction of the new edges

############################# Extracting remaining and new Interior Edge Information ##########################

  i=0;

  for line in range(1,Collection["INumEdgesLayer"+str(j)]+1):
    line = Collection['IEdgesLayer'+str(j)][line];
    i=i+1;
    edge = line;

    Collection['IEdgeIndicesLayer'+str(j)][edge[0],edge[1]] = [i];
    Collection['IEdgeIndicesLayer'+str(j)][edge[1],edge[0]] = [Collection['INumEdgesLayer'+str(j)]+i+Numnewsegments];

    Collection['IEdgesLayer'+str(j)][Collection['INumEdgesLayer'+str(j)]+i+Numnewsegments]=[edge[1],edge[0],edge[2]];

    Collection['IEdgeListLayer'+str(j)].append([edge[0],edge[1]]);
    Collection['IEdgeListLayer'+str(j)].append([edge[1],edge[0]]);

    Nout = Collection['INodeOutsLayer'+str(j)][edge[0]];
    Nout.append(edge[1]);
    Collection['INodeOutsLayer'+str(j)][edge[0]]= Nout;

    Nout = Collection['INodeOutsLayer'+str(j)][edge[1]];
    Nout.append(edge[0]);
    Collection['INodeOutsLayer'+str(j)][edge[1]]= Nout;

    Nin = Collection['INodeInsLayer'+str(j)][edge[1]];
    Nin.append(edge[0]);
    Collection['INodeInsLayer'+str(j)][edge[1]]= Nin;

    Nin = Collection['INodeInsLayer'+str(j)][edge[0]];
    Nin.append(edge[1]);
    Collection['INodeInsLayer'+str(j)][edge[0]]= Nin;

  Collection['INoPrintEdgeIndLayer'+str(j)] = [];
  count =0;
  for cedge in csegments:
    cedge = [int(cedge[0]),int(cedge[1])];

    a = Collection['INodeIndicesLayer'+str(j)][cvertices['vertices'][cedge[0]][0],cvertices['vertices'][cedge[0]][1],float(LayerGap*(j-1))][0];
    b = Collection['INodeIndicesLayer'+str(j)][cvertices['vertices'][cedge[1]][0],cvertices['vertices'][cedge[1]][1],float(LayerGap*(j-1))][0];
    cedge = sorted([a,b]);

    if cedge not in Collection['IEdgeListLayer'+str(j)]:
      count = count +1;
      i = i+1;
      Collection['IEdgeListLayer'+str(j)].append([cedge[0],cedge[1]]);
      Collection['IEdgeIndicesLayer'+str(j)][cedge[0],cedge[1]] = [i];
      Len = abs(scipy.spatial.distance.pdist([array(Collection['NodesLayer'+str(j)][cedge[0]]),array(Collection['NodesLayer'+str(j)][cedge[1]])]));
      Collection['IEdgesLayer'+str(j)][i]=[cedge[0],cedge[1],BlueEdgePenalty*Len];
      Collection['IEdgeListLayer'+str(j)].append([cedge[1],cedge[0]]);
      Collection['IEdgeIndicesLayer'+str(j)][cedge[1],cedge[0]] = [Collection['INumEdgesLayer'+str(j)]+i+Numnewsegments];
      Collection['IEdgesLayer'+str(j)][Collection['INumEdgesLayer'+str(j)]+i+Numnewsegments]=[cedge[1],cedge[0],BlueEdgePenalty*Len];
      Collection['INoPrintEdgeIndLayer'+str(j)].append(i);
      Collection['INoPrintEdgeIndLayer'+str(j)].append(Collection['INumEdgesLayer'+str(j)]+i+Numnewsegments);

      Nout = Collection['INodeOutsLayer'+str(j)][cedge[0]];
      Nout.append(cedge[1]);
      Collection['INodeOutsLayer'+str(j)][cedge[0]]= Nout;

      Nout = Collection['INodeOutsLayer'+str(j)][cedge[1]];
      Nout.append(cedge[0]);
      Collection['INodeOutsLayer'+str(j)][cedge[1]]= Nout;

      Nin = Collection['INodeInsLayer'+str(j)][cedge[1]];
      Nin.append(cedge[0]);
      Collection['INodeInsLayer'+str(j)][cedge[1]]= Nin;

      Nin = Collection['INodeInsLayer'+str(j)][cedge[0]];
      Nin.append(cedge[1]);
      Collection['INodeInsLayer'+str(j)][cedge[0]]= Nin;


  Collection['INumEdgesLayer'+str(j)] = len(Collection['IEdgesLayer'+str(j)]);   ############# This is all of the edges in the interior


###############################################################################################
########### End Interior Information Grab ############################
###################################################################################





#################################################################################################################################################






#######################################################################################################
###                                     Cplex  Loop                                                 ###
#######################################################################################################



###########################################################################
####     Main Problem Setup: Objective Statement maximize edgedists    ####
###########################################################################
####   max u^T*z with u=edgedists,edgedists, NodePosition  ####
###########################################################################


for k in range(1,NumLayers+1):                                                #### Set the Dictionaries and Lists to use for the current layer and pass (boundary or interior)
  for Boundary in [1,0]:
    if Boundary == 1:                                                         #### Setup for Boundary pass
      NumNodes = Collection['BNumNodesLayer'+str(k)]
      NumEdges = Collection['BNumEdgesLayer'+str(k)]
      MaxEdges = NumEdges;
      Nodes = Collection["BNodesLayer{0}".format(k)];
      NodeOuts = Collection["BNodeOutsLayer{0}".format(k)];
      NodeIns = Collection["BNodeInsLayer{0}".format(k)];
      NodeIndices = Collection["BNodeIndicesLayer{0}".format(k)];
      NodeList = Collection["BNodeListLayer{0}".format(k)];
      Edges = Collection["BEdgesLayer{0}".format(k)];
      EdgeIndices  = Collection["BEdgeIndicesLayer{0}".format(k)];
      EdgeList = Collection["BEdgeListLayer{0}".format(k)];
      NoPrintEdgeInds = Collection['BNoPrintEdgeIndLayer'+str(k)];

    if Boundary == 0:                                                          #### Setup for Interior pass
      NumNodes = Collection['INumNodesLayer'+str(k)]
      NumEdges = Collection['INumEdgesLayer'+str(k)]
      MaxEdges = NumEdges;
      Nodes = Collection["INodesLayer{0}".format(k)];
      NodeOuts = Collection["INodeOutsLayer{0}".format(k)];
      NodeIns = Collection["INodeInsLayer{0}".format(k)];
      NodeIndices = Collection["INodeIndicesLayer{0}".format(k)];
      NodeList = Collection["INodeListLayer{0}".format(k)];
      Edges = Collection["IEdgesLayer{0}".format(k)];
      EdgeIndices  = Collection["IEdgeIndicesLayer{0}".format(k)];
      EdgeList = Collection["IEdgeListLayer{0}".format(k)];
      NoPrintEdgeInds = Collection['INoPrintEdgeIndLayer'+str(k)];

    Overlaps = [];      #### Initiate lists pertaining to each pass for each layer
    EndNodes = [];
    StartNodes = [];



  ######### Potential Nodes  ########
  ###### Define the list of potential starting and ending nodes for each pass on each layer

    if k ==1:
      if Boundary ==1:
        NumStartNodes = int(NumNodes/2.0)
        for i in range(0,NumStartNodes):
          StartNodes.append(i+1);

        for i in range(NumStartNodes,NumNodes):
          EndNodes.append(i+1);

      else:
        LastNode = Collection['BNodesLayer'+str(k)][list(Uj).index(max(Uj))+1];
        LastNode = [LastNode[0],LastNode[1],float(LayerGap*(k-1))];
        for j in Nodes:
          Node1=[Nodes[j][0],Nodes[j][1],0.0];
          if scipy.spatial.distance.pdist([array(Node1),array(LastNode)])[0] <= PotentialRadius:
            StartNodes.append(j);

        for i in range(0,NumNodes):
          EndNodes.append(i+1);


  ######### Potential Nodes and Overlap Edges ########
    if k >= 2:

      if k == 2:
        if Boundary == 1:
          LastNode = Collection['INodesLayer'+str(k-1)][list(Uj).index(max(Uj))+1];
          LastNode = [LastNode[0],LastNode[1],float(LayerGap*(k-2))];
        else:
          LastNode = Collection['BNodesLayer'+str(k)][list(Uj).index(max(Uj))+1];
          LastNode = [LastNode[0],LastNode[1],float(LayerGap*(k-1))];

        for j in Nodes:
          Node1=[Nodes[j][0],Nodes[j][1],float(LayerGap*(k-1))];
          if scipy.spatial.distance.pdist([array(Node1),array(LastNode)])[0] <= PotentialRadius:
            StartNodes.append(j);

      else:
        if Boundary == 1:
          LastNode = [PBStartNode[0],PBStartNode[1],float(LayerGap*(k-1))];
        else:
          LastNode = [PIStartNode[0],PIStartNode[1],float(LayerGap*(k-1))];

        if LastNode in NodeList:
          StartNodes.append(NodeIndices[LastNode[0],LastNode[1],LastNode[2]][0]);


        else:
          for j in Nodes:
            Node1=[Nodes[j][0],Nodes[j][1],float(LayerGap*(k-1))];
            if scipy.spatial.distance.pdist([array(Node1),array(LastNode)])[0] <= PotentialRadius:
              StartNodes.append(j);

      if Boundary == 1:
        PBEdgePath = BEdgePath;
        for i in range(0,len(BEdgePath)):
          if BEdgePath[i] != 0:
            PrevEdge = Collection['BEdgesLayer'+str(k-1)][i+1];
            coord1= Collection['BNodesLayer'+str(k-1)][PrevEdge[0]];
            coord2= Collection['BNodesLayer'+str(k-1)][PrevEdge[1]];

            if [coord1[0],coord1[1],LayerGap*(k-1)] in Collection['BNodeListLayer'+str(k)] and [coord2[0],coord2[1],LayerGap*(k-1)] in Collection['BNodeListLayer'+str(k)]:
              Ind1 = Collection['BNodeIndicesLayer'+str(k)][coord1[0],coord1[1],LayerGap*(k-1)][0]; ###### For demo z coordinates are staged to be LayerGap apart
              Ind2 = Collection['BNodeIndicesLayer'+str(k)][coord2[0],coord2[1],LayerGap*(k-1)][0]; ###### For demo z coordinates are staged to be LayerGap apart

              if [Ind1,Ind2] in Collection['BEdgeListLayer'+str(k)]:
                EInd = Collection['BEdgeIndicesLayer'+str(k)][Ind1,Ind2][0];
                if EInd not in Overlaps:
                  Overlaps.append(EInd);
                  Collection['BEdgesLayer'+str(k)][EInd][2]=Collection['BEdgesLayer'+str(k)][EInd][2]*OverlapIncentive;

        for i in range(0,len(IEdgePath)):
          if IEdgePath[i] != 0:
            PrevEdge = Collection['IEdgesLayer'+str(k-1)][i+1];
            coord1= Collection['INodesLayer'+str(k-1)][PrevEdge[0]];
            coord2= Collection['INodesLayer'+str(k-1)][PrevEdge[1]];

            if [coord1[0],coord1[1],LayerGap*(k-1)] in Collection['BNodeListLayer'+str(k)] and [coord2[0],coord2[1],LayerGap*(k-1)] in Collection['BNodeListLayer'+str(k)]:
              Ind1 = Collection['BNodeIndicesLayer'+str(k)][coord1[0],coord1[1],LayerGap*(k-1)][0]; ###### For demo z coordinates are staged to be LayerGap apart
              Ind2 = Collection['BNodeIndicesLayer'+str(k)][coord2[0],coord2[1],LayerGap*(k-1)][0]; ###### For demo z coordinates are staged to be LayerGap apart

              if [Ind1,Ind2] in Collection['BEdgeListLayer'+str(k)]:
                EInd = Collection['BEdgeIndicesLayer'+str(k)][Ind1,Ind2][0];
                if EInd not in Overlaps:
                  Overlaps.append(EInd);
                  Collection['BEdgesLayer'+str(k)][EInd][2]=Collection['BEdgesLayer'+str(k)][EInd][2]*OverlapIncentive;

      else:
        for i in range(0,len(PBEdgePath)):
          if PBEdgePath[i] != 0:
            PrevEdge = Collection['BEdgesLayer'+str(k-1)][i+1];
            coord1= Collection['BNodesLayer'+str(k-1)][PrevEdge[0]];
            coord2= Collection['BNodesLayer'+str(k-1)][PrevEdge[1]];

            if [coord1[0],coord1[1],LayerGap*(k-1)] in Collection['INodeListLayer'+str(k)] and [coord2[0],coord2[1],LayerGap*(k-1)] in Collection['INodeListLayer'+str(k)]:
              Ind1 = Collection['INodeIndicesLayer'+str(k)][coord1[0],coord1[1],LayerGap*(k-1)][0]; ###### For demo z coordinates are staged to be LayerGap apart
              Ind2 = Collection['INodeIndicesLayer'+str(k)][coord2[0],coord2[1],LayerGap*(k-1)][0]; ###### For demo z coordinates are staged to be LayerGap apart

              if [Ind1,Ind2] in Collection['IEdgeListLayer'+str(k)]:
                EInd = Collection['IEdgeIndicesLayer'+str(k)][Ind1,Ind2][0];
                if EInd not in Overlaps:
                  Overlaps.append(EInd);
                  Collection['IEdgesLayer'+str(k)][EInd][2]=Collection['IEdgesLayer'+str(k)][EInd][2]*OverlapIncentive;

        for i in range(0,len(IEdgePath)):
          if IEdgePath[i] != 0:
            PrevEdge = Collection['IEdgesLayer'+str(k-1)][i+1];
            coord1= Collection['INodesLayer'+str(k-1)][PrevEdge[0]];
            coord2= Collection['INodesLayer'+str(k-1)][PrevEdge[1]];

            if [coord1[0],coord1[1],LayerGap*(k-1)] in Collection['INodeListLayer'+str(k)] and [coord2[0],coord2[1],LayerGap*(k-1)] in Collection['INodeListLayer'+str(k)]:
              Ind1 = Collection['INodeIndicesLayer'+str(k)][coord1[0],coord1[1],LayerGap*(k-1)][0]; ###### For demo z coordinates are staged to be LayerGap apart
              Ind2 = Collection['INodeIndicesLayer'+str(k)][coord2[0],coord2[1],LayerGap*(k-1)][0]; ###### For demo z coordinates are staged to be LayerGap apart

              if [Ind1,Ind2] in Collection['IEdgeListLayer'+str(k)]:
                EInd = Collection['IEdgeIndicesLayer'+str(k)][Ind1,Ind2][0];
                if EInd not in Overlaps:
                  Overlaps.append(EInd);
                  Collection['IEdgesLayer'+str(k)][EInd][2]=Collection['IEdgesLayer'+str(k)][EInd][2]*OverlapIncentive;

      for i in range(0,NumNodes):
        EndNodes.append(i+1);

  ##### Anticpation Edge Weights######

      if k < NumLayers:
        if Boundary == 1:
          for i in Collection['BEdgesLayer'+str(k+1)]:
            coord1= Collection['BNodesLayer'+str(k+1)][Collection['BEdgesLayer'+str(k+1)][i][0]];
            coord2= Collection['BNodesLayer'+str(k+1)][Collection['BEdgesLayer'+str(k+1)][i][1]];

            if [coord1[0],coord1[1],LayerGap*(k-1)] in Collection['BNodeListLayer'+str(k)] and [coord2[0],coord2[1],LayerGap*(k-1)] in Collection['BNodeListLayer'+str(k)]:
              Ind1 = Collection['BNodeIndicesLayer'+str(k)][coord1[0],coord1[1],LayerGap*(k-1)][0];
              Ind2 = Collection['BNodeIndicesLayer'+str(k)][coord2[0],coord2[1],LayerGap*(k-1)][0];

              if [Ind1,Ind2] in Collection['BEdgeListLayer'+str(k)]:
                EInd = Collection['BEdgeIndicesLayer'+str(k)][Ind1,Ind2][0];
                Collection['BEdgesLayer'+str(k)][EInd][2]=Collection['BEdgesLayer'+str(k)][EInd][2]*AnticipateIncentive;

          for i in Collection['IEdgesLayer'+str(k+1)]:
            coord1= Collection['INodesLayer'+str(k+1)][Collection['IEdgesLayer'+str(k+1)][i][0]];
            coord2= Collection['INodesLayer'+str(k+1)][Collection['IEdgesLayer'+str(k+1)][i][1]];

            if [coord1[0],coord1[1],LayerGap*(k-1)] in Collection['BNodeListLayer'+str(k)] and [coord2[0],coord2[1],LayerGap*(k-1)] in Collection['BNodeListLayer'+str(k)]:
              Ind1 = Collection['BNodeIndicesLayer'+str(k)][coord1[0],coord1[1],LayerGap*(k-1)][0];
              Ind2 = Collection['BNodeIndicesLayer'+str(k)][coord2[0],coord2[1],LayerGap*(k-1)][0];

              if [Ind1,Ind2] in Collection['BEdgeListLayer'+str(k)]:
                EInd = Collection['BEdgeIndicesLayer'+str(k)][Ind1,Ind2][0];
                Collection['BEdgesLayer'+str(k)][EInd][2]=Collection['BEdgesLayer'+str(k)][EInd][2]*AnticipateIncentive;

        if Boundary == 0:
          for i in Collection['IEdgesLayer'+str(k+1)]:
            coord1= Collection['INodesLayer'+str(k+1)][Collection['IEdgesLayer'+str(k+1)][i][0]];
            coord2= Collection['INodesLayer'+str(k+1)][Collection['IEdgesLayer'+str(k+1)][i][1]];

            if [coord1[0],coord1[1],LayerGap*(k-1)] in Collection['INodeListLayer'+str(k)] and [coord2[0],coord2[1],LayerGap*(k-1)] in Collection['INodeListLayer'+str(k)]:
              Ind1 = Collection['INodeIndicesLayer'+str(k)][coord1[0],coord1[1],LayerGap*(k-1)][0];
              Ind2 = Collection['INodeIndicesLayer'+str(k)][coord2[0],coord2[1],LayerGap*(k-1)][0];

              if [Ind1,Ind2] in Collection['IEdgeListLayer'+str(k)]:
                EInd = Collection['IEdgeIndicesLayer'+str(k)][Ind1,Ind2][0];
                Collection['IEdgesLayer'+str(k)][EInd][2]=Collection['IEdgesLayer'+str(k)][EInd][2]*AnticipateIncentive;

          for i in Collection['BEdgesLayer'+str(k+1)]:
            coord1= Collection['BNodesLayer'+str(k+1)][Collection['BEdgesLayer'+str(k+1)][i][0]];
            coord2= Collection['BNodesLayer'+str(k+1)][Collection['BEdgesLayer'+str(k+1)][i][1]];

            if [coord1[0],coord1[1],LayerGap*(k-1)] in Collection['INodeListLayer'+str(k)] and [coord2[0],coord2[1],LayerGap*(k-1)] in Collection['INodeListLayer'+str(k)]:
              Ind1 = Collection['INodeIndicesLayer'+str(k)][coord1[0],coord1[1],LayerGap*(k-1)][0];
              Ind2 = Collection['INodeIndicesLayer'+str(k)][coord2[0],coord2[1],LayerGap*(k-1)][0];

              if [Ind1,Ind2] in Collection['IEdgeListLayer'+str(k)]:
                EInd = Collection['IEdgeIndicesLayer'+str(k)][Ind1,Ind2][0];
                Collection['IEdgesLayer'+str(k)][EInd][2]=Collection['IEdgesLayer'+str(k)][EInd][2]*AnticipateIncentive;

    ####### Here begins the tuning of the cplex parameters timelimit etc..
    prob = cplex.Cplex();
    # prob.parameters.timelimit.set(timelim);
#     prob.parameters.mip.limits.solutions.set(Solim);
    prob.parameters.mip.tolerances.integrality.set(0);
    prob.objective.set_sense(prob.objective.sense.maximize);
    if Boundary == 0:
      # prob.parameters.mip.tolerances.absmipgap.set(2*Gaplim);
      prob.parameters.mip.tolerances.mipgap.set(InteriorGaplim);
    if Boundary == 1:
      if k==1:
        # prob.parameters.mip.tolerances.absmipgap.set(2*Gaplim);
        prob.parameters.mip.tolerances.mipgap.set(BoundaryGap);
      else:
        # prob.parameters.mip.tolerances.absmipgap.set(2*Gaplim);
        prob.parameters.mip.tolerances.mipgap.set(BoundaryGap);



#################################### Here begins the defining of the LP objects, edges, node positions, binaries etc..

    ## Edgedists ##
    edgedists = zeros(NumEdges);
    for i in range(0,NumEdges):  #### The NumEdges should capture all of the edges in either the boundary or interior
      edgedists[i]=Edges[i+1][2];

    ## Node Position ##

    nodes = ones(NumNodes)*NodeWeight;

    ## Binaries for Potential Starting and Ending Nodes

    bnodes = zeros(len(StartNodes)+len(EndNodes));

    ## Problem Objective ##

    my_obj = edgedists;
    prob.variables.add(obj = my_obj, types = [prob.variables.type.binary]*len(my_obj));
    my_obj = nodes;
    prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]*len(my_obj));
    my_obj = bnodes;
    prob.variables.add(obj = my_obj, types = [prob.variables.type.binary]*len(my_obj));

    ####################################################################################
    ####       Main Constraints: Fist Node, Last Node, Ingoing, Outgoing            ####
    ####################################################################################
    #### NodeOuts keeps track of all nodes attached to edges leaving a node         ####
    #### NodeIns keeps track of all nodes attached to edges coming into a  node     ####
    ####################################################################################

    constraintrowindex = 0;

    ### Outgoing Restraint (2a) ####### Can only leave a node once

    my_rhs = ones(NumNodes);
    prob.linear_constraints.add(rhs= my_rhs);

    for i in range(0,NumNodes):
      x = int(constraintrowindex);
      constraintrowindex = constraintrowindex+1
      v=1;
      for j in range(0,len(NodeOuts[i+1])):
        y=int(EdgeIndices[i+1, NodeOuts[i+1][j]][0]-1);
        prob.linear_constraints.set_coefficients(x, y, v);
      prob.linear_constraints.set_senses(x, "L");

    ### Outgoing Restraint (2b) ###### Cannot leave the last node and can only leave a node once

    my_rhs = ones(len(EndNodes));
    prob.linear_constraints.add(rhs= my_rhs);


    count = 0;
    for j in EndNodes:
      x= constraintrowindex;
      y = int(NumEdges+NumNodes+len(StartNodes)+count);
      count = count+1;
      prob.linear_constraints.set_coefficients(x, y, 1);
      prob.linear_constraints.set_senses(x, "L");
      for i in range(0,len(NodeOuts[j])):
        y=int(EdgeIndices[j, NodeOuts[j][i]][0]-1);
        prob.linear_constraints.set_coefficients(x, y, 1);
      constraintrowindex = constraintrowindex+1;

    ### Incoming Restraint (incoming <= outgoing)  (3) ###

    my_rhs = zeros(NumNodes);
    prob.linear_constraints.add(rhs= my_rhs);
    count = 0;
    for j in range(1,NumNodes+1):
      if j not in EndNodes:
        x = int(constraintrowindex);
        constraintrowindex = constraintrowindex+1;
        v=1;
        prob.linear_constraints.set_senses(x, "L");
        for i in range(0,len(NodeIns[j])):
          y=int(EdgeIndices[NodeIns[j][i],j][0]-1);
          z=int(EdgeIndices[j,NodeOuts[j][i]][0]-1);
          prob.linear_constraints.set_coefficients(x, y, v);
          prob.linear_constraints.set_coefficients(x, z, -v);
      else:
        y = int(NumEdges+NumNodes+len(StartNodes)+count)
        count = count +1;
        x = int(constraintrowindex);
        constraintrowindex = constraintrowindex+1;
        prob.linear_constraints.set_coefficients(x, y, -1);
        v=1;
        prob.linear_constraints.set_senses(x, "L");
        for i in range(0,len(NodeIns[j])):
          y=int(EdgeIndices[NodeIns[j][i],j][0]-1);
          z=int(EdgeIndices[j,NodeOuts[j][i]][0]-1);
          prob.linear_constraints.set_coefficients(x, y, v);
          prob.linear_constraints.set_coefficients(x, z, -v);

    ##### Starting Node Restraint (4a) ####### Must leave first node ## sum_j(x_sj)-b_s => 0

    my_rhs = zeros(len(StartNodes));
    prob.linear_constraints.add(rhs= my_rhs);

    x= constraintrowindex;
    v=1;
    count = 0;

    for j in StartNodes:
      x= constraintrowindex;
      y = int(NumEdges+NumNodes+count);
      count = count+1;
      prob.linear_constraints.set_coefficients(x, y, -1);
      prob.linear_constraints.set_senses(x, "G");
      for i in range(0,len(NodeOuts[j])):
        y=int(EdgeIndices[j, NodeOuts[j][i]][0]-1);
        prob.linear_constraints.set_coefficients(x, y, 1);

      constraintrowindex = constraintrowindex+1;

    ##### Starting Node Restraint (4b) ##### Cannot enter first node ## sum_j(x_js)+b_s <= 1

    my_rhs = ones(len(StartNodes));
    prob.linear_constraints.add(rhs= my_rhs);

    x= constraintrowindex;
    count = 0;
    for j in StartNodes:
      x= constraintrowindex;
      y = int(NumEdges+NumNodes+count);
      count = count+1;
      prob.linear_constraints.set_coefficients(x, y, 1);
      prob.linear_constraints.set_senses(x, "L");
      for i in range(0,len(NodeIns[j])):
        y=int(EdgeIndices[NodeIns[j][i],j][0]-1);
        prob.linear_constraints.set_coefficients(x, y, 1);

      constraintrowindex = constraintrowindex+1;

    ##### Ending Node Restraint (5) ###### Must enter last node

    my_rhs = zeros(len(EndNodes));
    prob.linear_constraints.add(rhs= my_rhs);
    count = 0;

    for j in EndNodes:
      x= constraintrowindex;
      y = int(NumEdges+NumNodes+len(StartNodes)+count);
      count = count+1;
      prob.linear_constraints.set_coefficients(x, y, -1);
      prob.linear_constraints.set_senses(x, "G");
      for i in range(0,len(NodeIns[j])):
        y=int(EdgeIndices[NodeIns[j][i],j][0]-1);
        prob.linear_constraints.set_coefficients(x, y, 1);

      constraintrowindex = constraintrowindex+1;

    ###### Unique Start and End Restraint (6a)  ## There can be only one first node ##

    my_rhs = ones(1);
    prob.linear_constraints.add(rhs= my_rhs);

    for i in range(0,len(StartNodes)):
      y = NumEdges+NumNodes+i;
      prob.linear_constraints.set_coefficients(constraintrowindex, y, 1);

    prob.linear_constraints.set_senses(constraintrowindex, "L");
    constraintrowindex = constraintrowindex+1;

    ###### Unique Start and End Restraint (6b)  #### There can be only one last node ##

    my_rhs = ones(1);
    prob.linear_constraints.add(rhs= my_rhs);

    for i in range(0,len(EndNodes)):
      y = NumEdges+NumNodes+len(StartNodes)+i;
      prob.linear_constraints.set_coefficients(constraintrowindex, y, 1);

    prob.linear_constraints.set_senses(constraintrowindex, "L");
    constraintrowindex = constraintrowindex+1;

    ## Subtour Constraint (7) ##

    for i in range(0,NumNodes):
      Ui = int(NumEdges+i);
      for j in range(0,len(NodeOuts[i+1])):
        x = int(constraintrowindex);
        constraintrowindex = constraintrowindex+1;
        my_rhs = [NumNodes-2];
        prob.linear_constraints.add(rhs= my_rhs);
        Uk = int(NumEdges+NodeOuts[i+1][j]-1);
        xij = int(EdgeIndices[i+1,NodeOuts[i+1][j]][0]-1);
        prob.linear_constraints.set_coefficients(x, Ui, 1);
        prob.linear_constraints.set_coefficients(x, Uk, -1);
        prob.linear_constraints.set_coefficients(x, xij, NumNodes-1);
        prob.linear_constraints.set_senses(x, "L");


    #### Node Position Restraint (8a) #### U_j's at least 1 if used #### Sum of outgoing xij - Ui <= 0 ##

    my_rhs = zeros(NumNodes);
    prob.linear_constraints.add(rhs= my_rhs);

    for j in range(0,NumNodes):
      x = int(constraintrowindex);
      constraintrowindex = constraintrowindex+1
      v=1;
      z=NumEdges+j;
      prob.linear_constraints.set_coefficients(x, z, -v);
      for i in range(0,len(NodeOuts[j+1])):
        y=int(EdgeIndices[j+1,NodeOuts[j+1][i]][0]-1);
        prob.linear_constraints.set_coefficients(x, y, v);
      prob.linear_constraints.set_senses(x, "L");


    #### Node Position Restraint (8b) ###### Consec U_j's are 1 apart if edge connects them #### In conjunction with Subtour add U_j - U_i + (NumNodes-1)x_ij<= (NumNodes-1) + 1 ##

    for i in range(0,NumNodes):
      Ui = int(NumEdges+i);
      for j in range(0,len(NodeOuts[i+1])):
        x = int(constraintrowindex);
        constraintrowindex = constraintrowindex+1;
        my_rhs = [NumNodes];
        prob.linear_constraints.add(rhs= my_rhs);
        Uk = int(NumEdges+NodeOuts[i+1][j]-1);
        xij = int(EdgeIndices[i+1,NodeOuts[i+1][j]][0]-1);
        prob.linear_constraints.set_coefficients(x, Ui, -1);
        prob.linear_constraints.set_coefficients(x, Uk, 1);
        prob.linear_constraints.set_coefficients(x, xij, NumNodes-1);
        prob.linear_constraints.set_senses(x, "L");


    #### Starting Node position Restraint (9) #### U_s + NumNodes*b_s <= NumNodes

    my_rhs = ones(len(StartNodes))*(NumNodes+1);
    prob.linear_constraints.add(rhs= my_rhs)
    count = 0;

    for j in StartNodes:
      count = count+1;
      x= constraintrowindex;
      y= NumEdges+j-1;
      z= NumEdges+NumNodes+count-1;
      prob.linear_constraints.set_coefficients(x, y, 1);
      prob.linear_constraints.set_coefficients(x, z, NumNodes);
      prob.linear_constraints.set_senses(x, "L");
      constraintrowindex = constraintrowindex+1;

    ##### Thickness Constraint (10) ########
############## Write a check at the end of the program that check if this is satisfied ########
    Beadpairs = dict();
    bpairs = 0;
    if k >=0:
      for j in range(0,NumNodes):
        MarkedEdges = [];
        # my_rhs = ones(1);
        # prob.linear_constraints.add(rhs= my_rhs);
        # x = int(constraintrowindex);
        # prob.linear_constraints.set_senses(x,"G");
        # constraintrowindex = constraintrowindex+1;
        for l in range(0,len(NodeOuts[j+1])):
          y=int(EdgeIndices[j+1,NodeOuts[j+1][l]][0]-1);
        #   prob.linear_constraints.set_coefficients(x, y, 1);
          MarkedEdges.append(y);
          y=int(EdgeIndices[NodeOuts[j+1][l],j+1][0]-1);
        #   prob.linear_constraints.set_coefficients(x, y, 1);
          MarkedEdges.append(y);
        for i in range(0,NumNodes):
          if scipy.spatial.distance.pdist([array(Nodes[i+1]),array(Nodes[j+1])])[0] <= VtxRad:
            for l in range(0,len(NodeOuts[i+1])):
              y1=int(EdgeIndices[i+1, NodeOuts[i+1][l]][0]-1);
              if y1 not in MarkedEdges:
                # prob.linear_constraints.set_coefficients(x, y1, 1);
              y2=int(EdgeIndices[NodeOuts[i+1][l],i+1][0]-1);
              if y2 not in MarkedEdges:
                # prob.linear_constraints.set_coefficients(x, y2, 1);
                bpairs += 1;
                Beadpairs[bpairs] = [j,i,y1,y2];

    ###### Cannot touch a node within a Beadswidth of a used node unless in  consecutive manner sum_k(xik)+sum_k(xjk)-xij-xji <= 1

    for i in range(1,bpairs+1):
        my_rhs = ones(1);
        prob.linear_constraints.add(rhs= my_rhs);
        x = int(constraintrowindex);
        prob.linear_constraints.set_senses(x,"L");
        constraintrowindex = constraintrowindex+1;
        for l in range(0,len(NodeOuts[Beadpairs[i][0]+1])):
            if NodeOuts[Beadpairs[i][0]+1][l] != Beadpairs[i][1]:
                y=int(EdgeIndices[Beadpairs[i][0]+1, NodeOuts[Beadpairs[i][0]+1][l]][0]-1);
                prob.linear_constraints.set_coefficients(x, y, 1);
        for l in range(0,len(NodeOuts[Beadpairs[i][1]+1])):
            if NodeOuts[Beadpairs[i][1]+1][l] != Beadpairs[i][0]:
                y=int(EdgeIndices[Beadpairs[i][1]+1, NodeOuts[Beadpairs[i][1]+1][l]][0]-1);
                prob.linear_constraints.set_coefficients(x, y, 1);



#####################################################################################################################################
####      Additional Constraints:  Optional constraints, Band, Curvature, Boundary Edge Print, Duplicate Edge Restraint          ####
#####################################################################################################################################


    ######## Must print all original boundary edges (14) ##########
    if Boundary == 1:
      my_rhs = ones(Collection["BONumEdgesLayer"+str(k)]/2.0);  ##### Number of constraints will be half of the number of original boundary edges as there is 1 constraint per 4 edges (2 pairs of duplicates)
      prob.linear_constraints.add(rhs= my_rhs);

      for j in range(0,int(Collection["BONumEdgesLayer"+str(k)]/2.0)):  ##### Set the constraints for each variable
        x= constraintrowindex;
        y1 = int(2*j);
        y2 = int(2*j+1);
        y3 = int(NumEdges/2.0+2*j);
        y4 = int(NumEdges/2.0+2*j+1);

        prob.linear_constraints.set_coefficients(x, y1, 1);
        prob.linear_constraints.set_coefficients(x, y2, 1);
        prob.linear_constraints.set_coefficients(x, y3, 1);
        prob.linear_constraints.set_coefficients(x, y4, 1);

        prob.linear_constraints.set_senses(x, "E");

        constraintrowindex = constraintrowindex+1;


    ##### Annular Constraints (15) ########
    # if Boundary == 0:
    #   center = [9.0,9.0,0.0];
    #   if k >=0:
    #     # for r in [2.5]:
    #     for r in [2, 4, 6, 8, 10]:
    #       my_rhs = ones(1);
    #       prob.linear_constraints.add(rhs= my_rhs);
    #
    #       for j in range(0,NumEdges):
    #         Nodea = int(Edges[j+1][0]);
    #         Nodeb = int(Edges[j+1][1]);
    #
    #         Nodea = Nodes[Nodea][0],Nodes[Nodea][1],0.0
    #         Nodeb = Nodes[Nodeb][0],Nodes[Nodeb][1],0.0
    #
    #         if scipy.spatial.distance.pdist([Nodea,array(center)])[0] > r and scipy.spatial.distance.pdist([Nodeb,array(center)])[0] < r:
    #           prob.linear_constraints.set_coefficients(constraintrowindex, j, 1);
    #
    #         if scipy.spatial.distance.pdist([Nodea,array(center)])[0] < r and scipy.spatial.distance.pdist([Nodeb,array(center)])[0] > r:
    #           prob.linear_constraints.set_coefficients(constraintrowindex, j, 1);
    #
    #       prob.linear_constraints.set_senses(constraintrowindex, "L");
    #       constraintrowindex = constraintrowindex+1;

    # if Boundary == 0:
    #   center = 2*[4.5,24.5,0.0];
    #   if k >=0:
    #     # for r in [2.5]:
    #     for r in [2, 3, 4, 5, 6, 7, 8]:
    #       my_rhs = ones(1);
    #       prob.linear_constraints.add(rhs= my_rhs);
    #
    #       for j in range(0,NumEdges):
    #         Nodea = int(Edges[j+1][0]);
    #         Nodeb = int(Edges[j+1][1]);
    #
    #         Nodea = Nodes[Nodea][0],Nodes[Nodea][1],0.0
    #         Nodeb = Nodes[Nodeb][0],Nodes[Nodeb][1],0.0
    #
    #         if scipy.spatial.distance.pdist([Nodea,array(center)])[0] > r and scipy.spatial.distance.pdist([Nodeb,array(center)])[0] < r:
    #           prob.linear_constraints.set_coefficients(constraintrowindex, j, 1);
    #
    #         if scipy.spatial.distance.pdist([Nodea,array(center)])[0] < r and scipy.spatial.distance.pdist([Nodeb,array(center)])[0] > r:
    #           prob.linear_constraints.set_coefficients(constraintrowindex, j, 1);
    #
    #       prob.linear_constraints.set_senses(constraintrowindex, "L");
    #       constraintrowindex = constraintrowindex+1;

    ##### Boxed Annular Constraints (15a) ########
    #
    # if Boundary == 0:
    #   center = [9.0,9.0,0.0];
    #   if k >=0:
    #     # for r in [2.5]:
    #     for r in [2, 4, 6, 8, 10]:
    #       my_rhs = ones(1);
    #       prob.linear_constraints.add(rhs= my_rhs);
    #
    #       for j in range(0,NumEdges):
    #         Nodea = int(Edges[j+1][0]);
    #         Nodeb = int(Edges[j+1][1]);
    #
    #         Nodea = Nodes[Nodea][0],Nodes[Nodea][1],0.0
    #         Nodeb = Nodes[Nodeb][0],Nodes[Nodeb][1],0.0
    #
    #         if abs(Nodeb[0]-center[0]) < r and abs(Nodeb[1]-center[1]) < r:
    #             if abs(Nodea[0]-center[0]) > r or abs(Nodea[1]-center[1]) > r:
    #                 prob.linear_constraints.set_coefficients(constraintrowindex, j, 1);
    #
    #
    #         elif abs(Nodea[0]-center[0]) < r and abs(Nodea[1]-center[1]) < r:
    #             if  abs(Nodeb[0]-center[0]) > r or abs(Nodeb[1]-center[1]) > r:
    #                 prob.linear_constraints.set_coefficients(constraintrowindex, j, 1);
    #
    #       prob.linear_constraints.set_senses(constraintrowindex, "L");
    #       constraintrowindex = constraintrowindex+1;
    #
    #
    # if Boundary == 0:
    #   center = [9.0,49.0,0.0];
    #   if k >=0:
    #     # for r in [2.5]:
    #     for r in [2, 4, 6, 8, 10]:
    #       my_rhs = ones(1);
    #       prob.linear_constraints.add(rhs= my_rhs);
    #
    #       for j in range(0,NumEdges):
    #         Nodea = int(Edges[j+1][0]);
    #         Nodeb = int(Edges[j+1][1]);
    #
    #         Nodea = Nodes[Nodea][0],Nodes[Nodea][1],0.0
    #         Nodeb = Nodes[Nodeb][0],Nodes[Nodeb][1],0.0
    #
    #         if abs(Nodeb[0]-center[0]) < r and abs(Nodeb[1]-center[1]) < r:
    #             if abs(Nodea[0]-center[0]) > r or abs(Nodea[1]-center[1]) > r:
    #                 prob.linear_constraints.set_coefficients(constraintrowindex, j, 1);
    #
    #
    #         elif abs(Nodea[0]-center[0]) < r and abs(Nodea[1]-center[1]) < r:
    #             if  abs(Nodeb[0]-center[0]) > r or abs(Nodeb[1]-center[1]) > r:
    #                 prob.linear_constraints.set_coefficients(constraintrowindex, j, 1);
    #
    #       prob.linear_constraints.set_senses(constraintrowindex, "L");
    #       constraintrowindex = constraintrowindex+1;
    # #
    #



        ##### Band Constraints (16) ########

    # if Boundary == 0:   ### Horizontal Bands
    #   if k >=0:
    #     for r in range(1,60,2):
    #       my_rhs = ones(1);
    #       prob.linear_constraints.add(rhs= my_rhs);
    #
    #       for j in range(0,NumEdges):
    #         Nodea = int(Edges[j+1][0]);
    #         Nodeb = int(Edges[j+1][1]);
    #
    #         Nodeay = Nodes[Nodea][1];
    #         Nodeby = Nodes[Nodeb][1];
    #
    #         if Nodeay > r and Nodeby < r:
    #           prob.linear_constraints.set_coefficients(constraintrowindex, j, 1);
    #
    #         if Nodeay < r and Nodeby > r:
    #           prob.linear_constraints.set_coefficients(constraintrowindex, j, 1);
    #
    #       prob.linear_constraints.set_senses(constraintrowindex, "L");
    #       constraintrowindex = constraintrowindex+1;

    if Boundary == 0:   ### Slant Bands
      if k >=0:
        for r in range(1,60,2):
          my_rhs = 2*ones(1);
          prob.linear_constraints.add(rhs= my_rhs);

          for j in range(0,NumEdges):
            Nodea = Nodes[int(Edges[j+1][0])];
            Nodeb = Nodes[int(Edges[j+1][1])];

            if Nodea[0]+Nodea[1] < r and Nodeb[0]+Nodeb[1] > r:
              prob.linear_constraints.set_coefficients(constraintrowindex, j, 1);

            if Nodea[0]+Nodea[1] > r and Nodeb[0]+Nodeb[1] < r:
              prob.linear_constraints.set_coefficients(constraintrowindex, j, 1);

          prob.linear_constraints.set_senses(constraintrowindex, "L");
          constraintrowindex = constraintrowindex+1;

    ###########################################################################################
    ####           Solve, calculate edgepath: print out solution           ####
    ###########################################################################################

    t0=time.clock();

    prob.solve();

    CPlexTime = time.clock()-t0;

    z = prob.solution.get_values();

    print z;

    #  EdgePath  ##
    if Boundary == 1:
      BEdgePath=[];
      for i in range(0,NumEdges):
        if i < NumEdges/2.0:
          BEdgePath.append(int(z[i]));
        if i >= NumEdges/2.0:
          BEdgePath[i-NumEdges/2]=BEdgePath[i-NumEdges/2]+int(z[i]);
    else:
      IEdgePath = [];
      for i in range(0,NumEdges):
        if i < NumEdges/2.0:
          IEdgePath.append(int(z[i]));
        if i >= NumEdges/2.0:
          IEdgePath[i-NumEdges/2]=IEdgePath[i-NumEdges/2]+int(z[i]);

    if Boundary == 1:
      EdgePath = BEdgePath;
    else:
      EdgePath = IEdgePath;

    print " "

    print "EdgePath"
    print EdgePath;
    tt = time.clock()-tti;

    ###################################################################################################################
    ########              Write Path Solution to Medit readable file and print solution            ###################
    ###################################################################################################################

    f = file('Layer'+str(k)+str(Boundary)+'Orig.mesh','w')
    sys.stdout = f

    print "MeshVersionFormatted 1";

    print "Dimension";
    print "3"
    print "Vertices"
    print NumNodes;

    # Set of mesh vertices
    for i in range(1,NumNodes+1):
      print Nodes[i][0], Nodes[i][1], Nodes[i][2],1;


    print "# Set of Edges";
    print "Edges";
    print NumEdges;

    for i in range(0,NumEdges/2):
      if Boundary == 1:
        if i+1 <= Collection["BONumEdgesLayer"+str(k)]:
          print Edges[i+1][0], Edges[i+1][1], 1;
        else:
          print Edges[i+1][0], Edges[i+1][1], 7;
      else:
        if i+1 <= Collection["IONumEdgesLayer"+str(k)]:
          print Edges[i+1][0], Edges[i+1][1], 1;
        else:
          print Edges[i+1][0], Edges[i+1][1], 7;


    # print 0;
    sys.stdout = orig_stdout
    f.close()


    NumSolEdges = 0;
    UsedNodes = [];
    print "Solution Edge list is ";
    for i in range(0,len(EdgePath)):
      if EdgePath[i] != 0:
        print Edges[i+1][0], Edges[i+1][1], "Edge",
        NumSolEdges = NumSolEdges +1;
        if Edges[i+1][0] not in UsedNodes:
          UsedNodes.append(Edges[i+1][0]);
        if Edges[i+1][1] not in UsedNodes:
          UsedNodes.append(Edges[i+1][1]);


    Uj=zeros(NumNodes);
    print "Used Node Uj's are ";
    for i in range(0,NumNodes):
      Uj[i]=z[NumEdges+i];

    print Uj;
    print "Ordered Used Ujs is", sorted(Uj);
    print "Number of Edges in Solution Path is", NumSolEdges;
    print "Maximum Number of Edges possible is", MaxEdges/2.0;

    NFractional = 0;

    for i in range(0,len(z)):
      if z[i]-round(z[i]) != 0:
        NFractional = NFractional+1;

    if NFractional >= 1:
      print "Fractional. Number of Fractional entries is", NFractional;
    else:
      print "Integral";

    print "CPlex Solver Processing Time is ", CPlexTime;

    print "Total Processing Time is ", tt;

    PathNodesList = [];
    IntNodeList = [];
    PathEdgeConnection = {}; ############### stores next node in path and whether to print the edge
    PESNodes = [];

    for i in range(0,NumEdges):
      if z[i] == 1:

        if Edges[i+1][2] == -1:
          if Nodes[Edges[i+1][0]][0] == Nodes[Edges[i+1][1]][0] and Nodes[Edges[i+1][0]][1] == Nodes[Edges[i+1][1]][1]:  #### Edges are directed thus first node of edge Edges[i+1][0] always comes first, blue edges are identified by -1 and if vertical do not print or move along edge 0
            PathEdgeConnection[Edges[i+1][0]] = [Edges[i+1][1],0,0];  #### [NextNodeOnPath, Don'tPrint, Don'tMove]
          else:
            PathEdgeConnection[Edges[i+1][0]] = [Edges[i+1][1],0,1];  #### [NextNodeOnPath, Don'tPrint, Move]

        else:
          PathEdgeConnection[Edges[i+1][0]] = [Edges[i+1][1],1,1];   #### [NextNodeOnPath, Print, Move]


        if Edges[i+1][0] not in PathNodesList: #### Add nodes on solution path to corresponding list
          PathNodesList.append(Edges[i+1][0]);
        else:
          IntNodeList.append(Edges[i+1][0])   #### Every node touched twice will not be a starting or ending node and hence in the interior of the edge path

        if Edges[i+1][1] not in PathNodesList:
          PathNodesList.append(Edges[i+1][1]);
        else:
          IntNodeList.append(Edges[i+1][1])

    for i in PathNodesList:  #####  Find Starting and Ending Nodes and add to list
      if i not in IntNodeList:
        PESNodes.append(i);

    print len(z), NumEdges+NumNodes, NumEdges+PESNodes[0]-1;
    if z[NumEdges+PESNodes[0]-1] < z[NumEdges+PESNodes[1]-1]:  ### determin which node is the starting node
      Startnode = PESNodes[0];
      EndNode = PESNodes[1];
    else:
      Startnode = PESNodes[1];
      EndNode = PESNodes[0];


    LastNode = Nodes[EndNode];   ##### Set Lastnode and currentnode to the starting node
    currentnode = Startnode;

    if Boundary == 1:
      PBStartNode = Nodes[Startnode];

    else:
      PIStartNode = Nodes[Startnode];

    ## Write out ##

    if k >= 1:
      f = file('Layer'+str(k)+str(Boundary)+'Sol.mesh','w')
      sys.stdout = f

      print "MeshVersionFormatted 1";

      print "Dimension";
      print "3"
      print "Vertices"
      print NumNodes;

      # Set of mesh vertices
      for i in range(1,NumNodes+1):
        print Nodes[i][0], Nodes[i][1], Nodes[i][2],1;


      print "# Set of Edges";
      print "Edges";
      print count_nonzero(EdgePath);

      for i in range(0,NumEdges/2):
        if EdgePath[i] != 0:
          if i+1 in NoPrintEdgeInds:
            print Edges[i+1][0], Edges[i+1][1], 7;
          else:
            print Edges[i+1][0], Edges[i+1][1], 1;

      # print 0;
      sys.stdout = orig_stdout
      f.close()



      f = file('Layer'+str(k)+str(Boundary)+'Orig.mesh','w')
      sys.stdout = f

      print "MeshVersionFormatted 1";

      print "Dimension";
      print "3"
      print "Vertices"
      print NumNodes;

      # Set of mesh vertices
      for i in range(1,NumNodes+1):
        print Nodes[i][0], Nodes[i][1], Nodes[i][2],1;


      print "# Set of Edges";
      print "Edges";
      print NumEdges;

      for i in range(0,NumEdges/2):
        if Boundary == 1:
          if i+1 <= Collection["BONumEdgesLayer"+str(k)]:
            print Edges[i+1][0], Edges[i+1][1], 1;
          else:
            print Edges[i+1][0], Edges[i+1][1], 7;
        else:
          if i+1 <= Collection["IONumEdgesLayer"+str(k)]:
            print Edges[i+1][0], Edges[i+1][1], 1;
          else:
            print Edges[i+1][0], Edges[i+1][1], 7;


      # print 0;
      sys.stdout = orig_stdout
      f.close()

      ########################################################################################################################
  ####### Vertex Path Out, should be modified for printer #######


      if k == 2 and Boundary == 1:
        f = file('VertexPath.txt','w')
        sys.stdout = f
        print "# x-coord, y-coord, z-coord, print when leaving node yes(1)/no(0)"

        for i in range(0,len(PathNodesList)-1):
          if PathEdgeConnection[currentnode][2] == 1:
            print Nodes[currentnode][0], Nodes[currentnode][1], float((k-2)*2), PathEdgeConnection[currentnode][1];
            currentnode = PathEdgeConnection[currentnode][0];
          else:
            currentnode = PathEdgeConnection[currentnode][0];

        print LastNode[0], LastNode[1], float((k-2)*2), 0;

        sys.stdout = orig_stdout
        f.close()

      else:
        f = file('VertexPath.txt','a')
        sys.stdout = f

        for i in range(0,len(PathNodesList)-1):
          if PathEdgeConnection[currentnode][2] == 1:
            print Nodes[currentnode][0], Nodes[currentnode][1], float((k-2)*2), PathEdgeConnection[currentnode][1];
            currentnode = PathEdgeConnection[currentnode][0];
          else:
            currentnode = PathEdgeConnection[currentnode][0];

        print LastNode[0], LastNode[1], float((k-2)*2), 0;

        sys.stdout = orig_stdout
        f.close()

#     m = count_nonzero(Uj);
#     ListUj = list(Uj);
#     for i in range(0,m):
#       print Nodes[ListUj.index(i+1)+1][0], Nodes[ListUj.index(i+1)+1][1], Nodes[ListUj.index(i+1)+1][2];



########################################################################################################################
####### Vertex Path Out, should be modified for printer #######

#     f = file('VertexPath.txt','a')
#     sys.stdout = f
#     m = count_nonzero(Uj);
#     ListUj = list(Uj);
#     for i in range(0,m):
#       print Nodes[ListUj.index(i+1)+1][0], Nodes[ListUj.index(i+1)+1][1], Nodes[ListUj.index(i+1)+1][2];

#     sys.stdout = orig_stdout
#     f.close()



########################################################################################################################
######################## Currently Unused and Extras ########################################

    ##############################################################
    ##############   90 deg turns only ###########################
    ##############################################################



#     for j in Edges:
#       StartNode = Nodes[Edges[j][0]];
#       EndNode = Nodes[Edges[j][1]];
#       Diff = substract_lists(StartNode,EndNode);

#       if Diff[0] != 0 and Diff[1] != 0 and Edges[j] not in csegments:
#         prob.linear_constraints.add(rhs= zeros(1));
#         prob.linear_constraints.set_coefficients(constraintrowindex, j-1, 1);
#         prob.linear_constraints.set_senses(constraintrowindex, "L");
#         constraintrowindex = constraintrowindex+1;

  ######################### Curvature Restraints #################

  # e = 1;
  # ETC.append(e);
  # CE.append(e);

  # StartCurvRes = constraintrowindex;
  # var = 0;
  # for i in NTEL[Edges[e][0]-1]:
  #   if i not in CE:
  #     PairEdgeCount = PairEdgeCount+1;
  #     if i not in ETC:
  #       ETC.append(i);
  #     if Edges[i][0] == Edges[e][0]:
  #       if CosAngleBtwConEdges(Nodes[Edges[i][1]],Nodes[Edges[e][0]],Nodes[Edges[e][1]],Edges[i][2],Edges[e][2]) > CosAngLim:
  #         prob.linear_constraints.add(rhs= ones(1));
  #         prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #         prob.linear_constraints.set_senses(constraintrowindex, "L");
  #         constraintrowindex = constraintrowindex+1;
  #       elif CosAngleBtwConEdges(Nodes[Edges[i][1]],Nodes[Edges[e][0]],Nodes[Edges[e][1]],Edges[i][2],Edges[e][2]) > FortyLim:
  #         my_obj = ones(1)*(-10.0);
  #         prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #         my_rhs = ones(1);
  #         prob.linear_constraints.add(rhs= my_rhs);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #         prob.linear_constraints.set_senses(constraintrowindex, "L");
  #         constraintrowindex = constraintrowindex+1;
  #         var = var+1;
  #       elif CosAngleBtwConEdges(Nodes[Edges[i][1]],Nodes[Edges[e][0]],Nodes[Edges[e][1]],Edges[i][2],Edges[e][2]) > FiftyLim:
  #         my_obj = ones(1)*(-5.0);
  #         prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #         my_rhs = ones(1);
  #         prob.linear_constraints.add(rhs= my_rhs);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #         prob.linear_constraints.set_senses(constraintrowindex, "L");
  #         constraintrowindex = constraintrowindex+1;
  #         var = var+1;

  #     else:
  #       if CosAngleBtwConEdges(Nodes[Edges[i][0]],Nodes[Edges[e][0]],Nodes[Edges[e][1]],Edges[i][2],Edges[e][2]) > CosAngLim:
  #         prob.linear_constraints.add(rhs= ones(1));
  #         prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #         prob.linear_constraints.set_senses(constraintrowindex, "L");
  #         constraintrowindex = constraintrowindex+1;
  #       elif CosAngleBtwConEdges(Nodes[Edges[i][0]],Nodes[Edges[e][0]],Nodes[Edges[e][1]],Edges[i][2],Edges[e][2]) > FortyLim:
  #         my_obj = ones(1)*(-10.0);
  #         prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #         my_rhs = ones(1);
  #         prob.linear_constraints.add(rhs= my_rhs);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #         prob.linear_constraints.set_senses(constraintrowindex, "L");
  #         constraintrowindex = constraintrowindex+1;
  #         var = var+1;
  #       elif CosAngleBtwConEdges(Nodes[Edges[i][0]],Nodes[Edges[e][0]],Nodes[Edges[e][1]],Edges[i][2],Edges[e][2]) > FiftyLim:
  #         my_obj = ones(1)*(-5.0);
  #         prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #         my_rhs = ones(1);
  #         prob.linear_constraints.add(rhs= my_rhs);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #         prob.linear_constraints.set_senses(constraintrowindex, "L");
  #         constraintrowindex = constraintrowindex+1;
  #         var = var+1;

  # for i in NTEL[Edges[e][1]-1]:
  #   if i not in CE:
  #     PairEdgeCount = PairEdgeCount+1;
  #     if i not in ETC:
  #       ETC.append(i);
  #     if Edges[i][0] == Edges[e][1]:
  #       if CosAngleBtwConEdges(Nodes[Edges[i][1]],Nodes[Edges[e][1]],Nodes[Edges[e][0]],Edges[i][2],Edges[e][2]) > CosAngLim:
  #         prob.linear_constraints.add(rhs= ones(1));
  #         prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #         prob.linear_constraints.set_senses(constraintrowindex, "L");
  #         constraintrowindex = constraintrowindex+1;
  #       elif CosAngleBtwConEdges(Nodes[Edges[i][1]],Nodes[Edges[e][1]],Nodes[Edges[e][0]],Edges[i][2],Edges[e][2])  > FortyLim:
  #         my_obj = ones(1)*(-10.0);
  #         prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #         my_rhs = ones(1);
  #         prob.linear_constraints.add(rhs= my_rhs);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #         prob.linear_constraints.set_senses(constraintrowindex, "L");
  #         constraintrowindex = constraintrowindex+1;
  #         var = var+1;
  #       elif CosAngleBtwConEdges(Nodes[Edges[i][1]],Nodes[Edges[e][1]],Nodes[Edges[e][0]],Edges[i][2],Edges[e][2])  > FiftyLim:
  #         my_obj = ones(1)*(-5.0);
  #         prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #         my_rhs = ones(1);
  #         prob.linear_constraints.add(rhs= my_rhs);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #         prob.linear_constraints.set_senses(constraintrowindex, "L");
  #         constraintrowindex = constraintrowindex+1;
  #         var = var+1;
  #     else:
  #       if CosAngleBtwConEdges(Nodes[Edges[i][0]],Nodes[Edges[e][1]],Nodes[Edges[e][0]],Edges[i][2],Edges[e][2]) > CosAngLim:
  #         prob.linear_constraints.add(rhs= ones(1));
  #         prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #         prob.linear_constraints.set_senses(constraintrowindex, "L");
  #         constraintrowindex = constraintrowindex+1;
  #       elif CosAngleBtwConEdges(Nodes[Edges[i][0]],Nodes[Edges[e][1]],Nodes[Edges[e][0]],Edges[i][2],Edges[e][2]) > FortyLim:
  #         my_obj = ones(1)*(-10.0);
  #         prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #         my_rhs = ones(1);
  #         prob.linear_constraints.add(rhs= my_rhs);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #         prob.linear_constraints.set_senses(constraintrowindex, "L");
  #         constraintrowindex = constraintrowindex+1;
  #         var = var+1;
  #       elif CosAngleBtwConEdges(Nodes[Edges[i][0]],Nodes[Edges[e][1]],Nodes[Edges[e][0]],Edges[i][2],Edges[e][2]) > FiftyLim:
  #         my_obj = ones(1)*(-5.0);
  #         prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #         my_rhs = ones(1);
  #         prob.linear_constraints.add(rhs= my_rhs);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #         prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #         prob.linear_constraints.set_senses(constraintrowindex, "L");
  #         constraintrowindex = constraintrowindex+1;
  #         var = var+1;

  # for i in range(1,NumEdges):
  #   e=ETC[i];
  #   CE.append(e);
  #   for i in NTEL[Edges[e][0]-1]:
  #     PairEdgeCount = PairEdgeCount+1;
  #     if i not in CE:
  #       if i not in ETC:
  #         ETC.append(i);
  #       if Edges[i][0] == Edges[e][0]:
  #         if CosAngleBtwConEdges(Nodes[Edges[i][1]],Nodes[Edges[e][0]],Nodes[Edges[e][1]],Edges[i][2],Edges[e][2]) > CosAngLim:
  #           prob.linear_constraints.add(rhs= ones(1));
  #           prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #           prob.linear_constraints.set_senses(constraintrowindex, "L");
  #           constraintrowindex = constraintrowindex+1;
  #         elif CosAngleBtwConEdges(Nodes[Edges[i][1]],Nodes[Edges[e][0]],Nodes[Edges[e][1]],Edges[i][2],Edges[e][2]) > FortyLim:
  #           my_obj = ones(1)*(-10.0);
  #           prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #           my_rhs = ones(1);
  #           prob.linear_constraints.add(rhs= my_rhs);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #           prob.linear_constraints.set_senses(constraintrowindex, "L");
  #           constraintrowindex = constraintrowindex+1;
  #           var = var+1;
  #         elif CosAngleBtwConEdges(Nodes[Edges[i][1]],Nodes[Edges[e][0]],Nodes[Edges[e][1]],Edges[i][2],Edges[e][2]) > FiftyLim:
  #           my_obj = ones(1)*(-5.0);
  #           prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #           my_rhs = ones(1);
  #           prob.linear_constraints.add(rhs= my_rhs);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #           prob.linear_constraints.set_senses(constraintrowindex, "L");
  #           constraintrowindex = constraintrowindex+1;
  #           var = var+1;

  #       else:
  #         if CosAngleBtwConEdges(Nodes[Edges[i][0]],Nodes[Edges[e][0]],Nodes[Edges[e][1]],Edges[i][2],Edges[e][2]) > CosAngLim:
  #           prob.linear_constraints.add(rhs= ones(1));
  #           prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #           prob.linear_constraints.set_senses(constraintrowindex, "L");
  #           constraintrowindex = constraintrowindex+1;
  #         elif CosAngleBtwConEdges(Nodes[Edges[i][0]],Nodes[Edges[e][0]],Nodes[Edges[e][1]],Edges[i][2],Edges[e][2]) > FortyLim:
  #           my_obj = ones(1)*(-10.0);
  #           prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #           my_rhs = ones(1);
  #           prob.linear_constraints.add(rhs= my_rhs);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #           prob.linear_constraints.set_senses(constraintrowindex, "L");
  #           constraintrowindex = constraintrowindex+1;
  #           var = var+1;
  #         elif CosAngleBtwConEdges(Nodes[Edges[i][0]],Nodes[Edges[e][0]],Nodes[Edges[e][1]],Edges[i][2],Edges[e][2]) > FiftyLim:
  #           my_obj = ones(1)*(-5.0);
  #           prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #           my_rhs = ones(1);
  #           prob.linear_constraints.add(rhs= my_rhs);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #           prob.linear_constraints.set_senses(constraintrowindex, "L");
  #           constraintrowindex = constraintrowindex+1;
  #           var = var+1;

  #   for i in NTEL[Edges[e][1]-1]:
  #     if i not in CE:
  #       PairEdgeCount = PairEdgeCount+1;
  #       if i not in ETC:
  #         ETC.append(i);
  #       if Edges[i][0] == Edges[e][1]:
  #         if CosAngleBtwConEdges(Nodes[Edges[i][1]],Nodes[Edges[e][1]],Nodes[Edges[e][0]],Edges[i][2],Edges[e][2]) > CosAngLim:
  #           prob.linear_constraints.add(rhs= ones(1));
  #           prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #           prob.linear_constraints.set_senses(constraintrowindex, "L");
  #           constraintrowindex = constraintrowindex+1;
  #         elif CosAngleBtwConEdges(Nodes[Edges[i][1]],Nodes[Edges[e][1]],Nodes[Edges[e][0]],Edges[i][2],Edges[e][2]) > FortyLim:
  #           my_obj = ones(1)*(-10.0);
  #           prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #           my_rhs = ones(1);
  #           prob.linear_constraints.add(rhs= my_rhs);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #           prob.linear_constraints.set_senses(constraintrowindex, "L");
  #           constraintrowindex = constraintrowindex+1;
  #           var = var+1;
  #         elif CosAngleBtwConEdges(Nodes[Edges[i][1]],Nodes[Edges[e][1]],Nodes[Edges[e][0]],Edges[i][2],Edges[e][2]) > FiftyLim:
  #           my_obj = ones(1)*(-5.0);
  #           prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #           my_rhs = ones(1);
  #           prob.linear_constraints.add(rhs= my_rhs);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #           prob.linear_constraints.set_senses(constraintrowindex, "L");
  #           constraintrowindex = constraintrowindex+1;
  #           var = var+1;
  #       else:
  #         if CosAngleBtwConEdges(Nodes[Edges[i][0]],Nodes[Edges[e][1]],Nodes[Edges[e][0]],Edges[i][2],Edges[e][2]) > CosAngLim:
  #           prob.linear_constraints.add(rhs= ones(1));
  #           prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #           prob.linear_constraints.set_senses(constraintrowindex, "L");
  #           constraintrowindex = constraintrowindex+1;
  #         elif CosAngleBtwConEdges(Nodes[Edges[i][0]],Nodes[Edges[e][1]],Nodes[Edges[e][0]],Edges[i][2],Edges[e][2]) > FortyLim:
  #           my_obj = ones(1)*(-10.0);
  #           prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #           my_rhs = ones(1);
  #           prob.linear_constraints.add(rhs= my_rhs);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #           prob.linear_constraints.set_senses(constraintrowindex, "L");
  #           constraintrowindex = constraintrowindex+1;
  #           var = var+1;
  #         elif CosAngleBtwConEdges(Nodes[Edges[i][0]],Nodes[Edges[e][1]],Nodes[Edges[e][0]],Edges[i][2],Edges[e][2]) > FiftyLim:
  #           my_obj = ones(1)*(-5.0);
  #           prob.variables.add(obj = my_obj, types = [prob.variables.type.integer]);
  #           my_rhs = ones(1);
  #           prob.linear_constraints.add(rhs= my_rhs);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, 2*NumEdges+NumNodes+var, -1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, i-1, 1);
  #           prob.linear_constraints.set_coefficients(constraintrowindex, e-1, 1);
  #           prob.linear_constraints.set_senses(constraintrowindex, "L");
  #           constraintrowindex = constraintrowindex+1;
  #           var = var+1;

  # EndCurvRes = constraintrowindex;

  # print count;
  # quit();
