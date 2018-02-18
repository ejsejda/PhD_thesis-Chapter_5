###############################################################################################
#####                        THIS SCRIPT WAS USED IN CHAPTER 5                              ###
####         TO CREATE DOMAINS BIGRAM NETWORK AND CALCULATING THE NETWORK PROPERTIES        ###
####      AS WELL AS TO DETECT LOUVAIN COMMUNITIES WITHIN THE LARGEST CONNECTED COMPONENT   ###
###############################################################################################

#!/usr/bin/python
import sys
import community
from community import *
import networkx as nx
import re
from scipy.stats.stats import pearsonr
import scipy
import matplotlib.pyplot as plt
from colour import Color

print "starting script.."

workDir1 = "working_directory"

####  Input files  
infile1 = "%s/PHI_FGSG_List" % (workDir1)
infile2 = "%s/pfamDomCount_inFG_withUniquePfamNo.txt" % (workDir1)
infile3 = "%s/hetero_bigrams_tabDel.txt" % (workDir1)
fh1 = open (infile1, "r")
fh2 = open (infile2, "r")
fh3 = open (infile3, "r")

###  output files
outfile1 ="%s/connectedComponentsInfo.txt" % (workDir1)
outfile2 = "%s/pfam_domains_to_community_detection.txt" % (workDir1)
outfile3 = "%s/nodes_clustering_coefficient.txt" % (workDir1)
outfile4 = "%s/nodes_degree.txt" % (workDir1)
outfile5 = "%s/nodes_degree_centrality.txt" % (workDir1)
fh_out1 = open (outfile1, "w")
fh_out2 = open (outfile2, "w")
fh_out3 = open (outfile3, "w")
fh_out4 = open (outfile4, "w")
fh_out5 = open (outfile5, "w")

######################################
fg_pfam = dict()
duf_nodes = dict()
domSet = set()	
			
#declaration of dictionaries of all possible mutants phenotypes in FG proteome as per PHI-base
phiFg_InVr = dict()#increased virulence
phiFg_Lethal = dict()#lethal
phiFg_LossPath = dict()#loss_of_pathogenicity
phiFg_RedVir = dict()#reduced_virulence
phiFg_UnaffPath = dict()#unaffected pathogenicity 
phiFg_MixOut = dict()#mixed outcome
phiFg_UnaffVir = dict()#Unaffected virulence

pfamInVr = set()#pfam domains only associated with increase virulence mutant proteins
pfamLetal = set()#pfam domains only associated with lethal mutant proteins
pfamLossPath = set()#pfam domains only associated with loss pathogenicity mutant proteins
pfamRedVir = set()#pfam domains only associated with reduced virulence mutant proteins
pfamUnaffPath = set()#pfam domains only associated with unaffected pathogenicity mutant proteins
pfamMixOut = set()#pfam domains only associated with mixed outcome mutant proteins
pfamUnaffVir = set()#pfam domains only associated with unaffected virulence mutant proteins
#######################################################################################
def getPhenoTyDict(a,d,k,v):
	if phenoType.startswith(a):
		d.setdefault(k,v)
     	return len(d)
########################################################################################
def pfamWithPhenotype(d,s,t):
	c =0
	for k,v in fg_pfam.items():
		if d.has_key(k):
			c=c+1
			print t, c, k, v
			for e in v:
				if e not in s:
					s.add(e)
			
			#print t, "pfam domains:"
	return ( t+ " pfam domains:", s)
##########################################################################
####            MAIN PROGRAM      ########################################
##########################################################################
for line in fh2.xreadlines():
	fields = line.split()
	fgId= fields[0].strip()
	pfamDoms = fields[3:]

	for dom in pfamDoms:
		fg_pfam.setdefault(fgId,set())
		if dom not in fg_pfam[fgId]:
			fg_pfam[fgId].add(dom)


for line in fh1.xreadlines():
	fields = line.split('\t')
	phiBaseId = fields[0].strip()
	fgId_inPhi = fields[1].strip()
	phenoType = fields[3].strip()
	
	lethal = getPhenoTyDict('Lethal', phiFg_Lethal, fgId_inPhi,phiBaseId)
	incVir = getPhenoTyDict('Increased virulence', phiFg_InVr, fgId_inPhi,phiBaseId)
	LossPath = getPhenoTyDict('Loss of pathogenicity', phiFg_LossPath, fgId_inPhi,phiBaseId)
	RedVir = getPhenoTyDict('Reduced virulence', phiFg_RedVir, fgId_inPhi,phiBaseId)
	UnaffPath = getPhenoTyDict('Unaffected pathogenicity', phiFg_UnaffPath, fgId_inPhi,phiBaseId)
	UnaffVir = getPhenoTyDict('Unaffected virulence', phiFg_UnaffVir, fgId_inPhi,phiBaseId)
	MixOut = getPhenoTyDict('mixed outcome', phiFg_MixOut, fgId_inPhi,phiBaseId)

print "Unique No of Lethal phenotypes in FG- source PHI-base:", lethal
print "Unique No of Increase virulence phenotypes in FG- source PHI-base:", incVir
print "Unique No of Increase virulence phenotypes in FG- source PHI-base:", LossPath
print "Unique No of Reduced virulence phenotypes in FG- source PHI-base:", RedVir
print "Unique No of Unaffected pathogenicity phenotypes in FG- source PHI-base:",UnaffPath
print "Unique No of Unaffected virulence phenotypes in FG- source PHI-base:",UnaffVir
print "Unique No of mixed outcome phenotypes in FG- source PHI-base:", MixOut

print pfamWithPhenotype(phiFg_InVr,pfamInVr, "Increased virulence")
print pfamWithPhenotype(phiFg_Lethal,pfamLetal, "Lethal")
#Lethal= pfamWithPhenotype(phiFg_InVr,pfamInVr, "Lethal")
#fhout.write(' '.join(map(str,Lethal)))
print pfamWithPhenotype(phiFg_LossPath,pfamLossPath, "Loss of pathogenicity")
print pfamWithPhenotype(phiFg_RedVir,pfamRedVir, "Reduced virulence")
print pfamWithPhenotype(phiFg_UnaffPath,pfamUnaffPath, "Unaffected pathogenicity")
print pfamWithPhenotype(phiFg_UnaffVir,pfamUnaffVir, "Unaffected virulence")
print pfamWithPhenotype(phiFg_MixOut,pfamMixOut, "mixed outcome")

#################################################################
###   Building domains bigram network
#################################################################
G = nx.Graph()
for line in fh3.xreadlines():
	fields = line.split('\t')
	node1 = fields[0].strip()
	node2 = fields[1].strip()
	G.add_edge(node1,node2, weight = float(fields[4].strip()))
fh3.close()

###################################################################
###   Calculating network properties  #############################
###################################################################
print "nodes:", G.number_of_nodes()
print "edges:", G.number_of_edges()

### connected components
network_CC = nx.number_connected_components(G)
print "transitivity",nx.transitivity(G)

# connected components statistics including number of nodes, edges  and number of DUFs per connected component
print "number of connected components:", network_CC
i = 1
for k in nx.connected_component_subgraphs(G):
	#print i, "CC has", k.number_of_nodes(),"Nodes and ", k.number_of_edges(), "Edges"
	fh_out1.write("%s\t%s\t%s\t%s" %(i, k.number_of_nodes(), k.number_of_edges(), len(duf_nodes)))
	i = i+1


### clustering coefficient
clusters = nx.clustering(G)
print "average aclustering coefficient of the network:" , nx.average_clustering(G)	
for k,v in clusters.items():
	c = 1
	#print c, k, v
	fh_out3.write("%s\t%s\t%s\n" %(c,k,v))
	c = c+1
fh_out3.close()

### node degree
nodeDegree = sorted(nx.degree(G).values(), reverse = True)
degree_dict = dict(nx.degree(G))
for k, v in degree_dict.items():
	c = 1
	#print c, k, v
	fh_out4.write("%s\t%s\t%s\n" % (c, k, v))
	c = c+1

### degree centrality
degree_centr_dict = dict(nx.degree_centrality(G))
for k, v in degree_centr_dict.items():
	c = 1
	#print c, k, v
	fh_out5.write("%s\t%s\t%s\n" % (c, k, v))
	c = c+1
fh_out5.close()

##########################################################
###  assigning number of DUF per connected component   ###
number = nx.number_connected_components(G)

for n in range(0,(number)):
	H = nx.connected_component_subgraphs(G)[n]
	for k in H.nodes():
		nodesNo = len(H.nodes())
		if k.startswith('DUF'):
			duf_nodes.setdefault(n, list())
			duf_nodes[n].append(k)		
	if duf_nodes.has_key(n):
		#print n, nodesNo,len(duf_nodes[n])
		fh_out1.write("%s\t%s\t%s\t%s\t%s\n" %(n, nodesNo, H.nodes(), len(duf_nodes[n]), duf_nodes[n]))
	else:
		#print n, nodesNo, 0
		fh_out1.write("%s\t%s\t%s\t%s\n" %(n, nodesNo, H.nodes(), 0))	 
fh_out1.close()
###############################################
####     louvian community detection       ####
whole_part = community.best_partition(G)
print "Modularity of a partition of the whole network:",modularity(whole_part, G)
H = nx.connected_component_subgraphs(G)[0]
part = community.best_partition(H)
print "Partitions: ", float(len(list(part.keys())))
print "Modularity of a partition of the largest CC:",modularity(part, H)

pfamDUF = set()
duf_InVr = set()
duf_Letal = set()
duf_LossPath = set()
duf_RedVir = set()
duf_UnaffPath = set()
duf_MixOut = set()
duf_UnaffVir = set()

for n in H.nodes():
	if n.startswith('DUF'):
		pfamDUF.add(n)

duf_InVr= pfamDUF.intersection(pfamInVr)
duf_Letal = pfamDUF.intersection(pfamLetal)
duf_LossPath = pfamDUF.intersection(pfamLossPath)
duf_RedVir = pfamDUF.intersection(pfamRedVir)
duf_UnaffPath = pfamDUF.intersection(pfamUnaffPath)
duf_MixOut = pfamDUF.intersection(pfamMixOut)
duf_UnaffVir = pfamDUF.intersection(pfamUnaffVir)
###############################################################################################
###   assigning pfam domain to community detected within first largest connected component  ###
for h in H.nodes():
	print h,(part[h])
	fh_out2.write("%s\t%s\n" %(h, part[h]))

############################################
####   drawing graph with communities   ####

size = float(len(set(part.values())))
pos = nx.graphviz_layout(H, prog = "neato")
count = 0
duf_UnaffPath_RedVir = set()
duf_UnaffPath_RedVir = duf_UnaffPath.intersection(duf_RedVir)

for com  in set(part.values()):
	count = count+1	
	list_nodes = list()
	duf_nodes = list()
	colorList = list()
	#list_nodes = [nodes for nodes in part.keys() if part[nodes] == com]
	for nodes in part.keys():
		if part[nodes] == com:
			list_nodes.append(nodes)
		
			#DUFs with associated phenotype decorated with diffrent colours. 
			#If no mutant phenotype associatedwith DUF the node is decorated in red
			if nodes in pfamDUF:
				
				if nodes in duf_UnaffPath:
					if nodes in duf_UnaffPath_RedVir:
						colorList.append('saddlebrown')
					else:
						colorList.append('cyan')
				elif nodes in duf_RedVir:
					colorList.append('deeppink')
				elif nodes in duf_Letal:
					colorList.append('chartreuse')
				else:
					colorList.append('red')			
			
			else:
				colorList.append(str(com/size))	#different shade of grey represent different communities detected
									
	nx.draw_networkx_nodes(H , pos, list_nodes, with_labels = False, node_size = 150, node_color = colorList, alpha = 1.0)
nx.draw_networkx_edges(H, pos, alpha=0.5)

##  ploting the largest connected component with communities detected and DUF depicted in diffrent colour.
##  DUFs with associated phenotype decorated with diffrent colours.
#Thus, if pfam node belong to unaffected pathogenicity and reduce virulence mutant protein node is depicted in brown. 
#If pfam is only part of unaffected pathogenicity mutant protein node is decorated in cyan.
#If pfam is only part of reduced virulence mutant protein node is decorated in deep pink.
#If pfam is only part of lethal mutant protein node is decorated in green-yellow (chartreuse).
#If no mutant phenotype associatedwith DUF the node is decorated in red
plt.show()
plt.savefig('/home/ela/Project/pfam(sinceJune2013)/output/fg_bigramsNetwork/bigramOrderNoMatter/heteroBigramsNetwork/partitions.png')		

