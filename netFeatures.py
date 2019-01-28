import sys
import getopt
from scipy.stats.stats import pearsonr
import networkx as nx
import numpy as np
from sklearn.preprocessing import normalize

option, args = getopt.getopt ( sys.argv[1:], 'n:i:r:p:o:s:' )

CONV_THRESHOLD = 0.000001

sInNet = ''
sInProfile = ''
#sNetOutFile = 
#sStatOutFile = 

dctProfile = {}

dctRho = {}
dctP = {}
dctNodeCorr = {}

g = nx.Graph()
g_sub = nx.Graph()
g_weighted = nx.Graph()

for op, p in option:
	if op == '-n':
		sInNet = p
	elif op == '-i':
		sInProfile = p
	elif op == '-r':
		nR = float ( p )
	elif op == '-p':
		nPval = float ( p )
	elif op == '-o':
		sNodeOutFile = p + '.node.corr.txt'
		sNetOutFile = p + '.weighted.ppi'
		sStatOutFile = p + '_corr.stats.txt'
		sNetPOutFile = p + '_netProp.txt'
	elif op == '-s':
		lstSamples = p.split ( ':' )
	else:
		print "Unknown option", op
		print "-n <PPI network> -i <profile> -r <row> -p <p-value> -o <OutPut Prefix> -s <# of Samples (treated:control)>"
		exit()

'''
lstSampleVector = np.ones ( ( 1, int ( lstSamples[0] ) ) )
print ( lstSampleVector )
lstSampleVector = np.append ( lstSampleVector, -1 * np.ones ( ( 1, int ( lstSamples[0] ) ) ) )
print ( -1 * np.ones ( ( 1, int ( lstSamples[0] ) ) ) )
print ( lstSampleVector )
'''
lstSampleVector = np.full ( ( 1, int ( lstSamples[0] ) ), 1 )
lstSampleVector = np.append ( lstSampleVector, np.full ( ( 1, int ( lstSamples[1] ) ), -1 ) )

sInFile = open ( sInProfile, 'r' )
sOutFile = open ( sNodeOutFile, 'w' )

sOutLine = ""

sReadLine = sInFile.readline()
for sReadLine in sInFile.readlines():
	sReadLine = sReadLine.replace ( '\n', '' )
	sToken = sReadLine.split ( '\t' )
	dctProfile[sToken[0]] = sToken[1:]
	c, p = pearsonr ( lstSampleVector.astype ( np.float ) , np.array ( sToken[1:] ).astype ( np.float ) )
	if str ( c ) == "nan":
		c = 0
	sOutLine += str ( sToken[0] ) + "\t" + str ( c ) + "\t" + str ( p ) + "\n"
	dctNodeCorr[sToken[0]] = abs(float(c))

sOutFile.write ( sOutLine )
sOutFile.close()

#exit()

#print pearsonr( np.array( dctProfile['1'] ).astype(np.float), np.array ( dctProfile['10'] ).astype(np.float) )
sInFile = open ( sInNet, 'r' )
for sReadLine in sInFile.readlines():
	sReadLine = sReadLine.replace ( '\n', '' )
	sToken = sReadLine.split ( '\t' )
	g.add_edge ( sToken[0], sToken[1] )

sOutFile = open ( sNetOutFile, 'w' )
"""
g.add_edge(1,2)
g.add_edge(1,9)
g.add_edge(2,10)
g.add_edge(9,13)
g.add_edge(1,14)
g.add_edge(2,13)
g.add_edge(2,16)
"""
sOutLine = ""

edge_list = []

for edge in g.edges():
#	lstA = dctProfile[str(edge[0])]
#	lstB = dctProfile[str(edge[1])]
	if edge[0] in dctProfile and edge[1] in dctProfile:
		c, p = pearsonr( np.array( dctProfile[str(edge[0])] ).astype(np.float), np.array ( dctProfile[str(edge[1])] ).astype(np.float) )
#		if c > nR and p < nPval:
		if str(c) == "nan":
			c = 0
		sOutLine += str ( edge[0] ) + "\t" + str ( edge[1] ) + "\t" + str ( c ) + "\t" + str ( p ) + "\n"
		try:
			dctRho[ str ( edge[0] ) ].append( float ( c ) )
		except KeyError:
			dctRho[ str ( edge[0] ) ] = []
			dctRho[ str ( edge[0] ) ].append( float ( c ) )

		try:
			dctRho[ str ( edge[1] ) ].append( float ( c ) )
		except KeyError:
			dctRho[ str ( edge[1] ) ] = []
			dctRho[ str ( edge[1] ) ].append( float ( c ) )

		try:			
			dctP[ str ( edge[0] ) ].append( float ( p ) )
		except KeyError:
			dctP[ str ( edge[0] ) ] = []
			dctP[ str ( edge[0] ) ].append( float ( p ) )

		try:
			dctP[ str ( edge[1] ) ].append( float ( p ) )
		except KeyError:
			dctP[ str ( edge[1] ) ] = []
			dctP[ str ( edge[1] ) ].append( float ( p ) )
			

#		dctRho[ str ( edge[0] ) ].append( float ( c ) )
#		dctRho[ str ( edge[1] ) ].append( float ( c ) )
#		dctP[ str ( edge[0] ) ].append( float ( p ) )
#		dctP[ str ( edge[1] ) ].append( float ( p ) )

		if c > nR and p < nPval:
			g_sub.add_edge ( edge[0], edge[1] )
		edge_list.append ( ( edge[0], edge[1], abs ( float ( c ) ) ) )
#		sLineOut = str ( edge[0] ) + "\t" + str ( edge[1] ) + "\t" + str ( c ) + "\t" + str ( p )
#		print sLineOut
#		sLineOut = str(edge[0]) + "\t" + str(edge[1]) + "\t" + str(c) + "\t" + str(p)
#		print sLineOut
#	print edge[0], edge[1]

sOutFile.write ( sOutLine )
sOutFile.close()

g_weighted.add_weighted_edges_from ( edge_list )

sOutFile = open ( sStatOutFile, 'w' )

sOutLine = ""

dgr = g.degree()
dgr_sub = g_sub.degree()


for gene in dctRho.keys():
	d = dgr[gene]
	d_sub = 0
	if gene in dgr_sub:
		d_sub = dgr_sub[gene]

	ratio = float( d_sub ) / float ( d )
	mean_r = np.mean( np.array ( dctRho[ str(gene) ] ) )
	std_r = np.std( np.array ( dctRho[ str(gene) ] ) )
	mean_p = np.mean( np.array ( dctP[ str(gene) ] ) )
	std_p = np.mean( np.array ( dctP[ str(gene) ] ) )

	sOutLine += str ( gene ) + "\t" + str ( d_sub ) + "\t" + str ( ratio ) + "\t" + str ( mean_r ) + "\t" + str ( std_r ) + "\t" + str ( mean_p ) + "\t" + str ( std_p ) + "\n"
#	print sLineOut

sOutFile.write ( sOutLine )
sOutFile.close()


################## Network Propagation ###########################
def normalize_cols ( matrix ):
	return normalize ( matrix, norm = 'l1', axis = 0 )

def generate_rank_list ( p_t ):
	gene_probs = zip ( g_weighted.nodes(), p_t.tolist() )
	for s in sorted ( gene_probs, key = lambda x: x[1], reverse = True ):
		yield s[0], s[1]

restart_prob = 0.7
original_graph = g_weighted
og_not_normalized = nx.to_numpy_matrix(original_graph)
og_matrix = normalize_cols(og_not_normalized )
p_0 = np.array([0.0] * g_weighted.number_of_nodes())

nodes = g_weighted.nodes()

nTotal = 0
for n in nodes:
	nTotal += dctNodeCorr[n]

for n in nodes:
	index = g_weighted.nodes().index( n )
	p_0[index] = dctNodeCorr[n] / float ( nTotal )
	
diff_norm = 1
p_t = np.copy(p_0)

while ( diff_norm > CONV_THRESHOLD ):
	epsilon = np.squeeze ( np.asarray ( np.dot ( og_matrix, p_t ) ) )
	no_restart = epsilon * ( 1 - restart_prob )
	restart = p_0 * restart_prob
	p_t_1 = np.add ( no_restart, restart )

	diff_norm = np.linalg.norm ( np.subtract ( p_t_1, p_t ), 1 )

	p_t = p_t_1

sOutLine = ""

for n, prob in generate_rank_list(p_t):
	print ( '{}\t{:.10f}'.format(n, prob) )
	sOutLine += str ( '{}\t{:.10f}'.format(n, prob) ) + "\n"

sOutFile = open ( sNetPOutFile, 'w' )
sOutFile.write ( sOutLine )
sOutFile.close()
