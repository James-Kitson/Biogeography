# -*- coding: utf-8 -*-
"""
This script is written in Python 3 and takes the parsed BEAST output from 
ggtree to write a sensibly formatted table of probabilities for locations on nodes.
This can then be used to put pie chards on nodes in ggtree.
"""

#% reset

pie_dict = {}

#create a list of islands we have included
my_locations = ['Aldabra','Anjouan','Europa','Grande_Comore','Grande_Glorieuse','Juan_de_Nova','Madagascar',\
'Mauritius','Moheli','Reunion','Rodrigues','Seychelles']

#read in the file and ignore the headers
pies = open('./BEAST_geo.txt').readlines()[1:]
#for each line do a series of replaces and splits to remove the structure placed in there by R
for line in pies:
    line = line.replace('\n','')
    line = line.replace("\\",'')
    line = line.replace('"','')
    line = line.replace('c(','')
    line = line.replace(')','')
    #the result is a list of node, node_locations and probabilities for each location
    line = line.split('\t')
    
    #for each line we'll split this list into separate items
    node = line[0]
    node_locations = line[1].replace(' ','').split(',')
    probs = line[2].replace(' ','').split(',')
    
    node_probs = {}
    #for each island in our complete list of locations
    for j in my_locations:
        #we check if each location in present on our node
        if j in node_locations:
            #iterate across nodes present on node
            for index, island in enumerate(node_locations):
                #if the node_location and complete location list match
                if island == j:
                    #append the node location and probabilities to a dictionary ({island:prob}..n)
                    node_probs[island] = probs[index]
        #drop outside that loop to prevent checking for all combinations
        else:
            #append the absent node locations and zero probabilities to the same dictionary above ({island:prob}..n)
            node_probs[j] = str(0)
    ##append the dictionaries of islands and probabilities to an overall dictionary with node as the key ({node:{island:prob}..n})
    pie_dict[node] = node_probs

#import pands so we can work with dataframe style structures
import pandas as pd

#make the dict of dicts into a dataframe
pie_dataframe = pd.DataFrame(pie_dict)

#it's the wrong way round (islands = rows, nodes = columns) so we'll transpose it
pie_dataframe = pie_dataframe.transpose()

#rename the dataframe index column to node
pie_dataframe = pie_dataframe.reindex(pie_dataframe.index.rename('node'))

#write the dataframe to a csv readable by R and ggtree
pie_dataframe.to_csv('BEAST_geo_parsed.csv')





