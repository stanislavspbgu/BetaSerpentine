#!/usr/bin/python

import os, shutil, sys, getopt, datetime

"""import functions from folder SerpenrinesFunctions"""
from SerpentinesFunctions.import_ArchCandy_data_v4 import *
from SerpentinesFunctions.find_Serpentines_v7 import *
from SerpentinesFunctions.print_serpentines_v8 import *
from SerpentinesFunctions.consensus_verification import *

def Log(Loglist, info):
    """Append info to LogFile and display info on the screen"""
    print info[:-1]
    Loglist.append(info)


def BetaSerpentine(FileName, OutputFolder, SortByConsensus = True, InputThreshold = 0, OutputThreshold = 0.2):
    """
    Main function of BetaSerpentines algorithm.
    :param FileName: path to the file with informatioin about beta-arches predicted by ArchCandy;
    :param InputThreshold: threshold for beta-Arches score, beta-arches with lower score are exlcuded;
    :param OutputThreshold: treshold for output beta-serpentines.
    :return:
    """

    version = 'SERPENTINES 1.0'
    Loglist = []
    """Create begin of LogFile (Version, value of InputThreshold and time of start)"""
    Log(Loglist, version + '\n')
    Log(Loglist, '\n')
    Log(Loglist, 'InputThreshold: ' + str(InputThreshold) + '\n')
    Log(Loglist, '\n')
    Log(Loglist, str(datetime.datetime.now().time())[:-7] + ' Start' + '\n')
    Log(Loglist, '\n')

    """ Import protein sequence and beta-arches """
    tmpList = inputArches(FileName)
    ListOfStrings = tmpList[1]
    Sequence = tmpList[0].strip()
    ProtName = tmpList[2]

    Log(Loglist, 'Filename: ' + FileName + '\n')
    Log(Loglist, 'Protein Name: ' + ProtName + '\n')
    Log(Loglist, 'Sequence: ' + Sequence + '\n')
    Log(Loglist, '\n')


    table_of_Arches = rm_ident_arches2(ListOfStrings) # import arches and remove identical ones

    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)
    os.chdir(OutputFolder)

    if not os.path.exists(ProtName +'_Inp'+ str(InputThreshold)):
        os.makedirs(ProtName +'_Inp'+ str(InputThreshold))

    os.chdir(ProtName +'_Inp'+ str(InputThreshold))

    Arches = ArchesImport(table_of_Arches, InputThreshold)  # creating a work set of arches
    MaxNArches = ArchesDestributions(Arches, ProtName)

    Log(Loglist, 'Total number of beta-arches: ' + str(len(Arches)) + '\n')

    Log(Loglist, 'Maximal number of beta-arches per AA: ' + str(MaxNArches) + '\n')

    """ Set of filters"""

    if MaxNArches > 80:
        Log(Loglist, '\n')
        Log(Loglist,  "Too much arches in the protein. Please increase threshold or try to analyse part of the protein sequence." + '\n')
        Log(Loglist, 'Maximal number of beta-arches per residue: ' + str(MaxNArches) + ', only 80 are acceptable.' + '\n')
        return
    if len(Sequence) < 19:
        Log(Loglist, '\n')
        Log(Loglist, "Input sequence is below the minimal allowed length. Please use sequence longer that 18 amino acids." + '\n')
        return
    if len(Arches) == 0:
        Log(Loglist, '\n')
        Log(Loglist, "There are no arches. You can try to decrease beta-arch threshold" + '\n')
        return

    """ BetaSerpentine execution"""
    matrix = create_compare_matrix(Arches, P=False)
    graph = graph_from_matrix(matrix)
    Cgraph = graph[0]
    Ngraph = graph[1]
    Serpentines_Res = AllSerpentines(Cgraph, Ngraph, Arches, Sequence, OutputThreshold)

    Serpentines = Serpentines_Res[0]
    NStrands = Serpentines_Res[1]
    NArcs = Serpentines_Res[2]
    Consensus = verify_consensus(NStrands, NArcs, Sequence)

    Log(Loglist, 'Sequence       : ' + str(Consensus[2]) + '\n')
    Log(Loglist, 'Structure      : ' + str(Consensus[1]) + '\n')
    Log(Loglist, 'p-value < 0.001: ' + str(Consensus[0]) + '\n')

    Log(Loglist, str(datetime.datetime.now().time())[:-7] + ' Total number of serpentines: ' + str(len(Serpentines))+'\n')

    """ Export beta-serpentines"""
    if SortByConsensus:
        ExportSerpintines(Serpentines, Consensus[3], FileName = ProtName + '_serpentines.txt', sep = ' ', SortByConsensus = True)
    else:
        ExportSerpintines(Serpentines, Consensus[3], FileName = ProtName + '_serpentines.txt', sep = ' ', SortByConsensus = False)

    ExportSerpintinesTableShort(Serpentines, FileName = ProtName + '_serpentines_table_short.txt')

    Log(Loglist, str(datetime.datetime.now().time())[:-7] + ' End' +'\n')

    LogFile = open(ProtName + '_log.txt', 'wb')
    for string in Loglist: LogFile.write(string)
    LogFile.close()

    return

error_message = '\nERROR, please check parameters:\n- Filename - path to the file with list of b-arches (required);\n- OutputFolder - path to folder for output files;\n- SortByConsensus - True or False (optional, default value True);\n- InputThreshold - from 0.0 to 1.0 (optional, default value 0.0);\n- OutputThreshold - from 0.o to 1.0 (optional, default value 0.2).\nExample: python Serpentines_CL.py test_arches.txt True 0 0.2\n'

if len(sys.argv) == 6: BetaSerpentine(str(sys.argv[1]), str(sys.argv[2]), sys.argv[3], float(sys.argv[4]), float(sys.argv[5]))
elif len(sys.argv) == 3: BetaSerpentine(str(sys.argv[1]), str(sys.argv[2]))
else: print error_message

