# SET OF FUNCTIONS TO PRINT SERPENTINES

from SerpentinesFunctions.consensus_verification import *
import operator

# script to visualise Serpintines


def PrintArc(arch, F=True, sep=' '):
    """
    Returns list of strings which display structure of Arc when printed.
    arch - amino acids sequence of Arc
    F - orientation of Arc (see Examples)
    sep - symbol to fill gaps.

    Input
    PrintArc('GTY', F = True, sep = ' ')

    Output:
    'g '
    '  '
    ' t'
    '
    '  '
    'y '

    Input
    PrintArc('GPYG', F = False, sep = ' ')

    Output:
    ' g'
    'p '
    '  '
    'y '
    '  '
    ' g'
    """
    L = len(arch)  # length of Arc
    strings = []  # empty list for output
    # dictionary to transform uppercase letters to lowercase
    dict = {'G': 'g', 'P': 'p', 'A': 'a', 'V': 'v', 'L': 'l', 'I': 'i', 'M': 'm', 'C': 'c', 'F': 'f', 'Y': 'y', 'W': 'w',
            'H': 'h', 'K': 'k', 'R': 'r', 'Q': 'q', 'N': 'n', 'E': 'e', 'D': 'd', 'S': 's', 'T': 't'}
    # transform uppercase letters to lowercase
    seq = ''  # empty string to collect lowercase letters
    for i in arch:
        seq = seq + dict[i]
        
    # collect strings to print structure of Arcs of different types
    if L == 3:
        strings.append(seq[0] + sep)
        strings.append(sep + sep)
        strings.append(sep + seq[1])
        strings.append(sep + sep)
        strings.append(sep + sep)
        strings.append(seq[2] + sep)

    elif L == 4:
        strings.append(seq[0] + sep)
        strings.append(sep + seq[1])
        strings.append(2*sep)
        strings.append(sep + seq[2])
        strings.append(sep + sep)
        strings.append(seq[3] + sep)

    elif L == 5:
        strings.append(seq[0] + sep)
        strings.append(sep + seq[1])
        strings.append(sep + seq[2])
        strings.append(2*sep)
        strings.append(sep + seq[3])
        strings.append(seq[4] + sep)

    elif L == 6:
        strings.append(seq[0] + sep)
        strings.append(sep + seq[1])
        strings.append(sep + seq[2])
        strings.append(sep + seq[3])
        strings.append(sep + seq[4])
        strings.append(seq[5] + sep)

    output = strings

    # flip strings
    if F == False:
        stringsR = []
        for i in strings:
            stringsR.append(i[::-1])
        output = stringsR        
    return output

    
def PrintStrand(strand, F=True, toArch=True, sep=' '):
    """
    Function returns list of two strings which display structure of Strand.
    strand - amino acids sequence of Strand
    F - direction of Strand (left to right or right to left)
    toArch - is strand before Arc
    sep - symbol to fill gaps

    Input
    PrintStrand('NQGNNQQNYQ', F = True, toArch = True, sep = ' ')

    Output:
    'N G N Q Y '
    ' Q N Q N Q'

    Input
    PrintStrand('NQGNNQQNYQ', F = False, toArch = True, sep = ' ')

    Output:
    'Q N Q N Q '
    ' Y Q N G N'

    Input
    PrintStrand('NQGNNQQNYQ', F = False, toArch = False, sep = ' ')

    Output:
    ' Y Q N G N'
    'Q N Q N Q '
    """
    seq = strand
    Line1 = ''  # output empty Lines
    Line2 = ''
    # fill Lines according to the orientation of strand (arguments F and toArc)
    if F == True and toArch == True:
        # counter to discriminate odd and even amino acids in strand
        index = 0

        # fill Lines according to orientation of strand (arguments F and toArc)
        for i in seq:
            if index % 2 == 0:
                Line2 = Line2 + i
                Line1 = Line1 + sep
            elif index % 2 != 0:
                Line2 = Line2 + sep
                Line1 = Line1 + i
            index = index + 1
        return [Line1,Line2]   

    # fill Lines according to orientation of strand (arguments F and toArc)
    elif F == False and toArch == False:
        index = 0  # fill Lines according to orientation of strand (arguments F and toArc)
        for i in seq:
            if index % 2 == 0:
                Line1 = Line1 + i
                Line2 = Line2 + sep
            elif index % 2 != 0:
                Line1 = Line1 + sep
                Line2 = Line2 + i
            index = index + 1
        return [Line1[::-1],Line2[::-1]]  # flip sequences

    # fill Lines according to orientation of strand (arguments F and toArc)
    elif F == True and toArch == False:
        index = 0  # fill Lines according to orientation of strand (arguments F and toArc)
        for i in seq:
            if index % 2 == 0:
                Line1 = Line1 + i
                Line2 = Line2 + sep
            elif index % 2 != 0:
                Line1 = Line1 + sep
                Line2 = Line2 + i
            index = index + 1
        return [Line1, Line2]


def splitSeqByStructure(structure):
    """
    Function returns dictionary of two lists: arcs sequences and strands sequences.
    Input: protein sequence, string describing structure,
    sep(symbol to fill gaps).

    Input.
    Protein sequence: 'NQGNNQQNYQQYSQNGNQQQGNNRYQGYQAYNAQAQPAGGYYQNYQGYSGYQQ'
    string describing structure: 'sSsSsSsSsSsSsSsAAAAAsSsSsSsSsSsSsAAAAAsSsSsSsSsSsSsSs'

    Output.
    'strands':{'NQGNNQQNYQQYSQN', 'GNNRYQGYQAYNA', 'GGYYQNYQGYSGYQQ'}
    'arches': {'GNQQQ', 'QAQPA'}
    """
    strands = []
    arches = []
    for i in range(len(structure)):
        if i % 2 == 0:
            strand = ''
            for aa in structure[i]: strand += aa
            strands.append(strand)
        else:
            arc = ''
            for aa in structure[i]: arc += aa
            arches.append(arc)

    data = {'arches': arches, 'strands': strands}  # dictionary with arches and strands
    return data


def insertString(Row,Col,ListOfStrings,Matrix):
    """
    Function takes a matrix and replaces elements in it with symbols from input strings,
    coordinates of replacements(left upper) are defined as Row and Col.
    Function returns a matrix with replaced elements.

    Input.
    Row = 2
    Col = 2
    ListOfStrings = ['CAG','GA']
    Matrix:
    _ _ _ _
    _ _ _ _
    _ _ _ _
    _ _ _ _


    Output.
    _ _ _ _
    _ C A G
    _ G A _
    _ _ _ _
    """
    for i in range(len(ListOfStrings)):
        index = 0
        for j in ListOfStrings[i]:
            Matrix[Row + i][Col + index] = j
            index = index + 1
            
    # for i in Matrix: print i
    return Matrix


def PrintSerpentine(structure, sep=' '):
    """
    Function returns a list of strings to print visual illustration of Serpentine.
    structure - string describing a structure.

    Input
    [['T', 'A', 'N'], ['N', 'T', 'C', 'A'], ['T', 'Q'], ['L', 'A', 'N'], ['F', 'L', 'V']]

    Output
    [' t c /', 'n a V ',' N  T  L ', 'A Q F ',' T l n ','/ a']
    """
    L = 0
    for l in structure: L += len(l)
    # create a dictionary of strands and lists
    data = splitSeqByStructure(structure)
    # create a set of coordinates
    NStrand = len(data['strands'])

    # module to collect coordinates of Strands and Arcs
    StrandLengths = []
    for i in data['strands']:
        StrandLengths.append(len(i))

    Cols = [int(L/2)]  # create a list for Col coordinates for Strands and set first one
    Rows = []  # an empty list for Row coordinates for Strands, collect Row and Col coordinates.
    # Row coordinates correspond to each forth line. Number of coordinates is equivalent to number of Strands
    # Col coordinates is calculated for s for each Strand.
    for row in range(0,NStrand):
        Rows = Rows + [row*4]
        if row < (NStrand - 1) and row % 2 == 0:
            d = (StrandLengths[row] - StrandLengths[row + 1])
        elif row < (NStrand - 1) and row % 2 != 0:
            d = 0

        Cols.append(Cols[row] + d)

    ArchCols = []  # create a list for Col coordinates for Arcs
    for i in range(len(data['arches'])):
        if i % 2 == 0: ArchCols.append(Cols[i] + StrandLengths[i])
        elif i % 2 != 0: ArchCols.append(Cols[i] - 2)

    # create empty matrix
    Matrix = [[sep for x in range(L * 2)] for x in range(NStrand * 8)]

    # print strands into matrix
    count = 0
    for i in data['strands']:
        if count == 0:            
            tmpList = PrintStrand(i, F = True, toArch = True)
            insertString(Rows[count],Cols[count],tmpList,Matrix)
        else:
            if count % 2 != 0:
                tmpList = PrintStrand(i, F = False, toArch = False)
                insertString(Rows[count],Cols[count],tmpList,Matrix)

            elif count % 2 == 0:
                tmpList = PrintStrand(i, F = True, toArch = False)
                insertString(Rows[count],Cols[count],tmpList,Matrix)
                
        count = count + 1

    # print arches into matrix
    Archcount = 0
    for i in data['arches']:
        if Archcount % 2 == 0:
            tmpList = PrintArc(i, F=True)
            insertString(Rows[Archcount],ArchCols[Archcount],tmpList,Matrix)
            
        elif Archcount % 2 != 0:
            tmpList = PrintArc(i, F=False)
            insertString(Rows[Archcount],ArchCols[Archcount],tmpList,Matrix)

        Archcount = Archcount + 1

    # add 'start' and 'end' symbols('/','\')
    Matrix[Rows[0]][Cols[0]-1] = '/' # start
    if len(data['strands']) % 2 != 0:
        Matrix[Rows[-1] + 1][Cols[-1] + StrandLengths[-1]] = '/'
    else:
        Matrix[Rows[-1] + 1][Cols[-2] - 1] = '\\'

    # create strings from matrix rows
    output = []
    for i in Matrix:
        string = ''
        for j in i:
            string = string + str(j)
        output.append(string)

    # transpose matrix and remove empty strings
    output_transpose = []
    for i in (range(L)[::-1]):
        string = ''
        for j in range(len(output)):
            n = output[j][i]
            string = string + n
        if string != sep*len(string): output_transpose.append(string)

    return output_transpose


def StructureToString(structure):
    """
    Function transform the list describing structure of serpentine into string.

    Input [['Q','Q','Q','Q','Q'],['G','G','G'],['Q','Q','Q','Q','Q']]

    Output 'QQQQQgggQQQQQ'
    """
    string = ''
    for i in range(len(structure)):
        if i % 2 == 0:
            for aa in structure[i]: string += aa
        else:
            for aa in structure[i]: string += aa.lower()
    return string




def ExportSerpintines(Serpentines, Consensus, FileName = 'SerpentinesStructure.txt', sep = ' ', SortByConsensus = True):

    """ MODIFY!!!
    Export serpentines structures and general information to .txt file.
    serpentines - dictionary of serpentines
    FileName - name of output File.

    Input
    {'6-29': [[0.479, 0.502], [['T', 'A', 'N'], ['N', 'T', 'C', 'A', 'T', 'Q'], ['L', 'A'], ['N', 'F', 'L'], ['V', 'Q', 'Q']],
     [3, 19], 1.0, [0, 0], [[3, 4, 5], [6, 7, 8, 9, 10, 11], [12, 13], [14, 15, 16], [17, 18, 19]], ['6 Res Arc 1', 'PPL'],
     0.9951159997383305, 0.99999594146826, [6, 29], 0.4905, 0.4881024168844599, 1, [0, 0]],
     ....}

     Output
    Arches:6-29
    Path[6, 29]
    ArchScores = [0.479, 0.502] Arches Types = ['6 Res Arc 1', 'PPL'] Positions: 3-19
    CompactScore1 = 0.995 CompactScore2 = 1.0 MAScore = 0.49 TwistScore = 1.0
    BSScore = 0.488
    Structure string: TANntcatqLAnflVQQ-3-19

     tcat    /
    n    q  Q
     N  L    Q
    A    A  V
     T  n    l
    /     f
    """

    OutputFile = open(FileName, 'w') # Create empty file

    if SortByConsensus:
        # sort Serpentines with Perimeter Score
        tmpDict = []
        for i in Serpentines:
            #print 'Structure ', Serpentines[i][5]
            #print 'Consensus ', Consensus
            Index = CompareToConsensus(Consensus, Serpentines[i][5])
            tup = (i, round(Index,3),round(Serpentines[i][11],3))
            tmpDict.append(tup)

        tmpTouple = sorted(tmpDict, key = operator.itemgetter(1,2), reverse=True) #(key=lambda item: item[1,2]) #sorted(student_tuples, key=itemgetter(1,2))
        #print 'Touple ', tmpTouple

    else:
        # sort Serpentines with Perimeter Score
        tmpDict = {}
        for i in Serpentines: tmpDict.update({i:Serpentines[i][11]})
        tmpTouple = list(tmpDict.items())
        tmpTouple.sort(key=lambda item: item[1], reverse = True)

    # cycle to print informtion about all serpentines in dictionary
    for j in tmpTouple:
        i = j[0]
        # create a list of strings corresponding to printed serpentine
        StringList = PrintSerpentine(Serpentines[i][1])

        # count SumScore
        sumscore = 0
        for j in Serpentines[i][0]: sumscore = sumscore + j

        # write information about Serpentine to file
        OutputFile.write('Arches:' + str(i) + '\n')

        OutputFile.write('Path' + str(Serpentines[i][9]) + '\n' +
                         'ArchScores = ' + str(Serpentines[i][0])+ ' Arches Types = ' + str(Serpentines[i][6]) +
                         ' Positions: ' + str(Serpentines[i][2][0]) + '-' + str(Serpentines[i][2][1]) + '\n' +
                         'CompactScore1 = ' + str(round(Serpentines[i][7],3)) +
                         ' CompactScore2 = ' + str(round(Serpentines[i][3],3)) +
                         ' MAScore = ' + str(round(Serpentines[i][10],3)) +
                         ' TwistScore = ' + str(round(Serpentines[i][8],3)) + '\n' +
                         'BSScore = ' + str(round(Serpentines[i][11],3)) + '\n' +
                         'Structure string: ' + StructureToString(Serpentines[i][1])+ '-' +str (Serpentines[i][2][0]) + '-' + str(Serpentines[i][2][1]) + '\n')
        OutputFile.write(' \n') # write empty line
        # write Serpentine 2D structure
        for j in StringList: OutputFile.write(j + '\n')
        OutputFile.write(' \n')



def ExportSerpintinesTableShort(Serpentines, FileName='SerpentinesTableShort.txt'):
    """
    Write sort table with information of serpentines

    Input
    {'6-29': [[0.479, 0.502], [['T', 'A', 'N'], ['N', 'T', 'C', 'A', 'T', 'Q'], ['L', 'A'], ['N', 'F', 'L'], ['V', 'Q', 'Q']],
     [3, 19], 1.0, [0, 0], [[3, 4, 5], [6, 7, 8, 9, 10, 11], [12, 13], [14, 15, 16], [17, 18, 19]], ['6 Res Arc 1', 'PPL'],
     0.9951159997383305, 0.99999594146826, [6, 29], 0.4905, 0.4881024168844599, 1, [0, 0]],
     ....}

    Output
    Name	BSScore	MAScore	CompactScore1	CompactScore2	TwistScore	Start	Structure	End
    6-29	 0.488	 0.49	 0.995	            1.0	            1.0	    3	TANntcatqLAnflVQQ	19
    5-29	 0.505	 0.505	 1.0	            1.0	            1.0	    1	QNTANntcaTQLAnflVQQHS	21
    13-29	 0.55	 0.552	 0.995	            1.0	            1.0	    4	ANNtcatqLAnflVQQ	19
    1-29	 0.563	 0.566	 0.995	            1.0	            1.0	    5	NNTcatqLAnflVQQ	19

    """
    File = open(FileName, 'wb')
    Head = 'Name' + '\t' \
          + 'BSScore' + '\t' \
          + 'MAScore' + '\t' \
          + 'CompactScore1' + '\t'\
          + 'CompactScore2' + '\t' \
          + 'TwistScore' +'\t' \
          + 'Start' + '\t' \
          + 'Structure' + '\t' \
          + 'End' + '\t' \
          + '\n'

    File.write(Head)
    for i in Serpentines:
        string = str(i) + '\t' \
                 + str(round(Serpentines[i][11],3)) +'\t' \
                 + str(round(Serpentines[i][10],3)) + '\t' \
                 + str(round(Serpentines[i][7],3)) + '\t' \
                 + str(round(Serpentines[i][3],3)) + '\t' \
                 + str(round(Serpentines[i][8],3)) + '\t'\
                 + str(Serpentines[i][2][0]) + '\t' \
                 + StructureToString(Serpentines[i][1]) + '\t'\
                 + str(Serpentines[i][2][1]) + '\t' \
                 + '\n'
        File.write(string)
    File.close()

