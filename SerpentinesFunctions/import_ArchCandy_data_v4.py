def importSequence(Filename):
    File = open(Filename, 'rb')
    Sequence = ''
    for line in File:
        line = line.strip()
        if line[0] != '>':
            sequence += line
    return Sequence


def findCoord(header):
    """
    Function return coordinate of each element from list of description in header
    """
    coords = [0]
    for i in ['Sequence', 'Score', 'Arc Type', 'Start', 'Stop', 'Length']:
        coords.append(header.index(i))
    return coords

def delSpaces(string):
    """
    Delete free spaces in given string, returns string without free space.
    """
    while string[0] == ' ': string = string[1:]
    while string[-1] == ' ': string = string[:-1]

    return string


def inputArches(FileName):
    """
    Import information about each of arch data from given ArchCandy file
    Return protein sequence, protein name and list of list with information about arch
    """
    File = open(FileName)
    output = []
    ProteinSequence = ''
    startImport = False
    isProteinSeq = False
    nArch = 0
    for line in File:
        if 'Protein sequence :' in line:
            ProteinSequence = line.replace('Protein sequence : ', '')

        if 'Gene name' in line:
            ProtName = line.replace('Gene name : ', '')

        if 'CANDIDATES DETAILS' in line:
            startImport = True

        if startImport and line[0] != '#' and len(line) > 1:
            line = line.strip('\n')

            if line[:6] == 'Number':
                coords = findCoord(line)
            else:
                nArch += 1
                Name = str(nArch)
                Score = delSpaces(line[coords[2]:coords[3]])
                ArchType = delSpaces(line[coords[3]:coords[4]])
                Start = delSpaces(line[coords[4]:coords[5]])
                End = delSpaces(line[coords[5]:coords[6]])
                output.append(Name + '\t' + Score + '\t' + ArchType + '\t' + Start + '\t' + End)

    File.close()
    return [ProteinSequence, output, ProtName[:-1]]



def structure(Type, start, end):
    """
    Function returns string defining structure of Arch.
    Type - arch type from Archcandy, start and end - coordinate of start and end in protein

    Output definitions:
    A - amino acid in Arc region
    S - amino acid in strand region (outside of Arch)
    s - amino acid in strand region (inside of Arch).

    Example: structure('BEPL',18)
    sSsSsSsAAAAsSsSsSs
    """
    length = end - start + 1
    dict = {'5 Res Arc': [((length - 5)/2), 5], '6 Res Arc 1': [((length - 5)/2),6],
            '6 Res Arc 2': [((length - 5)/2), 6], 'BEPL': [(((length - 4)/2)), 4],
            'GBEB': [(((length - 4)/2)), 4],
            'GBPL': [(((length - 4)/2)), 4], 'PPL': [(((length - 3)/2)), 3]}

    strandL = dict[Type][0]
    arcL = dict[Type][1]
    output = [range(start, start + strandL), range(start + strandL, start + strandL + arcL), range(start + strandL +
                                                                                                   arcL, end+1)]
    return output


def ArchesImport(table_of_arches, InputThreshold):
    """
    Function returns dictionary of Arches containing:
    arch structure, its score and Arch number in ArchCandy output.

    Example.
    Input table_of_arches:
    01	 0.644	 5 Res Arc	05-35
    02	 0.577	 5 Res Arc	05-31
    ...

    Output Dictionary:
    "01" : 0.644, [[5,6,7,...19],[20,21,22,23,24],[25,...34,35]], '01'
    "02" : 0.577, [[5,6,7,...17],[21,22,23,24,25],[26,...30,31]], '02'
    ...

    """
    Arches = {}  # empty dictionary for Arches

    import re

    for Line in table_of_arches:  # import information about Arches and process it
        Line = Line.strip('\n')  # delete end of line
        List = Line.split('\t')

        score = List[1]
        Type = List[2]
        Start = List[3]
        End = List[4]

        Structure = structure(Type, int(Start), int(End))
        L = 2 * len(Structure[0])
        LScore = 1 - (0.0003462*(2*L-7-45)**2)
        Score = float(score) / LScore  # correct score for length score
        Score = [round(float(Score), 3)]

        ArchName = List[0]

        list = [Score, Structure, ArchName, Type, Start, End]

        if list[0][0] >= InputThreshold:
            Arches.update({int(ArchName): list})  # assemble information into dictionary

    return Arches


def rm_ident_arches2(ListOfString):
    """
    Function removes identical arches in list of arches
    If two or more arches have the same coordinates, only one with bigger score is taken for further analysis

    Input
    1 ['1', '0.504', 'GBPL', '1', '18']
    2 ['2', '0.454', '5 Res Arc', '1', '19']
    3 ['3', '0.346', '6 Res Arc 1', '1', '20']
    4 ['4', '0.346', '6 Res Arc 2', '1', '20']
    5 ['5', '0.328', 'GBPL', '1', '14']
    ...

    Output
    1	0.504	GBPL	1	18
    2	0.454	5 Res Arc	1	19
    3	0.346	6 Res Arc 1	1	20
    5	0.328	GBPL	1	14
    ...
    """
    Arches = {}
    for line in ListOfString:
        line = line.strip('\n')
        Arch_info = line.split('\t')
        Arches.update({int(Arch_info[0]): Arch_info})

    output = []
    removed = []
    for i in range(1,len(Arches)):
        Arch1 = Arches[i]
        Arch2 = Arches[i+1]
        positions1 = Arch1[-1]
        positions2 = Arch2[-1]
        string = ''
        if positions1 == positions2 and len(Arch1[2]) == len(Arch2[2]) and i not in removed:
            if float(Arch1[1]) >= float(Arch2[1]):
                for j in Arch1:
                    string += j + '\t'
                output.append(string[:-1] + '\n')
                removed.append(i+1)
        elif i not in removed:
            for j in Arch1:
                string += j + '\t'
            output.append(string[:-1] + '\n')

    if len(Arches) not in removed:
        string = ''
        for j in Arches[len(Arches)]:
            string += j + '\t'
        output.append(string[:-1] + '\n')
    return output


def ArchesDestributions(Arches, FileName):
    """
    Count number of arch per 1 amino acid in protein and write it in the output file

    Input
    1 [[0.63], [[1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15, 16, 17, 18]], '1', 'GBPL', '1', '18']
    2 [[0.567], [[1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12], [13, 14, 15, 16, 17, 18, 19]], '2', '5 Res Arc', '1', '19']
    3 [[0.432], [[1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13], [14, 15, 16, 17, 18, 19, 20]], '3',
                                                                                            '6 Res Arc 1', '1', '20']


    Output in text file
    X	NArches
    0	6
    1	12
    2	12
    3	13
    4	17
    5	19
    """
    keys = Arches.keys()
    y = [0] * int(Arches[keys[0]][-1])

    for k in keys:
        Start = int(Arches[k][-2])
        End = int(Arches[k][-1])
        if End >= len(y):
            delta = End - len(y) + 1
            y = y + [0] * delta
        for p in range(Start-1, End-1):
            y[p] += 1

    '''File = open(FileName + '_ArchesDistrib.txt', 'wb')
    File.write('X\tNArches\n')
    for i in range(len(y)):
        File.write(str(i) + '\t' + str(y[i]) + '\n')
    File.close()'''

    return max(y)

def RewriteFile(FileName, Header, Replace1th=True):
    """
    Create new header for given protein file
    """
    OldFile = open(FileName, 'rb')
    Strings = []
    for line in OldFile:
        Strings.append(line)
    OldFile.close()
    if Replace1th: Strings[0] = Header
    else: Strings = [Header] + Strings

    NewFile = open(FileName, 'wb')
    for string in Strings:
        NewFile.write(string)
    NewFile.close()

def TakeFirstLine(Filename):
    """
    Gives the first line on file
    """
    File = open(Filename, 'rb')
    count = 0
    for line in File:
        if count == 0:
            string = line
            count = 1
        else: break
    return string

def CorrectFastaForAC(FileName):
    """
    Create fasta-file for ArchCandy if data not in UniProt fasta-format data

    Input
    >protein name
    QNTANNTCATQLANFLVQQHS

    Output
    >sp|000000|YFAP OS=Your Favorite Species GN=YFAP PE=1 SV=1
    QNTANNTCATQLANFLVQQHS
    """
    default_header = '>sp|000000|YFAP OS=Your Favorite Species GN=YFAP PE=1 SV=1' + '\n'
    FirstString = TakeFirstLine(FileName)
    if '>' not in FirstString:
        RewriteFile(FileName, default_header, Replace1th=False)
    if FirstString.count('|') < 2 or FirstString.count(' ') < 1:
        RewriteFile(FileName, default_header)
