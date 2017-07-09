BetaSerpentine RELEASE NOTES
====================

BetaSerpentine Version 1.0 2017/07/01

by Bondarev S.A., Bondareva O.V., Zhouravleva G.A. and Kajava A.V.

Contacts: andrey.kajava@crbm.cnrs.fr, s.bondarev@spbu.ru or stanislavspbgu@gmail.com

Numerous studies suggest that amyloid forming regions of more than 30 residues evoke structural arrangements consisting of several adjacent b-arches called b-serpentines or superpleated b-structure (Kajava A.V. et al., 2010, FASEB J). The BetaSerpentine program searches for possible b-serpentines within a protein and output them in order of preference starting from the highest score. 

N.B. BetaSerpentine program is not aimed at the detection of the amyloidogenic regions. If you want, first, to predict the amyloidogenicity of your protein, use ArchCandy. 

====================

The BetaSerpentine is written in Python v.2.7.6 and provided as a set of scripts, no additional packages are required. The program takes specific text files containing information about b-arches, predicted by ArchCandy, as an input. The BetaSerpentines is provided with a several examples of such files for well-known amyloidogenic proteins (folder "TestAmyloidExamples"). However for analisys of other proteins you also need the ArchCandy program to predict b-arches. OR you can use an online version of the program on dali.crbm,cnrs.fr, where you will need only protein sequence for the analisys. An we reccomend this online tool for common usage.

The minimal length of protein sequence is 18 amino acids. In highly amyloidogenic and low complexity sequences the program can find thousands of β-serpentines . Analysis of proteins with very big numbers of β-arches (for example, Rnq1) with default parameters will be interrupted because it will take several hours. To analyse them you should increase the β-arches threshold or try to analyse a part of the protein.

====================
Beta-serpentine is cross-platfom programm, it was tested on Windows and on Linux OS.

TO RUN the BetaSerpentine: 
run the main script (BetaSerpentine.py) using python v.2 in the command line and specify next parameters:
- Filename - name of the file with list of b-arches (required);
- OutputFolder - path to folder for output files (required);
- SortByConsensus - True or False - if True output arrangements will be sorted by similarity to consensus, if False - only by score (optional, default value True);
- InputThreshold - from 0.0 to 1.0 - a treshold for b-arches which will be analised (optional, default value 0.0);
- OutputThreshold - from 0.0 to 1.0 - a treshold for b-serpentines which will be presented in the output (optional, default value 0.2).

	python BetaSerpentine.py Filename OutputFolder SortByConsensus InputThreshold OutputThreshold

INPUT FILE: the file with b-arches, predicted by ArchCandy. Several examples are provided.

EXAMPLE RUN: 

	python BetaSerpentine.py TestAmyloidExamples/test_arches.txt TestAmyloidExamples/ True 0 0.2

or the same with default parameters

	python BetaSerpentine.py TestAmyloidExamples/test_arches.txt TestAmyloidExamples/

====================

INTERPRETATION OF THE OUTPUT:

BetaSerpentine after each run produces 3 files.
The Log file("ProteinName_log.txt") contaiins general information: timing, number of input b-arches, number of b-serpentines etc.

#and probability of amino acid be included in serpentine

The file "ProteinName_serpentines.txt" contains schemes of all serpentines, their location and corresponding scores. In this file b-serpentines are ranged according as specified: by similarity to consensus or only by score. "Structure" parameters represent consensus structure for an analysed sequence. "S", "A", and "-" mean β-strand, β-arc, and undefined structure, respectively. 
Individual scores of b-serpentines its localisation and structure are also written as a table in the file "ProteinName_table_short.txt".

====================

Please see the LICENSE file for the license terms for the software. It is
basically free for academic users, but a license fee applies to commercial
users. 

====================

THE PUBLICATION OF RESEARCH USING BetaSerpentine MUST INCLUDE AN APPROPRIATE
CITATION TO THE METHOD:

Bondarev S.A., Bondareva O.V., Zhouravleva G.A. and Kajava A.V.
“BetaSerpentine”: a bioinformatics tool for reconstruction of amyloid structures. in prep.


