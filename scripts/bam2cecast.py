#!/usr/bin/env python3
"""Author:  Cesare de Filippo
Contact: casare_filippo(at)eva.mpg.de
Generate the input file for cecast.
"""

##################################################
# Libraries required for the parsing arguments.
import re
import os
import copy
import random
from itertools import repeat
from signal import signal, SIGPIPE, SIG_DFL
import pysam
from os.path import exists
import argparse
import textwrap as _textwrap


class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        return _textwrap.wrap(text, 90)

##################################################


parser = argparse.ArgumentParser(
    description="This is a script to return the lineage assignment table(s) or to sample allele(s) at the provided informative sites.'\n",
    formatter_class=LineWrapRawTextHelpFormatter)

parser.add_argument('input_file', nargs='?')
parser.add_argument('-s', dest='sites',
                    help='File with informative sites in bed format with three extra columns: (1) ancestral and (2) derived alleles, (3) lineage name. It must be sorted, bgziped and tabixed.', required=True)
parser.add_argument('-b', dest='BlockSize',
                    help='The size of the blocks in Mb or the name of file with the coordinates of the blocks.', type=str)

parser.add_argument('-l', dest='minLength',
                    help='Minimum sequence length to be considered. Default is 20 bp.', default=20, type=int)

parser.add_argument('-L', dest='maxLength',
                    help='Maximum sequence length to be considered. Default is 100 bp.', default=100, type=int)

parser.add_argument('-q', dest='BaseQual',
                    help='Minimum base quality (BQ) value. Default is 30.',   type=int, default=30)

parser.add_argument('-m', dest='MapQual',
                    help='Minimum mapping quality (MQ) value. Default is 1.', type=int, default=1)

parser.add_argument('-d', dest='deam',
                    help="Consider putatively deaminated sequences at the last i,j (5',3') terminal positions.", type=str, default="0,0", required=False)

parser.add_argument('-D', dest='DoubleStrand',
                    help="The library is double stranded. Single stranded is the default.", default=False, required=False, action='store_true')

parser.add_argument('-I', dest='Indels',
                    help="Remove insertions and deletions (indels).", default=False, required=False, action='store_true')

parser.add_argument('-S', dest='StrandFilter',
                    help="Apply the strand orientation filter. This works only for single-stranded libraries.", default=False, required=False, action='store_true')

parser.add_argument('-F', dest='StrandFilterRelax',
                    help="Apply the strand orientation filter with a relaxed criteria, which is to apply it only at transitions sites. This works only for single-stranded libraries.", default=False, required=False, action='store_true')

parser.add_argument('-T', dest='Transversions',
                    help="Use only transversions from the informative sites (-s option).", default=False, required=False, action='store_true')

parser.add_argument('-O', dest='OutputAlleles',
                    help="Print the sampled allele(s) for each of the covered site.", default=False, required=False, action='store_true')

parser.add_argument('-C', dest='OutputCoverage',
                    help="Print also the coverage for each site.", default=False, required=False, action='store_true')

parser.add_argument('-R', dest='RandomSampling',
                    help="Random sample one read when multiple reads are spanning an informative site.", default=False, required=False, action='store_true')

parser.add_argument('-E', dest='ExcludeOtherAlleles',
                    help="Eclude the other alleles, i.e. only sample from the two given alleles.", default=False, required=False, action='store_true')


args = parser.parse_args()

##################################################
# Other libraries required


# To avoid the 'BrokenPipeError: [Errno 32] Broken pipe'
signal(SIGPIPE, SIG_DFL)

##################################################
# Define the parameters

input_file = args.input_file
if not exists(input_file):
    print("Input.bam file doesn't exist.")
    exit()

infosites = args.sites
if not exists(infosites):
    print("The table with informative sites doesn't exist.")
    exit()

# The blocks created with BlockSize (-b option)
BlockSize = args.BlockSize

# The lengths of the chromosomes
chrom_length = {
    '1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276,
    '5': 180915260, '6': 171115067, '7': 159138663, '8': 146364022,
    '9': 141213431, '10': 135534747, '11': 135006516, '12': 133851895,
    '13': 115169878, '14': 107349540, '15': 102531392, '16': 90354754,
    '17': 81195210, '18': 78077248, '19': 59128983, '20': 63025520,
    '21': 48129895, '22': 51304566, 'X': 155270560, 'Y': 59373566, 'MT': 16569}

# To be used to call the dictionary and avoid the random sorting of it for printing the outputs.
Chromosomes = [str(k) for k in range(1, 23)] + ['X', 'Y', 'MT']
blocks = []

if BlockSize is not None:
    OutputTableBlocks = True
    if os.path.isfile(BlockSize):
        with open(BlockSize, "r") as bins:
            for line in bins:
                line = line.rstrip().split('\t')
                block_name = str(line[0])+':'+str(line[1])+'-'+str(line[2])
                blocks += [block_name]
    else:
        try:
            BlockSize = round(float(BlockSize)*1e6)
            # It produces the bins of all blocks across chromosomes.
            for k in Chromosomes:
                # The start coordinates
                b1 = [b+1 for b in range(0, chrom_length[k], BlockSize)]
                # The end coordinates
                b2 = [b1[b]-1 for b in range(1, len(b1))] + [chrom_length[k]]
                for b in range(0, len(b1)):
                    block_name = k+':'+str(b1[b])+'-'+str(b2[b])
                    blocks += [block_name]
        except ValueError:
            print("ERROR: '-b {}'. '-b' should be a number specifying the block size or the path to a file with coordinates of the blocks.".format(args.BlockSize))
            exit()
else:
    OutputTableBlocks = False
    for k in Chromosomes:
        block_name = k+':'+'1'+'-'+str(chrom_length[k])
        blocks += [block_name]


##################
# Other otions: ##
##################

OutputAlleles = args.OutputAlleles
OutputCoverage = args.OutputCoverage
ExcludeOtherAlleles = args.ExcludeOtherAlleles
BQ_cutoff = args.BaseQual
MQ_cutoff = args.MapQual
minLength = args.minLength
maxLength = args.maxLength
rm_Indels = args.Indels
DoubleStrand = args.DoubleStrand
StrandFilter = args.StrandFilter
StrandFilterRelax = args.StrandFilterRelax
if StrandFilter:  # make sure that only the most strict is applied
    StrandFilterRelax = False
RandomSampling = args.RandomSampling

Transitions = ['CT', 'TC', 'GA', 'AG']

# -d: Deamination filter
terminal_deam = args.deam.split(',')
deam_filter_skip = True
if (terminal_deam[0] != '0' or terminal_deam[1] != '0'):
    deam_filter_skip = False

# Check the consistency of the -d option.
try:
    int(terminal_deam[0])
    int(terminal_deam[1])
except ValueError:
    print("-d {} incosistent. It must be two integers separeted by comma as follow:\n-d 1,1".format(args.deam))
    exit()


#########################
# Define some functions #
#########################

def alnseq(refseq, myseq, CIGAR, basequalities, start):
    """Align the referece sequence (refseq) and the sequence of interest (myseq), and return the positions and basequelities."""
    # a1 = refseq; a2 = myseq.
    a1 = a2 = ''
    # i1 and i2 will be the ith element modified for indels for a1 and a2 respectively.
    i1 = i2 = 0
    # The positions and the basequalities (bq) to be 'adjusted' in case of indels.
    positions = []
    bq = []
    pos = int(start)-1
    for (cigar_op, cigar_len) in CIGAR:
        if cigar_op == 0:
            a1 += refseq[i1:(cigar_len+i1)]
            a2 += myseq[i2:(cigar_len+i2)]
            bq += basequalities[i2:(cigar_len+i2)]
            positions += list(range(pos+1, pos+cigar_len+1))
            i1 += cigar_len
            i2 += cigar_len
        if cigar_op == 1:  # It's an insertion relative to the reference
            a2 += myseq[i2:(cigar_len+i2)]
            bq += basequalities[i2:(cigar_len+i2)]
            positions += list(range(pos+1, pos+cigar_len+1))
            i2 += cigar_len
            a1 += '-' * cigar_len
        if cigar_op == 2:  # It's a deletion relative to the reference
            a2 += '-' * cigar_len
            # 90 is an arbitrary base quality for the reference sequence.
            bq += [90] * cigar_len
            positions += [pos+j for j in range(1, cigar_len+1)]
            a1 += refseq[i1:(cigar_len+i1)]
            i1 += cigar_len
        pos = positions[-1]
    # must return refseq, myseq, pos, bq
    return (a1, a2, positions, bq)


def is_deaminated(refseq, myseq, terminal_positions=['1', '1'], isReverse=False, isDoubleStranded=False):
    """Determine whether a sequence is deaminated conditioned on the last 'i,j' terminal positions"""
    refseq = refseq.upper()
    output = False
    nt = 'CT'
    if (isDoubleStranded == False):
        if (isReverse):
            nt = 'GA'
        for p in range(0, int(terminal_positions[0])):
            if (refseq[p]+myseq[p] in nt):
                output = True
        if (output == False):
            for p in range(1, int(terminal_positions[1])+1):
                if(refseq[-p]+myseq[-p] in nt):
                    output = True
    else:
        for p in range(0, int(terminal_positions[0])):
            if (refseq[p]+myseq[p] in 'CT'):
                output = True
        if (output == False):
            for p in range(1, int(terminal_positions[1])+1):
                if(refseq[-p]+myseq[-p] in 'GT'):
                    output = True
    return(output)


def is_StrandFilter(Reference, Modified, isReverse=False, Relax=False):
    """ Whether the allele of the read passed the strand-orientation filter. """
    pass_SF = True
    if Relax:
        if Reference+Modified in ['CT', 'TC'] and isReverse == False:
            pass_SF = False
        elif Reference+Modified in ['GA', 'AG'] and isReverse == True:
            pass_SF = False
    else:
        if 'C' in Reference+Modified and isReverse == False:
            pass_SF = False
        elif 'G' in Reference+Modified and isReverse == True:
            pass_SF = False
    return (pass_SF)


def which_allele(Allele, Ancestral, Derived):
    """ Determine the state (Ancestral, Derived or Other) """
    ret_type = 'oth'
    Allele = Allele.upper()
    # Reference and Modified (i.e. alternative) are any types of alleles.
    if Derived.upper() == Allele:
        ret_type = 'der'
    elif Ancestral.upper() == Allele:
        ret_type = 'anc'
    return (ret_type)


# Functions to print the table and to sum over them in case no blocks are required.
def print_table(table, label="Lineage", alleles=['a', 'd', 'o'], rm_empty_lines=True):
    if rm_empty_lines:  # removes lineages without counts in the test ('0')
        table = rmlineage(table)
    if len(table) > 0:
        n = len(table[list(table.keys())[0]])
        human = ['X', 'Dai', 'Fre', 'Han', 'Man', 'Mbu',
                 'Pap', 'San', 'Sar', 'Yor', 'Kar'][:n]
        header = ''
        for j in human:
            for i in alleles:
                header += "\t"+j+"_"+i
        print('#'+label+header)
        for lin in sorted(table.keys()):
            line = lin
            for i in table[lin].values():
                line += '\t'+'\t'.join(str(s) for s in i)
            print(line)


def rmlineage(tab):
    """ Remove lineage without any count. """
    tab_out = dict()
    for lin in tab:
        if tab[lin]['0'].count(0) < 3:
            tab_out[lin] = tab[lin]
    return(tab_out)


def merge_table(tables, tables1=None, blocks=None):
    """ Merge for each lineage the table of the test with the other samples (tables1) """
    if blocks is None:
        blocks = sorted(tables.keys())
    tab_out = dict()
    for block in blocks:
        tab_block = dict()
        for lin in tables[block].keys():
            tab = tables[block][lin]
            tmp = {"0": list(tab.values())}
            # merge the tables
            if tables1 is not None:
                tab1 = copy.deepcopy(tables1[block][lin])
                for i in range(0, int(len(tab1['anc']))):
                    tmp[str(i+1)] = [tab1['anc'][i],
                                     tab1['der'][i], tab1['oth'][i]]
            tab_block[lin] = tmp
        tab_out[block] = tab_block
    return(tab_out)


def sum_tables(tables, blocks=None):
    if blocks is None:
        blocks = tables.keys()
    tab = dict()
    for block in blocks:
        for lin in tables[block].keys():
            tmp = copy.deepcopy(tables[block][lin])
            # merge the tables
            if lin not in tab:
                tab[lin] = tmp
            else:
                tmp0 = copy.deepcopy(tab[lin])
                for i in tmp:  # it assumes tmp0 and tmp1 have equal structure
                    for j in range(0, 2):
                        tmp0[i][j] += tmp[i][j]
                tab[lin] = tmp0
    return(tab)


def main():
    """ The main function. """
    with pysam.AlignmentFile(input_file, "rb", check_sq=False) as samfile, pysam.TabixFile(infosites) as tabixfile:
        output_table = dict()
        output_table_others = dict()
        # the template for the lineage assignment (or classes of sites) table.
        for block in blocks:
            tmp_tab = dict()
            tmp_tab_others = dict()
            chrom = block.split(":")[0]
            coordinates = block.split(":")[1].split('-')
            (start, end) = (int(coordinates[0]), int(coordinates[1]))
            output_table[block] = dict()
            output_table_others[block] = dict()
            infodict = {}
            try:
                for s in tabixfile.fetch(chrom, start, end, parser=pysam.asTuple()):
                    siteclass = s[5]  # or simply the Lineage
                    passFilter = True
                    if StrandFilter:
                        passFilter = not s[3] + s[4] in ['CG', 'GC']
                    if args.Transversions:
                        passFilter = not s[3] + s[4] in Transitions
                    if passFilter:
                        infodict[str(s[1])] = {
                            "anc": s[3], "der": s[4], "Lineage": siteclass, "pileup": '', "samples": s[6:]}
                    if siteclass not in tmp_tab:
                        tmp_tab[siteclass] = {"anc": 0, "der": 0, "oth": 0}
                        if len(s[6:]) > 0:
                            ac_tmp = {}
                            for a in ['anc', 'der', 'oth']:
                                ac_tmp[a] = list(repeat(0, len(s[6:])))
                            tmp_tab_others[siteclass] = copy.deepcopy(ac_tmp)
                # To prevent tmp_* objects to be overwritten by output_table
                output_table[block] = copy.deepcopy(tmp_tab)
                output_table_others[block] = copy.deepcopy(tmp_tab_others)
            except ValueError:
                del output_table[block]
                del output_table_others[block]
                pass
            for read in samfile.fetch(chrom, start, end):
                Cigar = read.cigarstring
                if (rm_Indels):
                    if 'I' in Cigar or 'D' in Cigar:
                        continue
                pos = read.get_reference_positions(full_length=False)
                L = len(pos)
                # Filter out softclip, hardclip, for MapQuality cutoff and sequence length
                if 'S' not in Cigar and 'H' not in Cigar and read.mapping_quality >= MQ_cutoff and L >= minLength and L <= maxLength:
                    posInTable = []
                    for p in pos:
                        if(str(p) in infodict.keys()):
                            posInTable.append(p)
                    if len(posInTable) > 0:
                        refseq = read.get_reference_sequence()
                        myseq = read.query_sequence
                        bq = read.query_qualities
                        if len(myseq) != len(refseq):
                            (refseq, myseq, pos, bq) = alnseq(
                                refseq, myseq, CIGAR=read.cigartuples, basequalities=bq, start=pos[0])
                        # Maybe the is_deaminated can be nested to speed up
                        if deam_filter_skip or is_deaminated(refseq, myseq, terminal_deam, read.is_reverse, isDoubleStranded=DoubleStrand):
                            for p in posInTable:
                                passStrandFilter = True
                                tab = infodict[str(p)]
                                Allele = myseq[pos.index(p)]
                                BQ = bq[pos.index(p)]
                                if StrandFilter or StrandFilterRelax:
                                    passStrandFilter = is_StrandFilter(
                                        Reference=tab['anc'], Modified=tab['der'], isReverse=read.is_reverse, Relax=StrandFilterRelax)
                                if passStrandFilter and BQ >= BQ_cutoff:
                                    if read.is_reverse:
                                        Allele = Allele.lower()
                                    infodict[str(p)]['pileup'] += Allele
            for i in sorted(infodict.keys(), key=int):
                pileup = infodict[i]['pileup']
                if len(pileup) > 0:
                    (Ancestral, Derived, Lineage) = (
                        infodict[i]['anc'], infodict[i]['der'], infodict[i]['Lineage'])
                    A = Ancestral.upper()
                    a = Ancestral.lower()
                    D = Derived.upper()
                    d = Derived.lower()
                    if ExcludeOtherAlleles:
                        pileup = ''.join(re.findall(
                            r"["+A+","+D+","+a+","+d+"]", pileup))
                    if len(pileup) > 0:
                        if RandomSampling:
                            pileup = random.sample(pileup, 1)[0]
                        for r in pileup.upper():
                            ret_type = which_allele(r, Ancestral, Derived)
                            output_table[block][Lineage][ret_type] += 1
                        if len(output_table_others) > 0:
                            alleles = infodict[i]['samples']
                            for j in range(0, len(alleles)):
                                ret_type = which_allele(
                                    alleles[j], Ancestral, Derived)
                                output_table_others[block][Lineage][ret_type][j] += 1
                        if (OutputAlleles):
                            if (OutputCoverage):
                                print(chrom, int(i)+1, pileup,
                                      len(pileup), sep='\t')
                            else:
                                print(chrom, int(i)+1, pileup, sep='\t')
    if not OutputAlleles:
        output = merge_table(output_table, output_table_others)
        if OutputTableBlocks:
            for block in blocks:
                chrom = block.split(":")[0]
                coordinates = block.split(":")[1].split('-')
                o = output.get(block, 0)
                if o != 0:
                    print_table(o, label=block)
        else:
            print_table(sum_tables(output))


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        exit()
    exit()
