#!/usr/bin/env python3

# Author: Cesare de Filippo

# Script to output chr, pos, allele(s) from vcf files.
# If the vcf contains more than one individual all of them will be considered.
# But this will not work with -C (Counts) options. Or at least not implemented yet.
# The script can be run with compressed and uncompressed files, as well as called after piping:
# For example:
# vcf2table.py input.vcf.gz
# zcat input.vcf.gz | vcf2table.py


######################
# Libraries requires
import argparse
import os
import re
import sys
import random
import pysam
import gzip
import mimetypes
import copy
from itertools import repeat
from os.path import exists
from collections import OrderedDict


parser = argparse.ArgumentParser(prog='vcf2table.py', description="Script to output chr, pos, allele(s) from vcf files'\n",
                                 usage='\n<%(prog)s [options] input.vcf> or <cat input.vcf | %(prog)s [options]>')

parser.add_argument('filename', nargs='?',
                    help='input.vcf (can also be compressed).')

parser.add_argument('-s', dest='sites',
                    help='File with informative sites in bed format with three extra columns: (1) ancestral and (2) derived alleles, (3) lineage name. It must be sorted, bgziped and tabixed.', default=None, required=False)

parser.add_argument('-b', dest='BlockSize',
                    help='The size of the blocks in Mb or the name of file with the coordinates of the blocks. This will work only with -s (sites) and -C (Counts) options.', type=str)

parser.add_argument('-k', dest='chromosome',
                    help='The chromosome name.', type=str)

parser.add_argument('-A', dest='Ambiguity',
                    help='Print ambiguity code at heterozygous positions instead of a random allele. This will print ? if the alleles are Indels or not A,C,G,T.', default=False, required=False, action='store_true')

parser.add_argument('-B', dest='Bed',
                    help='Output in bed format.', default=False, required=False, action='store_true')

parser.add_argument('-C', dest='Count',
                    help='Return a table with counts of ancestral, derived, and other alleles for the sites classes defined in the 5th column of the informative sites (-s).', default=False, required=False, action='store_true')

parser.add_argument('-D', dest='Diploid',
                    help='Print both alleles (diploid).', default=False, required=False, action='store_true')

parser.add_argument('-I', dest='Indels',
                    help='Remove indels, i.e. drop positions when one allele is longer than one character.', default=False, required=False, action='store_true')

parser.add_argument('-P', dest='Paste',
                    help='Add or paste the allele to the sites sepcified with -s.', default=False, required=False, action='store_true')

parser.add_argument('-R', dest='PrintRef',
                    help='Print also the reference allele.', default=False, required=False, action='store_true')

parser.add_argument('-T', dest='Transversions',
                    help="Use only transversions from the informative sites (-s option).", default=False, required=False, action='store_true')


args = parser.parse_args()

infosites = args.sites
if (infosites is not None):
    if not exists(infosites):
        print("The table with informative sites doesn't exist.")
        exit()


# I put also AA, CC, GG, TT to avoid some if conditions
ambiguity_code = {'AC': "M", 'CA': "M", 'AG': "R", 'GA': "R", 'AT': "W", 'TA': "W", 'CG': "S",
                  'GC': "S", 'CT': "Y", 'TC': "Y", 'GT': "K", 'TG': "K", 'AA': "A", 'CC': "C", 'GG': "G", 'TT': "T"}


# The blocks created with BlockSize (-b option)
BlockSize = args.BlockSize

# The lengths of the chromosomes
chrom_length = {
    '1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276,
    '5': 180915260, '6': 171115067, '7': 159138663, '8': 146364022,
    '9': 141213431, '10': 135534747, '11': 135006516, '12': 133851895,
    '13': 115169878, '14': 107349540, '15': 102531392, '16': 90354754,
    '17': 81195210, '18': 78077248, '19': 59128983, '20': 63025520,
    '21': 48129895, '22': 51304566, 'X': 155270560}


if args.chromosome is None:
    args.chromosome = chrom_length.keys()
elif args.chromosome not in chrom_length:
    print("Please provide the correct chromosome name (-k).")
    exit()


#########################
# Define some functions #
#########################


def sort_list(mylist):
    pattern = r'(\d+):(\d+)-(\d+)'

    def sort_key(item):
        match = re.match(pattern, item)
        if match:
            chrom = int(match.group(1))
            start = int(match.group(2))
            end = int(match.group(3))
            return (chrom, start, end)
        else:
            return (float('inf'),)
    sorted_list = sorted(mylist, key=sort_key)
    return sorted_list


# To be used to call the dictionary and avoid the random sorting of it for printing the outputs.
def make_blocks(BlockSize, chromosome):
    blocks = []
    if BlockSize is not None:
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
                for k in chromosome:
                    # The start coordinates
                    b1 = [b+1 for b in range(0, chrom_length[k], BlockSize)]
                    # The end coordinates
                    b2 = [
                        b1[b]-1 for b in range(1, len(b1))] + [chrom_length[k]]
                    for b in range(0, len(b1)):
                        block_name = k+':'+str(b1[b])+'-'+str(b2[b])
                        blocks += [block_name]
            except ValueError:
                print("ERROR: '-b {}'. '-b' should be a number specifying the block size or the path to a file with coordinates of the blocks.".format(args.BlockSize))
                exit()
    else:
        for k in [chromosome]:
            block_name = str(k)+':'+'1'+'-'+str(chrom_length[k])
            blocks += [block_name]
    return(blocks)


def myprint(line, allele, bed=True, printRef=False):
    if (bed):
        out = [line[0], int(line[1])-1, line[1]]
    else:
        out = line[0:2]
    if printRef:
        out.append(line[3])
    print('\t'.join(map(str, out)), allele, sep="\t")


def gt2allele(line, diploid=False, Ambiguity=False, Indels=False):
    """ Main function to process the line in the vcf in ordet to sample one allele, or print other forms of vcf genotypes."""
    for i in range(9, len(line)):
        gt = line[i].split(":")[0]
        if (gt != './.'):  # Missing data will not be printed as output
            alleles = [line[3]]
            alleles.extend(line[4].split(","))
            a1 = alleles[int(gt[0])]
            a2 = alleles[int(gt[2])]
            if (Indels):
                if(len(a1) > 1 or len(a2) > 1):
                    break
            if (diploid):
                out_a = a1+'/'+a2
            else:
                if (Ambiguity):
                    aa = a1+a2
                    if (aa in ambiguity_code.keys()):
                        out_a = ambiguity_code[aa]
                    else:
                        out_a = 'N'
                else:
                    out_a = random.sample([a1, a2], 1)[0]
            return(out_a)


def which_allele(allele, Ancestral, Derived):
    """ Get the type of allele """
    a_type = 'oth'
    if allele == Derived:
        a_type = 'der'
    if allele == Ancestral:
        a_type = 'anc'
    return(a_type)


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


def rmlineage(tab):
    """ Remove lineage without any count. """
    tab_out = dict()
    for lin in tab:
        if tab[lin]['0'].count(0) < 3:
            tab_out[lin] = tab[lin]
    return(tab_out)


def print_table(table, label="Lineage", alleles=['a', 'd', 'o'], rm_empty_lines=True):
    """ Print the table and sum over them in case no blocks are required. """
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


def read_sites(infosites, blocks, Transversions=False):
    """ Read the informative sites per chromosome and return a dictionary and as list """
    infodict = dict()
    out_table = dict()
    out_table_others = dict()
    for b in blocks:
        (chr, start_pos, end_pos) = proc_block(b)
        out_table[b] = dict()
        out_table_others[b] = dict()
        with pysam.TabixFile(infosites) as tabixfile:
            for s in tabixfile.fetch(chr, int(start_pos), int(end_pos), parser=pysam.asTuple()):
                siteclass = s[5]  # or simply the Lineage
                passFilter = True
                if Transversions:
                    passFilter = not s[3] + s[4] in ['CT', 'TC', 'GA', 'AG']
                if passFilter:
                    id = str(s[0]+"_"+s[2])
                    infodict[id] = {
                                 "anc": s[3], "der": s[4], "Lineage": siteclass, "test": '', "samples": s[6:]}
                if siteclass not in out_table[b]:
                    out_table[b][siteclass] = {"anc": 0, "der": 0, "oth": 0}
                    if len(s[6:]) > 0:
                        ac_tmp = {}
                        for a in ['anc', 'der', 'oth']:
                            ac_tmp[a] = list(repeat(0, len(s[6:])))
                        out_table_others[b][siteclass] = copy.deepcopy(ac_tmp)
    return(infodict, out_table, out_table_others)


def proc_block(block):
    block = block.split(':')
    chr = block[0]
    (start, end) = block[1].split('-')
    return(chr, start, end)


def process_file(file):
    if(infosites is not None):
        blocks = sort_list(make_blocks(BlockSize, args.chromosome))
        blocks_copy = copy.deepcopy(blocks)
        # Create a set to store the completed blocks
        completed_blocks = set()
        # The informative sites and the output tables
        (infodict, out_table, out_table_others) = read_sites(
            infosites, blocks, Transversions=args.Transversions)
        for line in file:
            if not line.startswith('#'):
                line = line.rstrip().split('\t')
                position = int(line[1])
                chrom = line[0]
                id_position = line[0]+'_'+line[1]
                if id_position not in infodict.keys():
                    continue
                # Iterate over the blocks
                for block in blocks:
                    block_chrom, start_end = block.split(':')
                    start, end = map(int, start_end.split('-'))
                    # Check if the chromosome and position fall within the block
                    if chrom == block_chrom and start <= position <= end:
                        allele = gt2allele(
                            line, args.Diploid, args.Ambiguity, args.Indels)
                        infosite = infodict[id_position]
                        del infodict[id_position]
                        if (allele is not None):
                            if (args.Count):
                                Lineage = infosite['Lineage']
                                a_type = which_allele(
                                    allele, Ancestral=infosite['anc'], Derived=infosite['der'])
                                out_table[block][Lineage][a_type] += 1
                                other_samples = infosite['samples']
                                if len(other_samples) > 0:
                                    for i in range(0, len(other_samples)):
                                        a = other_samples[i]
                                        if a in ['A', 'C', 'G', 'T']:
                                            a_type = which_allele(
                                                a, infosite['anc'], infosite['der'])
                                            out_table_others[block][Lineage][a_type][i] += 1
                            elif args.Paste:
                                print(line[0], position-1, position,
                                      * list(infosite.values()), allele, sep="\t")
                            else:
                                myprint(line, allele, args.Bed, args.PrintRef)
                        break
                    elif chrom == block_chrom and position > end:
                        # Add the completed block to the set
                        completed_blocks.add(block)
                        blocks.remove(block)
                        break
        # Print the table with lineage assignment
        if args.Count:
            output = merge_table(out_table, out_table_others)
            for block in blocks_copy:
                o = output.get(block, 0)
                if o != 0:
                    print_table(o, label=block)
    else:
        for line in file:
            if not line.startswith('#'):
                line = line.rstrip().split('\t')
                allele = gt2allele(line, diploid=args.Diploid,
                                   Ambiguity=args.Ambiguity, Indels=args.Indels)
                if allele is not None:
                    myprint(line, allele=allele, bed=args.Bed,
                            printRef=args.PrintRef)


try:
    if args.filename:
        file_input = args.filename
        ftype = mimetypes.guess_type(file_input)
        if 'gzip' in ftype:
            with gzip.open(file_input, 'rt') as file:
                process_file(file)
        else:
            with open(file_input, 'rt') as file:
                process_file(file)

    elif not sys.stdin.isatty():
        process_file(sys.stdin)

    else:
        parser.print_help()

except BrokenPipeError:
    pass
except KeyboardInterrupt:
    pass

exit()
