#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import numpy as np
import sys


def localize(item, ls):
    count = 0
    for i in ls:
        if i == item:
            break
        else:
            count += 1
    return count


def firstNonspace(ls):
    for i in ls:
        if i != "":
            break
    return i


def gc(seq):
    gc = 0
    for bp in seq:
        if bp == "C" or bp == "G":
            gc += 1
    return gc / len(seq)


def Dictparser(Dictionary):
    lowest = float(1000)
    for i in Dictionary:
        if float(Dictionary[i]) < float(lowest):
            lowest = Dictionary[i]
            key = i
    return [i, lowest]


def reverseComplement(seq):
    out = []
    for i in range(len(seq) - 1, -1, -1):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def Complement(seq):
    out = []
    for i in range(0, len(seq)):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def ribosome(seq):
    NTs = ['T', 'C', 'A', 'G']
    stopCodons = ['TAA', 'TAG', 'TGA']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1

    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    return protein


def SeqCoord(seq, start, end):
    return seq[start:end]


def howMany(ls, exclude):
    counter = 0
    for i in ls:
        if i != exclude:
            counter += 1
    return counter


def stabilityCounter(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = (str(0) + str(int))
        return (string)
    if len(str(int)) > 4:
        string = str(int)
        return (string)


def sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


def ave(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count / len(ls)


def derep(ls):
    outLS = []
    for i in ls:
        if i not in outLS:
            outLS.append(i)
    return outLS


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    # data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def GCcalc(seq):
    count = 0
    for i in seq:
        if i == "G" or i == "C":
            count += 1
    return count / len(seq)


def reject_outliers(data):
    m = 2
    u = np.mean(data)
    s = np.std(data)
    filtered = [e for e in data if (u - 2 * s < e < u + 2 * s)]
    return filtered


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x


def RemoveDuplicates(ls):
    empLS = []
    counter = 0
    for i in ls:
        if i not in empLS:
            empLS.append(i)
        else:
            pass
    return empLS


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length - 1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x) - 1]


def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls) - 1]:
        x = i
    return x


def pull(item, one, two):
    ls = []
    counter = 0
    for i in item:
        if counter == 0:
            if i != one:
                pass
            else:
                counter += 1
                ls.append(i)
        else:
            if i != two:
                ls.append(i)
            else:
                ls.append(i)
                counter = 0
    outstr = "".join(ls)
    return outstr


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def removeLS(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    return emptyList


def fasta(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def fasta2(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # if re.findall(r'locus_tag', header):
                #     header = header.split("locus_tag=")[1]
                #     header = header.split(";")[0]
                # else:
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # if re.findall(r'locus_tag', header):
                #     header = header.split("locus_tag=")[1]
                #     header = header.split(";")[0]
                # else:
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def allButTheFirst(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(1, length):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)]


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def filterRe(list, regex):
    ls1 = []
    ls2 = []
    for i in list:
        if re.findall(regex, i):
            ls1.append(i)
        else:
            ls2.append(i)
    return ls1, ls2


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls = filter(ls, [""])
    return ls


parser = argparse.ArgumentParser(
    prog="PseudoHunter.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************
              @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
             @@@@@@@.........#@@@@@@@@@@@@@(.........@@@@@@@
            @@@@@@@@..,,.........@@@@@@@.........*,..@@@@@@@@
            @@@@@(...@@@@@@@@@@@.........@@@@@@@@@@@...,@@@@@
            @@@@@...@@@.......@@@&..%..@@@@@......@@@...@@@@@
            @@@@@@..@@@ PSEUDO @@.,&%%..@@@ HUNTER @@..@@@@@@
            @@@@@@@.@@@.......@@*.&&%%%.%@@@......@@@.@@@@@@@
            @@@@@@@.&@@@@@@@@@@..&&&%%%%..@@@@@@@@@@/.@@@@@@@
            @@@@@@@@..%@@@@@/..@&&&&%%%%%@..(@@@@@(.,@@@@@@@@
            @@@@@@@@@UGA@@@@@@&&&&&&%%%%%%%@@@@@@UAA@@@@@@@@@
            @@@@@@UAG@@@@@@@@@@...&&%%%...@@@@@@@@@@UGA@@@@@@
            @@@@UAA@@@@@@@@@.................@@@@@@@@UAG@@@@@
            @@@@@..@@@@@@@.....................@@@@@@%..@@@@@
            @@@@@..@@@@@............@............@@@@@..@@@@@
            @@@@@.................@@@@@.................@@@@@
             @@@@@@@..........&@@@@@@@@@@@#..........@@@@@@@
              @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    Developed by Arkadiy Garber and John McCutcheon
    University of Montana, Biological Sciences
    Please send comments and inquiries to arkadiy.garber@mso.umt.edu
    ASCII art: https://manytools.org/hacker-tools/convert-images-to-ascii-art/ 

    You can provide this program with contigs via the \'-q and -r\'
    arguments, OR you can provide this program with pre-existing ORF-calls
    from Prodigal via the \'-n, -a, -rn, and -ra\' arguments. The ORF-calls
    must be from Prodigal, or formatted in such a way where the header name
    ends in an underscore and number (e.g. CP003881_1, CP003881_2, CP003881_3,...)

    The codeml control file is provided within the codeml program: paml4.8/codeml.ctl
    ************************************************************************
    '''))

# parser.add_argument('-blast', type=str, help="clu.tsv output from mmseqs", default="NA")

parser.add_argument('-q', type=str, help="input query genome (contigs in FASTA format). This is requried if you would "
                                         "like PseudoHunter to look for pseudogenes in intergenic regions.", default="NA")
parser.add_argument('-r', type=str,
                    help="input reference dataset (could be contigs from one genome or a set of closely-related genomes)",
                    default="NA")
parser.add_argument('-ref', type=str,
                    help="is the reference dataset one genome or multiple? (one/multiple) Default=one", default="one")

parser.add_argument('-n', type=str, help="ORFs in nucleic acid format", default="NA")
parser.add_argument('-a', type=str, help="ORFs in amino acid format", default="NA")
parser.add_argument('-rn', type=str, help="Reference ORFs in nucleic acid format", default="NA")
parser.add_argument('-ra', type=str, help="Reference ORFs in amino acid format", default="NA")

parser.add_argument('-gff', type=str,
                    help="gff file for the query genome. Provide this file if you would like PseudoHunter to "
                         "use your predicted ORFs, rather than starting from raw contigs.", default="NA")

parser.add_argument('-ctl', type=str, help="template control file for codeml", default="NA")

parser.add_argument('-out', type=str, help="name output directory", default="PseudoHunter_output")

parser.add_argument('-l', type=float,
                    help="minimum proportion of of target gene length that must be covered by alignment with query gene "
                         "for the query gene to be classified as \'intact\' (default = 0.75)",
                    default=0.75)

parser.add_argument('-L', type=float,
                    help="maximum run-on length, relative to reference, for a gene to be classified as \'intact\' (default = 1.25)",
                    default=1.25)

parser.add_argument('-d', type=float, help="maximum dN/dS value for gene too be considered \'intact\' (default = 0.3)",
                    default=0.3)
parser.add_argument('-M', type=float, help="maximum dS value for dN/dS calculation (default = 3)", default=3)
parser.add_argument('-m', type=float, help="minimum dS value for dN/dS calculation (default = 0.001)", default=0.001)
parser.add_argument('-t', type=int, help="number of threads to use for BLAST", default=1)
parser.add_argument('-e', type=str, help="e-value for BLAST search. Default = 1E-6", default=float(1E-6))
parser.add_argument('-s', type=str, help="search engine to use (blast/diamond). Default = blast", default="blast")
parser.add_argument('-delim', type=str, help="if you are poviding files with gene sequences, please provide a character "
                                             "that separates the gene number from the rest of the header name "
                                             "(ex: contig_1, that character is \'_\'). Default = \'_\'", default="_")
parser.add_argument('--skip', type=str,
                    help="By choosing this option, and providing pseudoHunter with the previously-created output directory, "
                         "you are choosing to skip time-consuming steps of this pipeline "
                         "(e.g. BLAST, Muscle, codeml), and would like to use the output files "
                         "created from a previous run to re-do the anlaysis, "
                         "perhaps with different parameters. All other arguments "
                         "(e.g. \'-a\', \'-n\', \'-q\', \'-r\', and particulary \'-out\') still need to be provided as before",
                    const=True, nargs="?")

args = parser.parse_args()

cwd = os.getcwd()

os.system("echo ${ctl} > ctl.txt")
file = open("ctl.txt")
for i in file:
    ctl = (i.rstrip())
os.system("rm ctl.txt")

if not re.findall(r'codeml', ctl):
    if args.ctl == "NA":
        print("\nUh-oh...looks like PseudoHunter cannot locate the necessary \ncodeml control file. Please double check that "
              "the installation\n with the \'setup.sh\' script worked without any errors. If issues \npersist, please do not hesistate to "
              "report an issue on GitHub: \nhttps://github.com/Arkadiy-Garber/PseudoHunter/issues\n")
        raise SystemExit
    else:
        ctl = args.ctl


mode = 0
if args.gff != "NA":
    if args.a == "NA" or args.n == "NA":
        print("\nLooks like you have provided a GFF file containing annotation data\n, but PseudoHunter did not find one "
              "or more files with query ORFs\n. These are necessary if you would like to incorporate your annotations. "
              "Please double-check your command, and make sure that \'-a\' and \'-n\' arguments are present")
        raise SystemExit
    else:
        faa = args.a
        fna = args.n
else:
    if args.a != "NA" or args.n != "NA":
        if args.q != "NA":
            print("\nLooks like you did not provide a GFF file, but have provided some ORF predictions.\n"
                  "For PseudoHunter to use ORF prediction, it needs an associated GFF file.\nSince that was not provided, "
                  "PseudoHunter can proceed with the provided contigs: " + args.q + "\n")
            choice = input("Would you like PseudoHunter to proceed with the provided contigs (y/n)? ")
            if choice == "n":
                print("Exiting...")
                raise SystemExit
            else:
                print("Alright! Moving on with " + args.q)
                mode = 1
                os.system("prodigal -i %s -a %s-proteins.faa -d %s-proteins.fna" % (args.q, args.q, args.q))
                faa = args.q + "-proteins.faa"
                fna = args.q + "-proteins.fna"

        else:
            print("\nLooks like you did not provide a GFF file, but have provided some ORF predictions.\n"
                  "For PseudoHunter to use ORF prediction, it needs an associated GFF file. \nPlease either "
                  "provide a GFF file via the \'-gff\' argument.\nOr provide contigs via the \'-q\' and \'-r\' arguments.\n")
            raise SystemExit
    else:
        if args.q != "NA":
            os.system("prodigal -i %s -a %s-proteins.faa -d %s-proteins.fna" % (args.q, args.q, args.q))
            faa = args.q + "-proteins.faa"
            fna = args.q + "-proteins.fna"
        else:
            print("PseudoHunter did not find your input files. "
                  "Please provide those via the \'-q\', or \'-a\' and \'-n\' arguments.")

counter = 0
if args.rn != "NA":
    refFna = args.rn
    counter += 1

if args.ra != "NA":
    refFaa = args.ra
    counter += 1

if counter < 2:
    # USER DID NOT PROVIDE ONE OR BOTH OF THE ORF FILES. SO WILL USE THE CONTIGS PROVIDED BY THE ARGS.R ARGUMENT

    if args.ref == "one":
        os.system("prodigal -i %s -a %s-proteins.faa -d %s-proteins.fna" % (args.r, args.r, args.r))
    else:
        os.system("prodigal -i %s -a %s-proteins.faa -d %s-proteins.fna -p meta" % (args.r, args.r, args.r))

    refFna = args.r + "-proteins.fna"
    refFaa = args.r + "-proteins.faa"


if not args.skip:
    print("Starting pipeline...")

    os.system("mkdir " + args.out)
    os.system('mkdir ' + args.out + "/dnds-analysis")

    if args.s == "blast":
        print("Running BLAST")
        os.system("makeblastdb -dbtype prot -in %s -out %s" % (refFaa, refFaa))
        # os.system("rm makeblastdb.out")
        os.system("blastp -query %s -db %s "
                  "-outfmt 6 -out %s/pseudogene.blast -evalue 1E-6 -num_threads %s -max_target_seqs 1" % (
                      faa, refFaa, args.out, args.t))

        os.system("rm %s.psq" % refFaa)
        os.system("rm %s.phr" % refFaa)
        os.system("rm %s.pin" % refFaa)

    elif args.s == "diamond":
        print("Running DIAMOND")
        os.system(
            "diamond makedb --in %s -d %s" % (args.ra, args.ra))
        # os.system("rm makedb.out")
        os.system("diamond blastp --db %s.dmnd --query %s --outfmt 6 --out %s/pseudogene.blast "
                  "--max-target-seqs 1 --evalue 1E-6 --threads %d" % (args.ra, args.a, args.out, args.t))

        # os.system("rm %s.dmnd" % args.ra)

    ####################################################################################################################
    faaDict = open(faa)
    faaDict = fasta2(faaDict)

    fnaDict = open(fna)
    fnaDict = fasta2(fnaDict)

    refFaaDict = open(refFaa)
    refFaaDict = fasta2(refFaaDict)

    refFnaDict = open(refFna)
    refFnaDict = fasta2(refFnaDict)

    prescreened = []
    alnLengthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    blast = open("%s/pseudogene.blast" % args.out)
    for i in blast:
        ls = i.rstrip().split("\t")
        if ls[0] not in prescreened:
            print(ls[0])
            print(ls[1])
            print("")
            alnLengthDict[ls[0]] = ls[3]

            outNUC = open(args.out + "/dnds-analysis/%s.faa.fna" % ls[0], "w")
            outNUC.write(">" + ls[1] + "\n")
            outNUC.write(refFnaDict[ls[1]] + "\n")
            outNUC.write(">" + ls[0] + "\n")
            outNUC.write(fnaDict[ls[0]] + "\n")
            outNUC.close()

            outAA = open(args.out + "/dnds-analysis/%s.faa" % ls[0], "w")
            outAA.write(">" + ls[1] + "\n")
            outAA.write(refFaaDict[ls[1]] + "\n")
            outAA.write(">" + ls[0] + "\n")
            outAA.write(faaDict[ls[0]] + "\n")
            outAA.close()

            if float(ls[3]) / len(refFaaDict) > args.l:
                prescreened.append(ls[0])

    # ALIGNING PROTEIN SEQUENCES AND CREATING A CODON ALIGNMENT
    print("aligning files...")
    DIR = args.out + "/dnds-analysis"
    os.system("for i in %s/*faa; do"
              " muscle -in $i -out $i.aligned.fa;"
              # " rm muscle.out;"
              " pal2nal.pl $i.aligned.fa $i.fna -output fasta > $i.codonalign.fa;"
              " done" % DIR)

    # BUILDING CONTROL FILES
    print("preparing for codeml analysis")
    DIR = args.out + "/dnds-analysis"
    codealign = os.listdir(DIR)
    count = 0
    for file in codealign:
        if re.findall(r'codonalign', file):
            count += 1
            clu = file.split(".faa")[0]
            setup = open(ctl)
            out = open("%s/%s.ctl" % (DIR, str(clu)), "w")

            for i in setup:
                if re.findall('seqfile', i):
                    out.write(
                        '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'seqfile' + ' ' + '= ' + args.out + '/dnds-analysis/' + file + ' ' + '*' + ' ' + 'sequence' + ' ' + 'data' + ' ' + 'filename\n')

                elif re.findall(r'outfile', i):
                    out.write(
                        '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'outfile' + ' ' + '=' + ' '
                        + args.out + '/dnds-analysis/mlcTree_' + str(
                            clu) + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' '
                        + '' + ' ' + '' + ' ' + '' + ' ' + '*' + ' ' + 'main' + ' ' + 'result' + ' ' + 'file' + ' ' + 'name\n')

                else:
                    out.write(i)
            out.close()
    print("")

    # RUNNING CODEML FOR DN/DS CALCULATION
    total = 0
    for i in codealign:
        if re.findall(r'codonalign', i):
            total += 1

    count = 0
    codealign = os.listdir(DIR)
    for file in codealign:
        if lastItem(file.split(".")) == "ctl":
            count += 1
            perc = (count / total) * 100
            sys.stdout.write("running codeml: %d%%   \r" % (perc))
            sys.stdout.flush()
            os.system("codeml %s/dnds-analysis/%s" % (args.out, file))
            # os.system("rm codeml.out")
    print("")

    # PARSING CODEML OUTPUT

    cwd = os.getcwd()
    DIR = args.out + "/dnds-analysis"

    if args.gff != "NA":
        gff = open(args.gff)
        gffDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'NA')))
        for i in gff:
            ls = i.rstrip().split("\t")
            if re.match(r'##FASTA', i):
                break
            else:
                if not re.match(r'#', i):
                    contig = ls[0]
                    orf = ls[8]
                    orf = orf.split(";")[0]
                    orf = orf.split("=")[1]

                    product = lastItem(ls[8].split(";")).split("=")[1]
                    product = replace(product, [","], ";")

                    gffDict[orf]["product"] = product
                    gffDict[orf]["contig"] = contig
                    gffDict[orf]["start"] = ls[3]
                    gffDict[orf]["end"] = ls[4]
                    gffDict[orf]["strand"] = ls[6]

    else:
        gffDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'NA')))
        faaDict1 = open(faa)
        faaDict1 = fasta(faaDict1)

        for i in faaDict1.keys():
            ls = i.split(" # ")
            start = ls[1]
            end = ls[2]
            contig = allButTheLast(ls[0], "_")
            gffDict[ls[0]]["contig"] = contig
            gffDict[ls[0]]["start"] = start
            gffDict[ls[0]]["end"] = end
            gffDict[ls[0]]["product"] = "NA"
            if ls[3] == "-1":
                gffDict[ls[0]]["strand"] = "-"
            else:
                gffDict[ls[0]]["strand"] = "+"

    if args.ra != "NA":
        faaRef = open(args.ra)
        faaRef = fasta2(faaRef)
    else:
        faaRef = open(args.r + "-proteins.faa")
        faaRef = fasta2(faaRef)

    if args.rn != "NA":
        fnaRef = open(args.rn)
        fnaRef = fasta2(fnaRef)
    else:
        fnaRef = open(args.r + "-proteins.fna")
        fnaRef = fasta2(fnaRef)

    if args.a != "NA" and mode == 0:
        faa = open(args.a)
        faa = fasta2(faa)

    else:
        faa = open(args.q + "-proteins.faa")
        faa = fasta2(faa)


    if args.n != "NA" and mode == 0:
        fna = open(args.n)
        fna = fasta2(fna)
    else:
        fna = open(args.q + "-proteins.fna")
        fna = fasta2(fna)

    alnLengthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    alnIdDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    blast = open(args.out + "/pseudogene.blast")
    for i in blast:
        ls = i.rstrip().split("\t")
        alnLengthDict[ls[0]][ls[1]] = ls[3]
        alnIdDict[ls[0]][ls[1]] = ls[2]

    print("summarizing codeml output")
    codealign = os.listdir(DIR)
    dndsDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in codealign:
        if re.findall(r'mlc', i):
            file = open(DIR + "/%s" % i, "r")
            for j in file:
                if re.search(r'#1', j):
                    ls = (j.rstrip().split(" "))
                    orf = ls[1]
                if re.search(r'#2', j):
                    ls = (j.rstrip().split(" "))
                    NODE = ls[1]
                line = j.rstrip()
            ls = line.split("  ")
            dS = remove(lastItem(ls), [" ", "=", "d", "S"])
            dN = remove(lastItem(ls[0:len(ls) - 1]), [" ", "=", "d", "N"])
            dndsDict[NODE]["orf"] = orf
            dndsDict[NODE]["dn"] = dN
            dndsDict[NODE]["ds"] = dS

    count = 0
    dndsList = []
    dndsDict2 = defaultdict(list)
    for i in sorted(dndsDict.keys()):
        count += 1
        if float(dndsDict[i]["dn"]) <= args.M and float(dndsDict[i]["ds"]) <= args.M and float(
                dndsDict[i]["ds"]) >= args.m and float(dndsDict[i]["dn"]) >= args.m:
            orf = dndsDict[i]["orf"]
            dn = dndsDict[i]["dn"]
            ds = dndsDict[i]["ds"]
            try:
                dnds = float(dndsDict[i]["dn"]) / float(dndsDict[i]["ds"])
            except ZeroDivisionError:
                dnds = "NA"
            dndsDict2[orf].append(i)
            if dnds != "NA":
                dndsList.append(dnds)

    print("preparing final output file: summary.csv")
    dsList = []
    dnList = []
    dndsList2 = []
    out = open(args.out + "/summary.csv", "w")
    out.write(
        "ORF_calls" + "," + "Ortholog" + "," + "Pseudogene" + "," + "contig" + "," + "start" + "," + "end" + "," + "strand" + "," + "geneLength" + "," + "AlignmentLength" + "," +
        "OrthologLength" + "," + "Identity" + "," + "Annotation" + "," + "Pseudogene_confidence" + "," +
        "NumberOfGeneFrags" + "," + "AlignmentLength/OrthologLength" + "," + "geneLength/OrthologLength" + "," + "dN" + "," + "dS" + "," + "dN/dS" + "," +
        "Translation" + "," + "Sequences" + "\n")

    for i in dndsDict2.keys():
        if len(dndsDict2[i]) > 1:
            dndsDict3 = defaultdict(list)
            for k in dndsDict2[i]:
                dndsDict3[allButTheLast(k, args.delim)].append(int(lastItem(k.split(args.delim))))
            for l in dndsDict3.keys():
                listOfLists = (cluster(dndsDict3[l], 2))
                for m in listOfLists:
                    try:
                        if len(m) > 1:
                            seq = ''
                            seq2 = ''
                            dnLS = []
                            dsLS = []
                            annotations = ''
                            ORFs = ''
                            TotalAlnLength = 0

                            idLS = []

                            for n in m:
                                originalN = stabilityCounter(n)
                                # originalN = n
                                ORF = (l + args.delim + str(originalN))

                                ORFs += ORF + "|"

                                annotation = gffDict[ORF]["product"]
                                annotations += annotation + "|"
                                seq += faa[ORF]
                                seq2 += fna[ORF]
                                dnLS.append(float(dndsDict[ORF]["dn"]))
                                dsLS.append(float(dndsDict[ORF]["ds"]))

                                alnLength = int(alnLengthDict[ORF][i])
                                TotalAlnLength += alnLength

                                identity = float(alnIdDict[ORF][i])
                                idLS.append(identity)

                            identity = ave(idLS)

                            dn = ave(dnLS)
                            ds = ave(dsLS)

                            try:
                                dnds = dn / ds
                            except ZeroDivisionError:
                                dnds = "NA"

                            fragments = len(m)
                            ORF = ORFs[0:len(ORFs) - 1]
                            annotation = annotations

                            strand = gffDict[ORF.split("|")[0]]["strand"]
                            contig = gffDict[ORF.split("|")[0]]["contig"]
                            start = gffDict[ORF.split("|")[0]]["start"]
                            end = gffDict[lastItem(ORF.split("|"))]["end"]

                        else:

                            fragments = 1
                            originalM0 = stabilityCounter(m[0])
                            # originalM0 = m[0]
                            ORF = (l + args.delim + str(originalM0))

                            contig = gffDict[ORF]["contig"]
                            start = gffDict[ORF]["start"]
                            end = gffDict[ORF]["end"]
                            strand = gffDict[ORF]["strand"]

                            annotation = gffDict[ORF]["product"]
                            seq = faa[ORF]
                            seq2 = fna[ORF]

                            TotalAlnLength = int(alnLengthDict[ORF][i])

                            identity = float(alnIdDict[ORF][i])

                            dn = dndsDict[ORF]["dn"]
                            ds = dndsDict[ORF]["ds"]

                            try:
                                dnds = float(dn) / float(ds)
                            except ZeroDivisionError:
                                dnds = "NA"

                    except (TypeError, ValueError):
                        if len(m) > 1:
                            seq = ''
                            seq2 = ''
                            dnLS = []
                            dsLS = []
                            ORFs = ''
                            TotalAlnLength = 0

                            idLS = []

                            for n in m:
                                # originalN = stabilityCounter(n)
                                originalN = n
                                ORF = (l + args.delim + str(originalN))

                                ORFs += ORF + "|"

                                annotation = gffDict[ORF]["product"]
                                annotations += annotation + "|"
                                seq += faa[ORF]
                                seq2 += fna[ORF]
                                dnLS.append(float(dndsDict[ORF]["dn"]))
                                dsLS.append(float(dndsDict[ORF]["ds"]))

                                alnLength = int(alnLengthDict[ORF][i])
                                TotalAlnLength += alnLength

                                identity = float(alnIdDict[ORF][i])
                                idLS.append(identity)

                            identity = ave(idLS)

                            dn = ave(dnLS)
                            ds = ave(dsLS)

                            try:
                                dnds = dn / ds
                            except ZeroDivisionError:
                                dnds = "NA"

                            fragments = len(m)
                            ORF = ORFs[0:len(ORFs) - 1]
                            annotation = annotations

                            contig = gffDict[ORF.split("|")[0]]["contig"]
                            strand = gffDict[ORF.split("|")[0]]["strand"]
                            start = gffDict[ORF.split("|")[0]]["start"]
                            end = gffDict[lastItem(ORF.split("|"))]["end"]

                        else:

                            fragments = 1
                            # originalM0 = stabilityCounter(m[0])
                            originalM0 = m[0]
                            ORF = (l + args.delim + str(originalM0))

                            strand = gffDict[ORF]["strand"]
                            contig = gffDict[ORF]["contig"]
                            start = gffDict[ORF]["start"]
                            end = gffDict[ORF]["end"]

                            annotation = gffDict[ORF]["product"]
                            seq = faa[ORF]
                            seq2 = fna[ORF]

                            TotalAlnLength = int(alnLengthDict[ORF][i])

                            identity = float(alnIdDict[ORF][i])

                            dn = dndsDict[ORF]["dn"]
                            ds = dndsDict[ORF]["ds"]
                            try:
                                dnds = float(dn) / float(ds)
                            except ZeroDivisionError:
                                dnds = "NA"

        else:
            fragments = 1
            ORF = dndsDict2[i][0]
            contig = gffDict[ORF]["contig"]
            strand = gffDict[ORF]["strand"]

            start = gffDict[ORF]["start"]
            end = gffDict[ORF]["end"]

            annotation = gffDict[ORF]["product"]
            seq = faa[ORF]
            seq2 = fna[ORF]
            dn = dndsDict[ORF]["dn"]
            ds = dndsDict[ORF]["ds"]

            try:
                dnds = float(dn) / float(ds)
            except ZeroDivisionError:
                dnds = "NA"

            TotalAlnLength = int(alnLengthDict[ORF][i])

            identity = float(alnIdDict[ORF][i])

        ratio = TotalAlnLength / len(faaRef[i])

        out.write(ORF + "," + i + ",")

        if dnds != "NA":
            if dnds > args.d or ratio < args.l or fragments > 1 or (len(seq) / len(faaRef[i])) > args.L:
                out.write("Y" + ",")
                prob = 1

                # GETTING MY STUPID PSEUDOGENE SCORE CALCULATED
                if dnds != "NA":
                    inflation = float(dnds) / ave(dndsList)
                else:
                    inflation = ave(dndsList)

                prob = prob * inflation
                prob = prob * float(fragments)
                if ratio < 1:
                    prob = prob / ratio
                else:
                    prob = prob * ratio

            else:
                out.write("N" + ",")

                # GETTING MY STUPID PSEUDOGENE SCORE CALCULATED
                prob = 1

                if dnds != "NA":
                    inflation = float(dnds) / ave(dndsList)
                else:
                    inflation = ave(dndsList)

                prob = prob * inflation
                prob = prob * float(fragments)
                if ratio < 1:
                    prob = prob / ratio
                else:
                    prob = prob * ratio
        else:
            if ratio < args.l or fragments > 1 or (len(seq) / len(faaRef[i])) > args.L:
                out.write("Y" + ",")
                prob = 1

                # GETTING MY STUPID PSEUDOGENE SCORE CALCULATED
                if dnds != "NA":
                    inflation = float(dnds) / ave(dndsList)
                else:
                    inflation = ave(dndsList)

                prob = prob * inflation
                prob = prob * float(fragments)
                if ratio < 1:
                    prob = prob / ratio
                else:
                    prob = prob * ratio

            else:
                out.write("N" + ",")

                # GETTING MY STUPID PSEUDOGENE SCORE CALCULATED
                prob = 1

                if dnds != "NA":
                    inflation = float(dnds) / ave(dndsList)
                else:
                    inflation = ave(dndsList)

                prob = prob * inflation
                prob = prob * float(fragments)
                if ratio < 1:
                    prob = prob / ratio
                else:
                    prob = prob * ratio

        # WRITING TO FILE

        out.write(
            str(contig) + "," + str(start) + "," + str(end) + "," + strand + "," + str(int(end) - int(start)) + "," + str(
                TotalAlnLength * 3) + "," +
            str(len(faaRef[i]) * 3) + "," + str(identity) + "," + str(annotation) + "," + str(prob) + "," +
            str(fragments) + "," + str(ratio) + "," + str(((int(end) - int(start))/3) / len(faaRef[i])) + "," + str(dn) + "," +
            str(ds) + "," + str(dnds) + "," + seq + "," + seq2 + "\n")

        dnList.append(dn)
        dsList.append(ds)

        if dnds != "NA":
            dndsList2.append(dnds)
    out.close()

    # INTERGENIC REGION ANALYSIS
    if args.q != "NA":
        print("Analyzing intergenic regions")
        contigs = open(args.q)
        contigs = fasta2(contigs)

        firstRow = ''
        count = 0
        summaryDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        summaryDict2 = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        summaryDict3 = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        summary = open(args.out + "/summary.csv")
        for i in summary:
            ls = i.rstrip().split(",")
            if ls[0] != "ORF_calls":
                count += 1
                contig = ls[3]
                start = int(ls[4])
                summaryDict[contig][start] = ls
                summaryDict2[ls[1]] = ls[0]

                summaryDict3[ls[0]]["ORF_calls"] = ls[0]
                summaryDict3[ls[0]]["Ortholog"] = ls[1]
                summaryDict3[ls[0]]["Pseudogene"] = ls[2]
                summaryDict3[ls[0]]["contig"] = ls[3]
                summaryDict3[ls[0]]["start"] = ls[4]
                summaryDict3[ls[0]]["end"] = ls[5]
                summaryDict3[ls[0]]["strand"] = ls[6]
                summaryDict3[ls[0]]["geneLength"] = ls[7]
                summaryDict3[ls[0]]["AlignmentLength"] = ls[8]
                summaryDict3[ls[0]]["OrthologLength"] = ls[9]
                summaryDict3[ls[0]]["Identity"] = ls[10]
                summaryDict3[ls[0]]["Annotation"] = ls[11]
                summaryDict3[ls[0]]["Pseudogene_confidence"] = ls[12]
                summaryDict3[ls[0]]["NumberOfGeneFrags"] = ls[13]
                summaryDict3[ls[0]]["AlignmentLength/OrthologLength"] = ls[14]
                summaryDict3[ls[0]]["geneLength/OrthologLength"] = ls[15]
                summaryDict3[ls[0]]["dN"] = ls[16]
                summaryDict3[ls[0]]["dS"] = ls[17]
                summaryDict3[ls[0]]["dN/dS"] = ls[18]
                summaryDict3[ls[0]]["Translation"] = ls[19]
                summaryDict3[ls[0]]["Sequences"] = ls[20]
            else:
                firstRow = i.rstrip()

        igDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        count = 0
        for i in summaryDict.keys():
            end = 0
            endORF = ""
            for j in (sorted(summaryDict[i])):
                ls = (summaryDict[i][j])
                contig = ls[3]
                start = int(ls[4])
                startOrf = ls[0]
                if start > end:
                    contigSeq = contigs[contig]
                    IG = contigSeq[end:start]
                    if len(IG) > 100:
                        count += 1
                        igDict[endORF + "!" + str(contig) + "-" + str(end) + "-" + str(start) + "!" + startOrf] = IG
                end = int(ls[5])
                endORF = ls[0]

            if end < len(contigSeq):
                IG = contigSeq[end:len(contigSeq)]
                if len(IG) > 100:
                    count += 1
                    igDict[endORF + "!" + str(contig) + "-" + str(end) + "-" + str(len(contigSeq)) + "!" + ""] = IG

        out = open(args.out + "/ig.fasta", "w")
        for i in igDict.keys():
            out.write(">" + i + "\n")
            out.write(igDict[i] + "\n")

        if args.rn != "NA":
            REF = args.rn
        else:
            REF = args.r + "-proteins.fna"

        os.system("makeblastdb -dbtype nucl -in %s -out %s" % (REF, REF))
        os.system(
            "blastn -query %s/ig.fasta -db %s -outfmt 6 -out %s/ig.blast -evalue 1E-6 -num_threads 4 -max_target_seqs 20" % (
            args.out, REF, args.out))
        os.system("rm %s.nhr" % REF)
        os.system("rm %s.nin" % REF)
        os.system("rm %s.nsq" % REF)

        igFasta = open("%s/ig.fasta" % args.out)
        igFasta = fasta2(igFasta)
        igFastaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in igFasta.keys():
            header = (i.split("!")[1])
            igFastaDict[header] = igFasta[i]

        refFna = open(REF)
        refFna = fasta2(refFna)

        out = open("%s/ig.csv" % args.out, "w")
        out.write(firstRow + "\n")
        blast = open("%s/ig.blast" % args.out)
        for i in blast:
            ls = i.rstrip().split("\t")

            if int(ls[8]) < int(ls[9]):
                strand = "+"
            else:
                strand = "-"

            start = int(ls[6])
            end = int(ls[7])
            if end - start > 100:

                contig = ls[0].split("!")[1].split("-")[0]

                ORF = "IG" + str(contig) + "_" + str(start + int(ls[0].split("!")[1].split("-")[1])) + "-" + str(
                    int(ls[0].split("!")[1].split("-")[1]) + end)

                percIdent = float(ls[2])
                alnLength = int(ls[3])

                FreeLivingOrtholog = (ls[1])

                ig = ls[0].split("!")[1]
                seq = igFastaDict[ig]
                HomologousRegion = seq[start - 1:end]

                lengthIntergenic = (int(ls[0].split("!")[1].split("-")[2]) - int(ls[0].split("!")[1].split("-")[1]))

                if len(summaryDict2[ls[1]]) > 0:
                    # THE ORTHOLOG IN FREE-LIVING RELATIVE HAD A PRIOR HIT TO A PREVIOSLY-PREDICTED ENDOSYMBIONT ORFS
                    PreExistingOrfHit = (summaryDict2[ls[1]])
                    hitLocation = (localize(PreExistingOrfHit, ls[0].split("!")))
                    if hitLocation > 2:
                        # THE PREVIOUSLY-PREDICTED ENDOSYMBIONT ORF THAT IS A
                        # BLAST-MATCH TO THE IG-HOMOLOGOUS ORTHOLOG IS NOT ADJACENT TO THE IG

                        out.write(
                            ORF + "," + FreeLivingOrtholog + "," + "Y" + "," + str(contig) + "," + str(start) + "," +
                            str(end) + "," + strand + "," + "NA" + "," + str(len(HomologousRegion)) + "," +
                            str(len(refFna[FreeLivingOrtholog])) + "," + str(
                                percIdent) + "," + "NA" + "," + "NA" + "," + "NA" +
                            "," + str(
                                alnLength / len(refFna[FreeLivingOrtholog])) + "," + "NA" + "," + "NA" + "," + "NA" +
                            "," + "NA" + "," + "NA" + "," + HomologousRegion + "\n")

                    else:
                        # THE PREVIOUSLY-PREDICTED ENDOSYMBIONT ORF THAT IS A
                        # BLAST-MATCH TO THE IG-HOMOLOGOUS ORTHOLOG IS ADJACENT TO THE IG, AND SHOULD BE COMBINED WITH IT.

                        NewGeneLength = int(summaryDict3[PreExistingOrfHit]["geneLength"]) + len(HomologousRegion)
                        NewAlignLength = int(summaryDict3[PreExistingOrfHit]["AlignmentLength"]) + len(HomologousRegion)
                        NewIdent = (float(summaryDict3[PreExistingOrfHit]["Identity"]) + float(percIdent)) / 2
                        Frags = int(summaryDict3[PreExistingOrfHit]["NumberOfGeneFrags"]) + 1
                        newAlnRatio = NewAlignLength / len(refFna[FreeLivingOrtholog])
                        newGeneRatio = NewGeneLength / len(refFna[FreeLivingOrtholog])

                        if hitLocation == 0:
                            contigStartPosition = summaryDict3[PreExistingOrfHit]["start"]
                            contigEndPosition = int(summaryDict3[PreExistingOrfHit]["end"]) + len(HomologousRegion)
                            newORF = PreExistingOrfHit + "|" + ORF
                            newNucSeq = summaryDict3[PreExistingOrfHit]["Sequences"] + HomologousRegion

                            out.write(
                                newORF + "," + FreeLivingOrtholog + "," + "Y" + "," + summaryDict3[PreExistingOrfHit][
                                    "contig"] + "," +
                                str(contigStartPosition) + "," + str(contigEndPosition) + "," + strand + "," + str(
                                    NewGeneLength) + "," +
                                str(NewAlignLength) + "," + str(len(refFna[FreeLivingOrtholog])) + "," + str(
                                    NewIdent) + "," +
                                summaryDict3[PreExistingOrfHit]["Annotation"] + "," + summaryDict3[PreExistingOrfHit][
                                    "Pseudogene_confidence"] +
                                "," + str(Frags) + "," + str(newAlnRatio) + "," + str(newGeneRatio) + "," +
                                summaryDict3[PreExistingOrfHit]["dN"] + "," +
                                summaryDict3[PreExistingOrfHit]["dS"] + "," + summaryDict3[PreExistingOrfHit][
                                    "dN/dS"] + "," + summaryDict3[PreExistingOrfHit]["Translation"] + "," +
                                newNucSeq + "\n")

                        else:
                            contigStartPosition = int(summaryDict3[PreExistingOrfHit]["start"]) - len(HomologousRegion)
                            contigEndPosition = int(summaryDict3[PreExistingOrfHit]["end"])
                            newORF = ORF + "|" + PreExistingOrfHit
                            newNucSeq = HomologousRegion + summaryDict3[PreExistingOrfHit]["Sequences"]

                            out.write(
                                newORF + "," + FreeLivingOrtholog + "," + "Y" + "," + summaryDict3[PreExistingOrfHit][
                                    "contig"] + "," +
                                str(contigStartPosition) + "," + str(contigEndPosition) + "," + strand + "," + str(
                                    NewGeneLength) + "," +
                                str(NewAlignLength) + "," + str(len(refFna[FreeLivingOrtholog])) + "," + str(
                                    NewIdent) + "," +
                                summaryDict3[PreExistingOrfHit]["Annotation"] + "," + summaryDict3[PreExistingOrfHit][
                                    "Pseudogene_confidence"] +
                                "," + str(Frags) + "," + str(newAlnRatio) + "," + str(newGeneRatio) + "," +
                                summaryDict3[PreExistingOrfHit]["dN"] + "," +
                                summaryDict3[PreExistingOrfHit]["dS"] + "," + summaryDict3[PreExistingOrfHit][
                                    "dN/dS"] + "," + summaryDict3[PreExistingOrfHit]["Translation"] + "," +
                                newNucSeq + '\n')

                else:
                    # THE ORTHOLOG IN FREE-LIVING RELATIVE HAD NO PRIOR HITS TO ANY OF THE PREDICTED ENDOSYMBIONT ORFS
                    out.write(
                        ORF + "," + FreeLivingOrtholog + "," + "Y" + "," + str(contig) + "," + str(start) + "," + str(
                            end) + "," + strand +
                        "," + "NA" + "," + str(len(HomologousRegion)) + "," + str(
                            len(refFna[FreeLivingOrtholog])) + "," +
                        str(percIdent) + "," + "NA" + "," + "NA" + "," + "NA" + "," +
                        str(alnLength / len(refFna[FreeLivingOrtholog])) + "," + "NA" + "," + "NA" + "," +
                        "NA" + "," + "NA" + "," + "NA" + "," + HomologousRegion + '\n')

        out.close()

        igSummaryDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        igSummary = open("%s/ig.csv" % args.out)
        for i in igSummary:
            ls = i.rstrip().split(",")
            if ls[0] != "ORF_calls":
                igSummaryDict[ls[0]] = i.rstrip()

        out = open("%s/summary-final.csv" % args.out, "w")
        originalSummary = open("%s/summary.csv" % args.out)
        for i in originalSummary:
            counter = 0
            ls = i.rstrip().split(",")
            if ls[0] != "ORF_calls":
                for j in igSummaryDict.keys():
                    if ls[0].split("|")[0] in j.split("|"):
                        counter += 1
                        out.write(igSummaryDict[j] + "\n")
                if counter == 0:
                    out.write(i.rstrip() + "\n")
                else:
                    counter = 0
            else:
                out.write(i.rstrip() + "\n")

        for i in igSummaryDict.keys():
            if not re.findall(r'\|', i):
                out.write(igSummaryDict[i] + "\n")
        out.close()

        os.system("rm %s/summary.csv" % args.out)

        count = 0
        end = 0
        out = open("%s/summary.csv" % args.out, "w")
        summaryExtend = open("%s/summary-final.csv" % args.out)
        for i in summaryExtend:
            ls = i.rstrip().split(",")
            if (re.findall(r'IG', ls[0]) and not re.findall(r'\|', ls[0])):
                start = lastItem(ls[0].split("_")).split("-")[0]
                if int(start) > int(end) - 10:
                    count += 1
                    out.write(i.rstrip() + "\n")
                end = lastItem(ls[0].split("_")).split("-")[1]
            else:
                count += 1
                out.write(i.rstrip() + "\n")
        out.close()

        os.system("rm %s/ig.blast %s/ig.csv %s/ig.fasta" % (args.out, args.out, args.out))
        os.system("rm %s/summary-final.csv" % args.out)

    ###################################################################################################################
    if len(dsList) > 0:
        print("")
        print("Identified " + str(count) + " orthologs in reference dataset")
        print("Average dN among orthologs: " + str(ave(dnList)))
        print("Average dS among orthologs: " + str(ave(dsList)))
        print("Average dN/dS among orthologs: " + str(ave(dndsList2)))
        print("")
        os.system("rm -f 2NG.t 2NG.dN 2NG.dS rst1 rst 2ML.t 2ML.dN 2ML.dS 4fold.nuc rub")
        print("Pipeline finished without any crashes. Thanks for using pseudoHunter!")
    else:
        print("")
        print("Identified " + str(count) + " orthologs in reference dataset")
        print("It looks like no orthologs were identified below the specified dS threshold of: " + str(args.M))
        print("Please try running again, with a higher value for -M")
        print(
            "Dont forget to add the \'-o %s\' flag, so that you don't need to wait for codeml to run again." % args.out)


else:
    cwd = os.getcwd()
    DIR = args.out + "/dnds-analysis"

    if args.gff != "NA":
        gff = open(args.gff)
        gffDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'NA')))
        for i in gff:
            ls = i.rstrip().split("\t")
            if re.match(r'##FASTA', i):
                break
            else:
                if not re.match(r'#', i):
                    contig = ls[0]
                    orf = ls[8]
                    # if not re.findall(r'Alias', ls[8]):
                    orf = orf.split(";")[0]
                    # else:
                    #     orf = orf.split(";")[1]
                    orf = orf.split("=")[1]

                    product = lastItem(ls[8].split(";")).split("=")[1]
                    product = replace(product, [","], ";")
                    gffDict[orf]["product"] = product
                    gffDict[orf]["contig"] = contig
                    gffDict[orf]["start"] = ls[3]
                    gffDict[orf]["end"] = ls[4]
                    gffDict[orf]["strand"] = ls[6]

    else:
        gffDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'NA')))

        faa = open(args.q + "-proteins.faa")
        faa = fasta(faa)

        for i in faa.keys():
            ls = i.split(" # ")
            start = ls[1]
            end = ls[2]
            contig = allButTheLast(ls[0], "_")
            gffDict[ls[0]]["contig"] = contig
            gffDict[ls[0]]["start"] = start
            gffDict[ls[0]]["end"] = end
            gffDict[ls[0]]["product"] = "NA"
            if ls[3] == "-1":
                gffDict[ls[0]]["strand"] = "-"
            else:
                gffDict[ls[0]]["strand"] = "+"

    if args.ra != "NA":
        faaRef = open(args.ra)
        faaRef = fasta2(faaRef)
    else:
        faaRef = open(args.r + "-proteins.faa")
        faaRef = fasta2(faaRef)

    if args.rn != "NA":
        fnaRef = open(args.rn)
        fnaRef = fasta2(fnaRef)
    else:
        fnaRef = open(args.r + "-proteins.fna")
        fnaRef = fasta2(fnaRef)

    if args.a != "NA" and mode == 0:
        faa = open(args.a)
        faa = fasta2(faa)

    else:
        faa = open(args.q + "-proteins.faa")
        faa = fasta2(faa)

    if args.n != "NA" and mode == 0:
        fna = open(args.n)
        fna = fasta2(fna)

    else:
        fna = open(args.q + "-proteins.fna")
        fna = fasta2(fna)

    print("\n\n")

    alnLengthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    alnIdDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    blast = open(args.out + "/pseudogene.blast")
    for i in blast:
        ls = i.rstrip().split("\t")
        alnLengthDict[ls[0]][ls[1]] = ls[3]
        alnIdDict[ls[0]][ls[1]] = ls[2]

    print("summarizing codeml output")
    codealign = os.listdir(DIR)
    dndsDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in codealign:
        if re.findall(r'mlc', i):
            file = open(DIR + "/%s" % i, "r")
            for j in file:
                if re.search(r'#1', j):
                    ls = (j.rstrip().split(" "))
                    orf = ls[1]
                if re.search(r'#2', j):
                    ls = (j.rstrip().split(" "))
                    NODE = ls[1]
                line = j.rstrip()
            ls = line.split("  ")
            dS = remove(lastItem(ls), [" ", "=", "d", "S"])
            dN = remove(lastItem(ls[0:len(ls) - 1]), [" ", "=", "d", "N"])
            dndsDict[NODE]["orf"] = orf
            dndsDict[NODE]["dn"] = dN
            dndsDict[NODE]["ds"] = dS

    count = 0
    dndsList = []
    dndsDict2 = defaultdict(list)
    for i in sorted(dndsDict.keys()):
        count += 1
        if float(dndsDict[i]["dn"]) <= args.M and float(dndsDict[i]["ds"]) <= args.M and float(
                dndsDict[i]["ds"]) >= args.m and float(dndsDict[i]["dn"]) >= args.m:
            orf = dndsDict[i]["orf"]
            dn = dndsDict[i]["dn"]
            ds = dndsDict[i]["ds"]
            try:
                dnds = float(dndsDict[i]["dn"]) / float(dndsDict[i]["ds"])
            except ZeroDivisionError:
                dnds = "NA"
            dndsDict2[orf].append(i)
            if dnds != "NA":
                dndsList.append(dnds)

    print("preparing final output file: summary.csv")
    dsList = []
    dnList = []
    dndsList2 = []
    out = open(args.out + "/summary.csv", "w")
    out.write(
        "ORF_calls" + "," + "Ortholog" + "," + "Pseudogene" + "," + "contig" + "," + "start" + "," + "end" + "," + "strand" + "," + "geneLength" + "," + "AlignmentLength" + "," +
        "OrthologLength" + "," + "Identity" + "," + "Annotation" + "," + "Pseudogene_confidence" + "," +
        "NumberOfGeneFrags" + "," + "AlignmentLength/OrthologLength" + "," + "geneLength/OrthologLength" + "," + "dN" + "," + "dS" + "," + "dN/dS" + "," +
        "Translation" + "," + "Sequences" + "\n")


    delimiter = ''
    for i in dndsDict2.keys():
        if len(dndsDict2[i]) > 1:
            dndsDict3 = defaultdict(list)
            for k in dndsDict2[i]:
                dndsDict3[allButTheLast(k, args.delim)].append(int(lastItem(k.split(args.delim))))
            for l in dndsDict3.keys():
                if l != "":
                    listOfLists = (cluster(dndsDict3[l], 2))
                    for m in listOfLists:
                        try:
                            if len(m) > 1:
                                seq = ''
                                seq2 = ''
                                dnLS = []
                                dsLS = []
                                annotations = ''
                                ORFs = ''
                                TotalAlnLength = 0

                                idLS = []

                                for n in m:
                                    originalN = stabilityCounter(n)
                                    originalN = n
                                    ORF = (l + args.delim + str(originalN))

                                    ORFs += ORF + "|"

                                    annotation = gffDict[ORF]["product"]
                                    annotations += annotation + "|"

                                    seq += faa[ORF]
                                    seq2 += fna[ORF]
                                    dnLS.append(float(dndsDict[ORF]["dn"]))
                                    dsLS.append(float(dndsDict[ORF]["ds"]))

                                    alnLength = int(alnLengthDict[ORF][i])
                                    TotalAlnLength += alnLength

                                    identity = float(alnIdDict[ORF][i])
                                    idLS.append(identity)

                                identity = ave(idLS)

                                dn = ave(dnLS)
                                ds = ave(dsLS)
                                try:
                                    dnds = dn / ds
                                except ZeroDivisionError:
                                    dnds = "NA"
                                fragments = len(m)
                                ORF = ORFs[0:len(ORFs) - 1]
                                annotation = annotations

                                strand = gffDict[ORF.split("|")[0]]["strand"]
                                contig = gffDict[ORF.split("|")[0]]["contig"]
                                start = gffDict[ORF.split("|")[0]]["start"]
                                end = gffDict[lastItem(ORF.split("|"))]["end"]

                            else:

                                fragments = 1
                                originalM0 = stabilityCounter(m[0])
                                if delimiter == ".":
                                    originalM0 = m[0]
                                    ORF = (l + args.delim + str(originalM0))
                                else:
                                    ORF = (l + args.delim + str(originalM0))

                                contig = gffDict[ORF]["contig"]
                                start = gffDict[ORF]["start"]
                                end = gffDict[ORF]["end"]
                                strand = gffDict[ORF]["strand"]

                                annotation = gffDict[ORF]["product"]
                                seq = faa[ORF]
                                seq2 = fna[ORF]

                                TotalAlnLength = int(alnLengthDict[ORF][i])

                                identity = float(alnIdDict[ORF][i])

                                dn = dndsDict[ORF]["dn"]
                                ds = dndsDict[ORF]["ds"]

                                try:
                                    dnds = float(dn) / float(ds)
                                except ZeroDivisionError:
                                    dnds = "NA"

                        except (TypeError, ValueError):
                            if len(m) > 1:
                                seq = ''
                                seq2 = ''
                                dnLS = []
                                dsLS = []
                                ORFs = ''
                                TotalAlnLength = 0

                                idLS = []

                                for n in m:
                                    # originalN = stabilityCounter(n)
                                    originalN = n
                                    ORF = (l + args.delim + str(originalN))

                                    ORFs += ORF + "|"

                                    annotation = gffDict[ORF]["product"]
                                    annotations += annotation + "|"
                                    seq += faa[ORF]
                                    seq2 += fna[ORF]
                                    dnLS.append(float(dndsDict[ORF]["dn"]))
                                    dsLS.append(float(dndsDict[ORF]["ds"]))

                                    alnLength = int(alnLengthDict[ORF][i])
                                    TotalAlnLength += alnLength

                                    identity = float(alnIdDict[ORF][i])
                                    idLS.append(identity)

                                identity = ave(idLS)

                                dn = ave(dnLS)
                                ds = ave(dsLS)

                                try:
                                    dnds = float(dn) / float(ds)
                                except ZeroDivisionError:
                                    dnds = "NA"

                                fragments = len(m)
                                ORF = ORFs[0:len(ORFs) - 1]
                                annotation = annotations

                                contig = gffDict[ORF.split("|")[0]]["contig"]
                                strand = gffDict[ORF.split("|")[0]]["strand"]
                                start = gffDict[ORF.split("|")[0]]["start"]
                                end = gffDict[lastItem(ORF.split("|"))]["end"]

                            else:

                                fragments = 1
                                originalM0 = stabilityCounter(m[0])

                                originalM0 = m[0]
                                ORF = (l + args.delim + str(originalM0))

                                strand = gffDict[ORF]["strand"]
                                contig = gffDict[ORF]["contig"]
                                start = gffDict[ORF]["start"]
                                end = gffDict[ORF]["end"]

                                annotation = gffDict[ORF]["product"]
                                seq = faa[ORF]
                                seq2 = fna[ORF]

                                TotalAlnLength = int(alnLengthDict[ORF][i])

                                identity = float(alnIdDict[ORF][i])

                                dn = dndsDict[ORF]["dn"]
                                ds = dndsDict[ORF]["ds"]

                                try:
                                    dnds = float(dn) / float(ds)
                                except ZeroDivisionError:
                                    dnds = "NA"

        else:
            fragments = 1
            ORF = dndsDict2[i][0]
            contig = gffDict[ORF]["contig"]
            strand = gffDict[ORF]["strand"]

            start = gffDict[ORF]["start"]
            end = gffDict[ORF]["end"]

            annotation = gffDict[ORF]["product"]
            seq = faa[ORF]
            seq2 = fna[ORF]
            dn = dndsDict[ORF]["dn"]
            ds = dndsDict[ORF]["ds"]

            try:
                dnds = float(dn) / float(ds)
            except ZeroDivisionError:
                dnds = "NA"

            TotalAlnLength = int(alnLengthDict[ORF][i])

            identity = float(alnIdDict[ORF][i])

        ratio = TotalAlnLength / len(faaRef[i])

        out.write(ORF + "," + i + ",")
        if dnds != "NA":
            if dnds > args.d or ratio < args.l or fragments > 1 or (len(seq) / len(faaRef[i])) > args.L:
                out.write("Y" + ",")
                prob = 1

                # GETTING MY STUPID PSEUDOGENE SCORE CALCULATED
                if dnds != "NA":
                    inflation = float(dnds) / ave(dndsList)
                else:
                    inflation = ave(dndsList)

                prob = prob * inflation
                prob = prob * float(fragments)
                if ratio < 1:
                    prob = prob / ratio
                else:
                    prob = prob * ratio

            else:
                out.write("N" + ",")

                # GETTING MY STUPID PSEUDOGENE SCORE CALCULATED
                prob = 1

                if dnds != "NA":
                    inflation = float(dnds) / ave(dndsList)
                else:
                    inflation = ave(dndsList)

                prob = prob * inflation
                prob = prob * float(fragments)
                if ratio < 1:
                    prob = prob / ratio
                else:
                    prob = prob * ratio
        else:
            if ratio < args.l or fragments > 1 or (len(seq) / len(faaRef[i])) > args.L:
                out.write("Y" + ",")
                prob = 1

                # GETTING MY STUPID PSEUDOGENE SCORE CALCULATED
                if dnds != "NA":
                    inflation = float(dnds) / ave(dndsList)
                else:
                    inflation = ave(dndsList)

                prob = prob * inflation
                prob = prob * float(fragments)
                if ratio < 1:
                    prob = prob / ratio
                else:
                    prob = prob * ratio

            else:
                out.write("N" + ",")

                # GETTING MY STUPID PSEUDOGENE SCORE CALCULATED
                prob = 1

                if dnds != "NA":
                    inflation = float(dnds) / ave(dndsList)
                else:
                    inflation = ave(dndsList)

                prob = prob * inflation
                prob = prob * float(fragments)
                if ratio < 1:
                    prob = prob / ratio
                else:
                    prob = prob * ratio

        # WRITING TO FILE

        # print(seq)
        # print(seq2)
        # print(strand)
        # print("")
        out.write(
            str(contig) + "," + str(start) + "," + str(end) + "," + str(strand) + "," + str(int(end) - int(start)) + "," + str(
                TotalAlnLength * 3) + "," +
            str(len(faaRef[i]) * 3) + "," + str(identity) + "," + str(annotation) + "," + str(prob) + "," +
            str(fragments) + "," + str(ratio) + "," + str(((int(end) - int(start))/3) / len(faaRef[i])) + "," + str(dn) + "," +
            str(ds) + "," + str(dnds) + "," + seq + "," + seq2 + "\n")

        dnList.append(dn)
        dsList.append(ds)
        if dnds != "NA":
            dndsList2.append(dnds)
    out.close()

    # INTERGENIC REGION ANALYSIS
    if args.q != "NA":
        print("Analyzing intergenic regions")
        contigs = open(args.q)
        contigs = fasta2(contigs)

        firstRow = ''
        count = 0
        summaryDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        summaryDict2 = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        summaryDict3 = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        summary = open(args.out + "/summary.csv")
        for i in summary:
            ls = i.rstrip().split(",")
            if ls[0] != "ORF_calls":
                count += 1
                contig = ls[3]
                start = int(ls[4])
                summaryDict[contig][start] = ls
                summaryDict2[ls[1]] = ls[0]

                summaryDict3[ls[0]]["ORF_calls"] = ls[0]
                summaryDict3[ls[0]]["Ortholog"] = ls[1]
                summaryDict3[ls[0]]["Pseudogene"] = ls[2]
                summaryDict3[ls[0]]["contig"] = ls[3]
                summaryDict3[ls[0]]["start"] = ls[4]
                summaryDict3[ls[0]]["end"] = ls[5]
                summaryDict3[ls[0]]["strand"] = ls[6]
                summaryDict3[ls[0]]["geneLength"] = ls[7]
                summaryDict3[ls[0]]["AlignmentLength"] = ls[8]
                summaryDict3[ls[0]]["OrthologLength"] = ls[9]
                summaryDict3[ls[0]]["Identity"] = ls[10]
                summaryDict3[ls[0]]["Annotation"] = ls[11]
                summaryDict3[ls[0]]["Pseudogene_confidence"] = ls[12]
                summaryDict3[ls[0]]["NumberOfGeneFrags"] = ls[13]
                summaryDict3[ls[0]]["AlignmentLength/OrthologLength"] = ls[14]
                summaryDict3[ls[0]]["geneLength/OrthologLength"] = ls[15]
                summaryDict3[ls[0]]["dN"] = ls[16]
                summaryDict3[ls[0]]["dS"] = ls[17]
                summaryDict3[ls[0]]["dN/dS"] = ls[18]
                summaryDict3[ls[0]]["Translation"] = ls[19]
                summaryDict3[ls[0]]["Sequences"] = ls[20]
            else:
                firstRow = i.rstrip()

        igDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        count = 0
        for i in summaryDict.keys():
            end = 0
            endORF = ""
            for j in (sorted(summaryDict[i])):
                ls = (summaryDict[i][j])
                contig = ls[3]
                start = int(ls[4])
                startOrf = ls[0]
                if start > end:
                    contigSeq = contigs[contig]
                    IG = contigSeq[end:start]
                    if len(IG) > 100:
                        count += 1
                        igDict[endORF + "!" + str(contig) + "-" + str(end) + "-" + str(start) + "!" + startOrf] = IG
                end = int(ls[5])
                endORF = ls[0]

            if end < len(contigSeq):
                IG = contigSeq[end:len(contigSeq)]
                if len(IG) > 100:
                    count += 1
                    igDict[endORF + "!" + str(contig) + "-" + str(end) + "-" + str(len(contigSeq)) + "!" + ""] = IG

        out = open(args.out + "/ig.fasta", "w")
        for i in igDict.keys():
            out.write(">" + i + "\n")
            out.write(igDict[i] + "\n")

        if args.rn != "NA":
            REF = args.rn
        else:
            REF = args.r + "-proteins.fna"

        os.system("makeblastdb -dbtype nucl -in %s -out %s" % (REF, REF))
        os.system("blastn -query %s/ig.fasta -db %s -outfmt 6 -out %s/ig.blast -evalue 1E-6 -num_threads 4 -max_target_seqs 20" % (args.out, REF, args.out))
        os.system("rm %s.nhr" % REF)
        os.system("rm %s.nin" % REF)
        os.system("rm %s.nsq" % REF)

        igFasta = open("%s/ig.fasta" % args.out)
        igFasta = fasta2(igFasta)
        igFastaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in igFasta.keys():
            header = (i.split("!")[1])
            igFastaDict[header] = igFasta[i]

        refFna = open(REF)
        refFna = fasta2(refFna)

        out = open("%s/ig.csv" % args.out, "w")
        out.write(firstRow + "\n")
        blast = open("%s/ig.blast" % args.out)
        for i in blast:
            ls = i.rstrip().split("\t")

            if int(ls[8]) < int(ls[9]):
                strand = "+"
            else:
                strand = "-"

            start = int(ls[6])
            end = int(ls[7])
            if end - start > 100:

                contig = ls[0].split("!")[1].split("-")[0]

                ORF = "IG" + str(contig) + "_" + str(start + int(ls[0].split("!")[1].split("-")[1])) + "-" + str(int(ls[0].split("!")[1].split("-")[1]) + end)

                percIdent = float(ls[2])
                alnLength = int(ls[3])

                FreeLivingOrtholog = (ls[1])

                ig = ls[0].split("!")[1]
                seq = igFastaDict[ig]
                HomologousRegion = seq[start - 1:end]

                lengthIntergenic = (int(ls[0].split("!")[1].split("-")[2]) - int(ls[0].split("!")[1].split("-")[1]))

                if len(summaryDict2[ls[1]]) > 0:
                    # THE ORTHOLOG IN FREE-LIVING RELATIVE HAD A PRIOR HIT TO A PREVIOSLY-PREDICTED ENDOSYMBIONT ORFS
                    PreExistingOrfHit = (summaryDict2[ls[1]])
                    hitLocation = (localize(PreExistingOrfHit, ls[0].split("!")))
                    if hitLocation > 2:
                        # THE PREVIOUSLY-PREDICTED ENDOSYMBIONT ORF THAT IS A
                        # BLAST-MATCH TO THE IG-HOMOLOGOUS ORTHOLOG IS NOT ADJACENT TO THE IG

                        out.write(ORF + "," + FreeLivingOrtholog + "," + "Y" + "," + str(contig) + "," + str(start) + "," +
                                  str(end) + "," + strand + "," + "NA" + "," + str(len(HomologousRegion)) + "," +
                                  str(len(refFna[FreeLivingOrtholog])) + "," + str(
                            percIdent) + "," + "NA" + "," + "NA" + "," + "NA" +
                                  "," + str(
                            alnLength / len(refFna[FreeLivingOrtholog])) + "," + "NA" + "," + "NA" + "," + "NA" +
                                  "," + "NA" + "," + "NA" + "," + HomologousRegion + "\n")

                    else:
                        # THE PREVIOUSLY-PREDICTED ENDOSYMBIONT ORF THAT IS A
                        # BLAST-MATCH TO THE IG-HOMOLOGOUS ORTHOLOG IS ADJACENT TO THE IG, AND SHOULD BE COMBINED WITH IT.

                        NewGeneLength = int(summaryDict3[PreExistingOrfHit]["geneLength"]) + len(HomologousRegion)
                        NewAlignLength = int(summaryDict3[PreExistingOrfHit]["AlignmentLength"]) + len(HomologousRegion)
                        NewIdent = (float(summaryDict3[PreExistingOrfHit]["Identity"]) + float(percIdent)) / 2
                        Frags = int(summaryDict3[PreExistingOrfHit]["NumberOfGeneFrags"]) + 1
                        newAlnRatio = NewAlignLength / len(refFna[FreeLivingOrtholog])
                        newGeneRatio = NewGeneLength / len(refFna[FreeLivingOrtholog])

                        if hitLocation == 0:
                            contigStartPosition = summaryDict3[PreExistingOrfHit]["start"]
                            contigEndPosition = int(summaryDict3[PreExistingOrfHit]["end"]) + len(HomologousRegion)
                            newORF = PreExistingOrfHit + "|" + ORF
                            newNucSeq = summaryDict3[PreExistingOrfHit]["Sequences"] + HomologousRegion

                            out.write(newORF + "," + FreeLivingOrtholog + "," + "Y" + "," + summaryDict3[PreExistingOrfHit]["contig"] + "," +
                                      str(contigStartPosition) + "," + str(contigEndPosition) + "," + strand + "," + str(NewGeneLength) + "," +
                                      str(NewAlignLength) + "," + str(len(refFna[FreeLivingOrtholog])) + "," + str(NewIdent) + "," +
                                      summaryDict3[PreExistingOrfHit]["Annotation"] + "," + summaryDict3[PreExistingOrfHit]["Pseudogene_confidence"] +
                                      "," + str(Frags) + "," + str(newAlnRatio) + "," + str(newGeneRatio) + "," +
                                      summaryDict3[PreExistingOrfHit]["dN"] + "," +
                                      summaryDict3[PreExistingOrfHit]["dS"] + "," + summaryDict3[PreExistingOrfHit]["dN/dS"] + "," + summaryDict3[PreExistingOrfHit]["Translation"] + "," +
                                      newNucSeq + "\n")

                        else:
                            contigStartPosition = int(summaryDict3[PreExistingOrfHit]["start"]) - len(HomologousRegion)
                            contigEndPosition = int(summaryDict3[PreExistingOrfHit]["end"])
                            newORF = ORF + "|" + PreExistingOrfHit
                            newNucSeq = HomologousRegion + summaryDict3[PreExistingOrfHit]["Sequences"]

                            out.write(newORF + "," + FreeLivingOrtholog + "," + "Y" + "," + summaryDict3[PreExistingOrfHit][
                                "contig"] + "," +
                                      str(contigStartPosition) + "," + str(contigEndPosition) + "," + strand + "," + str(NewGeneLength) + "," +
                                      str(NewAlignLength) + "," + str(len(refFna[FreeLivingOrtholog])) + "," + str(
                                NewIdent) + "," +
                                      summaryDict3[PreExistingOrfHit]["Annotation"] + "," + summaryDict3[PreExistingOrfHit][
                                          "Pseudogene_confidence"] +
                                      "," + str(Frags) + "," + str(newAlnRatio) + "," + str(newGeneRatio) + "," +
                                      summaryDict3[PreExistingOrfHit]["dN"] + "," +
                                      summaryDict3[PreExistingOrfHit]["dS"] + "," + summaryDict3[PreExistingOrfHit][
                                          "dN/dS"] + "," + summaryDict3[PreExistingOrfHit]["Translation"] + "," +
                                      newNucSeq + '\n')

                else:
                    # THE ORTHOLOG IN FREE-LIVING RELATIVE HAD NO PRIOR HITS TO ANY OF THE PREDICTED ENDOSYMBIONT ORFS
                    out.write(
                        ORF + "," + FreeLivingOrtholog + "," + "Y" + "," + str(contig) + "," + str(start) + "," + str(end) + "," + strand +
                        "," + "NA" + "," + str(len(HomologousRegion)) + "," + str(len(refFna[FreeLivingOrtholog])) + "," +
                        str(percIdent) + "," + "NA" + "," + "NA" + "," + "NA" + "," +
                        str(alnLength / len(refFna[FreeLivingOrtholog])) + "," + "NA" + "," + "NA" + "," +
                        "NA" + "," + "NA" + "," + "NA" + "," + HomologousRegion + '\n')

        out.close()

        igSummaryDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        igSummary = open("%s/ig.csv" % args.out)
        for i in igSummary:
            ls = i.rstrip().split(",")
            if ls[0] != "ORF_calls":
                igSummaryDict[ls[0]] = i.rstrip()

        out = open("%s/summary-final.csv" % args.out, "w")
        originalSummary = open("%s/summary.csv" % args.out)
        for i in originalSummary:
            counter = 0
            ls = i.rstrip().split(",")
            if ls[0] != "ORF_calls":
                for j in igSummaryDict.keys():
                    if ls[0].split("|")[0] in j.split("|"):
                        counter += 1
                        out.write(igSummaryDict[j] + "\n")
                if counter == 0:
                    out.write(i.rstrip() + "\n")
                else:
                    counter = 0
            else:
                out.write(i.rstrip() + "\n")

        for i in igSummaryDict.keys():
            if not re.findall(r'\|', i):
                out.write(igSummaryDict[i] + "\n")
        out.close()

        os.system("rm %s/summary.csv" % args.out)

        count = 0
        end = 0
        out = open("%s/summary.csv" % args.out, "w")
        summaryExtend = open("%s/summary-final.csv" % args.out)
        for i in summaryExtend:
            ls = i.rstrip().split(",")
            if (re.findall(r'IG', ls[0]) and not re.findall(r'\|', ls[0])):
                start = lastItem(ls[0].split("_")).split("-")[0]
                # start = (ls[0].split("_")[1].split("-")[0])
                if int(start) > int(end) - 10:
                    count += 1
                    out.write(i.rstrip() + "\n")
                end = lastItem(ls[0].split("_")).split("-")[1]
                # end = (ls[0].split("_")[1].split("-")[1])
            else:
                count += 1
                out.write(i.rstrip() + "\n")
        out.close()

        os.system("rm %s/ig.blast %s/ig.csv %s/ig.fasta" % (args.out, args.out, args.out))
        os.system("rm %s/summary-final.csv" % args.out)

    ##################################################################################################################
    if len(dsList) > 0:
        print("")
        print("Identified " + str(count) + " orthologs in reference dataset")
        print("Average dN among orthologs: " + str(ave(dnList)))
        print("Average dS among orthologs: " + str(ave(dsList)))
        print("Average dN/dS among orthologs: " + str(ave(dndsList2)))
        print("")

        os.system("rm -f 2NG.t 2NG.dN 2NG.dS rst1 rst 2ML.t 2ML.dN 2ML.dS 4fold.nuc rub")
        print("Pipeline finished without any crashes. Thanks for using pseudoHunter!")
    else:
        print("")
        print("Identified " + str(count) + " orthologs in reference dataset")
        print("It looks like no orthologs were identified below the specified dS threshold of: " + str(args.M))
        print("Please try running again, with a higher value for -M")
        print(
            "Dont forget to add the \'-o %s\' flag, so that you don't need to wait for codeml to run again." % args.out)

######################################################################################################################
