"""
This script counts alleles from mpileup output.
Note that it fails in irregular cases (e.g. insertions, deletions, indels, triallelic SNPs, or if there is no coverage at that SNP.
You must alter the mpileup output to a "," or something similar for it to work in the above cases.
Published with permission (written by Choongwon Jeong).
"""


#! /usr/bin/env python
import re, sys, os, gzip
arg = sys.argv

## Read in arguments
inFile = arg[1]
idlist = arg[2]
markerFile = arg[3]
outfile = arg[4]

def column(mat, i):
    return [row[i] for row in mat]


def parse_pileup(pileupstring):
    s1 = ""; read_end = 0
    for n in pileupstring:
        if n == "^" or n == "$":
            read_end = 1; continue
        elif read_end == 1:
            read_end = 0; continue
        else:
            s1 += n
    return s1


r1 = os.getcwd() + "/"; os.chdir(r1)

idlist2=[line.strip().split() for line in open(idlist, "r").readlines()]
ids = [val for val in idlist2[0][0].split(",")]
info = [line.strip().split() for line in open(markerFile, "r").readlines()[1:]]
info_header = [line.strip().split() for line in open(markerFile, "r").readlines()[0:1]][0]
tpileup = [line.strip().split("\t") for line in open(inFile + ".pileup", "r").readlines()]
nrec = max([len(tp) for tp in tpileup])

pileup = []
for tp in tpileup:
    if len(tp) == nrec:
        pileup.append(["*" if val == "" else val for val in tp])
    else:
        pileup.append(["*" if val == "" else val for val in tp])
        pileup[-1].extend(['*'] * (nrec - len(tp)))


pileup_pos = [val for val in column(pileup, 1)]
nind = (len(pileup[0]) - 3) / 3
nullvec = ["0,0", "0"] * nind

d1s = []
for v in info:
    vnum = [j for j,pp in enumerate(pileup_pos) if pp == v[2]]
    if len(vnum) == 0:
        d1s.append([val for val in v])
        d1s[-1].extend([val for val in nullvec])
        continue
    pvec = [val for val in pileup[vnum[0]]]; a_ref = pvec[2]
    pstr = [parse_pileup(pvec[4+j*3]) for j in range(nind)]
    td1 = [val for val in v]
    for pval in pstr:
        if pval == "*":
            td1.extend(["0,0", "0"])
        elif pval.count("-") > 0 or pval.count("+") > 0:
            "Indels are included in the pileup!"; break
        else:
            pval2 = ''.join([a_ref if val == "." or val == "," else val.upper() for val in pval])
            tnum1 = pval2.count(v[4])
            tnum2 = pval2.count(v[5])
            td1.extend([str(tnum1) + "," + str(tnum2), str(tnum1+tnum2)])
    d1s.append(td1)


hv = '\t'.join(info_header)
hv += "\t" + '\t'.join(sum([[v + ".AD", v + ".DP"] for v in ids], []))

F1 = open(outfile + ".allelecounts.txt", "w")
F1.writelines(hv + "\n")
F1.writelines('\n'.join(['\t'.join(d1) for d1 in d1s]) + "\n")
F1.close()
