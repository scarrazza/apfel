#! /usr/bin/env python

import glob, sys

def initvals():
    # the empty info selftructure, that for each value of nf will store
    # the values of alphas
    initvals_list = {"3" : [],
                     "4" : [],
                     "5" : [],
                     "6" : []}
    return initvals_list

def main(path = ''):
    path = path.strip()
    if path and not path.endswith('/'):
        path = path+'/'

    # dict.s of info structs
    infodict = {"lo" : initvals(),
                "nlo" : initvals()}

    nfs = [str(_) for _ in range(3,7)]
        
    for fname in glob.glob(path + '*.table'):
        fname = fname.split("/")[-1]
        ftype = fname.split("xtable_")[1].split("_nf")[0]
        nf    = fname.split("_nf")[1].split("_alphas")[0]
        alpha = fname.split("_alphas")[1].split(".table")[0]
        if nf in nfs:
            infodict[ftype][nf].append(alpha)

        for ftype in infodict:
            for nf in nfs:
                with open(path+"{0}_nf{1}.info".format(ftype,nf),'w') as f:
                    for alphas in sorted(infodict[ftype][nf]):
                        f.write('0.' + alphas[1:] + '\n')

if __name__ == '__main__':
    if len(sys.argv)>1:
        main(sys.argv[1])
    else:
        main()
