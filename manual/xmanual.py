#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2022-04-29 17:51:22 sander>

'''
xmanual: eXecute pdfLaTeX to create the KPP User Manual
Author: Rolf Sander, 2022-...
'''

import os, sys
assert sys.version_info >= (3, 6)
from glob import glob
import re

DONTEDIT = 'created automatically by xmanual.py, DO NOT EDIT!'
LATEXFILE = 'kpp_UserManual'

##############################################################################

def grep_i(searchstring, filewildcard, only_matching=False):
    # suffix _i indicates grep option '-i' for case insensitive
    results = []
    allfiles = glob(filewildcard)
    for onefile in allfiles:
        with open(onefile) as f:
            for line in f:
                result = re.search(searchstring, line, re.IGNORECASE)
                if (result):
                    if (only_matching):
                        results.append(result.group(0)) # return only match (grep -o)
                    else:
                        results.append(line.rstrip('\n')) # return complete line
    return results

##############################################################################

def clean():
    suffixes = ('log', 'toc', 'aux', 'bbl', 'blg')
    for suffix in suffixes:
        tmpfile = LATEXFILE + '.' + suffix
        if (os.path.isfile(tmpfile)):
            print(f'removing temporary file {tmpfile}')
            os.remove(tmpfile)
    os.remove('kpp_info.tex')

##############################################################################

def xmanual(cleanup=True):

    def run_pdflatex(filename):
        print(f'Running pdflatex {filename}...')
        os.system(f'pdflatex -halt-on-error {filename}.tex')

    # create kpp_info.tex which contains the kpp version number:
    current_version = re.sub('.*"(.*)".*', r'\1',
    grep_i('#define KPP_VERSION', '../src/gdata.h')[0])
    with open('kpp_info.tex', 'w') as infofile:
        print(f'% {DONTEDIT}', file=infofile)
        print(fr'\def\kppversion{{{current_version}}}', file=infofile)

    run_pdflatex(LATEXFILE)
    # BibTeX:
    # create *.bbl file:
    print('Running bibtex...')
    os.system(f'bibtex {LATEXFILE}')
    with open(LATEXFILE+'.blg') as f:
        for line in f:
            if ("Warning--I didn't find a database entry" in line):
                sys.exit(1)
    run_pdflatex(LATEXFILE)
    run_pdflatex(LATEXFILE)

    if (cleanup):
        clean()

##############################################################################

if __name__ == '__main__':

    if (len(sys.argv)>1):
        if (sys.argv[1]=='--keep'):
            xmanual(cleanup=False)
    else:
        xmanual()

##############################################################################
