#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os, signal, random, re, fnmatch, math, operator, requests, json
import time, locale, datetime


# nslides = 357
# filename = "20160121.exp4.12000.B3000.T0.10."
nslides = 463
filename = "20160121.exp3.12000.B3000.T0.10."
input_dir = "./movie/"

pre_latex = '''
\documentclass[10pt]{beamer}
\setbeamerfont{structure}{family=\\rmfamily}
\usepackage{amsthm}
\usepackage{graphicx}
\usepackage{graphics}
\usepackage{hyperref}
\\beamertemplatenavigationsymbolsempty
\setbeamertemplate{blocks}[rounded][shadow=true]
\setbeamertemplate{bibliography item}[text]
\setbeamertemplate{caption}[numbered]
\usetheme{default}
\usecolortheme{seahorse}
\mode<presentation>
{
   \setbeamercovered{transparent}
   \setbeamertemplate{items}[ball]
   \setbeamertemplate{theorems}[numbered]
   \setbeamertemplate{footline}[frame number]
}

\\begin{document}

'''

end_latex = '''
\end{document}

'''

std_frame = '''
\\begin{frame}
\includegraphics[scale=0.6]{'''

end_frame = '''}
\end{frame}

'''

fh = open(input_dir + filename + "tex", 'w')
fh.write(pre_latex)
for pi in xrange(1, nslides+1):
    fh.write("%s%s%02d.eps%s\n" % (std_frame, filename, pi, end_frame))

fh.write(end_latex)
fh.close()


