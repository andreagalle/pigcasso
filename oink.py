#!/usr/bin/env python2

# -*- coding: utf-8 -*-
"""
oink.py
"""

import os, sys, shutil

sys.dont_write_bytecode = True

import plotFields as pigcasso
import userPearls as pearls
import toolCase   as util

run_directory = "../../data/vtk/" ; run_version   = "a" ; filename_head = "PI" ; chk = False
res_directory = "../../data/fig/" 

if util.chk_dir(res_directory) == False: os.makedirs(res_directory)

fields_list = [] # put at least 'Uf_trz'

pigcasso.snapshots(run_directory,run_version,res_directory,filename_head,fields_list)

chk = pearls.throw(run_directory,run_version,res_directory)

if chk == False: print '\n--> no custom script provided by the user.'

