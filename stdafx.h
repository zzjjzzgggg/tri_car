#pragma once


#include <stdio.h>
#include "/home/jzzhao/git_project/netsnap/snap/Snap.h"

#define HEPTH_W1K


#ifdef HEPTH
const static TStr ROOT = "../hepth/";
const static TStr GFNm = "../cit-HepTh.gz";
const static int W = 1000;
#endif

#ifdef HEPTH_W1K
const static TStr ROOT = "../hepth/";
const static TStr GFNm = "../cit-HepTh_W1K.gz";
const static int W = 1000;
#endif


