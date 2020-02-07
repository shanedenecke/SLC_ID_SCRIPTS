/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 

    trimAl v1.2: a tool for automated alignment trimming in large-scale 
                 phylogenetics analyses 

    Copyright (C) 2009 Capella-Gutierrez S. and Gabaldon, T.
                       [scapella, tgabaldon]@crg.es

    This file is part of trimAl.

    trimAl is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl. If not, see <http://www.gnu.org/licenses/>.

 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#include <fstream>
#include <iostream>
#include <iomanip>

#include <stdlib.h>
#include <string.h>

#include "compareFiles.h"
#include "alignment.h"
#include "utils.h"

#define DNAType 1
#define RNAType 2
#define AAType  3

#define GAPPYOUT 1
#define STRICT   2

#define VERSION 1.2
#define REVISION 59

void menu(void);
void examples(void);

int main(int argc, char *argv[]){

  /* Parameters Control */
  bool appearErrors = false, complementary = false, colnumbering = false, nogaps = false, noallgaps = false, gappyout = false,
       strict = false, strictplus = false, automated1 = false, sgc = false, sgt = false, scc = false, sct = false, sfc = false,
       sft = false, sident = false;
  float conserve = -1, gapThreshold = -1, simThreshold = -1, comThreshold = -1, resOverlap = -1, seqOverlap = -1;
  int outformat = -1, prevType = -1, compareset = -1, stats = 0, select = -1, windowSize = -1, gapWindow = -1, simWindow = -1, conWindow = -1;

  /* Others varibles */
  ifstream compare;
  string nline;
  float *compareVect = NULL; 
  alignment *origAlig = NULL, *singleAlig = NULL;
  alignment **compAlig  = NULL;
  similarityMatrix *similMatrix = NULL;
  int i = 1, j, prev, lng, num = 0, maxAminos = 0, numfiles = 0, referFile = 0, *delColumns = NULL;
  char c, *infile = NULL, *outfile = NULL, *outhtml = NULL, *matrix = NULL, **filesToCompare = NULL, *aux = NULL, *str = NULL, line[256];

  /* ------------------------------------------------------------------------------------------------------ */

  /* Exec: TrimAl - Shows the menu. */

  /* ------------------------------------------------------------------------------------------------------ */
  if(argc == 1) { 
    menu();
    return 0;
  }

  /* ------------------------------------------------------------------------------------------------------ */

  /*                                        Help and Version Menu                                           */

  /* ------------------------------------------------------------------------------------------------------ */
  if(!strcmp(argv[i], "-h") && (i+1 == argc)) {
    menu(); examples();
    return 0;
  }

  if(!strcmp(argv[i], "--version") && (i+1 == argc)) {
    cout << endl << "trimAl " << VERSION << "rev" << REVISION << endl << endl;
    return 0;
  }

  /***** ***** ***** ***** ***** ***** ***** Parameters Processing ***** ***** ***** ***** ***** ***** *****/
  origAlig = new alignment;

  while(i < argc) {

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                Input and Output files and format output                                */

   /* Option -in ------------------------------------------------------------------------------------------- */
    if(!strcmp(argv[i], "-in") && (i+1 != argc) && (infile == NULL)) {

      if((sfc) || (sft) || (comThreshold != -1)) {
        cerr << endl << "ERROR: Not allowed in combination of file comparision." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if(compareset == -1) {
        lng = strlen(argv[++i]);
        infile = new char[lng + 1];
        strcpy(infile, argv[i]);

        if(!origAlig -> loadAlignment(infile)) {
          cerr << endl << "ERROR: Alignment not loaded: \"" << infile << "\" Check the file's content." << endl << endl;
          appearErrors = true;
        }
      }

      else {
        cerr << endl << "ERROR: Option \"" << argv[i] << "\" not valid. A reference file exists with alignments to compare." << endl << endl;
        appearErrors = true;
        i++;
      } 
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -out ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-out")) && (i+1 != argc) && (outfile == NULL)) {
      lng = strlen(argv[++i]);
      outfile = new char[lng + 1];
      strcpy(outfile, argv[i]);
    }

   /* Option -htmlout -------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-htmlout")) && (i+1 != argc) && (outhtml == NULL)) {
      lng = strlen(argv[++i]);
      outhtml = new char[lng + 1];
      strcpy(outhtml, argv[i]);
    }

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                           Output File format                                           */

   /* Option -clustal -------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-clustal") && (outformat == -1))
      outformat = 1;

   /* Option -fasta -------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-fasta") && (outformat == -1))
      outformat = 8;

   /* Option -nbrf ------------------------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-nbrf") && (outformat == -1))
      outformat = 3;

   /* Option -nexus ------------------------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-nexus") && (outformat == -1))
      outformat = 17;

   /* Option -mega ------------------------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-mega") && (outformat == -1))
      outformat = 21;

   /* Option -phylip3.2 ------------------------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-phylip3.2") && (outformat == -1))
      outformat = 11;

   /* Option -phylip ------------------------------------------------------------------------------------ */
    else if(!strcmp(argv[i], "-phylip") && (outformat == -1))
      outformat = 12;

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                         Similarity Matrix File                                         */

   /* Option -matrix --------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-matrix") && (i+1 != argc) && (matrix == NULL)) {
      lng = strlen(argv[++i]);
      matrix = new char[lng + 1];
      strcpy(matrix, argv[i]);
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /*                                   File with a alignments' set to compare                               */

   /* Option -compareset ----------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-compareset") && (i+1 != argc) && (compareset == -1)) {

      if(infile == NULL) {
        compare.open(argv[++i], ifstream::in);
        if(!compare) {
          cerr << endl << "ERROR: Check the reference file with the alignments to compare." << endl << endl;
          appearErrors = true;
        }

        while(compare.getline(line, 256)) numfiles++;
        compare.close();

        compareset = i;
      }

      else {
        cerr << endl << "ERROR: Option \"" << argv[i] << "\" not valid. A reference file exists with alignments to compare." << endl << endl;
        appearErrors = true;
        i++;
      }
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /*                                           Manual Method Values                                         */

   /* Option -gt, gapthreshold ----------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-gapthreshold") || !strcmp(argv[i], "-gt")) && (i+1 != argc) && (gapThreshold == -1)) {

      if(delColumns != NULL) {
        cerr << endl << "ERROR: Not allowed in combination of other manual methods." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
        cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        appearErrors = true;
        i++;
      }

      else {
        if(utils::isNumber(argv[i+1])) {
          gapThreshold = 1 - atof(argv[++i]);
          if((gapThreshold < 0) || (gapThreshold > 1)) {
            cerr << endl << "ERROR: The gap threshold value should be between 0 and 1." << endl << endl;
            appearErrors = true;
            i++;
          }
        }
        else {
          cerr << endl << "ERROR: The gap threshold value should be a positive real number." << endl << endl;
          appearErrors = true;
          i++;
        }
      }
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -st -simthreshold ----------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-simthreshold") || !strcmp(argv[i], "-st")) && (i+1 != argc) && (simThreshold == -1)) {

      if(delColumns != NULL) {
        cerr << endl << "ERROR: Not allowed in combination of other manual methods." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
        cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        appearErrors = true;
        i++;
      }

      else {
        if(utils::isNumber(argv[i+1])) {
          simThreshold = atof(argv[++i]);
          if((simThreshold < 0) || (simThreshold > 1)) {
            cerr << endl << "ERROR: The similarity threshold value should be between 0 and 1." << endl << endl;
            appearErrors = true;
            i++;
          }
        }
        else {
          cerr << endl << "ERROR: The similarity threshold value should be a positive real number." << endl << endl;
          appearErrors = true;
          i++;
        }
      }
    }
   /* ------------------------------------------------------------------------------------------------------ */


   /* Option -ct -conthreshold ----------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-conthreshold") || !strcmp(argv[i], "-ct")) && (i+1 != argc) && (comThreshold == -1)) {

      if(delColumns != NULL) {
        cerr << endl << "ERROR: Not allowed in combination of other manual methods." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
        cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        appearErrors = true;
        i++;
      }

      else if(infile != NULL) {
        cerr << endl << "ERROR: Not allowed in combination with -in option." << endl << endl;
        appearErrors = true;
        i++;
      }

      else {
        if(utils::isNumber(argv[i+1])) {
          comThreshold = atof(argv[++i]);
          if((comThreshold < 0) || (comThreshold > 1)) {
            cerr << endl << "ERROR: The consistency threshold value should be between 0 and 1." << endl << endl;
            appearErrors = true;
            i++;
          }
        }
        else {
          cerr << endl << "ERROR: The consistency threshold value should be a positive real number." << endl << endl;
          appearErrors = true;
          i++;
        }
      }
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -cons ----------------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-cons")) && (i+1 != argc) && (conserve == -1)) {

      if(delColumns != NULL) {
        cerr << endl << "ERROR: Not allowed in combination of other manual methods." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus)  || (automated1)) {
        cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        appearErrors = true;
        i++;
      }

      else {
        if(utils::isNumber(argv[i+1])) {
          conserve = atof(argv[++i]);
          if((conserve < 0) || (conserve > 100)) {
            cerr << endl << "ERROR: The minimal positions value should be between 0 and 100." << endl << endl;
            appearErrors = true;
            i++;
          }
       }
        else {
          cerr << endl << "ERROR: The minimal positions value should be a positive real number." << endl << endl;
          appearErrors = true;
          i++;
        }
      }
    }

   /* Option -coloverlap ----------------------------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-resoverlap")) && (i+1 != argc) && (resOverlap == -1)) {

      if(delColumns != NULL) {
        cerr << endl << "ERROR: Not allowed in combination of other methods." << endl << endl;
        appearErrors = true;
        i++;
      }

      else {
        if(utils::isNumber(argv[i+1])) {
          resOverlap = atof(argv[++i]);
          if((resOverlap < 0) || (resOverlap > 1)) {
            cerr << endl << "ERROR: The residue overlap value should be between 0 and 1." << endl << endl;
            appearErrors = true;
            i++;
          }
       }
        else {
          cerr << endl << "ERROR: The residue overlap value should be a positive real number." << endl << endl;
          appearErrors = true;
          i++;
        }
      }
    }

   /* Option -seqoverlap ----------------------------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-seqoverlap")) && (i+1 != argc) && (seqOverlap == -1)) {

      if(delColumns != NULL) {
        cerr << endl << "ERROR: Not allowed in combination of other methods." << endl << endl;
        appearErrors = true;
        i++;
      }

      else {
        if(utils::isNumber(argv[i+1])) {
          seqOverlap = atof(argv[++i]);
          if((seqOverlap < 0) || (seqOverlap > 100)) {
            cerr << endl << "ERROR: The sequences overlap value should be between 0 and 100." << endl << endl;
            appearErrors = true;
            i++;
          }
       }
        else {
          cerr << endl << "ERROR: The minimal positions value should be a positive real number." << endl << endl;
          appearErrors = true;
          i++;
        }
      }
    }

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                           Windows Size Values                                           */

   /* Option -w -------------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-w") && (i+1 != argc) && (windowSize == -1)){

      if((gapWindow != -1) || (simWindow != -1) || (conWindow != -1)) {
        cerr << endl << "ERROR: Not allowed in combination with this specific window value." << endl << endl;
        appearErrors = true;
      }

      else if(delColumns != NULL) {
        cerr << endl << "ERROR: It's imposible to use this windows size in combination of selection method." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
        cerr << endl << "ERROR: Not allowed in combination of automatic methods." << endl << endl;
        appearErrors = true;
        i++;
      }

      else {
        if(utils::isNumber(argv[i+1])) {
          windowSize = atoi(argv[++i]);
          if(windowSize <= 0){
            cerr << endl << "ERROR: The window value should be a positive integer number." << endl << endl;
            appearErrors = true;
            i++;
          }
        }
        else {
          cerr << endl << "ERROR: The window value should be a number." << endl << endl;
          appearErrors = true;
          i++;
        }
      }
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -gw -------------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-gw") && (i+1 != argc) && (gapWindow == -1)){

      if(windowSize != -1) {
        cerr << endl << "ERROR: Not allowed in combination of general window value." << endl << endl;
        appearErrors = true;
      }

      else if(delColumns != NULL) {
        cerr << endl << "ERROR: It's imposible to use this windows size in combination of selection method." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
        cerr << endl << "ERROR: Not allowed in combination of automatic methods." << endl << endl;
        appearErrors = true;
        i++;
      }

      else {
        if(utils::isNumber(argv[i+1])) {
          gapWindow = atoi(argv[++i]);
          if(gapWindow <= 0){
            cerr << endl << "ERROR: The window value should be a positive integer number." << endl << endl;
            appearErrors = true;
            i++;
          }
        }
        else {
          cerr << endl << "ERROR: The window value should be a number." << endl << endl;
          appearErrors = true;
          i++;
        }
      }
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -sw -------------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-sw") && (i+1 != argc) && (simWindow == -1)){

      if(windowSize != -1) {
        cerr << endl << "ERROR: Not allowed in combination of general window value." << endl << endl;
        appearErrors = true;
      }

      else if(delColumns != NULL) {
        cerr << endl << "ERROR: It's imposible to use this windows size in combination of selection method." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
        cerr << endl << "ERROR: Not allowed in combination of automatic methods." << endl << endl;
        appearErrors = true;
        i++;
      }

      else {
        if(utils::isNumber(argv[i+1])) {
          simWindow = atoi(argv[++i]);
          if(simWindow <= 0){
            cerr << endl << "ERROR: The window value should be a positive integer number." << endl << endl;
            appearErrors = true;
            i++;
          }
        }
        else {
          cerr << endl << "ERROR: The window value should be a number." << endl << endl;
          appearErrors = true;
          i++;
        }
      }
    }

   /* Option -cw -------------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-cw") && (i+1 != argc) && (conWindow == -1)){

      if(windowSize != -1) {
        cerr << endl << "ERROR: Not allowed in combination of general window value." << endl << endl;
        appearErrors = true;
      }

      else if(delColumns != NULL) {
        cerr << endl << "ERROR: It's imposible to use this windows size in combination of selection method." << endl << endl;
        appearErrors = true;
        i++;
      }

      else {
        if(utils::isNumber(argv[i+1])) {
          conWindow = atoi(argv[++i]);
          if(conWindow <= 0){
            cerr << endl << "ERROR: The window value should be a positive integer number." << endl << endl;
            appearErrors = true;
            i++;
          }
        }
        else {
          cerr << endl << "ERROR: The window value should be a number." << endl << endl;
          appearErrors = true;
          i++;
        }
      }
    }
   /* ------------------------------------------------------------------------------------------------------ */
   /* ------------------------------------------------------------------------------------------------------ */

   /*                                            Automatics Method                                           */

   /* Option -nogaps --------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-nogaps") && (!nogaps)) {

      if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
        cerr << endl << "ERROR: Not allowed in combination of window values." << endl << endl;
        appearErrors = true;
      }

      else if((noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
        cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || 
        (comThreshold != -1) || (delColumns != NULL)) {
        cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        appearErrors = true;
        i++;
      }

      else
        nogaps = true;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -noallgaps --------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-noallgaps") && (!noallgaps)) {

      if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
        cerr << endl << "ERROR: Not allowed in combination of window values." << endl << endl;
        appearErrors = true;
      }

      else if((nogaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
        cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || 
        (comThreshold != -1) || (delColumns != NULL)) {
        cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        appearErrors = true;
        i++;
      }

      else
        noallgaps = true;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -gappyout ------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-gappyout") && (!strict)) {

      if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
        cerr << endl << "ERROR: Not allowed in combination of window values." << endl << endl;
        appearErrors = true;
      }

      else if((nogaps) || (noallgaps) || (strict) || (strictplus) || (automated1)) {
        cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || 
        (comThreshold != -1) || (delColumns != NULL)) {
        cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        appearErrors = true;
        i++;
      }

      else
        gappyout = true;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -strict --------------------------------------------------------------------------------------- */
    else if(!strcmp(argv[i], "-strict") && (!strict)) {

      if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
        cerr << endl << "ERROR: Not allowed in combination of window values." << endl << endl;
        appearErrors = true;
      }

      else if((nogaps) || (noallgaps) || (gappyout) || (strictplus) || (automated1)) {
        cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || 
        (comThreshold != -1) || (delColumns != NULL)) {
        cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        appearErrors = true;
        i++;
      }

      else
        strict = true;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -strictplus ----------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-strictplus")) && (!strictplus)) {

      if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
        cerr << endl << "ERROR: Not allowed in combination with this window value." << endl << endl;
        appearErrors = true;
      }

      else if((nogaps) || (noallgaps) || (gappyout) || (strict) || (automated1)) {
        cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || 
        (comThreshold != -1) || (delColumns != NULL)) {
        cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        appearErrors = true;
        i++;
      }

      else
        strictplus = true;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -automated1 ----------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-automated1")) && (!automated1)) {

      if((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
        cerr << endl << "ERROR: Not allowed in combination with this window value." << endl << endl;
        appearErrors = true;
      }

      else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || 
        (comThreshold != -1) || (delColumns != NULL)) {
        cerr << endl << "ERROR: Combinations between automatic methods are not allowed." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || (comThreshold != -1) || (delColumns != NULL)) {
        cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed" << endl << endl;
        appearErrors = true;
        i++;
      }

      else
        automated1 = true;
    }

   /* ------------------------------------------------------------------------------------------------------ */

   /*                                               Statistics                                               */

   /* Option -sgc ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-sgc")) && (!sgc)) {
      sgc = true;
      stats--;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -sgt ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-sgt")) && (!sgt)) {
      sgt = true;
      stats--;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -scc ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-scc")) && (!scc)) {
      scc = true;
      stats--;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -sct ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-sct")) && (!sct)) {
      sct = true;
      stats--;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -sident --------------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-sident")) && (!sident)) {
      sident = true;
      stats--;
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -sfc ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-sfc")) && (!sfc)) {

      if(infile != NULL) {
        cerr << endl << "ERROR: Not allowed in combination with -in option." << endl << endl;
        appearErrors = true;
        i++;
      }

      else {
        sfc = true;
        stats--;
      }
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -sft ------------------------------------------------------------------------------------------ */
    else if((!strcmp(argv[i], "-sft")) && (!sft)) {

      if(infile != NULL) {
        cerr << endl << "ERROR: Not allowed in combination with -in option." << endl << endl;
        appearErrors = true;
        i++;
      }

      else {
        sft = true;
        stats--;
      }
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /*                                            Others parameters                                           */

   /* ------------------------------------------------------------------------------------------------------ */

   /* Option -complementary -------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-complementary")) && (complementary == false)) {
      complementary = true;
    }

   /* Option -colnumbering ------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-colnumbering")) && (colnumbering == false)) {
      colnumbering = true;
    }

   /* Option -select -------------------------------------------------------------------------------- */
    else if((!strcmp(argv[i], "-select")) && (select == -1) && ((i+3) < argc) && (!strcmp(argv[++i], "{")) && (!strcmp(argv[i+2], "}"))) {

      if((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
        cerr << endl << "ERROR: Combinations between automatic and manual methods are not allowed." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || (comThreshold != -1)) {
        cerr << endl << "ERROR: Not allowed in combination of other manual methods." << endl << endl;
        appearErrors = true;
        i++;
      }

      else if((windowSize != -1) || (gapWindow != -1)|| (simWindow != -1)) {
        cerr << endl << "ERROR: It's imposible to use this windows size in combination of selection method." << endl << endl;
        appearErrors = true;
        i++;
      }

      else {
        aux = new char[strlen(argv[++i]) + 1];
        aux[strlen(argv[i])] = '\0';

        strcpy(aux, argv[i]);
        str = strtok(aux, ",");

        num = 0;
        while(str != NULL) {
          num += 2;
          str = strtok(NULL, ",");
        }

        delete [] aux;
        delColumns = new int[num];
        num = 0;

        str = strtok(argv[i], ",");
        while(str != NULL) {

          for(j = 0, prev = 0; ((j < (int) strlen(str)) && (!appearErrors)); j++)
            if(str[j] == '-') 
              prev = j + 1;
            else if(!isdigit(str[j]))
              appearErrors = true;

          if((appearErrors) || (prev == 1)) {
            cerr << endl << "ERROR: This option only accepts integer numbers between 0 and the alignment's length - 1." << endl << endl;
            appearErrors = true;
            break;
          }

          if(prev) {
            aux = new char[prev];
            aux[(prev - 1)] = '\0';
            strncpy(aux, str, prev - 1);
            if((num > 0) && (atoi(aux) < delColumns[num - 1]))
              appearErrors = true;
            else
              delColumns[num++] = atoi(aux);
            delete [] aux;
          }

          j = (j - prev) + 1;
          aux = new char[j];
          aux[(j - 1)] = '\0';
          strncpy(aux, &str[prev], j - 1);
          if((num > 0) && (atoi(aux) < delColumns[num -1]))
          appearErrors = true;
          else if(prev)
            delColumns[num++] = atoi(aux);
          else {
            delColumns[num++] = atoi(aux);
            delColumns[num++] = atoi(aux);
          }
          delete [] aux;

          str = strtok(NULL, ",");

          if(appearErrors) {
            cerr << endl << "ERROR: The number of columns to remove should be consecutive." << endl << endl;
            break;
          }
        }
        i++;
      }
    }
   /* ------------------------------------------------------------------------------------------------------ */

   /*                                          Not Valids Parameters                                         */

   /* ------------------------------------------------------------------------------------------------------ */
    else {
      cerr << endl << "ERROR: Parameter \"" << argv[i] << "\" not valid." << endl << endl;
      appearErrors = true;
    }
   /* ------------------------------------------------------------------------------------------------------ */
    i++;

    if(appearErrors)
      break;

  }
   /* ------------------------------------------------------------------------------------------------------ */

   /*                                       Postprocessing Parameters                                        */

   /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors) && (infile != NULL)) {

    if(((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1) ||
      (gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || (delColumns != NULL) ||
      (resOverlap != -1) || (seqOverlap != -1) || (stats < 0)) &&
      (!origAlig -> isFileAligned())) {
        cerr << endl << "ERROR: The sequences in the input alignment should be aligned in order to use trimming method." << endl << endl;
        appearErrors = true;
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */
  if((!appearErrors) && (windowSize != -1) && (compareset != -1))
    cerr << "INFO: Try with specific comparison file window value. parameter -cw." << endl << endl;
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */
  if((matrix != NULL) && (!appearErrors)) {
    if((!strict) && (!strictplus) && (!automated1) && (simThreshold == -1.0) && (!scc) && (!sct)) {
      cerr << endl << "ERROR: The Similarity Matrix can only be used with methods that use this matrix." << endl << endl;
      appearErrors = true;
    }

    if((gapWindow != -1) ||((compareset == -1) && (conWindow != -1))) {
      cerr << endl << "ERROR: The Similarity Matrix can only be used with general/similarity windows size." << endl << endl;
      appearErrors = true;
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */
  if((complementary) && (!appearErrors))
    if((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1)
      && (gapThreshold == -1) && (conserve == -1) && (simThreshold == -1) && (delColumns == NULL)) {
      cerr << endl << "ERROR: This parameter can only be used with either an automatic or a cleanning method." << endl << endl;
      appearErrors = true;
    }
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */
  if((colnumbering) && (!appearErrors)) {
    if((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1)
      && (gapThreshold == -1) && (conserve == -1) && (simThreshold == -1) && (delColumns == NULL)) {
      cerr << endl << "ERROR: This parameter can only be used with either a trimming method." << endl << endl;
      appearErrors = true;
    }
    else if(stats < 0) {
      cerr << endl << "ERROR: This parameter is not valid when statistics' parameters are defined." << endl << endl;
      appearErrors = true;
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */
  if((outhtml != NULL) && (outfile != NULL) && (!appearErrors)) {
    if(!strcmp(outhtml, outfile)) {
      cerr << endl << "ERROR: The output and html files should not be the same." << endl << endl;
      appearErrors = true;
    }
  }

  /* ------------------------------------------------------------------------------------------------------ */

  if((outhtml != NULL) && (!appearErrors)) {
   if((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1) &&
      (gapThreshold == -1) && (conserve == -1) && (simThreshold == -1) && (comThreshold == -1) &&
      (delColumns == NULL) && (resOverlap == -1) && (seqOverlap == -1)) {
      cerr << endl << "ERROR: This parameter can only be used with either a trimming method." << endl << endl;
      appearErrors = true;
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */

  if((outhtml != NULL) && (!appearErrors)) {
   if(((gapThreshold != -1) || (simThreshold != -1)) && (comThreshold != -1)) {
      cerr << endl << "ERROR: Impossible to generate the HTML file using two consecutive trimming methods." << endl << endl;
      appearErrors = true;
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */


  /* ------------------------------------------------------------------------------------------------------ */
  if(((resOverlap != -1) || (seqOverlap != -1)) && (!appearErrors)) {

    if((resOverlap != -1) && (seqOverlap == -1)) {
      cerr << endl << "ERROR: The sequence overlap value should be defined." << endl << endl;
      appearErrors = true;
    }

    else if((resOverlap == -1) && (seqOverlap != -1)) {
      cerr << endl << "ERROR: The residue overlap value should be defined." << endl << endl;
      appearErrors = true;
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */
  if((stats < 0) && (!appearErrors)) {
    stats--;

    if(((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)
      || (gapThreshold != -1) || (conserve != -1) || (simThreshold != -1)) && (outfile == NULL)) {
      cerr << endl << "ERROR: An output file should be defined in order to get the alignment's statistics." << endl << endl;
      appearErrors = true;
    }
  }
  /* ------------------------------------------------------------------------------------------------------ */

  /* ------------------------------------------------------------------------------------------------------ */
  if((comThreshold != -1) && (conserve != -1) && (!appearErrors)) {

    if((gapThreshold != -1) || (simThreshold != -1)) {
      cerr << endl << "ERROR: Combinations among thresholds are not allowed." << endl << endl;
      appearErrors = true;
    }
  }
  /* **** ***** ***** ***** ***** ***** **** End of Parameters Processing **** ***** ***** ***** ***** ***** **** */

  /* **** ***** ***** ***** ***** ***** ***** Files Comparison Methods ***** ***** ***** ***** ***** ***** **** */
  if((compareset != -1) && (!appearErrors)) {

    compAlig = new alignment*[numfiles];
    filesToCompare = new char*[numfiles];

    /* -------------------------------------------------------------------- */
    compare.open(argv[compareset], ifstream::in);

    for(i = 0; (i < numfiles)  && (!appearErrors); i++) {

      /* -------------------------------------------------------------------- */
      for(nline.clear(), compare.read(&c, 1); (c != '\n') && ((!compare.eof())); compare.read(&c, 1))
        nline += c;

      filesToCompare[i] = new char [nline.size() + 1];
      strcpy(filesToCompare[i], nline.c_str());
      /* -------------------------------------------------------------------- */

      /* -------------------------------------------------------------------- */
      compAlig[i] = new alignment;
      if(!compAlig[i] -> loadAlignment(filesToCompare[i])) {
        cerr << endl << "Alignment not loaded: \"" << filesToCompare[i] << "\" Check the file's content." << endl << endl;
        appearErrors = true;
      }

      else {
        if(!compAlig[i] -> isFileAligned()) {
          cerr << endl << "ERROR: The sequences in the input alignment should be aligned in order to use this method." << endl << endl;
          appearErrors = true;
        } else {
          compAlig[i] -> sequenMatrix();

          if(compAlig[i] -> getNumAminos() > maxAminos)
            maxAminos = compAlig[i] -> getNumAminos();

          if((compAlig[i] -> getTypeAlignment() != prevType) && (prevType != -1)) {
            cerr << endl << "ERROR: The alignments' datatypes are different. Check your dataset." << endl << endl;
            appearErrors = true;
          } else
            prevType = compAlig[i] -> getTypeAlignment();
        }
      }
    }
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    if(!appearErrors) {

      compareVect = new float[maxAminos];
      if((stats >= 0) && (outfile != NULL))
        referFile = compareFiles::algorithm(compAlig, filesToCompare, compareVect, numfiles, true);
      else
        referFile = compareFiles::algorithm(compAlig, filesToCompare, compareVect, numfiles, false);

      if(windowSize != -1)
        compareFiles::applyWindow(compAlig[referFile] -> getNumAminos(), windowSize, compareVect);
      else if(conWindow != -1)
        compareFiles::applyWindow(compAlig[referFile] -> getNumAminos(), conWindow, compareVect);

      origAlig -> loadAlignment(filesToCompare[referFile]);

    }
    /* -------------------------------------------------------------------- */

    /* -------------------------------------------------------------------- */
    for(i = 0; i < numfiles; i++)
      compAlig[i] -> destroySequenMatrix();

    for(i = 0; i < numfiles; i++) {
      delete compAlig[i];
      delete[] filesToCompare[i];
    }
    /* -------------------------------------------------------------------- */
  }

  /* **** ***** ***** ***** ***** ***** **** Errors Control **** ***** ***** ***** ***** ***** **** */
  if(appearErrors) {

    delete singleAlig;
    delete origAlig;
    delete[] compAlig;

    delete similMatrix;
    delete []delColumns;

    delete[] filesToCompare;
    delete[] compareVect;

    delete[] outfile;
    delete[] outhtml;

    delete[] infile;
    delete[] matrix;

    return -1;
  }
  /* **** ***** ***** ***** ***** ***** ** End Errors Control ** ***** ***** ***** ***** ***** **** */

  /* -------------------------------------------------------------------- */
  if(conserve == -1)
    conserve  = 0;
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(windowSize != -1) {
    gapWindow = windowSize;
    simWindow = windowSize;
  }
  else {
    if(gapWindow == -1)
      gapWindow = 0;
    if(simWindow == -1)
      simWindow = 0;
  }
  origAlig -> setWindowsSize(gapWindow, simWindow);

  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(outformat != -1)
    origAlig -> setOutputFormat(outformat);
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if((strict) || (strictplus) || (automated1) || (simThreshold != -1.0) || (scc == 1) || (sct == 1)) {
    similMatrix = new similarityMatrix();

    if(matrix != NULL)
      similMatrix -> loadSimMatrix(matrix);

    else {
      if((origAlig -> getTypeAlignment()) == AAType)
        similMatrix -> defaultAASimMatrix();
      else
        similMatrix -> defaultNTSimMatrix();
    }

    if(!origAlig -> setSimilarityMatrix(similMatrix)) {
      cerr << endl << "ERROR: It's imposible to proccess the Similarity Matrix." << endl << endl;
      return -1;
    }
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(sgc) {
    origAlig -> printStatisticsGapsColumns();
    stats++;
    if(stats < -1)
      cout << endl;
  }

  if(sgt) {
    origAlig -> printStatisticsGapsTotal();
    stats++;
    if(stats < -1)
      cout << endl;
  }

  if(scc) {
    origAlig -> printStatisticsConservationColumns();
    stats++;
    if(stats < -1)
      cout << endl;
  }

  if(sct) {
    origAlig -> printStatisticsConservationTotal();
    stats++;
    if(stats < -1)
      cout << endl;
  }

  if(sident) {
    origAlig -> printSeqIdentity();
    stats++;
    if(stats < -1)
      cout << endl;
  }

  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(compareset != -1) {
    if(sfc)
      compareFiles::printStatisticsFileColumns(origAlig -> getNumAminos(), compareVect);
    if(sft)
      compareFiles::printStatisticsFileAcl(origAlig -> getNumAminos(), compareVect);
  }
  /* -------------------------------------------------------------------- */
  /* -------------------------------------------------------------------- */
  if(nogaps)
    singleAlig = origAlig -> cleanGaps(0, 0, complementary);

  else if(noallgaps)
    singleAlig = origAlig -> cleanNoAllGaps(complementary);

  else if(gappyout)
    singleAlig = origAlig -> clean2ndSlope(complementary);

  else if(strict)
    singleAlig = origAlig -> cleanCombMethods(complementary, false);

  else if(strictplus)
    singleAlig = origAlig -> cleanCombMethods(complementary, true);


  else if(automated1) {
    if(origAlig -> calculateSeqIdentity() == GAPPYOUT)
      singleAlig = origAlig -> clean2ndSlope(complementary);
    else
      singleAlig = origAlig -> cleanCombMethods(complementary, false);
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(comThreshold != -1)
    singleAlig = origAlig -> cleanCompareFile(comThreshold, conserve, compareVect, complementary);
 /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if((resOverlap != -1) && (seqOverlap != -1)) {
    singleAlig = origAlig -> cleanSpuriousSeq(resOverlap, (seqOverlap/100));
    singleAlig = singleAlig -> cleanNoAllGaps(complementary);
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(simThreshold != -1.0) {
    if(gapThreshold != -1.0)
      singleAlig = origAlig -> clean(conserve, gapThreshold, simThreshold, complementary);
    else
      singleAlig = origAlig -> cleanConservation(conserve, simThreshold, complementary);
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  else if(gapThreshold != -1.0)
    singleAlig = origAlig -> cleanGaps(conserve, gapThreshold, complementary);
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(delColumns != NULL) {

    if(delColumns[(num - 1)] >= origAlig -> getNumAminos()) {
      cerr << endl << "ERROR: This option only accepts integer numbers between 0 and the alignment's length - 1." << endl << endl;
      appearErrors = true;
    }
    else
      singleAlig = origAlig -> removeColumns(delColumns, num, complementary);
  }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if(singleAlig == NULL) {
    singleAlig = origAlig;
    origAlig = NULL;
  }

  if((outfile != NULL) && (!appearErrors)) {
    if(!singleAlig -> saveAlignment(outfile)) {
      cerr << endl << "ERROR: It's imposible to generate the output file." << endl << endl;
      appearErrors = true;
    }
  }
  else if((stats >= 0) && (!appearErrors))
    singleAlig -> printAlignment();
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if((colnumbering) && (!appearErrors))
    singleAlig -> printCorrespondence();
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  if((outhtml != NULL) && (!appearErrors))
    if(!origAlig -> alignmentSummaryHTML(outhtml, singleAlig -> getNumAminos(), singleAlig -> getNumSpecies(),
                                     singleAlig -> getCorrespResidues(), singleAlig -> getCorrespSequences())) {
      cerr << endl << "ERROR: It's imposible to generate the HTML output file." << endl << endl;
      appearErrors = true;
    }
  /* -------------------------------------------------------------------- */

  /* -------------------------------------------------------------------- */
  delete singleAlig;
  delete[] compAlig;

  delete similMatrix;
  delete []delColumns;

  delete[] filesToCompare;
  delete[] compareVect;

  delete[] outfile;
  delete[] outhtml;

  delete[] infile;
  delete[] matrix;
  /* -------------------------------------------------------------------- */

  return 0;
}

void menu(void) {

  cout << endl;
  cout << "trimAl " << VERSION << "rev" << REVISION << ". Copyright (C) 2009. Salvador Capella-Gutierrez and "
       << "Toni Gabaldón." << endl << endl;

  cout << "trimAl webpage: http://trimal.cgenomics.org" << endl << endl;

  cout << "This program is free software: you can redistribute it and/or modify " << endl
       << "it under the terms of the GNU General Public License as published by " << endl
       << "the Free Software Foundation, the last available version." << endl << endl;

  cout << "Please cite:\tSalvador Capella-Gutierrez, Jose M. Silla-Martinez and" << endl
       << "            \tToni Gabaldon. trimAl: a tool for automated alignment " << endl
       << "            \ttrimming (2009)." << endl << endl;

  cout << "Basic usage" << endl
       << "\ttrimal -in <inputfile> -out <outputfile> -(other options)." << endl << endl;

  cout << "Common options (for a complete list please see the User Guide or visit http://trimal.cgenomics.org):" << endl << endl;
  cout << "    -h                       " << "Print this information and show some examples." << endl;
  cout << "    --version                " << "Print the trimAl version." << endl << endl;

  cout << "    -in <inputfile>          " << "Input file in several formats (clustal, fasta, NBRF/PIR, nexus, phylip3.2, phylip)." << endl;
  cout << "    -compareset <inputfile>  " << "Input list of paths for the files containing the alignments to compare." << endl;
  cout << "    -matrix <inpufile>       " << "Input file for user-defined similarity matrix (default is Blosum62)." << endl << endl;

  cout << "    -out <outputfile>        " << "Output alignment in the same input format (default stdout). (default input format)" << endl;
  cout << "    -htmlout <outputfile>    " << "Get a summary of trimal's work in an HTML file." << endl << endl;

  cout << "    -clustal                 " << "Output file in CLUSTAL format" << endl;
  cout << "    -fasta                   " << "Output file in FASTA format" << endl;
  cout << "    -nbrf                    " << "Output file in NBRF/PIR format" << endl;
  cout << "    -nexus                   " << "Output file in NEXUS format" << endl;
  cout << "    -mega                    " << "Output file in MEGA format" << endl;
  cout << "    -phylip3.2               " << "Output file in PHYLIP3.2 format" << endl;
  cout << "    -phylip                  " << "Output file in PHYLIP/PHYLIP4 format" << endl << endl;

  cout << "    -complementary           " << "Get the complementary alignment." << endl;
  cout << "    -colnumbering            " << "Get the relationship between the columns in the old and new alignment." << endl << endl;

  cout << "    -select { n,l,m-k }      " << "Selection of columns to be removed from the alignment. (see User Guide)." << endl;
  cout << "    -gt -gapthreshold <n>    " << "1 - (fraction of sequences with a gap allowed)." << endl;
  cout << "    -st -simthreshold <n>    " << "Minimum average similarity allowed." << endl;
  cout << "    -ct -conthreshold <n>    " << "Minimum consistency value allowed." << endl;
  cout << "    -cons <n>                " << "Minimum percentage of the positions in the original alignment to conserve." << endl << endl;

  cout << "    -nogaps                  " << "Remove all positions with gaps in the alignment." << endl;
  cout << "    -noallgaps               " << "Remove columns composed only by gaps." << endl << endl;

  cout << "    -gappyout                " << "Use automated selection on \"gappyout\" mode. This method only uses "
                                          << "information based on gaps' distribution. (see User Guide)." << endl;
  cout << "    -strict                  " << "Use automated selection on \"strict\" mode. (see User Guide)." << endl;
  cout << "    -strictplus              " << "Use automated selection on \"strictplus\" mode. (see User Guide)."  << endl;
  cout << "                             " << "(Optimized for Neighbour Joining phylogenetic tree reconstruction)."<< endl << endl;

  cout << "    -automated1              " << "Use a heuristic selection of the automatic method based on similarity statistics. "
                                          << "(see User Guide)." << endl;
  cout << "                             " << "(Optimized for Maximum Likelihood phylogenetic tree reconstruction)."<< endl << endl;

  cout << "    -resoverlap              " << "Minimum overlap of a positions with other positions in the column to be considered a "
                                          << "\"good position\". (see User Guide)." << endl;
  cout << "    -seqoverlap              " << "Minimum percentage of \"good positions\" that a sequence must have in order to be conserved. "
                                          << "(see User Guide)." << endl << endl;

  cout << "    -w <n>                   " << "(half) Window size, score of position i is the average of the window (i - n) to (i + n)."
                                          << endl;
  cout << "    -gw <n>                  " << "(half) Window size only applies to statistics/methods based on Gaps." << endl;
  cout << "    -sw <n>                  " << "(half) Window size only applies to statistics/methods based on Similarity." << endl;
  cout << "    -cw <n>                  " << "(half) Window size only applies to statistics/methods based on Consistency." << endl << endl;

  cout << "    -sgc                     " << "Print gap percentage count for columns in the input alignment." << endl;
  cout << "    -sgt                     " << "Print accumulated gap percentage count." << endl;
  cout << "    -scc                     " << "Print conservation values for columns in the input alignment." << endl;
  cout << "    -sct                     " << "Print accumulated conservation values count." << endl;  
  cout << "    -sfc                     " << "Print compare values for columns in the selected alignment from compare files method."
                                          << endl;
  cout << "    -sft                     " << "Print accumulated compare values count for the selected alignment from compare files method."
                                          << endl;
  cout << "    -sident                  " << "Print identity statistics for all sequences in the alignemnt. (see User Guide)." 
                                          << endl << endl;

}

void examples(void) {

  cout << "Some Examples:" << endl << endl;

  cout << "1) Removes all positions in the alignment with gaps in 10% or more of" << endl
       << "   the sequences, unless this leaves less than 60%. In such case, print" << endl
       << "   the 60% best (with less gaps) positions." << endl << endl;

  cout << "   trimal -in <inputfile> -out <outputfile> -gt 0.9 -cons 60" << endl << endl;

  cout << "2) As above but, the gap percentage is averaged over a window starting" << endl
       << "   3 positions before and ending 3 positions after the column." << endl << endl;

  cout << "   trimal -in <inputfile> -out <outputfile> -gt 0.9 -cons 60 -w 3" << endl << endl;

  cout << "3) Uses an automatic method to decide optimal thresholds, based in the gap percentage" << endl
       << "   count over the whole alignment. (see User Guide for details)." << endl << endl;

  cout << "   trimal -in <inputfile> -out <outputfile> -gappyout" << endl << endl;

  cout << "4) Uses automatics methods to decide optimal thresholds, based on the combination " << endl
       << "   of strict method and similarity values. (see User Guide for details)." << endl << endl;

  cout << "   trimal -in <inputfile> -out <outputfile> -strictplus" << endl << endl;

  cout << "5) Uses an heuristic to decide the optimal method to trimming the alignment. " << endl
       << "   (see User Guide for details)." << endl << endl;

  cout << "   trimal -in <inputfile> -out <outputfile> -automated1" << endl << endl;

  cout << "6) Uses residue and sequences overlap thresholds to delete some sequences from the " << endl
       << "   alignemnt. (see User Guide for details)." << endl << endl;

  cout << "   trimal -in <inputfile> -out <outputfile> -resoverlap 0.8 -seqoverlap 75" << endl << endl;

  cout << "7) Selection of columns to be deleted from the alignment. The selection can " << endl 
       << "   be a column's number or a column's number interval." << endl << endl;

  cout << "   trimal -in <inputfile> -out <outputfile> -select { 2,3,10,45-60,68,70-78 }" << endl << endl;

  cout << "8) Get the complementary alignment from the alignment previously trimmed." << endl << endl; 

  cout << "   trimal -in <inputfile> -out <outputfile> -select { 2,3,45-60 } -complementary" << endl << endl;
}

