/*
 *  YASS 1.15
 *  Copyright (C) 2004-2016
 *  the YASS team
 *  Laurent Noe, Gregory Kucherov, Mikhail Roytberg, 
 *  Steven Corroy, Antoine De Monte, Christophe Valmir.
 *
 *  laurent.noe|<A>|univ-lille1.fr
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the CeCILL License as published by
 *  the CEA-CNRS-INRIA; either version 2 of the License, or (at your
 *  option) any later version, and the GNU General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 *
 *  This software contains code derived from the GNU libavl library.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __DISPLAY_H_
#define __DISPLAY_H_

#include "global_var.h"
#include "tuple.h"
#include "threads.h"

void Display_Alignements(MA * first_MA);
void Display_Stats();
void Display_Params();
void Display_Progress(long int pos, long int maxpos, Feature * feature);
void DisplayHistoScore(Feature ** feature /* feature[2] */);
void DisplaySubstitutionMatrix();
#endif
