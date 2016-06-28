//
//  oslo.h
//  Oslo Model
//
//  Created by Henry O'Hagan on 17/03/2014.
//  Copyright (c) 2014 Henry O'Hagan. All rights reserved.
//
/*--------------------------------------------------REFERENCES-----------------------------------------------------*/
/*
 http://www.cplusplus.com/reference/ctime/
 http://www.cplusplus.com/reference/unordered_map/unordered_map/
 http://www.cplusplus.com/reference/map/map/map/
 */

#ifndef Oslo_Model_oslo_h
#define Oslo_Model_oslo_h

#include "OsloPile.h"
#include <ctime>
#include <unordered_map>
#include <map>

void heightvstime(std::ofstream&);
void hvt_mov_av(std::ofstream&,std::ofstream&, int);
void avh_stdv_hp(std::ofstream&,std::ofstream&,std::ofstream&,std::ofstream&);
void avalanch_pdf(std::ofstream&);
void avalanch_pdfb(std::ofstream&);
void avalanch_pdfbl(std::ofstream&);
void avalanche_mom(std::ofstream&);
void dropsize(std::ofstream&);
void modropsize(std::ofstream&);



std::vector<oslopile> construct(int, double, int );
void printtime(time_t);



#endif
