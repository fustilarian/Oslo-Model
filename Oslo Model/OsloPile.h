//
//  OsloPile.h
//  Oslo Model
//
//  Created by Henry O'Hagan on 16/03/2014.
//  Copyright (c) 2014 Henry O'Hagan. All rights reserved.
//
/*------------------------------------------------------NOTES------------------------------------------------------*/
/*
A class that constructs a ricepile and evolves it according to the oslo model, storing relevant variables as grains
 are added.
*/
#ifndef Oslo_Model_OsloPile_h
#define Oslo_Model_OsloPile_h
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <gsl/gsl_rng.h>
gsl_rng * r;
const gsl_rng_type * T;

class oslopile
{
private:
    std::vector<std::vector<int>> state; //vector containing all the local slopes and thresholds
    double p; //probability of z_th=1
    int height; //height of pile
    int time; //time, number of grains added
    int dropsize; //number of grains dropping out during previous relaxation
    int avsize; //avalanch size during previous relaxation
    
    void relax(int i)
    {
        rethreshold(i);
        unsigned long l=state.size();
        if(l==1)
        {
            state[0][0]--;
            height--;
        }else
        {
            if(i==0)
            {
                state[0][0]-=2;
                state[1][0]++;
                height--;
            }else if (i==l-1)
            {
                state[l-1][0]--;
                state[l-2][0]++;
            }else
            {
                state[i][0]-=2;
                state[i+1][0]++;
                state[i-1][0]++;
            }
        }

    }//relax site i, reassign the threshold slope of i.
    
    void rethreshold(unsigned long i)
    {
        if(gsl_rng_uniform(r)<p)
        {state[i][1]=1;}
        else
        {state[i][1]=2;}
    }//With probability p assign the threshold slope to 1 and 1-p to 2.
    
    bool over(unsigned long i)
    {
        return state[i][0]>state[i][1];
    }//returns true if site i is unstable
    
public:
/*------------------------------------------------CONSTRUCTORS-----------------------------------------------------*/

    oslopile()
    {}
    
    oslopile(unsigned long l, double prob)
    {
        p=prob;
        height=0;
        time=0;
        dropsize=0;
        avsize=0;
        std::vector<int> vec(2);
        state.assign(l,vec);
        for(unsigned long i=0 ; i<l ; i++)
        {
            state[i][0]=0;
            rethreshold(i);
        }
    }
/*Initialise system with local slopes set to zero. With probability p assign the threshold slope to 1 and 1-p to
2.
*/
    
/*---------------------------------------------------ACCESSORS-----------------------------------------------------*/
    int get_localslope(unsigned long i)
    {
        return state[i][0];
    }
    //returns local slope
    int get_slopethreshold(unsigned long i)
    {
        return state[i][1];
    }
    //returns threshold slope
    int get_height()
    {
        return height;
    }
    //returns height
    int get_time()
    {
        return time;
    }//returns the time (total number of grains added)
    int get_size()
    {
        return int(state.size());
    }//returns system size
    double get_thprobability()
    {
        return p;
    }//returns probability of reassigning a probability to z_th=1
    int get_avalanche_size()
    {
        return avsize;
    }//returns the size of the previous avalanche
    int get_dropsize()
    {
        return dropsize;
    }//returns the dropsize associated with the previous grain added

    
/*---------------------------------------------------MODIFIERS-----------------------------------------------------*/
    void addgrain()
    {
        time++;
        height++;
        state[0][0]++;
        int i=0, j, l=0; //i number of relaxations, j checks whether i has changed
        do
        {
            j=i;
            for(int k=0; k<state.size() ; k++)
            {
                if(over(k))
                {
                    relax(k);
                    i++;
                    if(k==state.size()-1)
                    {
                        l++;//l counts the drop size
                    }
                }
            }//for all sites, relax if site is over threshold
        }while(j!=i);//while at least one site is over z_th, check that all subsequent sites aren't over z_th
        std::vector<int> ivec(2);
        avsize=i;
        dropsize=l;
    }
    //add a grain to site 0, allow system to relax, return avalanch and drop size
    
    void evolve_till_recurrent()
    {
        do
        {
            addgrain();
        }while(double(height)<double(state.size())*1.728);
    }//adds grains to the ricepile until it reaches the set of recurrent configurations
    /*
     this works for p=0.5, the relation h<L*1.728 will work for sizes up to 1024 (because of the 3SF acurracy, the configs
     visited could become skewed if the state that the pile ends up in has very low probability to occur)
     */

};






#endif


