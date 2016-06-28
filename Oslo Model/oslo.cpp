//
//  oslo.cpp
//  Oslo Model
//
//  Created by Henry O'Hagan on 16/03/2014.
//  Copyright (c) 2014 Henry O'Hagan. All rights reserved.
//
/*-----------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------NOTES_&_REFERENCES-----------------------------------------------*/
/*
 http://www.cplusplus.com/reference/ctime/
 http://www.cplusplus.com/reference/unordered_map/unordered_map/
 http://www.cplusplus.com/reference/map/map/map/

 Use of gsl will require custom header/library search paths, 
 use the command gsl-config --cflags --libs-without-cblas in
 the shell in order to find the necessary information.
 GSL REFERENCE MANUAL:
 http://www.gnu.org/software/gsl/manual/gsl-ref.html
*/


#include "oslo.h"
using namespace std;

/*-----------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------1D_OSLO_MODEL--------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------------------------*/
int main(int argc, const char * argv[])
{
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T); //select rng, the default is "Mersenne twister" (mt19937)
    
    cout<<"RNG :"<<'\t'<<gsl_rng_name(r)<<endl; //this is just the RNG algorithm used
    
    /*--uncomment each function to output data for each respective task--*/
    
    //----TASK 2 FILES & FUNCTIONS----
    //ofstream outfile001("hvt.dat");//height versus time file
    //ofstream outfile002("hvtm.dat");//height versus running mean time file
    //ofstream outfile003("hvtmc.dat");//"  "  "  "  " with data collapse (rescaled axes using scaling ansatz)
    //heightvstime(outfile001);
    //hvt_mov_av(outfile002, outfile003, 100);
    //----TASK 3 FILES & FUNCTIONS----
    //ofstream outfile011("mvl.dat");//mean(height) vs system size
    //ofstream outfile012("svl.dat");//stddev(height) vs system size
    //ofstream outfile013("pvh.dat");//probability of a given height vs height for different system sizes
    //ofstream outfile014("pvhc.dat");//"  "  "  "  "  "   "   "   "    "   " with data collapsed
    //avh_stdv_hp(outfile011, outfile012, outfile013, outfile014);
    //----TASK 4 FILES & FUNCTIONS----
    //ofstream outfile021("pvs.dat");//probability of an avalanch vs avalanch size for different iterations
    //ofstream outfile022("pvsb.dat");//" " with data binning
    //avalanch_pdf(outfile021);
    //avalanch_pdfb(outfile022);
    //ofstream outfile023("pvsbl.dat");//probability of an avalanch vs avalanch size w/ different system sizes, binned
    //avalanch_pdfbl(outfile023);
    //ofstream outfile024("movl.dat");//moments vs systemsize
    //avalanche_mom(outfile024);
    //----TASK 5 FILES & FUNCTIONS----
    //ofstream outfile031("pvd.dat");
    //dropsize(outfile031);
    ofstream outfile032("mdvl.dat");
    modropsize(outfile032);
    
    gsl_rng_free (r); //frees all the memory associated with the generator r.
    return 0;
}
/*-----------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------SIMULATIONS-----------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------------------------*/

/*---------------------------------------------------TASK_2A-------------------------------------------------------*/
void heightvstime(ofstream& fout)
{
    vector<oslopile> piles=construct(11,0.5,1);//creates 11 rice piles

    time_t timer, now;
    for(int k=1; k<=100;k++)//counts every 1% of iterations, outputs time for each percent completed
    {
        time(&timer);
        for(int j=0; j<8000;j++)
        {
            for(int i=0; i<11; i++)//adds grain, outputs height, time not needed since iterations=time
            {
                fout<<piles[i].get_height()<<'\t';
                piles[i].addgrain();
            }   fout<<endl;
        }
        time(&now);
        printtime(now);
        cout<<k<<"% done, "<<difftime(now, timer)<<" s/%"<<endl;
    }
}
/*---------------------------------------------------TASK_2B-------------------------------------------------------*/
void hvt_mov_av(ofstream& fout, ofstream& foutc, int avsize)
{
    vector<oslopile> piles=construct(11, 0.5,1);//creates 11 rice piles

    vector<vector<int>> height(11); //all the heights used for each average
    vector<double> m_height(11,0.0);// contains the mean average height
    for(int i=0; i<11 ; i++) height[i].resize(avsize);
    for(int j=0; j<avsize;j++)
    {
        for(int i=0; i<11; i++)
        {
            height[i][j]=piles[i].get_height();
            m_height[i]+=double(height[i][j]);
            piles[i].addgrain();
        }
    }
    
    for(int i=0; i<11; i++) m_height[i]/=double(avsize);
    
    time_t timer, now;
    for(int k=1; k<=100;k++)//counts every 1% of iterations, outputs time for each percent completed
    {
        time(&timer);
        for(int j=0; j<8000;j++)
        {
            fout<<piles[0].get_time()<<'\t';//time is outputted since it doesn't start at t=0
            for(int i=0; i<11; i++)//adds grain, outputs running mean height and collapsed heights
            {
                foutc<<double(piles[i].get_time())/double(piles[i].get_size()*piles[i].get_size())<<'\t';
                foutc<<m_height[i]/double(piles[i].get_size())<<'\t';
                fout<<m_height[i]<<'\t';
                piles[i].addgrain();
                m_height[i]+=double(piles[i].get_height()-height[i][j%avsize])/double(avsize);
                height[i][j%avsize]=piles[i].get_height();//update the heights stored
            }   fout<<endl; foutc<<endl;
        }
        time(&now);
        printtime(now);
        cout<<k<<"% done, "<<difftime(now, timer)<<" s/%"<<endl;
    }
}
/* Rather than recalculating the moving average at each time step, the heights at each time step required is
 remembered and then updated as the ricepiles are updated. This requires a couple of extra calculations i.e.
 (time mod (mean size) ), which is however not as severe as the number of operations from summing up all the
 heights.
 */

/*---------------------------------------------------TASK_3--------------------------------------------------------*/
void avh_stdv_hp(ofstream& foutm, ofstream& fouts, ofstream& foutp, ofstream& foutpc)
{
    time_t timer, now;
    vector<oslopile> piles=construct(10, 0.5,2); //initialise 11 piles with size 1,2 etc..
    for(int i=0; i<10; i++)
    {
        time(&timer);
        printtime(timer);
        cout<<"Creating "<<i<<"th recurrent pile"<<endl;
        piles[i].evolve_till_recurrent(); //evolves all piles 'til the recurrent configurations
    }
    vector<double> mean(10,0.0), meansquare(10,0.0); //variables for mean and mean square
    vector<vector<int>> heightfreq(10); //variables for height frequency
    for(int i=0; i<10; i++) heightfreq[i].assign(piles[i].get_size()+1, 0);
    
    for(int i=1; i<=100 ; i++)
    {
        time(&timer);
        for(int j=0; j<10000; j++)
        {
            for(int k=0; k<10; k++)
            {
                heightfreq[k][piles[k].get_height()-piles[k].get_size()]++;
                mean[k]+=piles[k].get_height();
                meansquare[k]+=piles[k].get_height()*piles[k].get_height();
                piles[k].addgrain();
            }
        }
        time(&now);
        printtime(now);
        cout<<"Cycling through pile configurations... "<<i<<"% done, "<<difftime(now, timer)<<" s/%"<<endl;
    }
    
    cout<<"Writing Results..."<<endl;
    vector<double> stddev(10);
    for(int i=0;i<10;i++)//write results
    {
        mean[i]/=1000000.0;
        meansquare[i]/=1000000.0;
        stddev[i]=sqrt(meansquare[i]-mean[i]*mean[i]);
        foutm<<mean[i]<<'\t'<<piles[i].get_size()<<endl;
        fouts<<stddev[i]<<'\t'<<piles[i].get_size()<<endl;
    }
    double p_hl, h_t;
    for(int i=0; i<heightfreq[9].size(); i++)
    {
        for(int j=0; j<10; j++)
        {
            if(i<heightfreq[9-j].size())
            {
                p_hl=double(heightfreq[9-j][i])/1000000.0;
                h_t=double(i+piles[9-j].get_size());
                foutp<<p_hl<<'\t'<<h_t<<'\t';
                foutpc<<pow(double(piles[9-j].get_size()),(0.245))*p_hl<<'\t'//*stddev[9-j]
                <<pow(double(piles[9-j].get_size()),(-0.245))*(h_t-1.7287*piles[9-j].get_size())<<'\t';//(h-mean)/stddev
            }
        }foutp<<endl; foutpc<<endl;
    }
    
}
/*The dynamics of L=1 is different, so it is excluded. Lower thresholds are more likely to be reassigned a
 new threshold as grains are added. In the L=1 case, the threshold assignment is completely random.
 -h_t~L and <h>~L so their difference doesn't need to be compensated for.
 */

/*---------------------------------------------------TASK_4A-------------------------------------------------------*/
/*----------AVALANCHE SIZE PROBABILITY MASS FUNCTION----------*/
/*
 For different numbers of iterations
 */
void avalanch_pdf(ofstream& fout)
{
    time_t timer, now;
    int orders=5;
    oslopile pile(128,0.5);
    pile.evolve_till_recurrent();
    vector<unordered_map<int, int>> sizefreq(orders);
    int n=100;
    for(int k=0; k<orders ; ++k)
    {
        for(int i=1; i<=100 ; i++)
        {

            time(&timer);
            for(int j=0; j<n; j++)
            {
                unordered_map<int, int>::const_iterator s=sizefreq[k].find(pile.get_avalanche_size());
                if(s==sizefreq[k].end())
                {
                    sizefreq[k].insert(make_pair<int, int>(pile.get_avalanche_size(),1) );
                }else
                {
                    sizefreq[k][pile.get_avalanche_size()]++;
                }
                pile.addgrain();
        }
        time(&now);
        printtime(now);
        cout<<"Cycling through "<<k+1<<"th pile configurations... "<<i<<"% done, "<<difftime(now, timer)<<" s/%"<<endl;
    }
        n*=10;
    }
    cout<<"writing results..."<<endl;
    vector<unordered_map<int, int>::iterator> s(orders);
    for(int i=0; i<orders ; ++i) s[i]=sizefreq[i].begin();
    do
    {
        double n=10000.0;
        for(int i=1; i<orders; i++) n*=10.0;
        for(int i=0; i<orders; i++)
        {
            if(s[orders-1-i]!=sizefreq[orders-1-i].end())
            {
                fout<<double(s[orders-1-i]->second)/n<<'\t'<<s[orders-1-i]->first<<'\t';
                ++s[orders-1-i];
                n/=10.0;
            }
        }fout<<endl;
    }while(s[orders-1]!=sizefreq[orders-1].end());
    
}
/*
Since it's difficult to order the sizes, an unordered map is used, which also removes redundant values where freq=0.
 the map uses binary search on the keys.
 Upper bound on avalanch size: sum_I (L-I)^2, for z_th 2 -> 1 for all sites.
 
*/

/*----------AVALANCHE SIZE PROBABILITY WITH BINNING----------*/
/*
 for different numbers of iterations
 */
void avalanch_pdfb(ofstream& fout)
{
    double base=1.2;
    int orders=5;
    time_t timer, now;
    oslopile pile(128,0.5);
    pile.evolve_till_recurrent();
    vector<map<int, int>> sizefreq(orders);
    int min=100, bin_number;
    for(int k=0; k<orders ; ++k)
    {
        for(int i=1; i<=100 ; i++)
        {
            time(&timer);
            for(int j=0; j<min; j++)
            {
                if(pile.get_avalanche_size()!=0)
                {
                    bin_number=int(log(double(pile.get_avalanche_size()))/log(base));
                }else
                {
                    bin_number=0;
                }
                map<int, int>::const_iterator s=sizefreq[k].find(bin_number);
                if(s==sizefreq[k].end())
                {
                    sizefreq[k].insert(make_pair(bin_number,1) );
                }else
                {
                    sizefreq[k][bin_number]++;
                }
                pile.addgrain();
            }
            time(&now);
            printtime(now);
            cout<<"Cycling through "<<k+1<<"th pile configurations... "<<i<<"% done, "<<difftime(now, timer)<<" s/%"<<endl;
        }
        min*=10;
    }
    cout<<"writing results..."<<endl;
    vector<map<int, int>::iterator> s(orders);
    double s_min=1.0, s_max=base, s_j, p_sj;
    for(int i=0; i<orders ; ++i) s[i]=sizefreq[i].begin();
    do
    {
        double n=10000.0;
        for(int i=1; i<orders; i++) n*=10.0;
        for(int i=0; i<orders; i++)
        {
            if(s[orders-i-1]!=sizefreq[orders-i-1].end())
            {
                for(int j=0; j<s[orders-i-1]->first; j++)
                {
                    s_max*=base;
                    s_min*=base;
                }
                s_min=floor(s_min);
                s_max=floor(s_max);
                if(s_min!=s_max)
                {
                    s_min++;//so that s_max(j) is not equal to s_min(j+1)
                }
                s_j=sqrt(s_min*s_max);
                p_sj=double(s[orders-i-1]->second)/(n*(s_max-s_min+1.0));
                fout<<s_j<<'\t'<<p_sj<<'\t';
                ++s[orders-i-1];
                n/=10.0;
                s_min=1.0;
                s_max=base;
            }
        }fout<<endl;
    }while(s[orders-1]!=sizefreq[orders-1].end());
    
}
/*
Map has slower access than unordered map, but the overhead insetion is justified for when the data is to be plotted
 vector<map<int,int>::iterator> s(4) is used because we need to iterate over the frequency map primary keys for each
 different pile i, in parallel, so that the file is outputted row by row.
 */
/*---------------------------------------------------TASK_4B-------------------------------------------------------*/
/*
 for different system sizes
 */
void avalanch_pdfbl(ofstream& fout)
{
    double base=1.2;
    time_t timer, now;
    vector<oslopile> piles=construct(4, 0.5, 128);
    for(int i=0; i<piles.size(); i++) piles[i].evolve_till_recurrent();
    vector<map<int,int>> sizefrequency(4);
    int bin_number;
    for(int i=0; i<4 ; i++)
    {
        for(int j=1; j<=100 ; ++j)
        {
            time(&timer);
            for(int k=0 ; k<1000000 ;++k)
            {
                if(piles[i].get_avalanche_size()!=0)
                {
                    bin_number=int(log(double(piles[i].get_avalanche_size()))/log(base));
                }else
                {
                    bin_number=0;
                }
                map<int,int>::const_iterator s=sizefrequency[i].find(bin_number);
                if(s==sizefrequency[i].end())
                {
                    sizefrequency[i].insert(make_pair(bin_number, 1));
                }else
                {
                    sizefrequency[i][bin_number]++;
                }
                piles[i].addgrain();
            }
            time(&now);
            printtime(now);
            cout<<"Cycling through "<<i<<"th pile configurations... "<<j<<"% done, "<<difftime(now, timer)<<" s/%"<<endl;
        }
    }
    cout<<"Writing results..."<<endl;
    vector<map<int,int>::iterator> s(4);
    map<unsigned long,int> orderedsize;
    for(int i=0; i<4; ++i)
    {
        s[i]=sizefrequency[i].begin();
        orderedsize.insert(make_pair(sizefrequency[i].size(), i));
    }
    double s_min=1.0, s_max=base, s_j, p_sj;
    vector<bool> break_condition(5,false);
    do
    {
        break_condition[4]=true;
        for(map<unsigned long, int>::reverse_iterator m=orderedsize.rbegin(); m!=orderedsize.rend(); m++)
        {
            if(s[m->second]!=sizefrequency[m->second].end())
            {
                for(int j=0; j<s[m->second]->first; j++)
                {
                    s_max*=base;
                    s_min*=base;
                }
                s_min=floor(s_min);
                s_max=floor(s_max);
                if(s_min!=s_max)
                {
                    s_min++;//so that s_max(j) is not equal to s_min(j+1)
                }
                s_j=sqrt(s_min*s_max);
                p_sj=double(s[m->second]->second)/(100000000.0*(s_max-s_min+1.0));
                fout<<s_j<<'\t'<<p_sj<<'\t';
                s_min=1.0;
                s_max=base;
                ++s[m->second];
            }else
            {
                break_condition[m->second]=true;
            }
            break_condition[4]=break_condition[4]&&break_condition[m->second];
        }fout<<endl;
    }while(!break_condition[4]);
}
/*
 map<unsigned long, int> orderedsize is used to sort the vector references i, in order of the size of the vectors
 within i. the primary key is their size and the associated value is the reference i. then iterate backwards over
 orderedsize so that the longest columns are outputted first, so that the columns are not mixed up.
 */
/*---------------------------------------------------TASK_4D-------------------------------------------------------*/
/*---------------------------------------------------MOMENTS-------------------------------------------------------*/
void avalanche_mom(ofstream& fout)
{
    time_t timer, now;
    vector<oslopile> piles=construct(7, 0.5, 16);
    for(int i=0; i<piles.size(); i++) piles[i].evolve_till_recurrent();
    vector<vector<vector<double>>> moment( 7 , vector<vector<double>>
                                         ( 6 , vector<double>
                                         (100, 0.0)          )        );
    double intermediate;
    for(int i=0; i<7; i++)//iterates over system sizes i=0-> L=128, i=1-> L=256 etc
    {
        for(int j=0; j<100; j++)// iterates over mean calculations
        {
            time(&timer);
            for(int k=0; k<10000; k++)//iterates over pile configurations
            {
                for(int l=0; l<6; l++)//iterates over moment number
                {
                    intermediate=double(piles[i].get_avalanche_size());
                    for(int m=0; m<l; m++)
                    {
                        intermediate*=double(piles[i].get_avalanche_size());
                    }
                    moment[i][l][j]+=intermediate;
                }
                piles[i].addgrain();
            }
            for(int l=0 ; l<6 ; l++)
            {
                moment[i][l][j]/=10000.0;
            }
            time(&now);
            printtime(now);
            cout<<"Cycling through "<<i<<"th pile configurations... "
            <<j+1<<"% done, "<<difftime(now, timer)<<" s/%"<<endl;
        }
    }
    cout<<"writing results..."<<endl;
    double meanofmean=0.0, meansquareofmean=0.0;
    for(int l=0 ; l<7 ; l++)//different system sizes
    {
        fout<<piles[l].get_size()<<'\t';
        for(int k=0; k<6 ; k++)//different moments
        {
            for(int j=0;j<100; j++)//means
            {
                meanofmean+=moment[l][k][j];
                meansquareofmean+=moment[l][k][j]*moment[l][k][j];
            }
            meanofmean/=100.0;
            meansquareofmean/=100.0;
            fout<<'\t'<<meanofmean<<'\t'<<sqrt((meansquareofmean-meanofmean*meanofmean)/100.0)<<'\t';
            meanofmean=0.0;
            meansquareofmean=0.0;
        }fout<<endl;
    }
}
/*---------------------------------------------------TASK_5--------------------------------------------------------*/
/*-------------------------------------------------DROP_SIZE-------------------------------------------------------*/

void dropsize(ofstream& fout)
{
    double base=1.1;
    time_t timer, now;
    vector<oslopile> piles=construct(4, 0.5, 32);
    for(int i=0; i<piles.size(); i++) piles[i].evolve_till_recurrent();
    vector<map<int,int>> sizefrequency(4);
    int bin_number, ite=400000;
    vector<int> total(4,0);
    for(int k=0; k<4; k++)
    {
    for(int i=1; i<=100 ; i++)
    {
        
        time(&timer);
        for(int j=0; j<ite; j++)
        {
            if(piles[k].get_dropsize()!=0)
            {
                total[k]++;
                bin_number=int(log(double(piles[k].get_dropsize()))/log(base));
                map<int, int>::const_iterator s=sizefrequency[k].find(bin_number);
                if(s==sizefrequency[k].end())
                {
                    sizefrequency[k].insert(make_pair(bin_number,1) );
                }else
                {
                    sizefrequency[k][bin_number]++;
                }
            }
            piles[k].addgrain();
        }
        time(&now);
        printtime(now);
        cout<<"Cycling through "<<k+1<<"th pile configurations... "
        <<i<<"% done, "<<difftime(now, timer)<<" s/%"<<endl;
    }
    }cout<<"writing results..."<<endl;
    vector<map<int, int>::iterator> s(4);
    map<double,int> orderedsize;
    vector<int> deg(4);
    for(int i=0; i<4 ; ++i)
    {
        s[i]=sizefrequency[i].begin();
        if(orderedsize.find(double(sizefrequency[i].size()))!=orderedsize.end())
        {
            orderedsize.insert(make_pair(double(sizefrequency[i].size()),i) );
        }else
        {
            orderedsize.insert(make_pair(double(sizefrequency[i].size())+0.5,i) );
        }
    }
    vector<bool> break_condition(5,false);
    double d_min=1.0, d_max=base, d_j, p_dj;
    do
    {
        break_condition[4]=true;
        for(map<double, int>::reverse_iterator m=orderedsize.rbegin(); m!=orderedsize.rend(); m++)
        {
            if(s[m->second]!=sizefrequency[m->second].end())
            {
                
                for(int j=0; j<s[m->second]->first; j++)
                {
                    d_max*=base;
                    d_min*=base;
                }
                d_min=floor(d_min);
                d_max=floor(d_max);
                if(d_min!=d_max)
                {
                    d_min++;//so that s_max(j) is not equal to s_min(j+1)
                }
                d_j=sqrt(d_min*d_max);
                p_dj=double(s[m->second]->second)/(double(total[m->second])*(d_max-d_min+1.0));
                fout<<d_j<<'\t'<<p_dj<<'\t';
                d_min=1.0;
                d_max=base;
                ++s[m->second];
            }else
            {
                break_condition[m->second]=true;
            }
            break_condition[4]=break_condition[m->second]&&break_condition[4];
        }fout<<endl;
    }while(!break_condition[4]);

}
/*----Moments OF DROPSIZE----*/
void modropsize(ofstream& fout)
{
    double inter;
    vector<oslopile> piles=construct(7, 0.5, 4);
    for(int i=0; i<7 ; i++) piles[i].evolve_till_recurrent();
    vector<vector<double>> mom(8,vector<double>(6,0.0));
    vector<int> total(7,0);
    time_t timer, now;
    for(int k=0; k<7; k++)
    {
            for(int i=1; i<=100 ; i++)
            {
            
                time(&timer);
                for(int j=0; j<100000; j++)
                {
                    if(piles[k].get_dropsize()!=0)
                    {
                        inter=piles[k].get_dropsize();
                        for(int l=0; l<6; l++)
                        {
                            mom[k][l]+=inter;
                            inter*=piles[k].get_dropsize();
                        }
                        total[k]++;
                    }
                    piles[k].addgrain();
                }
                time(&now);
                printtime(now);
                cout<<"Cycling through "<<k+1<<"th pile configurations... "
                <<i<<"% done, "<<difftime(now, timer)<<" s/%"<<endl;
            }
    }
    for(int k=0; k<7; k++)
    {
        fout<<piles[k].get_size()<<'\t';
        for(int l=0; l<6 ; l++)
        {
            fout<<mom[k][l]/double(total[k])<<'\t';
        }
        fout<<endl;
    }
}

/*-----------------------------------------------OTHER_GUBBINS-----------------------------------------------------*/
void printtime(time_t now)
{
    tm *ltm=localtime(&now);
    cout<<'['<< 1 + ltm->tm_hour << ":"<< 1 + ltm->tm_min <<":"<< 1 + ltm->tm_sec<<']'<<'\t';
}//prints time given by "time_t now"


vector<oslopile> construct(int logmaxsize, double p, int s)
{
    vector<oslopile> piles(logmaxsize);
    int size=s;
    for(int i=0; i<logmaxsize ; i++)
    {
        piles[i]=oslopile(size,p);
        size*=2;
    }
    return piles;
}
//creates a vector of ricepiles each with double the size of the previous up to logmaxsize
/*
void dropsize(ofstream& fout)
{
    time_t timer, now;
    vector<oslopile> piles=construct(4, 0.5, 128);
    for(int i=0; i<piles.size(); i++) piles[i].evolve_till_recurrent();
    vector<map<int,int>> sizefrequency(4);
    int ite=100000;
    vector<int> total(4,0);
    for(int k=0; k<4; k++)
    {
        for(int i=1; i<=100 ; i++)
        {
            
            time(&timer);
            for(int j=0; j<ite; j++)
            {
                if(piles[k].get_dropsize()!=0)
                {
                    map<int, int>::const_iterator s=sizefrequency[k].find(piles[k].get_dropsize());
                    if(s==sizefrequency[k].end())
                    {
                        sizefrequency[k].insert(make_pair(piles[k].get_dropsize(),1) );
                        total[k]++;
                    }else
                    {
                        sizefrequency[k][piles[k].get_dropsize()]++;
                        total[k]++;
                    }
                }
                piles[k].addgrain();
            }
            time(&now);
            printtime(now);
            cout<<"Cycling through "<<k+1<<"th pile configurations... "
            <<i<<"% done, "<<difftime(now, timer)<<" s/%"<<endl;
        }
    }cout<<"writing results..."<<endl;
    vector<map<int, int>::iterator> s(4);
    map<unsigned long,int> orderedsize;//sorts in order of size
    for(int i=0; i<4 ; ++i)
    {
        s[i]=sizefrequency[i].begin();
        orderedsize.insert(make_pair(sizefrequency[i].size(), i));
    }
    vector<bool> break_condition(5,false);
    do
    {
        break_condition[4]=true;
        for(map<unsigned long, int>::reverse_iterator m=orderedsize.rbegin(); m!=orderedsize.rend(); m++)
        {
            if(s[m->second]!=sizefrequency[m->second].end())
            {
                fout<<s[m->second]->first<<'\t'<<double(s[m->second]->second)/double(total[m->second])<<'\t';
                ++s[m->second];
            }else
            {
                break_condition[m->second]=true;
            }
            break_condition[4]=break_condition[m->second]&&break_condition[4];
        }fout<<endl;
    }while(!break_condition[4]);
    
}
*/

