#include <iostream>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include "event.h"
#include "e_list.h"
#include "random_tsaimh.h"

#define ALPHA 0.7	//The probability that another on-off period arrives
#define lambda 2
#define mu 6
#define sita 5
#define Tw_Shape 2
#define Tw_Rate 1
#define Tl_Shape 2
#define Tl_Rate 1

#define SP_EXPIRED 0	//Silent period expired
#define OF_ARRIVAL 1	//On-off period arrives
#define TD_EXPIRED 2	//Timer Td expired
#define LW_ARRIVAL 3	//UE handover from VoLTE to VoWiFi
#define WT_ARRIVAL 4	//UE moves put of VoWiFi and timer Td is triggered
#define WW_ARRIVAL 5	//UE handover from a VoWiFi cell to another VoWiFi cell during Td, which means Td is not expired.

using namespace std;

E_List List;

Event* Generate_arrival_event(int _type, double _time){
      Event* arrival_event;
      arrival_event = new Event; 
      arrival_event->setEventType(_type);
      arrival_event->setTimeStamp(_time);
      return arrival_event;
}

int main()
{
	int simuTimes;	//The number of simulation times
	int randomSeed = (int) time(NULL);
	int Nt;
	int Crh;	//The counting number of reduced handover with delay timer Td and standard rounds Nt

	double mean_Ta = 1.0/lambda;
	double mean_Ts = 1.0/mu;
	double mean_Td = 1.0/sita;
	double mean_Tw = Tw_Shape/Tw_Rate;
	double mean_Tl = Tl_Shape/Tl_Rate;
	double var_Tw = Tw_Shape/(Tw_Rate*Tw_Rate);
	double var_Tl = Tl_Shape/(Tl_Rate*Tl_Rate);

	Prob p(randomSeed++);
	double alpha;
	Expon Ta(randomSeed++, mean_Ta);
	Expon Ts(randomSeed++, mean_Ts);
	Expon Td(randomSeed++, mean_Td);
	Gamma Tw(randomSeed++, mean_Tw, var_Tw);
	Gamma Tl(randomSeed++, mean_Tl, var_Tl);

	Event* initE;
	Event* e;
	double User_Time = 0;	//The timeline of user's active-silent (on-off) actions
	double Mob_Time = 0;	//The timeline of mobility 

	printf("Please enter the following parameter:\n");
	printf("1. Simulation times: ");
	scanf("%d",&simuTimes);
	printf("\n2. The number of standard rounds - Nt: ");
	scanf("%d",&Nt);		//這邊之後要改放在for迴圈裡面，讓每一次的Nt值為random?還是應該要固定Nt值？

	for(int i=0; i<simuTimes; i++)
	{
		Crh=0;
		User_Time = 0;
		Mob_Time = 0;

		Mob_Time += Tw++;
		initE = Generate_arrival_event(LW_ARRIVAL, Mob_Time);
		List << *initE;
		User_Time += Ts++;
		initE = Generate_arrival_event(SP_EXPIRED, User_Time);
		List << *initE;

		alpha = p++;

		while(alpha <= ALPHA)
		{
			List >> e;
			alpha = p++;	//next alpha value
		}

		while(!List)
		{
			List >> e;
		}
	}

	printf("E[Nh(Nt, Td)] = %f\n", (double)Crh/simuTimes);	//E[Nh(Nt, Td)]=Crh/simuTimes
	return 0;
}
