//We use DTA when it's VoWiFi disconnected and is in silent period
//Silent period seems to be 40%-60% during one call session

//Td is fixed
#include <iostream>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <math.h>
#include "event.h"
#include "e_list.h"
#include "random_tsaimh.h"
#include <fstream>

/* Parameter definition */
#define ALPHA 0.8	//The probability that another on-off period arrives
#define lambda 0.05	//Ta is an exponential distribution with mean E[Ta] = 1/lambda
#define mu 0.05		//Ts is an exponential distribution with mean E[Ts] = 1/mu
#define Da 0.02		//In active periods, every Da milli-seconds comes a packet.
#define Ds 0.16		//In silent periods, every Ds milli-seconds comes a packet.

/* Event definition */
#define SP_EXPIRED 0	//Silent period expired
#define AP_EXPIRED 1	//Active period expired
#define TD_EXPIRED 2	//Timer Td expired
#define LW_ARRIVAL 3	//UE handover from VoLTE to VoWiFi
#define WT_ARRIVAL 4	//UE moves out of VoWiFi and timer Td is triggered
#define WW_ARRIVAL 5	//UE handover from a VoWiFi cell to another VoWiFi cell during Td, which means Tl expired before Td.

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
	double mean_Td;	//Td is an exponential distribution with mean E[Td] = 1/sita
	int randomSeed = (int) time(NULL);
	int Crh;		//The number of handover with DTA during whole simulation times;
	int Crh_0;		//The number of handover with standard procedure(Td = 0) during whole simulation times
	int sum_Crh;
	int sum_Crh_0;
	int Np;			//The number of packets during whole simulation times
	int Npd;		//The number of packet drops during whole simulation times
	int Npda;		//The number of packet drops during active periods in one td delay period
	int Npds;		//The number of packet drops during silent periods in one td delay period
	int total_Npda;	//The total packet drops during active periods
	int total_Npds;	//The total packet drops during silent periods
	double sum_PLR;	//The sum of packet loss rate through all simulations
	double sum_PLR_a;		//The sum of packet loss rate in active periods through all simulations.
	double sum_PLR_s;		//The sum of packet loss rate in silent periods through all simulations.
	bool inActivePeriod;	//The boolean value to see whether it's in active period or not.
	bool inTdPeriod;		//The boolean value to see whether it's in td period or not.
	bool isFirstE;			//The boolean value to see whether it's the first event after td timer triggers.
	bool on;

	/* Input of simulation times and E[td] */
	simuTimes = 1000000;
	double list_Td[]={1, 2, 5, 10, 15, 20, 25};

	/* Input finished */

	// double sita = 1.0/mean_Td;
	double mean_Ta = 1.0/lambda;
	double mean_Ts = 1.0/mu;
	double mean_Tw = 10;
	double mean_Tl = 1;
	double var_Tw = mean_Tw*mean_Tw;
	double var_Tl = mean_Tl*mean_Tl;

	double alpha;			//If alpha <= ALPHA, there is another on-off period following up
	Prob p(randomSeed++);	//The generator of alpha

	Expon Ta(randomSeed++, mean_Ta);			//Ta is the interval time during active period
	Expon Ts(randomSeed++, mean_Ts);			//Ts is the interval time during silent period
	Gamma Tw(randomSeed++, mean_Tw, var_Tw);	//Tw is the interval time during VoWiFi is connected
	Gamma Tl(randomSeed++, mean_Tl, var_Tl);	//Tl is the interval time during VoWiFi is disconnected

	Event* initE;
	Event* e;
	double preU_Time;	//The time of previous User event happens.
	double startTd_Time;//The time of Td timer starts.
	double User_Time;	//The timeline of user's active-silent (on-off) actions
	double Mob_Time;	//The timeline of mobility 
	double Td_Time;		//The time when td expired
	double sum_Rh;

	/* Mathematical definitions */
	// double E_Nh_0;			//E[Nh(0)]
	// double E_Nrh_td;		//E[Nrh(Td)]
	// double E_Nh_td;			//E[Nh(Td)]
	// double Rh_td;
	// double Laplace_fsL;		//The Laplace tranform function of fsL, which is the p.d.f of TsL
	// double eta_Tl;			//mean_Tl = 1/eta_Tl
	// double eta_Tw;			//mean_Tw = 1/eta_Tw
	// double base;

	/* Mathematical Analysis */
	// eta_Tw = (double)Tw_Rate/Tw_Shape;
	// eta_Tl = (double)Tl_Rate/Tl_Shape;
	// base = (Tl_Shape*eta_Tl)/(Tl_Shape*eta_Tl+sita);
	// Laplace_fsL = pow(base, Tl_Shape);
	// E_Nh_0 = 2*(lambda+mu)*eta_Tl*eta_Tw/(lambda*mu*(eta_Tl+eta_Tw)*(1-ALPHA));
	// E_Nrh_td = Laplace_fsL*eta_Tl*eta_Tw/(mu*(eta_Tl+eta_Tw)*(1-ALPHA));
	// E_Nh_td = E_Nh_0 - E_Nrh_td;
	// Rh_td = Laplace_fsL*lambda/(2*(lambda+mu));
	
	//Output file: The results will be written in the file
	ofstream fp;
	fp.open("td_WiFidisconnection_result_twd1.out");
	if(!fp)
	{
		printf("Fail to open the file!\n");
	}


	// Experiments for different E[Td]
	for(int a=0; a<7;a++)
	{
		mean_Td = list_Td[a];

		Td_Time = 0;

		sum_Crh = 0;
		sum_Crh_0 = 0;

		sum_Rh = 0;

		sum_PLR = 0;
		sum_PLR_a = 0;
		sum_PLR_s = 0;

		/* Simulation starts */
		for(int i=0; i<simuTimes; i++)
		{
			printf("\n** %d simulation **\n",i+1);

			/* Reset starts */
			User_Time = 0;
			Mob_Time = 0;
			preU_Time = 0;
			startTd_Time = 0;

			Crh_0 = 0;
			Crh = 0;

			Np = 0;
			Npd = 0;
			total_Npda = 0;
			total_Npds = 0;

			while(!List)
			{
				List >> e;	
			}
			/* Reset ends */

			/* Initialization starts */
			randomSeed++;
			//Randomly simulate that UE is going to be either VoWiFi connected or disconnected at the beginning
			if(randomSeed%2==1)
			{
				Mob_Time += Tl++;
				initE = Generate_arrival_event(LW_ARRIVAL, Mob_Time);	//Going to be VoWiFi connected
			}
			else
			{
				Mob_Time += Tw++;
				initE = Generate_arrival_event(WT_ARRIVAL, Mob_Time);	//Going to be VoWiFi disconnected
			}
			List << *initE;

			//Generate an AP_EXPIRED
			User_Time += Ta++;
			initE = Generate_arrival_event(AP_EXPIRED, User_Time);
			List << *initE;
			inActivePeriod = true;
			Np += floor(User_Time/Da);
			
			inTdPeriod = false;
			isFirstE = false;

			on = true;
			/* Initialization ends */

			/* Events handling starts */
			while(on)
			{
				List >> e;
				switch(e->getEventType())
				{
					case SP_EXPIRED:
						//printf("== SP_EXPIRED ==\n");
						if(inTdPeriod)
						{
							if(isFirstE)
							{
								isFirstE = false;
								Npds += floor((User_Time - startTd_Time)/Ds);
							}
							else
							{
								Npds += floor((User_Time - preU_Time)/Ds);
							}
						}
						alpha = p++;		//next alpha value
						if(alpha <= ALPHA)	//There is another on-off period following up
						{	
							preU_Time = User_Time;
							User_Time += Ta++;
							Np += floor((User_Time - preU_Time)/Da);
							e = Generate_arrival_event(AP_EXPIRED, User_Time);
							List << *e;
							inActivePeriod = true;
						}
						else				//There is no other on-off period following up. Call session ends.
						{	
							on = false;
						}
						break;

					case AP_EXPIRED:
						//printf("== AP_EXPIRED ==\n");
						if(inTdPeriod)
						{
							Npda += floor((User_Time - preU_Time)/Da);
						}
						preU_Time = User_Time;
						User_Time += Ts++;
						Np += floor((User_Time - preU_Time)/Ds);
						e = Generate_arrival_event(SP_EXPIRED, User_Time);
						List << *e;
						inActivePeriod = false;
						break;

					case LW_ARRIVAL:
						//printf("== LW_ARRIVAL ==\n");
						Crh++;
						Crh_0++;

						Mob_Time += Tw++;
						e = Generate_arrival_event(WT_ARRIVAL, Mob_Time);
						List << *e;				
						break;

					case WT_ARRIVAL:
						//printf("== WT_ARRIVAL ==\n");
						Crh_0++;
						if(!inActivePeriod)	//If it's in silent period, we can use DTA to reduce handover
						{
							//Td timer starts
							inTdPeriod = true;
							isFirstE = true;
							Npda = 0;
							Npds = 0;

							// Td_Time = Mob_Time + Td++;
							Td_Time = Mob_Time + mean_Td;
							startTd_Time = Mob_Time;
							Mob_Time += Tl++;
							if(Td_Time > Mob_Time)	//Tl expired first, UE handovers from one VoWiFi to another VoWiFi
							{
								e = Generate_arrival_event(WW_ARRIVAL, Mob_Time);
								List << *e;
							}
							else		//Td expired first, UE handovers from VoWiFi to VoLTE when Td expires
							{
								e = Generate_arrival_event(TD_EXPIRED, Td_Time);
								List << *e;
								e = Generate_arrival_event(LW_ARRIVAL, Mob_Time);
								List << *e;
							}
						}
						else			//If it's in active period, we don't use DTA. And so UE handovers from VoWiFi to VoLTE immediately
						{
							Crh++;
							Mob_Time += Tl++;
							e = Generate_arrival_event(LW_ARRIVAL, Mob_Time);
							List << *e;
						}
						break;

					case TD_EXPIRED:
						//printf("== TD_EXPIRED ==\n");
						//If td expired before tl, UE handovers from VoWiFi to VoLTE
						inTdPeriod = false;
						if(inActivePeriod)
						{
							Npda += floor((Td_Time - preU_Time)/Da);
						}
						else
						{
							if(isFirstE)
							{
								isFirstE = false;
								Npds += floor((Td_Time - startTd_Time)/Ds);
							}
							else
							{
								Npds += floor((Td_Time - preU_Time)/Ds);
							}
						}
						Npd += Npda + Npds;

						total_Npda += Npda;
						total_Npds += Npds;

						Crh++;
						break;

					case WW_ARRIVAL:
						//printf("== WW_ARRIVAL ==\n");
						inTdPeriod = false;
						if(inActivePeriod)
						{
							Npda += floor((Mob_Time - preU_Time)/Da);
						}
						else
						{
							if(isFirstE)
							{
								isFirstE = false;
								Npds += floor((Mob_Time - startTd_Time)/Ds);
							}
							else
							{
								Npds += floor((Mob_Time - preU_Time)/Ds);
							}
						}
						Npd += Npda + Npds;

						total_Npda += Npda;
						total_Npds += Npds;

						Crh++;
						Crh_0++;
						Mob_Time += Tw++;
						e = Generate_arrival_event(WT_ARRIVAL, Mob_Time);
						List << *e;
						break;

					default:
						break;			
				}
			}
			/* Events handling ends */
			
			
			if(Crh_0!=0)
			{
				sum_Crh += Crh;
				sum_Crh_0 += Crh_0;
				sum_Rh += ((double)(Crh_0-Crh)/Crh_0);
			}
			if(Np!=0)
			{
				sum_PLR += (double)Npd/Np;
				sum_PLR_a += (double)total_Npda/Np;
				sum_PLR_s += (double)total_Npds/Np;
			}
		}
		/* Simulation ends */


		// printf("\nE[Td] = %f", mean_Td);
		fp << "E[Td] = " << mean_Td << endl;
		// printf("\nMathematical Analysis\n");
		// printf("E[Nh(0)] = %f\n", E_Nh_0);
		// printf("E[Nrh(Td)] = %f\n", E_Nrh_td);
		// printf("E[Nh(Td)] = %f\n", E_Nh_td);
		// printf("Rh(Td) = %f\n", Rh_td);

		// printf("\nSimulation Analysis\n");
		// printf("E[Nh(0)] = %f\n", (double)sum_Crh_0/simuTimes);
		// printf("E[Nrh(Td)] = %f\n", ((double)(sum_Crh_0-sum_Crh))/simuTimes);
		// printf("E[Nh(Td)] = %f\n", (double)sum_Crh/simuTimes);
		// printf("Rh(Td) = %f\n", 1-(double)sum_Crh/sum_Crh_0);
		
		// the number of handover
		// fp << (double)sum_Crh_0/simuTimes << endl;
		// fp << (double)sum_Crh/simuTimes << endl;

		// Rh(td)
		fp << 1-(double)sum_Crh/sum_Crh_0 << endl;

		// printf("\nError Rate\n");
		// printf("E[Nh(0)] = %.6f%%\n", (((double)sum_Crh_0/simuTimes) - E_Nh_0)/E_Nh_0*100);
		// printf("E[Nrh(Td)] = %.6f%%\n", (((double)(sum_Crh_0-sum_Crh)/simuTimes) - E_Nrh_td)/E_Nrh_td*100);
		// printf("E[Nh(Td)] = %.6f%%\n", (((double)sum_Crh/simuTimes) - E_Nh_td)/E_Nh_td*100);
		// printf("Rh(Td) = %.6f%%\n", (abs(1-(double)sum_Crh/sum_Crh_0 - Rh_td)/Rh_td)*100);

		// printf("\nPacket Loss Rate\n");
		// printf("Packet loss = %.3f %%\n", sum_PLR/simuTimes*100);
		// printf("Packet loss during active period = %.3f %%\n", sum_PLR_a/simuTimes*100);
		// printf("Packet loss during silent period = %.3f %%\n", sum_PLR_s/simuTimes*100);
		
		// packet loss rate
		fp << sum_PLR/simuTimes << endl;
		fp << sum_PLR_a/simuTimes << endl;
		// fp << sum_PLR_s/simuTimes << endl;
	}
	fp.close();
	printf("\nSuccess!\n");
	return 0;
}
