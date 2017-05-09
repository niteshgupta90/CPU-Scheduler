#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "Scheduler.h"
#include "rngs.h"

#undef DBG_RESPONSE
#undef DBG_SYSTEM_STATE
#undef DBG_FUTURE_EVENT
#undef DBG_PRINT_RANDOM
#define DBG_VALIDATION
#undef DBG_REGENRATION
#undef DBG_SYSTEM_STATE
#undef DBG_CYCLETIME
#undef DBG_MVA
//#define DBG_REGENRATION

int num_cycles = 0;
int old_num_cycles = 0;
int RCN = 0;
double sum_of_cycle_t = 0.0;
double old_sum_of_cycle_t = 0.0;
double ThinkTime = 0.0;
int NumberThinkTimes = 0;
double SwInTime = 0.0;
int NumberSwIn = 0;
double CPUTime = 0.0;
int NumberCPU = 0;
double IO1Time = 0.0;
int NumberIO1 = 0;
double IO2Time = 0.0;
int NumberIO2 = 0;

int sim_Num = 0;
double lower_Out_R = 0.0;
double InRange_R = 0.0;
double upper_Out_R = 0.0;
double lowerConfidence_R = 0.0;
double upperConfidence_R = 0.0;
double Mean_R = 0.0;
int aboveMean_R = 0;
int belowMean_R = 0;

double lower_Out_A = 0.0;
double InRange_A = 0.0;
double upper_Out_A = 0.0;
double lowerConfidence_A = 0.0;
double upperConfidence_A = 0.0;
double Mean_A = 0.0;
int aboveMean_A = 0;
int belowMean_A = 0;

int q_size = 0;
bool halt;      /* End of simulation flag */
int nsys = 0;          /* Number of customers in system */
int nWorkStation = 0;
int nActive = 0;
int narr = 0;          /* Number of Arrivals */
int nReserve = 0;
int nSwapIN = 0;
int nCPU = 0;
int nIO1 = 0;
int nIO2 = 0;
int ncom_Reserve = 0;
int ncom_CPU = 0;
int ncom_IO1 = 0;
int ncom_IO2 = 0;
int ncom_SwapIN = 0;
int ncom = 0;          /* Number of Completions */
int node_number = 0;
int event_counter = 0; /* Number of the events processed by the simulkator*/
double End_time = SIMULATON_END;

double Busy = 0.0;       /* Busy time */
double stime = 0.0;       /* Clock of the simulator - Simulation time */
double departure_t = 0.0;   /* time of arrival read from trace */
double service_t = 0.0;   /* service time read from trace */
double process_t = 0.0;   /* service time read from trace */

double q_size1 = 0.0,q_size2 = 0.0,q_size3 = 0.0,q_size4 = 0.0,q_size5 = 0.0;
double Area_q1 = 0.0,Area_q2 = 0.0,Area_q3 = 0.0,Area_q4 = 0.0,Area_q5 = 0.0;
double Area_w1 = 0.0,Area_w2 = 0.0,Area_w3 = 0.0,Area_w4 = 0.0,Area_w5 = 0.0;
double Area_sys = 0.0;
double Busy_sys = 0.0;

double nRegeneration_R = 0.0;
double Sv_R = 0.0;
double Sa_R = 0.0;
double Sn_R = 0.0;
long double Saa_R = 0.0;
long double Svv_R = 0.0;
long double Snn_R = 0.0;
long double Sav_R = 0.0;
long double San_R = 0.0;
int ncomOld_R = 0;
double BetaOld_R = 0.0;
double wbar_R = 0.0;
double nbar_R = 0.0;
double width_nbar_R = 0.0;
double width_wbar_R = 0.0;
int n_cycle_WS[NUM_OF_MACHINES];
double R_avg = 0.0;
double Sum_R = 0.0;

double nRegeneration_A = 0.0;
double Sv_A = 0.0;
double Sa_A = 0.0;
double Sn_A = 0.0;
long double Saa_A = 0.0;
long double Svv_A = 0.0;
long double Snn_A = 0.0;
long double Sav_A = 0.0;
long double San_A = 0.0;
int ncomOld_A = 0;
double BetaOld_A = 0.0;
double wbar_A = 0.0;
double nbar_A = 0.0;
double width_nbar_A = 0.0;
double width_wbar_A = 0.0;

int ncom_Active = 0;          /* Number of Completions */
double Area_WS = 0.0;
double Busy_WS = 0.0;
int n_cycle_Active[NUM_OF_MACHINES];
double Sum_A = 0.0;
double Area_Active = 0.0;
int nActive_P = 0;

double sum = 0.0;
double V[]={1,2.5,250,16.25,6.25};
double S[]={0,0.210,0.0030,0.040,0.180};
double ST[]={D,LI,LI,LI,LI};
double Z_bar[]={5,5,5,5,5};

//double V[]={1,2.5,25,2.25};
//double S[]={0,0.210,0.030,0.040};
//double ST[]={D,LI,LI,LI};
//double Z_bar[]={5,5,5,5};

double n[M][N+1];
double w[M][N+1];
double X[M][N+1];
double U[M][N+1];

double Exponential(double m)
{
    return (-m * log(1.0 - Random()));
}

double Uniform(double a, double b)
{
    return (a + (b - a) * Random());
}

double HyperExponential()
{
    double A1,A2,Y,val;
    
    A1 = ALPHA;
    A2 = A1 + BETA;

    Y = Uniform(0,1);
    if(Y<=A1){
        val = Exponential(MU1);
    }
    else{
        val = Exponential(MU2);
    }
    
    return val;
}

double GetDeparture_Workstation(void)
{
    double val;
    SelectStream(1);
    val = Exponential(Z);
    ThinkTime = ThinkTime + val;
    NumberThinkTimes++;
#ifdef DBG_PRINT_RANDOM
    printf("Arrival value: %6.2f\n",val);
#endif
    return val;
}

double GetService_SwapIN(void)
{
    double val;
    SelectStream(2);
    val = Exponential(S2);
    SwInTime = SwInTime + val;
    NumberSwIn++;
#ifdef DBG_PRINT_RANDOM
    printf("Service value: %6.2f\n",val);
#endif
    return val;
}

double GetService_CPU(void)
{
    double val;
    SelectStream(3);
    val = Exponential(S3);
    CPUTime = CPUTime + val;
    NumberCPU++;
#ifdef DBG_PRINT_RANDOM
    printf("Service value: %6.2f\n",val);
#endif
    return val;
}

double GetService_IO1(void)
{
    double val;
    SelectStream(4);
    val = Exponential(S4_IO1);
    IO1Time = IO1Time + val;
    NumberIO1++;
#ifdef DBG_PRINT_RANDOM
    printf("Service value: %6.2f\n",val);
#endif
    return val;
}

double GetService_IO2(void)
{
    double val;
    SelectStream(5);
    val = Exponential(S4_IO2);
    IO2Time = IO2Time + val;
    NumberIO2++;
#ifdef DBG_PRINT_RANDOM
    printf("Service value: %6.2f\n",val);
#endif
    return val;
}

int possibleInterruption(){
    
    SelectStream(6);
    double p = Uniform(0,1);
#ifdef DBG_PRINT_RANDOM
    printf("Decision probability value: %6.2f\n",p);
#endif
    if(p<PROB_IO1){
        return IO1;
    }
    else if(p<PROB_IO2+PROB_IO1){
        return IO2;
    }
    else if(p<PROB_IO2+PROB_IO1+PROB_PROCEED_COMP){
        return PROCEED_COMP;
    }
    else{
        return PROCEED_CPU;
    }
}

bool isTrueCompeletion(){
    
    SelectStream(7);
    double p = Uniform(0,1);
#ifdef DBG_PRINT_RANDOM
    printf("Decision probability value: %6.2f\n",p);
#endif
    if(p<PROB_COMPLETEION){
        return true;
    }
    else
        return false;
}

void computeConfidenceInterval_Response(){
    double delta, z;
    double rhat;
    double numerator, factor;
    z= FACTOR;
    /* Confidence interval for average waiting time */
    rhat = Sa_R/Sn_R;
    wbar_R = rhat;
    numerator = Saa_R - 2.0*rhat*San_R + rhat*rhat*Snn_R;
    numerator = sqrt(numerator);
    factor = nRegeneration_R/(nRegeneration_R-1.0);
    factor = sqrt(factor);
    delta = factor*(numerator/Sn_R);
    width_wbar_R = z*delta;
    /* Confidence interval for average number of machines */
    rhat = Sa_R/Sv_R;
    nbar_R = rhat;
    numerator = Saa_R - 2.0*rhat*Sav_R + rhat*rhat*Svv_R;
    numerator = sqrt(numerator);
    factor = nRegeneration_R/(nRegeneration_R-1.0);
    factor = sqrt(factor);
    delta = factor*(numerator/Sv_R);
    width_nbar_R = z*delta;
}

void computeConfidenceInterval_Active(){
    double delta, z;
    double rhat;
    double numerator, factor;
    z= FACTOR;
    /* Confidence interval for average waiting time */
    rhat = Sa_A/Sn_A;
    wbar_A = rhat;
    numerator = Saa_A - 2.0*rhat*San_A + rhat*rhat*Snn_A;
    numerator = sqrt(numerator);
    factor = nRegeneration_A/(nRegeneration_A-1.0);
    factor = sqrt(factor);
    delta = factor*(numerator/Sn_A);
    width_wbar_A = z*delta;
    /* Confidence interval for average number of machines */
    rhat = Sa_A/Sv_A;
    nbar_A = rhat;
    numerator = Saa_A - 2.0*rhat*Sav_A + rhat*rhat*Svv_A;
    numerator = sqrt(numerator);
    factor = nRegeneration_A/(nRegeneration_A-1.0);
    factor = sqrt(factor);
    delta = factor*(numerator/Sv_A);
    width_nbar_A = z*delta;
}

int decideToStop_Response(int event_type){
    int limit;
    limit = MINIMUM_NUMBER_REG_CYCLES;
    if ((nRegeneration_R > limit)&&((width_wbar_R/wbar_R)< (PRECISION/2.0)))
        event_type = END;
    return event_type;
}

int decideToStop_Active(int event_type){
    int limit;
    limit = MINIMUM_NUMBER_REG_CYCLES;
    if ((nRegeneration_A > limit)&&((width_wbar_A/wbar_A)< (PRECISION/2.0)))
        event_type = END;
    return event_type;
}

void reset_measures_Response(){
    int i;
    Area_sys = 0;
    for( i = 0; i < NUM_OF_MACHINES; i++)
        n_cycle_WS[i] = -1;
}

void reset_measures_Active(){
    int i;
    Area_Active = 0;
    for( i = 0; i < NUM_OF_MACHINES; i++)
        n_cycle_Active[i] = -1;
}

bool regPoint_Response(int event_type) {
    bool Reg= false;
    int i;
    if ((event_type == DEPARTURE_FROM_IO2)&&(nIO2==NUM_OF_MACHINES)) {
        Reg = true;
        i = 0;
        while ((n_cycle_WS[i] >= N_CYCLE_MIN) && (i < NUM_OF_MACHINES))
            i++;
        if(i < NUM_OF_MACHINES)
            Reg = false;
    }
    return Reg;
}

bool regPoint_Active(int event_type) {
    bool Reg= false;
    int i;
    if ((event_type == DEPARTURE_FROM_IO2)&&(nIO2==NUM_OF_MACHINES)) {
        Reg = true;
        i = 0;
        while ((n_cycle_Active[i] >= N_CYCLE_MIN) && (i < NUM_OF_MACHINES))
            i++;
        if(i < NUM_OF_MACHINES)
            Reg = false;
    }
    return Reg;
}

void getRegenrationValues_Response()
{
    double RegCycleT, RegCycleN;
    nRegeneration_R++;
    RegCycleT = stime - BetaOld_R;
    BetaOld_R = stime;
    RegCycleN = ncom - ncomOld_R;
    ncomOld_R = ncom;
    Sa_R = Sa_R + Area_sys;
    Saa_R = Saa_R + Area_sys * Area_sys;
    Sv_R = Sv_R + RegCycleT;
    Svv_R = Svv_R + RegCycleT * RegCycleT;
    Sn_R = Sn_R + RegCycleN;
    Snn_R = Snn_R + RegCycleN * RegCycleN;
    Sav_R = Sav_R + Area_sys * RegCycleT;
    San_R = San_R + Area_sys * RegCycleN;
#ifdef DBG_REGENRATION
    printf("\nLength of Regeneration Time            = %10.2f",RegCycleT);
    printf("\nLength of Regeneration Number of Cycle = %10.2f",RegCycleN);
    printf("\nAverage Response from Area          = %10.2f",Area_sys/RegCycleN);
#endif
}

void getRegenrationValues_Active()
{
    double RegCycleT, RegCycleN;
    nRegeneration_A++;
    RegCycleT = stime - BetaOld_A;
    BetaOld_A = stime;
    RegCycleN = ncom_Active - ncomOld_A;
    ncomOld_A = ncom_Active;
    Sa_A = Sa_A + Area_Active;
    Saa_A = Saa_A + Area_Active * Area_Active;
    Sv_A = Sv_A + RegCycleT;
    Svv_A = Svv_A + RegCycleT * RegCycleT;
    Sn_A = Sn_A + RegCycleN;
    Snn_A = Snn_A + RegCycleN * RegCycleN;
    Sav_A = Sav_A + Area_Active * RegCycleT;
    San_A = San_A + Area_Active * RegCycleN;
}

double avgResponseTime(){
    return Sum_R/ncom;
}

double avgActiveTime(){
    return Sum_A/ncom_Active;
}

void validateResponse()
{
    double Response = 0.0;
    Response = getAverageResponse();
    if(Response<lowerConfidence_R){
        lower_Out_R++;
    }
    else if(Response<=upperConfidence_R){
        InRange_R++;
    }
    else{
        upper_Out_R++;
    }
    
    if(wbar_R>=Response)
        aboveMean_R++;
    else
        belowMean_R++;
}

void validateActive()
{
    double Active = 0.0;
    Active = getAverageActive();
    if(Active<lowerConfidence_A){
        lower_Out_A++;
    }
    else if(Active<=upperConfidence_A){
        InRange_A++;
    }
    else{
        upper_Out_A++;
    }
    
    if(wbar_A>=Active)
        aboveMean_A++;
    else
        belowMean_A++;
}

void setConfidenceParameters(){
    lowerConfidence_R = wbar_R-width_wbar_R;
    upperConfidence_R = wbar_R+width_wbar_R;
    Mean_R = wbar_R;
    lowerConfidence_A = wbar_A-width_wbar_A;
    upperConfidence_A = wbar_A+width_wbar_A;
    Mean_A = wbar_A;
}

int main(int argc, char* argv[]){
    PlantSeeds(DEFAULT);
    doMVA();
#ifdef DBG_VALIDATION
    for(sim_Num=0;sim_Num<NUM_OF_SIMULATION;sim_Num++){
        simulate();
        setConfidenceParameters();
        validateResponse();
        validateActive();
        print_results();
    }
    //report();
    destroy_list();
    destroy_queue();
    printValidationResults();
#else
    simulate();
    report();
#endif
    //report();
    //getchar();
    return 0;
}

void simulate() {
    /* Simulation core */
    initialize();
#ifdef DBG_SYSTEM_STATE
    printf("\nSystem State for first 30 events:");
#endif
    while (!(halt)){
        engine();
#ifdef DBG_SYSTEM_STATE
        if(event_counter <= 30){
           printSystemState();
        }
#endif
    }
}

void engine(void){
    int event_type;
    double  oldtime;
    double  interval;
    tree new_event;
    int q1,q2,q3,q4,q5;
    q1 = getQueue_size(0);
    q2 = getQueue_size(1);
    q3 = getQueue_size(2);
    q4 = getQueue_size(3);
    q5 = getQueue_size(4);
    
    
    /* Get the first event notice from Future Event List */
    new_event = event_pop();
    
    /* update clock */
    oldtime = stime;
    stime = new_event->event.occur_time;
    interval = stime - oldtime;
    
    /* Collect statistics */
    
    if(nReserve>0){
        Area_w1 += nReserve * interval;
        Area_q1 += (nReserve-1) * interval;
    }
    if(nCPU>0){
        Area_w3 += nCPU * interval;
        Area_q3 += (nCPU-1) * interval;
    }
    if(nSwapIN>0){
        Area_w2 += nSwapIN * interval;
        Area_q2 += (nSwapIN-1) * interval;
    }
    if(nIO1>0){
        Area_w4 += nIO1 * interval;
        Area_q4 += (nIO1-1) * interval;
    }
    if(nIO2>0){
        Area_w5 += nIO2 * interval;
        Area_q5 += (nIO2-1) * interval;
    }

    if (nsys > 0){
        Busy_sys = Busy_sys + interval;
        Area_sys = Area_sys + nsys * interval;
    }
    
    if (nActive > 0){
        Area_Active = Area_Active + nActive * interval;
    }
    
    if (nWorkStation > 0){
        Area_WS = Area_WS + nWorkStation * interval;
    }
    
    /* Identify and process current event */
    event_type = new_event->event.type;
   
    /* Identify and process regeneration point for Response Time */
    if (regPoint_Response(event_type)) {
        getRegenrationValues_Response();
        getRegenrationValues_Active();
#ifdef DBG2
        printRegenrationvalues();
        printFutureEventList();
        printQueue(0);
        printQueue(1);
        printQueue(2);
        printQueue(3);
        printQueue(4);
#endif
#ifdef DBG_RESPONSE
        printSystemValues();
#endif
        if (nRegeneration_R > 30 && NUM_OF_MACHINES>1 ) {
            computeConfidenceInterval_Response();
            computeConfidenceInterval_Active();
            event_type = decideToStop_Response(event_type);
            //printSystemState();
            //print_results();
            //printAllDepartures();
        }
        if (nRegeneration_R > 30 && NUM_OF_MACHINES==1 ) {
            event_type = END;
            printMeanvalues();
        }
        reset_measures_Response();
        reset_measures_Active();
    }
    
    switch(event_type){
        case ARRIVAL_AT_WORKSTATION :
            arrivalAtWorkStation(new_event);
            break;
        case ARRIVAL_AT_RESERVE :
            arrivalAtReserve(new_event);
            break;
        case ARRIVAL_AT_SWAP_IN :
            arrivalAtSwapIN(new_event);
            break;
        case ARRIVAL_AT_CPU :
            arrivalAtCPU(new_event);
            break;
        case ARRIVAL_AT_IO1 :
            arrivalAtIO1(new_event);
            break;
        case ARRIVAL_AT_IO2 :
            arrivalAtIO2(new_event);
            break;
        case ARRIVAL_AT_SWAPOUT :
            arrivalAtSwapOut(new_event);
            break;
        case DEPARTURE_FROM_WORKSTATION :
            departureFromWorkStation(new_event);
            break;
        case DEPARTURE_FROM_RESERVE :
            departureFromReserve(new_event);
            break;
        case DEPARTURE_FROM_SWAP_IN :
            departureFromSwapIN(new_event);
            break;
        case DEPARTURE_FROM_CPU :
            departureFromCPU(new_event);
            break;
        case DEPARTURE_FROM_IO1 :
            departureFromIO1(new_event);
            break;
        case DEPARTURE_FROM_IO2 :
            departureFromIO2(new_event);
            break;
        case DEPARTURE_FROM_SWAPOUT :
            departureFromSwapOut(new_event);
            break;
        case END : halt = true;
            break;
    }
    event_counter++;
}

void initialize(){
    struct node* first_notice;
    int i;
    /* Control Settings  */
    halt = false;

    /* Basic Statistic Measures  */
    nWorkStation    = 0;
    nsys   = 0;
    narr   = 0;
    ncom   = 0;
    Busy   = 0;
    stime  = 0;
    Area_q1 = 0;

    q_size          = 0;
    /* Basic  counters  */
    job_number      = 0;
    node_number     = 0;
    nSwapIN  = 0;
    q_size = 0;
    nActive = 0;
    nReserve = 0;
    nSwapIN = 0;
    nCPU = 0;
    nIO1 = 0;
    nIO2 = 0;
    ncom_Reserve = 0;
    ncom_CPU = 0;
    ncom_IO1 = 0;
    ncom_IO2 = 0;
    ncom_SwapIN = 0;
    ncom = 0;          /* Number of Completions */
    node_number = 0;
    event_counter = 0; /* Number of events processed by the simulator*/
    End_time = SIMULATON_END;

    Busy = 0.0;
    stime = 0.0;
    departure_t = 0.0;
    service_t = 0.0;
    process_t = 0.0;

    q_size1 = 0.0,q_size2 = 0.0,q_size3 = 0.0,q_size4 = 0.0,q_size5 = 0.0;
    Area_q1 = 0.0,Area_q2 = 0.0,Area_q3 = 0.0,Area_q4 = 0.0,Area_q5 = 0.0;
    Area_w1 = 0.0,Area_w2 = 0.0,Area_w3 = 0.0,Area_w4 = 0.0,Area_w5 = 0.0;
    Area_sys = 0.0;
    Busy_sys = 0.0;

    nRegeneration_R = 0.0;
    Sv_R = 0.0;
    Sa_R = 0.0;
    Sn_R = 0.0;
    Saa_R = 0.0;
    Svv_R = 0.0;
    Snn_R = 0.0;
    Sav_R = 0.0;
    San_R = 0.0;
    ncomOld_R = 0;
    BetaOld_R = 0.0;
    wbar_R = 0.0;
    nbar_R = 0.0;
    width_nbar_R = 0.0;
    width_wbar_R = 0.0;
    R_avg = 0.0;
    Sum_R = 0.0;

    nRegeneration_A = 0.0;
    Sv_A = 0.0;
    Sa_A = 0.0;
    Sn_A = 0.0;
    Saa_A = 0.0;
    Svv_A = 0.0;
    Snn_A = 0.0;
    Sav_A = 0.0;
    San_A = 0.0;
    ncomOld_A = 0;
    BetaOld_A = 0.0;
    wbar_A = 0.0;
    nbar_A = 0.0;
    width_nbar_A = 0.0;
    width_wbar_A = 0.0;
    ncom_Active = 0;
    Area_WS = 0.0;
    Busy_WS = 0.0;
    Sum_A = 0.0;
    Area_Active = 0.0;
    nActive_P = 0;
    sum = 0.0;
    /* Future Event List and additional structures*/
    root  = NULL;
    queue = NULL;
    last  = NULL;
    
    for(i=0;i<QUEUE_NUM;i++)
        q[i] = createQueue();
        
    /* Initialize Event notice of first arrival and Schedule first event */
    for(i=0;i<NUM_OF_MACHINES;i++){
        nWorkStation++;
        departure_t = GetDeparture_Workstation();
        first_notice = get_new_node();
        first_notice->event.departure_time = departure_t;
        first_notice->event.service_time = 0;
        first_notice->event.occur_time =  departure_t + stime;
        first_notice->event.start_time = first_notice->event.occur_time;
        
        first_notice->event.type = DEPARTURE_FROM_WORKSTATION;
        first_notice->next = NULL;
        first_notice->event.index = i;
        first_notice->event.processing = false;
        sprintf(first_notice->event.name, "P%d", (job_number++));
        
        recordCycle_Time[first_notice->event.index].id = first_notice->event.index; //set process id
        recordCycle_Time[first_notice->event.index].start = stime; //set cycle start time
        recordCycle_Time[first_notice->event.index].end = INFINITY;
        schedule(first_notice);
    }
}

void printCycleInformation(int processIndex){
    double difference;
    num_cycles++;
    difference = recordCycle_Time[processIndex].end-recordCycle_Time[processIndex].start;
    printf("\nCycle Nr. %d;  Process ID   = %d",num_cycles,processIndex);
    printf("\nStart Time   = %f",recordCycle_Time[processIndex].start);
    printf("\nEnd Time     = %f;   Cycle Time    = %f\n",recordCycle_Time[processIndex].end, difference );
    sum_of_cycle_t = sum_of_cycle_t + difference;
}

void arrivalAtWorkStation(struct node* node_event){
    nWorkStation++;
    nsys--;
    ncom++;
    n_cycle_WS[node_event->event.index]++;
    recordCycle_Time[node_event->event.index].end = stime;  //set cycle end time
#ifdef DBG_CYCLETIME
    printCycleInformation(node_event->event.index);  //print cycle information
#endif
    //Calculate Response Time
    node_event->event.time_R = stime - node_event->event.start_time;
    departure_t = GetDeparture_Workstation();
    Sum_R += node_event->event.time_R;
    node_event->event.departure_time = departure_t;
    node_event->event.occur_time =  departure_t + stime;
    node_event->event.start_time = node_event->event.occur_time;
    recordCycle_Time[node_event->event.index].start = stime; //set cycle start time
    node_event->event.processing = false;
    node_event->event.type = DEPARTURE_FROM_WORKSTATION;
    schedule(node_event);
}

void departureFromWorkStation(struct node* node_event){
    nsys++;
    nWorkStation--;
    node_event->event.type = ARRIVAL_AT_SWAP_IN;
    node_event->event.occur_time = stime;
    schedule(node_event);
}

void arrivalAtReserve(struct node* node_event)
{
    /* Update statistics */
    nReserve++;
    if (nActive<MPD)
    {
        node_event->event.type = DEPARTURE_FROM_RESERVE;
        node_event->event.occur_time = stime;
        schedule(node_event);
    }
    else {
        /* Process arrival at busy server */
        enqueue(node_event,0);
    }
}

void departureFromReserve(struct node* node_event)
{
    ncom_Reserve++;
    nActive++;
    nReserve--;
    node_event->event.start_time_Active = stime;
    node_event->event.type = ARRIVAL_AT_SWAP_IN;
    node_event->event.occur_time = stime;
    schedule(node_event);
}

void arrivalAtSwapIN(struct node* node_event){
    nSwapIN++;
    nActive++;

    if (nSwapIN == 1)
    {
        service_t = GetService_SwapIN();
        node_event->event.type = DEPARTURE_FROM_SWAP_IN;
        node_event->event.service_time = service_t;
        node_event->event.occur_time = stime + node_event->event.service_time;
        schedule(node_event);
    }
    else {
        /* Process arrival at busy server */
        enqueue(node_event,1);
    }
}

void departureFromSwapIN(struct node* node_event){
    struct node* next_job;
    ncom_SwapIN++;
    nSwapIN--;
    node_event->event.type = ARRIVAL_AT_CPU;
    node_event->event.start_time_Active = stime;
    node_event->event.occur_time = stime;
    schedule(node_event);
    
    q_size = getQueue_size(1);
    if (q_size > 0) {
        /* Process departure from a server with a queue*/
        next_job = dequeue(1);
        service_t = GetService_SwapIN();
        next_job->event.service_time = service_t;
        next_job->event.type = DEPARTURE_FROM_SWAP_IN;
        next_job->event.occur_time = stime + next_job->event.service_time;
        schedule(next_job);
    }
}

void arrivalAtCPU(struct node* node_event){
    nCPU++;
    if (nCPU == 1)
    {
        process_t = GetService_CPU();
        node_event->event.service_time = process_t;
        node_event->event.type = DEPARTURE_FROM_CPU;
        node_event->event.occur_time = stime + node_event->event.service_time;
        schedule(node_event);
    }
    else {
        /* Process arrival at busy server */
        enqueue(node_event,2);
    }
}

void departureFromCPU(struct node* node_event){
    int possibleIntrpt;
    ncom_CPU++;
    nCPU--;
    possibleIntrpt = possibleInterruption();
    {
        if(possibleIntrpt == IO1){
            node_event->event.type = ARRIVAL_AT_IO1;
        }
        else if(possibleIntrpt == IO2){
            node_event->event.type = ARRIVAL_AT_IO2;
        }
        else if(possibleIntrpt == PROCEED_COMP){
            node_event->event.type = ARRIVAL_AT_SWAPOUT;
        }
        else if(possibleIntrpt == PROCEED_CPU){
            node_event->event.type = ARRIVAL_AT_CPU;
        }
    }
    node_event->event.occur_time = stime;
    schedule(node_event);
    
    q_size = getQueue_size(2);
    if (q_size > 0) {
        struct node* next_job;
        /* Process departure from a server with a queue*/
        next_job = dequeue(2);
        process_t = GetService_CPU();
        next_job->event.service_time = process_t;
        next_job->event.type = DEPARTURE_FROM_CPU;
        next_job->event.occur_time = stime + node_event->event.service_time;
        schedule(next_job);
    }
}

void arrivalAtIO1(struct node* node_event){
    nIO1++;
    if (nIO1 == 1)
    {
        service_t = GetService_IO1();
        node_event->event.type = DEPARTURE_FROM_IO1;
        node_event->event.service_time = service_t;
        node_event->event.occur_time = stime + node_event->event.service_time;
        schedule(node_event);
    }
    else {
        /* Process arrival at busy server */
        enqueue(node_event,3);
    }
}

void departureFromIO1(struct node* node_event){
    ncom_IO1++;
    nIO1--;
    node_event->event.type = ARRIVAL_AT_CPU;
    node_event->event.occur_time = stime;
    schedule(node_event);
    
    q_size = getQueue_size(3);
    if (q_size > 0) {
        struct node* next_job;
        /* Process departure from a server with a queue*/
        next_job = dequeue(3);
        service_t = GetService_IO1();
        next_job->event.service_time = service_t;
        next_job->event.type = DEPARTURE_FROM_IO1;
        next_job->event.occur_time = stime + next_job->event.service_time;
        schedule(next_job);
    }
}

void arrivalAtIO2(struct node* node_event){
    nIO2++;
    if (nIO2 == 1)
    {
        service_t = GetService_IO2();
        node_event->event.type = DEPARTURE_FROM_IO2;
        node_event->event.service_time = service_t;
        node_event->event.occur_time = stime + node_event->event.service_time;
        schedule(node_event);
    }
    else {
        /* Process arrival at busy server */
        enqueue(node_event,4);
    }
}

void departureFromIO2(struct node* node_event){
    ncom_IO2++;
    nIO2--;
    node_event->event.type = ARRIVAL_AT_CPU;
    node_event->event.occur_time = stime;
    schedule(node_event);
    
    q_size = getQueue_size(4);
    if (q_size > 0) {
        struct node* next_job;
        /* Process departure from a server with a queue*/
        next_job = dequeue(4);
        service_t = GetService_IO2();
        next_job->event.service_time = service_t;
        next_job->event.type = DEPARTURE_FROM_IO2;
        next_job->event.occur_time = stime + next_job->event.service_time;
        schedule(next_job);
    }
}

void arrivalAtSwapOut(struct node* node_event){
    nActive--;
    ncom_Active++;
    n_cycle_Active[node_event->event.index]++;
    
    //Calculate Active Time
    node_event->event.time_A = stime - node_event->event.start_time_Active;
    Sum_A += node_event->event.time_A;
    
    node_event->event.type = DEPARTURE_FROM_SWAPOUT;
    node_event->event.occur_time = stime;
    schedule(node_event);
}

void departureFromSwapOut(struct node* node_event){
    if(isTrueCompeletion()==true){
        node_event->event.type = ARRIVAL_AT_WORKSTATION;
    }
    else{
        node_event->event.type = ARRIVAL_AT_SWAP_IN;
    }
    node_event->event.occur_time = stime;
    schedule(node_event);
}

void print_results(){
    printf("\nSimulation Results:");
    printf("\n---------------------------------");
    printf("\nNumber of considered events                 = %10d", event_counter);
    printf("\nLast Simulation Time                        = %10.4f", stime);
    printf("\nNumber of Regenration                       = %10d\n", (int)nRegeneration_R);
    printf("\n---------------------------------");
    printf("\nResults for Response Time:\n");
    printf("\nConfidence Interval with 90 pcnt confidence = %10.4f - %2.4f", wbar_R-width_wbar_R,wbar_R+width_wbar_R);
    printf("\nConfidence Range                            = %10.4f", 2*width_wbar_R);
    printf("\nMean Response Time with Regeneration Method = %10.4f", wbar_R);
    printf("\nMean Response Time with MVA                 = %10.4f", getAverageResponse());
    printf("\n\n");
    printf("\nResults for Active Time:\n");
    printf("\nConfidence Interval with 90 pcnt confidence = %10.4f - %2.4f", wbar_A-width_wbar_A,wbar_A+width_wbar_A);
    printf("\nConfidence Range                            = %10.4f", 2*width_wbar_A);
    //printf("\nNumber of Regenration                       = %10.4f", nRegeneration_A);
    printf("\nMean Active Time with Regeneration Method   = %10.4f", wbar_A);
    printf("\nMean Active Time with MVA                   = %10.4f", getAverageActive());

    printf("\n--------------------------------------\n\n");
}

void report(){
    print_results();
    destroy_list();
    destroy_queue();
}

double getVariance(double sum_square, double sum)
{
    double variance =0;
    variance = (sum_square / (narr-1)) - ((sum * sum) / (narr * (narr-1)));
    return variance;
}

// A utility function to create a new linked list node.
struct QNode* newNode(struct node* k)
{
    struct QNode *temp = (struct QNode*)malloc(sizeof(struct QNode));
    temp->key = k;
    temp->next = NULL;
    return temp;
}

// A utility function to create an empty queue
struct Queue *createQueue()
{
    struct Queue *queue = (struct Queue*)malloc(sizeof(struct Queue));
    queue->front = queue->rear = NULL;
    return queue;
}

// The function to add a key k to q
void enqueue(struct node* k, int num)
{
    // Create a new LL node
    struct QNode *temp = newNode(k);
    
    // If queue is empty, then new node is front and rear both
    if (q[num]->rear == NULL)
    {
        q[num]->front = q[num]->rear = temp;
        return;
    }
    
    // Add the new node at the end of queue and change rear
    q[num]->rear->next = temp;
    q[num]->rear = temp;
}

// Function to remove a key from given queue q
struct node *dequeue(int num)
{
    // If queue is empty, return NULL.
    if (q[num]->front == NULL)
        return NULL;
    
    // Store previous front and move front one node ahead
    struct QNode *temp = q[num]->front;
    q[num]->front = q[num]->front->next;
    
    // If front becomes NULL, then change rear also as NULL
    if (q[num]->front == NULL)
        q[num]->rear = NULL;
    return temp->key;
}

void destroy_queue()
{
    int i=0;
    for(i=0;i<QUEUE_NUM;i++)
        free(q[i]);
}

void destroy_list()
{
    free(head);
}

//function to pop event from list
struct node* event_pop()
{
    struct node* n;
    n = head;
    head = head->next;
    return n;
}

struct node* returnHead(){
    return head;
}

// function to insert a new_node in a list.
void schedule(struct node* new_node)
{
    struct node* current;
    /* Special case for the head end */
    if (head == NULL || head->event.occur_time >= new_node->event.occur_time)
    {
        new_node->next = head;
        head = new_node;
    }
    else
    {
        /* Locate the node before the point of insertion */
        current = head;
        while (current->next!=NULL &&
               current->next->event.occur_time < new_node->event.occur_time)
        {
            current = current->next;
        }
        new_node->next = current->next;
        current->next = new_node;
    }
}

// A function to create a new node
struct node* get_new_node()
{
    /* allocate node */
    struct node* new_node =
    (struct node*) malloc(sizeof(struct node));
    
    /* put in the data  */
    new_node->next =  NULL;
    
    return new_node;
}

//A function to get the queue size
double getQueue_size(int num){
    double i=0.0;
    struct QNode *temp = q[num]->front;
    while(temp!=NULL){
        i++;
        temp= temp->next;
    }
    return i;
}

void initializeMVA(){
    int i,k;
    for(k=0;k<=N;k++){
        for(i=0;i<M;i++){
            n[i][k]=0;
            w[i][k]=0;
            U[i][k]=0;
            X[i][k]=0;
        }
    }
}

double getAverageResponse()
{
    double R = 0.0;
    //R = V[1]*S[1] + V[2]*S[2] + V[3]*S[3]+ V[4]*S[4];
    //CPU is the bottleneck station
    double X2_MAX,X0;
    X2_MAX = 1/S[2];
    X0 = X2_MAX/V[2];
    //printf("\nCPU Throughput = %f",X2_MAX);
    //printf("\nOverall Throughput = %f",X0);
    R = ((NUM_OF_MACHINES/X[0][N]) - Z_bar[0]);
    //printf("\nResponse Time = %f\n",R);
    return R;
}

double getAverageActive()
{
    int i;
    double A = 0.0, Sum = 0.0;
    for(i=1;i<M;i++){
        Sum += n[i][20];
    }
    A = Sum/X[1][20];
    return A;
}

void doMVA(){
    initializeMVA();
    int i,k;
    //(* INPUT: [ {Vi}, {Si }, {STi} ] *)
    //(* OUTPUT: [{Xi}, {Ui }, {ni}, {wi}] *)
    //(*Initialization*)
    for (i = 0;i<M;i++){
        n[i][0] = 0.0;
    }
    //(*Compute the performance measures*)
    for (k=1;k<=N;k++){
        for (i=0;i<M;i++){
            if (ST[i] == D){
                w[i][k] = Z_bar[i];
            }
            else{
                w[i][k] = S[i] * (1 + n[i][k-1]);
            }
        }
        sum = 0.0;
        for (i=0;i<M;i++){
            sum = sum + V[i] * w[i][k];
        }
        X[0][k] = k/sum;
        for (i=0;i<M;i++){
            X[i][k] = V[i] * X[0][k];
            if (ST[i] == D){
                n[i][k] = Z_bar[i] * X[i][k];
                U[i][k] = n[i][k]/k;
            }
            else{
                U[i][k] = S[i] * X[i][k];
                n[i][k] = U[i][k] * (1 + n[i][k-1]);
            }
        }
    }
#ifdef DBG_MVA
    printMVAResult("Throughput",X);
    printMVAResult("Utilization",U);
    printMVAResult("Average Queue Length",n);
    printMVAResult("Average Waiting Time",w);
#endif
}

void printMVAResult(char name[], double val[M][N+1]){
    int i,k;
    printf("\n%s:\n\n",name);
    for(i=0;i<M;i++){
        printf("\t      Station %d\t",i);
    }
    printf("\n");
    for(k=N;k<=N;k++){
        printf("C%d\t",k);
        for(i=0;i<M;i++){
            printf("\t%f\t",val[i][k]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void printSystemState(){
    printf("\nEvent_Num\tstime \tRegenration\tnCompleted\tnActive");
    printf("\n  %d     \t%.4f  \t%d         \t%d        \t%d",event_counter,stime,(int)nRegeneration_R,ncom ,ncom_Active);
    
    printf("\nnWorkStation  nSystem  nActive  nReserve  nSwapIN   nCPU    nIO1   nIO2");
    printf("\n\t%d          %d       %d       %d            %d        %d       %d     %d",nWorkStation,nsys,nActive,nReserve ,nSwapIN,nCPU,nIO1,nIO2);
    
    printf("\n");
}

void printTopQueueEvent(int num){
    if (q[num]->front != NULL){
        printf("\nQueue Event\n");
        printf("\nEvent Type        : %10d",q[num]->front->key->event.type);
        printf("\nEvent Name        : %10s",q[num]->front->key->event.name);
        printf("\nEvent Occur Time  : %10.2f",q[num]->front->key->event.occur_time);
        printf("\nEvent Service Time: %10.2f",q[num]->front->key->event.service_time);
        printf("\nEvent Departure Time: %10.2f\n",q[num]->front->key->event.departure_time);
    }
}

void printSystemvalues(){
    double partial_sum_of_cycle_t;
    int partial_num_cycles;
    partial_sum_of_cycle_t = sum_of_cycle_t - old_sum_of_cycle_t;
    old_sum_of_cycle_t = sum_of_cycle_t;
    partial_num_cycles = num_cycles - old_num_cycles;
    old_num_cycles = num_cycles;
    printf("\n  ======================= ARE WE STOPPING HERE? ====================");
    printf("\n  Total Cycle Time = %f; Number of Cycles = %d",
           sum_of_cycle_t,num_cycles);
    printf("\n Avg Cycle time = %f", (sum_of_cycle_t/num_cycles));
    printf("\n Avg Think time = %f (%f, %d) N=%d, V=%f",
           (ThinkTime/NumberThinkTimes),ThinkTime,NumberThinkTimes,ncom,(1.0*ncom/ncom));
    printf("\n Avg SwapIn time = %f (%f, %d) N=%d, V=%f",
           (SwInTime/NumberSwIn),SwInTime,NumberSwIn,ncom_SwapIN,(1.0*ncom_SwapIN/ncom));
    printf("\n Avg CPU time = %f (%f, %d) N=%d, V=%f",
           (CPUTime/NumberCPU),CPUTime,NumberCPU,ncom_CPU,(1.0*ncom_CPU/ncom));
    printf("\n Avg IO1 time = %f (%f, %d) N=%d, V=%f",
           (IO1Time/NumberIO1),IO1Time,NumberIO1,ncom_IO1,(1.0*ncom_IO1/ncom));
    printf("\n Avg IO2 time = %f (%f, %d) N=%d, V=%f",
           (IO2Time/NumberIO2),IO2Time,NumberIO2,ncom_IO2,(1.0*ncom_IO2/ncom));
    printf("\n  Partial Total Cycle Time = %f; Partial Number of Cycles = %d",
           partial_sum_of_cycle_t,partial_num_cycles);
    printf("\n   Partial Average Cycle time = %f", (partial_sum_of_cycle_t/partial_num_cycles));
    RCN = ncom - ncomOld_R;
    printf("\n   AREA_SYS = %f; RCN = %d; AVERAGE = %f",Area_sys,RCN, (Area_sys/RCN));
    printf("\nAverage Response Time               = %10.2f", Sum_R/ncom);
}

void printEvent(struct node* node_event){
    printf("\nEvent\n");
    printf("\nEvent Type        : %10d",node_event->event.type);
    printf("\nEvent Name        : %10s",node_event->event.name);
    printf("\nEvent Occur Time  : %10.2f",node_event->event.occur_time);
    printf("\nEvent Service Time: %10.2f",node_event->event.service_time);
    printf("\nEvent Departure Time: %10.2f\n",node_event->event.departure_time);
}

void printQueue(int num){
    struct QNode *queue_event = q[num]->front;
    if (q[num]->front != NULL){
        printf("\nQueue List:\n");
        printf("\nType\t\tName\tOccur Time\tService Time\tDeparture Time\n");
        
        while(queue_event!=NULL){
            printf("%4d\t\t%s\t\t\t%.3f\t\t%.3f\t\t%.3f\n",queue_event->key->event.type,queue_event->key->event.name,queue_event->key->event.occur_time,queue_event->key->event.service_time,queue_event->key->event.departure_time);
            queue_event = queue_event->next;
        }
    }
}

void printFutureEventList(){
    struct node* node_event;
    node_event = returnHead();
    static int num = 1;
    printf("\nEvent List: %d\n",num++);
    printf("\nType\t\tName\tOccur Time\tService Time\tDeparture Time\n");
    
    while(node_event!=NULL){
        printf("%4d\t\t%s\t\t\t%.4f\t\t%.4f\t\t%.4f\n",node_event->event.type,node_event->event.name,node_event->event.occur_time,node_event->event.service_time,node_event->event.departure_time);
        
        node_event = node_event->next;
    }
}

void printMeanvalues(){
    printf("\nNumber of considered events                 = %10d", event_counter);
    printf("\nLast Simulation Time                        = %10.2f", stime);
    printf("\nMean Response time                          = %10.2f", Sum_R/ncom);
    printf("\nMean Active Time                            = %10.2f", Sum_A/ncom_Active);
    printf("\n");
}

void printValidationResults(){
    printf("\nValidtion Results:\n");
    printf("\nNumber of Confidence Interval higher than theoretical Response Time    = %d",(int)upper_Out_R);
    printf("\nNumber of Confidence Interval Covering theoretical Response Time       = %d",(int)InRange_R);
    printf("\nNumber of Confidence Interval lower than theoretical Response Time     = %d",(int)lower_Out_R);
    printf("\n");
    printf("\nNumber of Confidence Interval higher than theoretical Active Time      = %d",(int)upper_Out_A);
    printf("\nNumber of Confidence Interval Covering theoretical Active Time         = %d",(int)InRange_A);
    printf("\nNumber of Confidence Interval lower than theoretical Active Time       = %d",(int)lower_Out_A);
    printf("\n");
    printf("\nNumber of Average Response Time higher than theoretical Response Time  = %d",(int)aboveMean_R);
    printf("\nNumber of Average Response Time lower than theoretical Response Time   = %d",(int)belowMean_R);
    printf("\nNumber of Average Active Time higher than theoretical Active Time      = %d",(int)aboveMean_A);
    printf("\nNumber of Average Active Time lower than theoretical Active Time       = %d",(int)belowMean_A);
    printf("\n_______________________\n");
}

void printRegenrationvalues()
{
    printf("\nNumber of Regeneration_R:  = %10.2f", nRegeneration_R);
    printf("\nArea_sys:                  = %10.2f", Area_sys);
    printf("\nSa_R:                      = %10.2f", Sa_R);
    printf("\nSaa_R:                     = %10.2Lf", Saa_R);
    printf("\nSn_R:                      = %10.2f", Sn_R);
    printf("\nSnn_R:                     = %10.2Lf", Snn_R);
    printf("\nSan_R:                     = %10.2Lf", San_R);
    printf("\nSa_R:                      = %10.2f", Sa_R);
}

void printAllDepartures()
{
    printf("\nSwapIN Q:                  = %10.2f", (double)ncom_SwapIN/ncom);
    printf("\nCPU Q:                     = %10.2f", (double)ncom_CPU/ncom);
    printf("\nIO1 Q:                     = %10.2f", (double)ncom_IO1/ncom);
    printf("\nIO2 Q:                     = %10.2f\n", (double)ncom_IO2/ncom);
}

