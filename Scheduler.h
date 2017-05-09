#define ARRIVAL 100
#define DEPARTURE 101
#define REST 102
#define END 103

#define ARRIVAL_AT_WORKSTATION 1
#define ARRIVAL_AT_RESERVE 2
#define ARRIVAL_AT_SWAP_IN 3
#define ARRIVAL_AT_CPU 4
#define ARRIVAL_AT_IO1 5
#define ARRIVAL_AT_IO2 6
#define ARRIVAL_AT_SWAPOUT 17
#define DEPARTURE_FROM_WORKSTATION 8
#define DEPARTURE_FROM_RESERVE 9
#define DEPARTURE_FROM_SWAP_IN 10
#define DEPARTURE_FROM_CPU 11
#define DEPARTURE_FROM_IO1 12
#define DEPARTURE_FROM_IO2 13
#define DEPARTURE_FROM_SWAPOUT 14

#define NUM_OF_SIMULATION 100
#define NUM_OF_MACHINES 20
#define PROBABILITY 0.3
#define SIMULATON_END 2000000
#define DEFAULT    123456789L  /* initial seed, use 0 < DEFAULT < MODULUS   */
#define NUM_OF_RUN 50
#define QUEUE_NUM 5
#define Z 5

#define S2 0.210
#define S3 0.003
#define S4_IO1 0.040
#define S4_IO2 0.180
#define CPU_QUANTUM 0.030
#define MPD 100
#define ALPHA 0.8
#define BETA 0.2
#define MU1 0.015
#define MU2 0.075

#define PROB_IO1 0.065
#define PROB_IO2 0.025
#define PROB_PROCEED_COMP 0.01
#define PROB_PROCEED_CPU 0.90
#define PROB_SWAP_OUT 0.6
#define PROB_COMPLETEION 0.4

#define IO1 200
#define IO2 201
#define PROCEED_COMP 202
#define SWAP_OUT 203
#define COMPLETION 204
#define PROCEED_CPU 205

#define N_CYCLE_MIN 1
#define MINIMUM_NUMBER_REG_CYCLES 30
#define FACTOR 1.645
#define PRECISION 0.05

#define M   5
#define N   20
#define D   1.0
#define LI  2.0

typedef struct node* tree;
typedef struct node* list;

int node_number;
tree root;         /* Pointer to the root of the Future Event List implemented with a tree structure */
list queue;
list last;
int job_number;    /* (progressive) Job identification number */

typedef struct {
    int index;
    int type;
    char name[256];
    double service_time;
    double occur_time;
    double departure_time;
    bool processing;
    double start_time;
    double start_time_Active;
    double time_R;
    double time_A;
    //long double time_RT;
    //long double time_AT;
} event_notice;

struct node{
    event_notice event;
    struct node* next;
};

//queue Data Structure

//this link always point to first Link
struct node *head = NULL;
//this link always point to last Link
struct node *last = NULL;
//this link will point to current Link
struct node *current = NULL;

// An array of linked list (LL) node to store a queue entry
struct Queue *q[QUEUE_NUM];

struct QNode
{
    struct node* key;
    struct QNode *next;
};

// The queue, front stores the front node of LinkedList and rear stores ths
// last node of LinkedList
struct Queue
{
    struct QNode *front, *rear;
};

struct Cycle{
    int id;
    double start;
    double end;
};

struct Cycle recordCycle_Time[20];

struct QNode* newNode(struct node* k);
struct Queue *createQueue();
void enqueue(struct node* k, int num);
struct node *dequeue(int num);

void simulate(void);
void initialize(void);
void engine(void);
void arrival(struct node*);
void departure(struct node*);
void rest(struct node* node_event);
void report(void);
void print_results();
double GetArrival(void);
double GetService(void);
void schedule(struct node* newevent) ;

struct node* returnHead();
struct node* event_pop();
struct node* get_new_node();
void destroy_list();
void destroy_queue();
double getQueue_size(int num);
void printTopQueueEvent();
void printQueue();
void printEvent(struct node* node_event);
void printFutureEventList();

void arrivalAtWorkStation(struct node* node_event);
void arrivalAtReserve(struct node* node_event);
void arrivalAtSwapIN(struct node* node_event);
void arrivalAtCPU(struct node* node_event);
void arrivalAtIO1(struct node* node_event);
void arrivalAtIO2(struct node* node_event);
void arrivalAtSwapOut(struct node* node_event);

void departureFromWorkStation(struct node* node_event);
void departureFromReserve(struct node* node_event);
void departureFromSwapIN(struct node* node_event);
void departureFromCPU(struct node* node_event);
void departureFromIO1(struct node* node_event);
void departureFromIO2(struct node* node_event);
void departureFromSwapOut(struct node* node_event);

bool regPoint(int event_type);
void collectRegStatistics();
void doMVA();
void printMVAResult(char name[], double val[M][N+1]);
void initializeMVA();
void printValue(char* name, double Area_w, int nQueue, int nDeparture);
void printSystemState();
void printMeanvalues();
double getAverageResponse();
double getAverageActive();
void printRegenrationvalues();
void printAllDepartures();
void printSystemValues();
void setConfidenceParameters();
void validateActive();
void validateResponse();
void printValidationResults();
void printCycleInformation(int processIndex);
