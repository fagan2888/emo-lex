/*========================================================================
  PISA  (www.tik.ee.ethz.ch/pisa/)
 
  ========================================================================
  Computer Engineering (TIK)
  ETH Zurich
 
  ========================================================================
  PISALIB 
  
  Implementation file.

  The functions in this file provide a framework for the functions you
  need to implement. These functions will be called by other parts of
  PISALib
  
  file: selector_user.c
  author: William La Cava, lacava@upenn.edu

  
  ========================================================================
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "selector.h"
#include "selector_user.h"

/*--------------------| global variable definitions |-------------------*/

/*==== declared in selector_user.h used in other files as well ====*/


char *log_file = "lex_error.log"; /**** Changed for LEX. */

char paramfile[FILE_NAME_LENGTH]; /* file with local parameters */

bool eplex = false;     // whether to use epsilon lexicase (for continuous problems)
/*-------------------------| individual |-------------------------------*/
int set_objective_value(individual *ind, int index, double obj_value)
/* Sets the objective_value of an individual.
   
   post: It is checked if index is in the correct range. If it is the
         objective value for the objective number i in individual ind
         is set to obj_value.
         Return value is 0 if successful and 1 otherwise.
*/        
{
     if(ind == NULL || index < 0 || index >= dimension)
          return (1);
     else
     {
          /**********| added for LEX |**************/
          ind->objective_value[index] = obj_value;
          /**********| addition for LEX end |*******/
          
          return (0);
     }
}


individual *create_individual()
/* Allocates memory for a new individual and initializes values.

   post: returns a pointer to the allocated memory
         returns NULL if allocation failed
*/
{
     /**********| added for LEX |**************/
     double *obj_value;
     /**********| addition for LEX end |*******/

     individual *return_ind;

     return_ind = (individual *) malloc(sizeof(individual));
     if (return_ind == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "selector out of memory");
          return(NULL);
     }

     /**********| added for LEX |**************/

     obj_value = (double *) malloc(sizeof(double) * dimension);
     if (obj_value == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "selector out of memory");
          return(NULL);
     }

     return_ind->objective_value = obj_value;     

     /**********| addition for LEX end |*******/

     return (return_ind);
}


void free_individual(individual* ind)
/* Frees the memory for one indiviual.

   post: memory for ind is freed
*/
{
     /**********| added for LEX |**************/
     free(ind->objective_value);
     /**********| addition for LEX end |*******/
     
     free(ind);
}


/*-------------------------| statemachine functions |-------------------*/

int state1() 
/* Do what needs to be done in state 1.

   pre: The global variable 'paramfile' contains the name of the
        parameter file specified on the commandline.
        The global variable 'alpha' contains the number of indiviuals
        you need to read using 'read_ini()'.
        The global variable 'mu' contains the number of individuals
        you need to select.
        
   post: optionally read parameter specific for the module
         optionally do some initialization
         read_ini() called
         'mu' parents selected
         undesired individuals deleted from the global population
         write_sel() called
         write_arc() called
         return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file operation failed.
*/
{
     int result; /* stores return values of called functions */
     int *PISA_identities;
     int *result_identities; /* array of identities returned by read_ini */

     /**** Here you could call a function to read the local parameter
           values from the 'paramfile'. */

     /**********| added for LEX |**************/
     result = read_local_parameters();  
     if (result != 0)
     { 
          log_to_file(log_file, __FILE__, __LINE__,
                      "couldn't read local parameters");
          return (1);
     }
     /**********| addition for LEX end |*******/

     
     result_identities = (int *) malloc(alpha * sizeof (int));
     if (result_identities == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "selector out of memory");
          return (1);
     }
     
     result = read_ini(result_identities);   /* read ini file */
     if (result == 1)
          return (2); /* reading ini file failed */
    

     PISA_identities = (int *) malloc(mu * sizeof (int));
     if (result_identities == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "selector out of memory");
          return (1);
     }
     
     /**********| added for LEX |**************/

     result = select_ind(alpha, result_identities,
                         PISA_identities, dimension); /* changedddd */

     if (result != 0)
     {
          log_to_file(log_file, __FILE__, __LINE__, "selection failed");
          return (1);
     }
     
     /**********| addition for LEX end |*******/

     result = write_sel(PISA_identities);     /* write sel file */
     if(result != 0)
     {
          log_to_file(log_file, __FILE__, __LINE__, "failed write_sel()");
          return(1);
     }
     
     free(PISA_identities);
     free(result_identities);

     result = write_arc();     /* write arc file
                                  (individuals in global population) */
     if(result != 0)
     {
          log_to_file(log_file, __FILE__, __LINE__, "failed write_arc()");
          return(1);
     }
  
     return (0);   
}  


int state3()
/* Do what needs to be done in state 3.

   pre: The global variable 'lambda' contains the number of indiviuals
        you need to read using 'read_var()'.
        The global variable 'mu' contains the number of individuals
        you need to select and then write to the 'sel' file.
        
   post: read_var() called
         mu parents selected
         undesired individuals deleted from the global population
         write_sel() called
         write_arc() called
         return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     int result; /* stores return values of called functions */
     int *offspring_identities; /* array with IDs filled by read_var() */
     int *parent_identities; /* array with IDs of parents */

     offspring_identities = (int *) malloc(lambda * sizeof(int)); 
     if (offspring_identities == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "selector out of memory");
          return (1);
     }

     parent_identities = (int *) malloc(mu * sizeof(int)); 
     if (parent_identities == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "variator out of memory");
          return (1);
     }
     
     result = read_var(offspring_identities);
     if (result == 1) /* if some file reading error occurs, return 2 */
          return (2);

     /**********| added for LEX |**************/

     //result = select_ind(lambda, offspring_identities,
     //                    parent_identities, dimension);
     // for now, just set parent_identities to offspring, like normal lexicase
     assert(mu == lambda);
     //parent_identities = offspring_identities;
     memcpy(parent_identities,offspring_identities,lambda*sizeof(int));
     if (result != 0)
     {
          log_to_file(log_file, __FILE__, __LINE__, "selection failed");
          return (1);
     }
         
     /**********| addition for LEX end |*******/

     result = write_sel(parent_identities);
     if(result != 0)
     {
          log_to_file(log_file, __FILE__, __LINE__, "failed write_sel()");
          return(1);
     }

     free(parent_identities);
     free(offspring_identities);

     result = write_arc();
     
     if(result != 0)
     {
          log_to_file(log_file, __FILE__, __LINE__, "failed write_arc()");
          return (1);
     }

     return (0);   
}  


int state5()
/* Do what needs to be done in state 5.

   pre: state 5 means the variator has just terminated

   post: you probably don't need to do anything, just return 0
         return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.

*/
{
     return(0); 
}


int state6()
/* Do what needs to be done in state 6.

   pre: state 6 means that the selector has to terminate.

   post: free all memory
         return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     int current_id;
    
     current_id = get_first();
     while(current_id != -1) /* freeing memory */
     {
          remove_individual(current_id);
          current_id = get_next(current_id);
     }
     return (0);
}


int state9()
/* Do what needs to be done in state 9.

   pre: State 9 means that the variator is reset and ready to start
        again in state 0.

   post: You probably don't have to do anything, just return 0.
         return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     return (0);
}


int state10()
/* Do what needs to be done in state 10.

   pre: State 10 means that the selector need to reset and get ready
        to start again in state 1.

   post: Get ready to start again in state 1.
         return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
   /* freeing memory is done in selector.c */
   return (0);
}



int is_finished() 
/* Tests if ending criterion of your algorithm applies.

   post: return value == 1 if optimization should stop
         return value == 0 if optimization should continue

   remark: Normally the variator decides when to terminate. So
           just return 0.
*/
{
     return (0);
}


/**********| added for LEX |**************/


int read_local_parameters()
{
     
     FILE *fp;

     int result;
     char str[CFG_NAME_LENGTH];
     int seed, tmpep;
     
     /* reading parameter file with parameters for selection */
     fp = fopen(paramfile, "r"); 
     assert(fp != NULL);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "seed") == 0);
     result = fscanf(fp, "%d", &seed); /* fscanf() returns EOF if
                                          reading fails. */
     assert(result != EOF); /* no EOF, 'seed' correctly read */
     
     srand(seed); /* seeding random number generator */

     fscanf(fp, "%s", str);
     assert(strcmp(str, "epsilon") == 0);
     result = fscanf(fp, "%i", &tmpep); /* fscanf() returns EOF if
                                          reading fails. */
     assert(result != EOF); /* no EOF, 'seed' correctly read */
     
     eplex = tmpep; 
     //printf("eplex: %i\n",eplex);
     fclose(fp);
  
     /* do some other initialization steps... */
     return (0);
}


/* Implements LEX. Takes size individual from variation, updates global
   population and selects mu new individual for variation. */
int select_ind(int size, int *population, int *parents,
               int dimension)
{
     /* @params size: number of offspring produced by variation
      * @params population: population indices to choose from
      * @returns parents: indices of selected individuals
      * @params dimension: number of objectives
      */
     
     printf("select_ind...\n");
     printf("size:%i\n",size);
     printf("population size:%i\n",sizeof(population)/sizeof(population[0]));
     assert(dimension >= 0);
          // if continuous objectives, calculate epsilon
     double * epsilon = (double *) malloc(dimension * sizeof(double));
     if (eplex){
         printf("calculating epsilon\n"); 
         calculate_epsilon(population, size, dimension, epsilon);
     }
     else{
         printf("setting epsilon to zero\n"); 
         for (int i =0; i<dimension; ++i)
             epsilon[i] = 0;
     }
     
     printf("starting selection...\n");
     
     /* choose mu individuals by lexicase selection */
     //#pragma omp parallel for 
     for(int i = 0; i < mu; i++)
     {
         printf("selection %i\n",i);
         // cases, i.e. objectives, to consider each selection event
         int *cases = (int *) malloc(dimension * sizeof(int));
         for (int i = 0; i<dimension; ++i)
            cases[i] = i;                
         // shuffle objectives
         shuffle(cases, dimension);
         int pos = lex_choose(cases,epsilon,dimension);
        // if (pos == -1) /* Choosing failed. */
        //      return (1);
         parents[i] = pos;
         free(cases);
     }
     return (0);
}
void calculate_epsilon(int *ids, int popsize, int dimension, double *epsilon)
{
    // calculate median absolute deviation (MAD) of each objective across the population
    for (int i  = 0; i < dimension; ++i)
    {
        printf("calc epsilon for dimension %i...\n",i);
        double * objs = (double *) malloc(popsize * sizeof(double));
        for (int j = 0; j<popsize; ++j)
            objs[j] = get_objective_value(j,i);
        printf("mad..\n");
        epsilon[i] = mad(objs,popsize);
        free(objs);
    }

}
int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}
/// calculate median
double median(double *v, int size) 
{

    // instantiate a vector
    printf("\tinstantiate a vector\n");
    printf("\tsize:%i\n",size);
    double * x = (double *) malloc(size * sizeof(double));
    if (x == NULL)
    {
        printf("median out of memory\n");
        log_to_file(log_file, __FILE__, __LINE__, "median out of memory");
        return (-1);
    }
	// copy v to x 
    printf("\tcopy v to x\n");
    memcpy(x, v, sizeof(double)*size);
    printf("\tqsort\n");
 	qsort (x, size, sizeof(double), compare);   
    double answer; 
    // if evenly sized, return average of middle two elements
    if (size % 2 == 0) 
        answer= (x[size/2] + x[size/2-1]) / 2;
    // otherwise return middle element
    else
        answer= x[size/2];
    free(x);
    return answer;
}

/// median absolute deviation
double mad(double * x, int size) 
{
    // returns median absolute deviation (MAD)
    // get median of x
    printf("get median of x\n");
    double x_median = median(x,size);
    //calculate absolute deviation from median
    // instantiate a vector
    double * dev = (double *) malloc(size * sizeof(double));
    if (dev == NULL)
    {
        log_to_file(log_file, __FILE__, __LINE__, "median out of memory");
        return (-1);
    }
    printf("calculate deviation\n");    
    for (int i =0; i < size; ++i)
        dev[i] = abs(x[i] - x_median);
    // return median of the absolute deviation
    double answer = median(dev,size);
    free(dev);
    return answer;
}
/* Generate a random integer. */
int irand(int range)
{
     int j;
     j = (int) ((double)range * (double) rand() / (RAND_MAX+1.0));
     return (j);
}


/* Chooses one individual uniformly among those with the lowest counter.
   Returns ID of chosen individual.
   Returns -1 of choosing failed for any reason. */
int lex_choose(int *cases, double *epsilon, int dimension) 
{
     printf("lex_choose\n");
     int starting_pool_size = get_size();
     
     printf("%i: population semantics:\n---\n",__LINE__);
     for (int i = 0; i < starting_pool_size; ++i)
     {
         printf("ind %i:\t",i);
         for (int j = 0; j < dimension; ++j)
            printf("%i,",get_objective_value(i,j));    
         printf("\n");
     }
     printf("ep:\t");
     for (int i = 0; i<dimension; ++i)
         printf("%d,",epsilon[i]);
     printf("\n---\n");
     
     int *pool = (int *) malloc(starting_pool_size * sizeof(int));
     for (int i =0;i<starting_pool_size; ++i)
         pool[i] = i;
     printf("%i: initial pool: \t",__LINE__);
     for (int i = 0; i<starting_pool_size; ++i)
         printf("%i,",pool[i]);
     printf("\n");
     printf("%i: cases: \t",__LINE__);
     for (int i = 0; i<dimension; ++i)
         printf("%i,",cases[i]);
     printf("\n");  
     
     bool pass = true;
     int c = 0; 
     int pool_size = 0;
     while (pass) 
     { 
         double best_val = 0;
         int *sel_pool = NULL;
         //int *tmp_pool = NULL;
         printf("cases[c=%i]:%i\n",c,cases[c]);

         printf("size of sel_pool:%i\n",sizeof(sel_pool)/sizeof(sel_pool[0]));
         for (int i = 0; i < starting_pool_size; ++i)
         {
             if (i==0 || get_objective_value(pool[i],cases[c]) < best_val + epsilon[cases[c]])
             {
                 best_val = get_objective_value(pool[i],cases[c]);  // reset best val 
                 printf("\t%i: new best_val (%i,%i):%d\n",pool[i],cases[c],best_val);
                 // restart the winner pool
                 printf("\t%i: free sel_pool\n",__LINE__);
                 free(sel_pool);
                 //sel_pool = NULL;
                 //printf("malloc sel_pool\n");
                 sel_pool = (int *) malloc(sizeof(int));
                 if (sel_pool == NULL)
                 {
                    printf("out of memory line %i\n",__LINE__);
                    log_to_file(log_file, __FILE__, __LINE__, "selector out of memory");
                    return (-1);
                 }
                 pool_size = 1;
                 // push this individual into the pool
                 //printf("%i: pushing this individual into the pool\n",__LINE__);
                 //printf("size of tmp_pool:%i\n",sizeof(tmp_pool)/sizeof(tmp_pool[0]));
                 //printf("%i: size of sel_pool:%i\n",__LINE__,sizeof(sel_pool)/sizeof(sel_pool[0]));
                 //printf("sel_pool=tmp_pool\n");         
                 //sel_pool = tmp_pool;
                 printf("\t%i: sel_pool[0] = pool[%i] = %i\n",__LINE__,i,pool[i]);
                 sel_pool[0] = pool[0];
             }
             else if (get_objective_value(pool[i],cases[c]) == best_val + epsilon[cases[c]])
             {
                 // add pool[i] to winners
                 ++pool_size;
                 printf("\tadding pool[%i] to winners (best_val + ep=%d)\n",i,best_val+epsilon[cases[c]]);
                 int * tmp_pool = (int *) malloc(pool_size * sizeof(int));
                 

                 if (tmp_pool == NULL)
                 {
                    log_to_file(log_file, __FILE__, __LINE__, "selector out of memory");
                    return (-1);
                 }
                 for (int j = 0; j<pool_size-1; ++j)
                    tmp_pool[j] = sel_pool[j];
                 tmp_pool[pool_size-1] = pool[i];
                 // push this individual into the pool
                 //printf("\tsize of tmp_pool:%i\n",sizeof(tmp_pool)/sizeof(tmp_pool[0]));
                 free(sel_pool);
                 sel_pool = tmp_pool;
                 printf("\t%i: sel_pool[%i-1] = pool[%i] = %i\n",__LINE__,pool_size,i,pool[i]);
                 //sel_pool[pool_size-1] = pool[i];
             }
         }
         printf("%i: size of sel_pool:%i\n",__LINE__,sizeof(sel_pool)/sizeof(sel_pool[0]));
         printf("%i: pool_size: %i\n",__LINE__,pool_size);
         ++c;   // increment case            
         pass = (pool_size > 1 && c < dimension); // keep going if needed
         printf("%i: resetting pool for next case\n",__LINE__); 
         printf("%i: sel_pool:\t",__LINE__);
         for (int j = 0; j < pool_size; ++j)
             printf("%i,",sel_pool[j]);
         printf("\n");
         free(pool);
         pool = (int *) malloc ( pool_size * sizeof(int));
         memcpy(pool,sel_pool, pool_size * sizeof(int)); 
         printf("%i: new pool:\t",__LINE__);
         for (int j = 0; j < pool_size; ++j)
             printf("%i,",pool[j]);
         printf("\n");
     }
     int pick = irand(pool_size);
     int selection = pool[pick];
     
     free(pool);
     return selection;

     /// FEMO ///////////
    // int *ids_to_choose;
    // int size, min, current_id, pick_id, return_id;
    // 
    // ids_to_choose = (int *) malloc(get_size() * sizeof(int));
    // if (ids_to_choose == NULL)
    // {
    //      log_to_file(log_file, __FILE__, __LINE__, "selector out of memory");
    //      return (-1);
    // }  
    // 
    // current_id = get_first();
    // if(current_id == -1)
    //      return(-1);
    // 
    // min = get_counter(current_id);
    // size = 1;
    // ids_to_choose[0] = current_id;
    // current_id = get_next(current_id);
    // while(current_id != -1)
    // {
    //      if(min > get_counter(current_id))
    //      {
    //           size = 1;
    //           ids_to_choose[0] = current_id;
    //           min = get_counter(current_id);
    //      }
    //      else if(min == get_counter(current_id))
    //      {
    //           size++;
    //           ids_to_choose[size - 1] = current_id;
    //      }
    //      current_id = get_next(current_id);
    // }
    // 
    // pick_id = irand(size);
    // 
    // return_id = ids_to_choose[pick_id];

    // free(ids_to_choose);
    // 
    // return (return_id);
}


double get_objective_value(int id, int index)
{
     individual *temp;
     temp = get_individual(id);
     if(temp == NULL || index < 0 || index >= dimension)
          return(-1);
     return(temp->objective_value[index]);  
}
// shuffle contents of int array.
void shuffle(int *array, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}
/**********| addition for LEX end |*******/
