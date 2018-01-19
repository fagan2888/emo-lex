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

#include <algorithm>
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
     printf("read_local_parameters\n");
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
     printf("read_ini\n");
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
     printf("select_ind\n");
     result = select_ind(alpha, result_identities,
                         PISA_identities, dimension); /* changedddd */

     if (result != 0)
     {
          log_to_file(log_file, __FILE__, __LINE__, "selection failed");
          return (1);
     }
     
     /**********| addition for LEX end |*******/
     //printf("write_sel\n");
     result = write_sel(PISA_identities);     /* write sel file */
     if(result != 0)
     {
          log_to_file(log_file, __FILE__, __LINE__, "failed write_sel()");
          return(1);
     }
     
     free(PISA_identities);
     free(result_identities);
     //printf("write_arc\n");
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
     // get current archive ids
     std::vector<int> current_arc; 
     current_arc.push_back(get_first());
     for (unsigned i = 1; i< get_size(); ++i)
         current_arc.push_back(get_next(current_arc[i-1]));

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

     //printf("read_var\n");
     result = read_var(offspring_identities);
     if (result == 1) /* if some file reading error occurs, return 2 */
          return (2);

     /**********| added for LEX |**************/

     // clear old global population 
     for (unsigned i = 0; i < current_arc.size(); -++i)
         remove_individual(current_arc[i]);
    
     
     result = select_ind(lambda, offspring_identities,
                         parent_identities, dimension);
     
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

     // clear global population
     //printf("write_arc\n");
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
     result = fscanf(fp, "%e", &seed); /* fscanf() returns EOF if
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
     
     //printf("select_ind...\n");
     //printf("size:%i\n",size);
     assert(dimension >= 0);
     // if continuous objectives, calculate epsilon
     double * epsilon = (double *) malloc(dimension * sizeof(double));
     if (eplex){
         //printf("calculating epsilon\n"); 
         calculate_epsilon(population, size, dimension, epsilon);
     }
     else{
         //printf("setting epsilon to zero\n"); 
         for (int i =0; i<dimension; ++i)
             epsilon[i] = 0;
     }
     
     //printf("%i: population semantics:\n---\n",__LINE__);
     //for (int i = 0; i < get_size(); ++i)
     //{
     //    printf("ind %i:\t",population[i]);
     //    for (int j = 0; j < dimension; ++j)
     //       printf("%e,",get_objective_value(population[i],j));    
     //    printf("\n");
     //}
     //printf("---\nmin:\t");
     //std::vector<int> pop(mu);
     //std::iota(pop.begin(),pop.end(),0);
     //for (int i = 0; i<dimension; ++i)
     //    printf("%e,",min_obj(pop,i));
     //
     //printf("\n---\nep:\t");
     //for (int i = 0; i<dimension; ++i)
     //    printf("%e,",epsilon[i]);
     //printf("\n---\n");
     //printf("starting selection...\n");
     
     /* choose mu individuals by lexicase selection */
     #pragma omp parallel for   
     for(int i = 0; i < mu; i++)
     {
         //printf("selection %i\n",i);
         // cases, i.e. objectives, to consider each selection event
         int *cases = (int *) malloc(dimension * sizeof(int));
         for (int i = 0; i<dimension; ++i)
            cases[i] = i;                
         // shuffle objectives
         shuffle(cases, dimension);
         int pos = lex_choose(population, size, cases, dimension, epsilon);
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
        //printf("calc epsilon for dimension %i...\n",i);
        std::vector<double> objs;
        for (int j = 0; j<popsize; ++j)
            objs.push_back( get_objective_value(j,i) );
        //printf("epsilon[%i] = %e..\n",i,mad(objs));
        epsilon[i] = mad(objs);
    }

}

/// calculate median
double median(const std::vector<double>& v) 
{
    // instantiate a vector
    std::vector<double> x(v.size());
    x.assign(v.data(),v.data()+v.size());
    // middle element
    size_t n = x.size()/2;
    // sort nth element of array
    std::nth_element(x.begin(),x.begin()+n,x.end());
    // if evenly sized, return average of middle two elements
    if (x.size() % 2 == 0) {
        nth_element(x.begin(),x.begin()+n-1,x.end());
        return (x[n] + x[n-1]) / 2;
    }
    // otherwise return middle element
    else
        return x[n];
}

/// median absolute deviation
double mad(const std::vector<double>& x) 
{
    // returns median absolute deviation (MAD)
    // get median of x
    double x_median = median(x);
    //printf("x_median: %e\n",x_median);
    //calculate absolute deviation from median
    std::vector<double> dev;
    for (int i =0; i < x.size(); ++i)
        dev.push_back(fabs(x[i] - x_median));
    //printf("deviations: ");
    for (int i =0; i<dev.size(); ++i)
        //printf("%e ",dev[i]);
    //printf("\n");
    // return median of the absolute deviation
    return median(dev);
}
/* Generate a random integer. */
int irand(int range)
{
     int j;
     j = (int) ((double)range * (double) rand() / (RAND_MAX+1.0));
     return (j);
}

double min_obj(const std::vector<int>& p, int o)
{
    // returns minimum objective value among individuals p for objective o.
    double mino = get_objective_value(p[0],o);
    for (unsigned i = 1; i<p.size(); ++i)
    {
        if (get_objective_value(p[i],o) < mino)
            mino = get_objective_value(p[i],o);
    }
    return mino;
}

/* Chooses one individual uniformly among those with the lowest counter.
   Returns ID of chosen individual.
   Returns -1 of choosing failed for any reason. */
int lex_choose(int * pop, int starting_pool_size, int *cases, int dimension, double *epsilon) 
{  
     std::vector<int> pool;
     for (int i =0;i<starting_pool_size; ++i)
         pool.push_back(pop[i]);
     //printf("%i: initial pool: \t",__LINE__);
     //for (int i = 0; i<starting_pool_size; ++i)
     //    printf("%i,",pool[i]);
     //printf("\n");
     //printf("%i: cases: \t",__LINE__);
     //for (int i = 0; i<dimension; ++i)
     //    printf("%i,",cases[i]);
     //printf("\n");  
     
     bool pass = true;
     int c = 0; 
     int pool_size = 0;
     double p = 10000; //precision
     while (pass) 
     { 
         
         std::vector<int> sel_pool;

         //printf("pool: ");
         //for (unsigned i =0; i<pool.size();++i) 
         //    printf("%i(%e) ",pool[i],get_objective_value(pool[i],cases[c]));
         //printf("\n");
         double threshold = min_obj(pool,cases[c]) + epsilon[cases[c]];
         //printf("case %i threshold: %e + %e = %e\n",cases[c],min_obj(pool,cases[c]),epsilon[cases[c]],threshold);
         //printf("sel_pool: ");
         for (int i = 0; i < pool.size(); ++i)
         {
             double f = get_objective_value(pool[i],cases[c]);
                          
             if (f <= threshold)
             {
                 // add pool[i] to winners
                 //printf("%i (%e), ",pool[i],f);
                 sel_pool.push_back(pool[i]);
             }
         }
         //printf("\n");
         //printf("%i: size of sel_pool:%i\n",__LINE__,sizeof(sel_pool)/sizeof(sel_pool[0]));
         //printf("%i: pool_size: %i\n",__LINE__,pool_size);
         ++c;   // increment case            
         pass = (sel_pool.size() > 1 && c < dimension); // keep going if needed
         assert(sel_pool.size()>0); 
         pool = sel_pool;
     }
     //printf("final pool size: %i\n",pool.size());
     int pick = irand(pool.size());
     int selection = pool[pick];
     //printf("selecting %i (pick = %i)\n",selection, pick);  
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
