#include<stdio.h>
#include <stdlib.h> // For the memory allocation
#include<math.h> // fabs (float absolute) function
#include <stddef.h> // typedef macros (to have parameters in functions as the kernels)

// Parameters for the integration convergence
#define TOL 1.0e-5
#define JMAX 20

//  f: R -> R 
//  function to integraete
float f( float x){
    return 1. / (x * x);
}


// Implement the kernel change and the new variabels and limits automatically
// f(x)
typedef float (*IntegrandFunc)(float x, void *user_data);
/* x -> t: transforms (t_left, t_right) and evaluates g(t)  */
typedef float (*SubstitutionKernel)(float t, IntegrandFunc func, void *user_data);

// We change the integrate strcut we defined previuosly
struct integrate{
    int n_step;                // iteration step
    double x_left, x_right;    // interval borders x_left < x_right
    double sum;                // accumulated value of the integral

    // 
    IntegrandFunc integrand;   // float (*integrand)(float);
    SubstitutionKernel kernel; // Variable transformation kernel (optional, NULL for standard) 
    void *user_data
};

// We define different kernels 
float kernel_identity(float t, IntegrandFunc func, void *user_data) { 
    return func(t);  // x = t, dx = dt
}

float eval_step(struct integrate *t_struct, float x_eval) {  // Necesario para no usar wrappers y tener kernels genericos
    if (t_struct->kernel != NULL) {
        return t_struct->kernel(x_eval, t_struct->integrand, t_struct->user_data);
    }
    return t_struct->integrand(x_eval, t_struct->user_data);
}

/* Routine implementing the extended midpoint rule which has a structure as argument */
void midpoint( struct integrate *t){
    float x, len, s, del, ddel;  // x: Evaluating point
                            // len: Number of grid points at current level
                            // s: partial sum of the current level
                            // del: step size at current level
    int it, j;              // it: Grid partition (Powers of 2)

    (*t).n_step++;  // alternative: t -> n_step++
                    // Increases the level for the partition grid when function is called

    if ( t -> n_step == 1){ // First evaluation: extrema values
        (*t).sum = ((*t).x_right - (*t).x_left) * eval_step( t, 0.5 * ( (*t).x_right + (*t).x_left)  );
//     (*t).sum = ((*t).x_right - (*t).x_left) * (*t).integrand( 0.5 * ( (*t).x_right + (*t).x_left)  ) ;
    } else {
        for(it = 1, j = 1; j < (*t).n_step - 1 ; j++)   // Loop to build  the grid by increasing j
            it *= 3;                                   // Each time all it is calles, it increases a power of 2 with a bit-wise operation
        len = it;
        del = ((*t).x_right - (*t).x_left) / (3.0 * len);
        ddel = del + del;
        x = (*t).x_left + 0.5 * del;                    // Starts the grid by shiftinf it from the x_left  point
        s = 0.0;
        for(j = 0; j < it; j++) {     // Evaluation of the funtion at the points of the current grid (Level)
            s += eval_step(t, x);
            x += ddel;
            s += eval_step(t, x);
            x += del;
        }
        (*t).sum = (*t).sum / 3.0 + del * s;          // adding values to our structure
    }

}

/* Routine with a tolerance to stop the calculations */
void midpoint_tol(struct integrate *t){
    int j;
    float old = 0.0;    // previous step

    for(j=0; j <= JMAX; j++){
        midpoint( &(*t) );      // Performing 1 level of integration in the grid
        if( j >= 4 ){           // Forcing at least five iterations
            if( fabs( (*t).sum - old ) < TOL*fabs(old) || ((*t).sum == 0.0 && old == 0.0)){ // Checking convergence
                return ;        // Void functions cannot return any value, thus this is how we exit the function
            }
        }
        old = (*t).sum;
    }
    printf("Too many steps in routine trapzoid_tol");
}


// infinite integrals
float kernel_midinf(float t, IntegrandFunc func) {
    return func(1. / t) / (t * t);    // f(1/t) / t^2
}

float integrate_infinite_upper(const struct integrate *original) {
    struct integrate new_t = *original;

    new_t.n_step = 0;      // Reset integral; Por qué no funciona new_t->n_step ???
    new_t.sum = 0.0;
    
    /* Modify integration limits for x = 1/t */
    new_t.x_left = 0.0;          /* t at infinity = 0 */
    new_t.x_right = 1.0 / (*original).x_left;     /* t at lower bound = 1/a */
    
    /* Assign function pointer and substitution kernel */
    new_t.integrand =  (*original).integrand;     ;
    new_t.kernel = kernel_midinf;

    midpoint_tol( &new_t ); // Actual integration

    printf("--- Execution Results ---\n");
    printf("Transformed Integration Result : %.6f\n", new_t.sum);
    printf("Transformed Local Bounds       : t_left = %.2f, t_right = %.2f\n\n", new_t.x_left, new_t.x_right);

    return new_t.sum;
}

float integrate_infinite_lower(const struct integrate *original) {
    struct integrate new_t = *original;

    new_t.n_step = 0;      
    new_t.sum = 0.0;
    
    /* Modify integration limits for x = 1/t */
    new_t.x_right = 0.0;          /* t at -infinity = 0 */
    new_t.x_left = 1.0 / (*original).x_right;     /* t at lower bound = 1/a */
    
    /* Assign function pointer and substitution kernel */
    new_t.integrand =  (*original).integrand;     ;
    new_t.kernel = kernel_midinf;

    midpoint_tol( &new_t ); // Actual integration

    printf("--- Execution Results ---\n");
    printf("Transformed Integration Result : %.6f\n", new_t.sum);
    printf("Transformed Local Bounds       : t_left = %.2f, t_right = %.2f\n\n", new_t.x_left, new_t.x_right);

    return new_t.sum;
}

float integrate_infinite_both(const struct integrate *original) {
    struct integrate new_t = *original;

    new_t.n_step = 0;      
    new_t.sum = 0.0;
    
    /* Modify integration limits for x = 1/t */
    new_t.x_right = 1.0 / (*original).x_left;          /* t at -infinity = 0 */
    new_t.x_left = 1.0 / (*original).x_right;     /* t at lower bound = 1/a */
    
    /* Assign function pointer and substitution kernel */
    new_t.integrand =  (*original).integrand;     ;
    new_t.kernel = kernel_midinf;

    midpoint_tol( &new_t ); // Actual integration

    printf("--- Execution Results ---\n");
    printf("Transformed Integration Result : %.6f\n", new_t.sum);
    printf("Transformed Local Bounds       : t_left = %.2f, t_right = %.2f\n\n", new_t.x_left, new_t.x_right);

    return new_t.sum;
}


// Sqrt singularities Eca. 4.4.5 y 4.4.6
/* Struct to hold kernel-specific parameters */
typedef struct {
    float a;           /* singularity boundary */
    void *user_params; /* Outer user parameters for f(x) if any */
} SingularParams;

/* Substitution kernel for lower square-root singularity */
float kernel_sql(float t, IntegrandFunc func, void *user_data) {
    SingularParams *p = (SingularParams *)user_data;
    float x = p->a + (t * t);
    /* Returns 2 * t * f(a + t^2) */
    return 2.0f * t * func(x, p->user_params);
}
float kernel_squ(float t, IntegrandFunc func, void *user_data) {
    SingularParams *p = (SingularParams *)user_data;
    float x = p->a - (t * t);
    /* Returns 2 * t * f(a + t^2) */
    return 2.0f * t * func(x, p->user_params);
}

float integrate_sql(const struct integrate *original) {
    struct integrate new_t = *original;

    new_t.n_step = 0;      // Reset integral; Por qué no funciona new_t->n_step ???
    new_t.sum = 0.0;
    
    /* Modify integration limits for x = 1/t */
    new_t.x_left = 0.0;          /* t at infinity = 0 */
    new_t.x_right = sqrtf((*original).x_right - (*original).x_left);    
    new_t.integrand =  (*original).integrand;     
    
    /* Falaa un if para seleccionar a < b o al reves */
    SingularParams kernel_params = {
        .a = (float)original->x_left,
        .user_params = original->user_data
    };
    new_t.kernel = kernel_sql;



    

    midpoint_tol( &new_t ); // Actual integration

    printf("--- Execution Results ---\n");
    printf("Transformed Integration Result : %.6f\n", new_t.sum);
    printf("Transformed Local Bounds       : t_left = %.2f, t_right = %.2f\n\n", new_t.x_left, new_t.x_right);

    return new_t.sum;
}

float integrate_infinite_lower(const struct integrate *original) {
    struct integrate new_t = *original;

    new_t.n_step = 0;      
    new_t.sum = 0.0;
    
    /* Modify integration limits for x = 1/t */
    new_t.x_right = 0.0;          /* t at -infinity = 0 */
    new_t.x_left = 1.0 / (*original).x_right;     /* t at lower bound = 1/a */
    
    /* Assign function pointer and substitution kernel */
    new_t.integrand =  (*original).integrand;     ;
    new_t.kernel = kernel_midinf;

    midpoint_tol( &new_t ); // Actual integration

    printf("--- Execution Results ---\n");
    printf("Transformed Integration Result : %.6f\n", new_t.sum);
    printf("Transformed Local Bounds       : t_left = %.2f, t_right = %.2f\n\n", new_t.x_left, new_t.x_right);

    return new_t.sum;
}



int main(){
    struct  integrate t = {
        .n_step = 0,
        .x_left = 5.0,        // Lower bound a = 1.0 
        .x_right = 1e30,      // Abstract representation of upper bound 
        .sum = 0.0,
        .integrand = f,
        .kernel = NULL
    };

    float result; 

    printf("--- Before Integ ration ---\n");
    printf("Original Bounds : x_left = %.2f, x_right = %.2e\n", t.x_left, t.x_right);
    printf("Original Steps  : %d\n\n", t.n_step);

    result = integrate_infinite_upper(&t);

    t.x_right = -5.0;
    result += integrate_infinite_lower(&t);


    /*printf("--- After Integration (Original Preserved) ---\n");
    printf("Original Bounds : x_left = %.2f, x_right = %.2e\n", t.x_left, t.x_right);
    printf("Original Steps  : %d\n", t.n_step);*/
    return 0;
}



/*
int main(){
    int i;   // Iterative variable

    // Trepezoid method is roubust specilally for not so smooth integrands
    struct integrate t;  

    t.x_left = 0.0;
    t.x_right = 1.0;
    t.integrand = f;

// Integration with a predifined number of iterations. Not convergence guaranteed
    printf("Fixed iterations:\n");

    for(i = 1; i <= 5; i++ ){
        printf("\tIteration n = %d -> %.10f\n", t.n_step, t.sum );
        midpoint( &t );
    }

// We can reboot the pointer by forcing the sum and step to zero
// Integration with a given toleance a a maximum number of steps

    t.n_step = 0;
    t.sum = 0;

    printf("\nTOL = %.2e, JMAX = %d\n", TOL, JMAX);
    midpoint_tol( &t );
    printf("\tIteration n = %d -> %.10f\n", t.n_step, t.sum );

// Integration with Simpson
// In general it is most efficient since we need less evaluations. Useful for
// functions with  a finith 4th derivative.

    t.n_step = 0;
    t.sum = 0;

    printf("\nSimpson\nTOL = %.2e, JMAX = %d\n", TOL, JMAX);
    simp_midpoint_tol( &t );
    printf("\tIteration n = %d -> %.10f\n", t.n_step, t.sum );

      
    t.n_step = 0;
    t.sum = 0;
    romberg_tol( &t );
    printf("\tIteration n = %d -> %.10f\n", t.n_step, t.sum );
    
    return 0;
}
*/