#include<stdio.h>
#include<math.h> // abs function

#define TOL 1.0e-5
#define JMAX 20

//  f: R -> R 
float f( float x){
    return x * x;
}

/* C structures are like list's whose element can be of diferent data-tyepes.
 They cannot contain functions as elements, but they can be allocated with pointers. */
struct integrate
{
    int n_step;                // iteration step
    double x_left, x_right;    // interval borders x_left < x_right
    double sum;                // accumulated value of the integral

    float (*integrand)(float);
};

/* Routine implementing the extended trapezoidal rule which has a structure as argument */
void trapzoid( struct integrate *t){
    float x, len, s, del;  // x: Evaluating point
                            // len: Number of grid points at current level
                            // s: partial sum of the current level
                            // del: step size at current level
    int it, j;              // it: Grid partition (Powers of 2)

    (*t).n_step++;  // alternative: t -> n_step++
                    // Increases the level for the partition grid when function is called

    if ( t -> n_step == 1){ // First evaluation: extrema values
        ((*t).sum = (0.5) * ((*t).x_right - (*t).x_left) * ((*t).integrand( (*t).x_right ) + (*t).integrand( (*t).x_left) ));
    } else {
        for(it = 1, j = 1; j < (*t).n_step - 1 ; j++)   // Loop to build  the grid by increasing j
            it <<= 1;                                   // Each time all it is calles, it increases a power of 2 with a bit-wise operation
        len = it;
        del = ((*t).x_right - (*t).x_left) / len;
        x = (*t).x_left + 0.5 * del;                    // Starts the grid by shiftinf it from the x_left  point
        for(s = 0.0, j = 0; j < it; j++, x += del)      // Evaluation of the funtion at the points of the current grid (Level)
            s += (*t).integrand(x); 
        (*t).sum = 0.5 * ((*t).sum + del * s);          // adding values to our structure
    }

}

/* Routine with a tolerance to stop the calculations */
void trapzoid_tol(struct integrate *t){
    int j;
    float old = 0.0;    // previous step

    for(j=0; j <= JMAX; j++){
        trapzoid( &(*t) );      // Performing 1 level of integration in the grid
        if( j >= 4 ){           // Forcing at least five iterations
            if( fabs( (*t).sum - old ) < TOL*fabs(old) || ((*t).sum == 0.0 && old == 0.0)){ // Checking convergence
                return ;        // Void functions cannot return any value, thus this is how we exit the function
            }
        }
        old = (*t).sum;
    }
    printf("Too many steps in routine trapzoid_tol");
}

/* Simpson with a tolerance to stop the calculations */
void simpson_tol(struct integrate *t){

    float old_t = 0.0;    // previous step of trapeziod
    float old_s = 0.0;    // previous step of Simpson 
    float j, simp;        //  iterator, Simpson 

    for(j=0; j <= JMAX; j++){
        trapzoid( &(*t) );      // Performing 1 level of integration in the grid
        simp =  (4. * (*t).sum - old_t) / 3.;
        if( j >= 4 ){           // Forcing at least five iterations
            if( fabs( simp - old_s ) < TOL*fabs(old_s) || (simp == 0.0 && old_s == 0.0)){ // Checking convergence
                (*t).sum = simp;
                return ;        // Void functions cannot return any value, thus this is how we exit the function
            }
        }
        old_s = simp;
        old_t = (*t).sum;
    }
    printf("Too many steps in routine simpson_tol");
}

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
        trapzoid( &t );
    }

// We can reboot the pointer by forcing the sum and step to zero
// Integration with a given toleance a a maximum number of steps

    t.n_step = 0;
    t.sum = 0;

    printf("\nTOL = %.2e, JMAX = %d\n", TOL, JMAX);
    trapzoid_tol( &t );
    printf("\tIteration n = %d -> %.10f\n", t.n_step, t.sum );

// Integration with Simpson
// In general it is most efficient since we need less evaluations. Useful for
// functions with  a finith 4th derivative.

    t.n_step = 0;
    t.sum = 0;

    printf("\nSimpson\nTOL = %.2e, JMAX = %d\n", TOL, JMAX);
    simpson_tol( &t );
    printf("\tIteration n = %d -> %.10f\n", t.n_step, t.sum );

    return 0;
}