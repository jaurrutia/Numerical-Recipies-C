#include<stdio.h>

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

int main(){
    int i;

    struct integrate t;
    t.x_left = 0.0;
    t.x_right = 1.0;
    t.integrand = f;



    for(i = 1; i <= 10; i++ ){
        printf("Iteration n = %d -> %.10f\n", t.n_step, t.sum );
        trapzoid( &t );
    }
    return 0;
}