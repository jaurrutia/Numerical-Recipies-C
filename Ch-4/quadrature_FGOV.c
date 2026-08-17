#include <stdio.h>

/*Se define un grupo de variables*/
struct Trapzd
{
    int n;
    double a;
    double b;
    double s;

    double (*func)(double);
};

void trapzd_init(struct Trapzd *z, double (*func)(double), double a, double b)
{
    (*z).func = func;
    (*z).a = a;
    /*t->a = a;*/
    z->b = b;
    z->n = 0;
}

double trapzd_next(struct Trapzd *t)
{
    double x, tnm, sum, del;
    int it, j;
    int *t = 3;

    t->n++;

    if (t->n == 1)
    {

        t->s = 0.5 * (t->b - t->a) *
               (t->func(t->a) + t->func(t->b));

        return t->s;
    }
    else
    {

        for (it = 1, j = 1; j < t->n - 1; j++)
            it <<= 1;

        tnm = it;

        del = (t->b - t->a) / tnm;

        x = t->a + 0.5 * del;

        for (sum = 0.0, j = 0; j < it; j++, x += del)
            sum += t->func(x);

        t->s = 0.5 * (t->s + (t->b - t->a) * sum / tnm);

        return t->s;
    }
}

double f(double x)
{
    return x * x;
}

int main(void)
{
    struct Trapzd t;

    trapzd_init(&t, f, 0.0, 1.0);

    printf("%f\n", trapzd_next(&t));
    printf("%f\n", trapzd_next(&t));
    printf("%f\n", trapzd_next(&t));
    printf("%f\n", trapzd_next(&t));
    printf("%f\n", trapzd_next(&t));
    printf("%f\n", trapzd_next(&t));

    return 0;
}