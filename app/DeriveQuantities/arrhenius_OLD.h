int arrhenius_f(const gsl_vector * x, void *fit_data, gsl_vector * f);
int arrhenius_df(const gsl_vector * x, void *fit_data, gsl_matrix * J);
int arrhenius_fdf(const gsl_vector * x, void *fit_data, gsl_vector * f, gsl_matrix * J);
void fit_(size_t *n, double *t, double *y, double *a, double *b, double *c, double *err);
void print_state(size_t iter, gsl_multifit_fdfsolver * s);
