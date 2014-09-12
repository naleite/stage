//
//  fonction.h
//  Simulation
//
//  Created by NALEITE on 14-5-22.
//
//

 extern void generer_image_XY(int Nw, double Imoy_w, double DOP_w, int Nw_, double Imoy_w_,double DOP_w_, int ordre, int reals, double **image_X, double **image_Y);


 extern void generer_3image(int taille, int reals,double **image_X,double **image_Y, double **image_A_tot, double **image_A_ICEO, double **image_A_delta);

extern double d_min( double a, double b);

extern double d_min_1d(double *a, int dim);

extern void calcul_histogram(double *vecteur, int size, double min, double max, int pt_nb_bin, double *histo);

extern void calcul_histogram(double *vecteur, int size, double min, double max, int pt_nb_bin, double *histo);

extern void  calculer_Pd_Pfa(int nb_niv, double *tab_P, double *tab_Q, double *hist_P, double *hist_Q, int nb_pts);

extern void generer_image_XY(int Nw, double Imoy_w, double DOP_w, int Nw_, double Imoy_w_,double DOP_w_, int ordre, int reals, double **image_X, double **image_Y);

extern void generer_3image(int taille, int reals,double **image_X,double **image_Y, double **image_A_tot, double **image_A_ICEO, double **image_A_delta);

extern double *test_tot(int reals, int nbw, int nbw_, double **tab_tot);

extern double *test_delta(int reals, int nbw, int nbw_, double **tab_delta);

extern double log_Proba_p(double ro,double p);

extern double *test_LRT(double DOP_F, double DOP_w_, double DOP_w, int reals, int nbw, int nbw_, double **tab);
extern void generer_COR(double *h0_test,double *h1_test,int reals,char nom_fichier[]);

extern double *test_GLRT(int reals, int nbw, int nbw_, double **tab);

extern double log_Proba_up(double s,double u,double DOP);

extern double *test_LRT_s(double moy_h1_w,double moy_h0,double DOP_h1_w,double DOP_h0,int reals,int nbw,int nbw_,double **tab);

extern double *test_GLRT_s(int reals, int nbw, int nbw_, double **tab_tot);

extern double calculer_auc(double *vect_x,double *vect_y,int nb);

extern double *test_var_gausse(int reals, int nbw, int nbw_, double **tab_tot);

extern double *test_intensite(int reals, int nbw, int nbw_, double **tab_tot);

extern double *test_gausse_multi(int reals, int nbw,int nbw_, double **tab_test);

extern double *test_diff_logmoment(int reals, int nbw,int nbw_, double **tab_test);

extern double *test_quot_logmoment(int reals, int nbw,int nbw_, double **tab_test);
