/***************************************************************
Nom du fichier: /.../utils_fade.h
Descriptif    :	  Bibliothèque de fonctions diverses (nouvelles + reprogrammation d'anciennes fonctions phyti
Auteur        : J.Fade .
Date          : 03/08 
Nom du coordinateur informatique: F.Galland.
****************************************************************/

extern void d_carre_rect(int dim_y,int dim_x,double **input_array,double **out_array);

extern void d_racine_carree_2d(int nbre_l,int nbre_c,double **image);

extern double **d_2D_Lire_bin_image(char *nom, int *w, int *h);

extern int d_2D_Sauv_bin_image(char *image_file_name,double **image,int colonne,int ligne);

extern double bar_err_est( double *tab, int sqreals);

extern double bar_err_var( double *tab, int sqreals);

extern void saveVec2oct(double *input, int dim_Ligne, char* name_var, char* name_file, int erase_file);

extern void saveVecI2oct(int *input, int dim_Ligne, char* name_var, char* name_file, int erase_file);

extern void saveMat2oct(double **input, int dim_Ligne, int dim_Col, char* name_var, char* name_file, int erase_file);

extern void stat_image( double **image, double *variance, double *moy, int haut, int larg);

extern double d_round( double reel);

extern int mathieu(int age);

extern double d_moyenne(int card, double *tab);

extern double moyenne(int start_index,int end_index, double *tab);

extern double d_variance(int card, double *tab, double moyenne);

extern double d_variance_unbias(int card, double *tab, double moyenne);

extern float *f_alloue_1d(int size);

extern double *d_alloue_1d(int size);

extern double **d_alloue_2d(int h, int w);

extern void  **d_efface_2d(int h, int w, double **img);

extern double **d_annul_mat(int idim, int jdim,double **mat);

extern double **d_addmat(int idim, int jdim,double **sortie, double **in1, double **in2);

extern int compare_doubles(const void *a, const void *b);

extern double d_max( double a, double b);

extern double d_max_1d( double *a, int dim);

extern double variance(int start_index, int end_index, double *tab, double moyenne);

