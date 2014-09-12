/***************************************************************
Nom du fichier: /.../utils_fade.c
Descriptif    :	 Bibliothèque de fonctions diverses (nouvelles + reprogrammation d'anciennes fonctions phyti
mises à jour.
Auteur        : J.FADE .
Date          : 03/08
Nom du coordinateur informatique: F Galland .
****************************************************************/

//# include "lib_std.h"
# include <stdio.h>
# include <math.h>
# include <stdlib.h>

//#include "lib_math.h"
//#include "lib_phyti.h"
//# include "math_vect_gcc.h"
//# include "im_utils_gcc.h"
//# include "util_srf_gcc.h"
//# include "im_manip_gcc.h"
//# include "math_vect_2_gcc.h"

//# include "lib_quantique.h"


 /*=========================================================================*/
/*                                                                         */
/* 02 :        CARRE DE 1 TABLEAU 2dim (CAS Double)                         */

/*                                                                         */
/*=========================================================================*/
/* entree : dim_y,dim_x  : dimensions du tableau 2d (dim_y=hauteur)        */
/*          input_array  : tableau 2d                                      */
/* sortie : output_array : tableau 2d                                      */
/*=========================================================================*/

/*=========================================================================*/
/* auteur : Guerault */  
/* date   : 03/10/97 */  
/* modif  : 06/10/97 */  
/* Coordinateur : F. Goudail */  
/*=========================================================================*/

void d_carre_rect(int dim_y,int dim_x,double **input_array,double **out_array)
{
  int i, j;
  
  for (i = 0; i < dim_y ; i++){
    for (j = 0; j < dim_x ; j++){
      out_array[i][j] =input_array[i][j] * input_array[i][j];
    }/* next j */
  }/* next i */

}
// End of d_carre_rect


/*=========================================================================*/
/*                                                                         */
/* 01 :      MISE A LA RACINE CARREE DE CHAQUE PIXEL DU TABLEAU            */
/*                                                                         */
/*=========================================================================*/
/* descriptif : 						           */
/*	Cette fonction effectue la racine carree de chaque                 */
/*	pixel du tableau : le test pour s'assurer que les pixels sont      */
/*	positifs est realise : message d'erreur dans ce cas, et les        */
/*	pixels incrimines recoivent le niveau 255     !!!!!                */
/* formule    : image(i,j) = sqrt(image(i,j)                               */
/*=========================================================================*/
/* entree : nbre_l,nbre_c : dimensions du tableau                          */
/* sortie : image         : tableau 2d  	                           */
/*=========================================================================*/
/* auteur : Tramoni & Martin                                               */
/* date   : 03/10/97                                                       */
/* modif  : 06/10/97                                                       */
/* Coordinateur : T. Gaidon                                                */
/*=========================================================================*/

void d_racine_carree_2d(int nbre_l,int nbre_c,double **image)
{
int       i,j,a=0;   /* variables de compteur */
double    val;       /* variable temporaire */


for ( i = 0 ; i < nbre_l ; i++)
	for ( j = 0 ; j < nbre_c ; j++)
		{
		val = (double) image[i][j];
		if (val >= 0.0 )
			{
			val = sqrt(val);
			image[i][j] = (double) val;
			}
		else
			{
			image[i][j] = 255.0;
			a++;
			}
		}
if ( a > 0 )
   {
     fprintf (stderr," Des pixels a valeurs negatives se trouvaient dans l'image\n");
   }

}

/***********************************************end of f_racine_carree_2d()*/

double moyenne(int start_index, int end_index, double *tab)
{
	int i=0;
	double sum=0.0;
	
	for (i=start_index;i<end_index;i++)
	{
		sum+=tab[i];
	}
	return (sum/(end_index-start_index));
}


double d_moyenne(int card, double *tab)
{
 int i=0; 
 double sum=0.0;

	for (i=0;i<card;i++)
		{
		sum+=tab[i];
		}
return (sum/card);
}

double variance(int start_index, int end_index, double *tab, double moyenne)
{
	int i;
	double sum=0.0;
	int card=end_index-start_index;
	for (i=start_index;i<end_index;i++)
	{
		sum+=tab[i]*tab[i];
	}
	return (sum/card-moyenne*moyenne);
}


double d_variance(int card, double *tab, double moyenne)
{
 int i=0; 
 double sum=0.0;

	for (i=0;i<card;i++)
		{
		sum+=tab[i]*tab[i];
		}
return (sum/card-moyenne*moyenne);
}

double d_variance_unbias(int card, double *tab, double moyenne)
{
 int i=0; 
 double sum=0.0;

	for (i=0;i<card;i++)
		{
		sum+=tab[i]*tab[i];
		}
return ((sum-card*moyenne*moyenne)/(card-1));
}

float *f_alloue_1d(int size)
{
 float  *vector;

	vector = calloc (size, sizeof(float));
	if( vector == NULL )
		{
		     fprintf(stderr,"Allocation impossible");
		     exit(EXIT_FAILURE);
		}
return (vector);
}

double *d_alloue_1d(int size)
{
 double  *vector;

	vector = calloc (size, sizeof(double));
	if( vector == NULL )
		{
		     fprintf(stderr,"Allocation impossible");
		     exit(EXIT_FAILURE);
		}
return (vector);
}

double **d_alloue_2d(int h, int w)
{
 double  **tableau;
 int i=0;

	tableau = calloc (h, sizeof(double*));

	if( tableau == NULL )
		{
		     fprintf(stderr,"Allocation impossible");
		     exit(EXIT_FAILURE);
		}

	for (i=0;i<h;i++)
		{
		tableau[i]=calloc (w, sizeof(double));
		if( tableau[i] == NULL )
			{
			     fprintf(stderr,"Allocation impossible");
			     exit(EXIT_FAILURE);
			}
		}
return (tableau);
}

// Fonction de désallocation
void d_efface_2d(int h, int w, double **img)
{
 int i=0;

	for (i=0;i<h;i++)
		{
		free(img[i]);
		}
}


/*-------------caracteristiques de la fonction -------------------

nom de la fonction : **d_2D_Lire_bin_image
donnees globales   :
entree             : char *nom : nom du fichier filtre .
		     int * ligne,colonne : ptrs pour retour ligne colonne
retour             : NULL --->  Pb
		     double **data : le filtre sous forme double 2D .
appel              :
remarques          : Lit un filtre dans un fichier binaire .
		     Renvoie le nombre de lignes et de colonnes  .
--------------------------------------------------------------------------*/
double **d_2D_Lire_bin_image(char *nom, int *w, int *h)
{
  FILE *file;
  float *tmp;
  double **data,**d_alloue_2d();
  char 	chaine1[100], *chaine2[1];
  int i,j;
  
  
  if(nom==NULL) file = stdin;
  else
    {
      file = fopen(nom,"r");
      if(file==NULL) {printf("erreur a l'ouverture de |%s|\n",nom);return(NULL);}
    }
  //fscanf(file,"%d %d\n",h,w);
  // premiere ligne du fichier stockee dans chaine1
  fgets (chaine1, 100, file); 
  // le premier chiffre de la chaine est mis dans h,
  // le reste de la chaine de caracteres dans chaine2.
  *h = (int) strtol(chaine1, chaine2, 10); 
  // chaine 2 est recupere pour avoir la largeur
  *w = (int) strtol(*chaine2, NULL, 10);
  data=d_alloue_2d((*h),(*w));
  tmp = f_alloue_1d(*w);
  if(data==NULL) 
    {
      fprintf(stderr,"Fonction d_2D_Lire_raw_image(): Impossible d'allouer de la place pour ptr image");
      if(nom!=NULL) fclose(file);
      free(tmp);
      return(NULL);
    }
  
  for(i=0;i<(*h);i++)
    {
      fread( tmp, sizeof(float),(*w),file);
      
      for (j=0;j<(*w);j++)
	{
	  data[i][j] = (double) tmp[j];
	  // fprintf(stderr,"%g ", (float) data[i][j]); // VERIF
	}
      //fprintf(stderr,"\n"); //VERIF 
    }
  
  if(nom!=NULL) fclose(file);
  free(tmp);
  return(data);
} // end of d_2D_Lire_bin_image


/*-------------caracteristiques de la fonction -------------------

nom de la fonction : d_2D_Sauv_bin_image
donnees globales   :
entree             : char *nom : nom du fichier filtre .
		     int  ligne,colonne : nb ligne colonne
		     float **image : le image sous forme double 2D .

retour             : -1 Pb
		      0 Ok .
appel              :
remarques          : Ecrit un une image en format raw data.
--------------------------------------------------------------------------*/
int d_2D_Sauv_bin_image(char *image_file_name,double **image, int colonne, int ligne)
{
  FILE *file;
  float *tmp;
  int i,j;
  
  tmp = f_alloue_1d(colonne);
  
  if(image_file_name==NULL) file = stdout;
  else
    {
      file = fopen(image_file_name,"w");
      if(file==NULL) {printf("erreur a la creation de |%s|\n",image_file_name);free(tmp);return(-1);}
    }
  fprintf(file,"%d %d\n",ligne,colonne);
  
  for(i=0;i<ligne;i++)
    {
      for (j=0;j<colonne;j++)
	{
	  tmp[j]= (float) image[i][j]; 
	  // fprintf(stderr,"%g ", tmp[j]); VERIF
	}
      // fprintf(stderr,"\n"); VERIF
      fwrite( tmp, sizeof(float), colonne, file);
    }
  
  
  if(image_file_name!=NULL) fclose(file);
  free(tmp);
  return(1);
  
} // End of d_2D_Sauv_bin_image



double bar_err_est( double *tab, int sqreals)
{
  double *tab1,*tab2;
  int k,i;

  tab1= (double*) calloc( (unsigned) sqreals, sizeof(double));
  tab2= (double*) calloc( (unsigned) sqreals, sizeof(double));


  for (k=0;k<sqreals;k++)
    {
      for (i=0;i<sqreals;i++)
	{ tab1[i]=tab[i+k*sqreals];

	}
      tab2[k]=d_moyenne(sqreals,tab1);
    }
  
  return sqrt(d_variance(sqreals,tab2,d_moyenne(sqreals,tab2))/sqreals)/2;
  free(tab1);
  free(tab2);
}

double bar_err_var( double *tab, int sqreals)
{
  double *tab1,*tab2;
  int k,i;

  tab1= (double*) calloc( (unsigned) sqreals, sizeof(double));
  tab2= (double*) calloc( (unsigned) sqreals, sizeof(double));


  for (k=0;k<sqreals;k++)
    {
      for (i=0;i<sqreals;i++)
	{ tab1[i]=tab[i+k*sqreals];

	}
      tab2[k]=d_variance(sqreals,tab1,d_moyenne(sqreals,tab1));
    }
  
  return sqrt(d_variance(sqreals,tab2,d_moyenne(sqreals,tab2))/sqreals)/2;
  free(tab1);
  free(tab2);
}






// vecteur double
void saveVec2oct(double *input, int dim_Ligne, char* name_var, char* name_file, int erase_file)
{
 FILE *file ;
 int i ;

 if (erase_file == 1)
   {
     file=fopen(name_file,"w");
     fprintf(file, "# Created by Physics and Images Processing Group \n") ;
   }
 else
   file=fopen(name_file,"a");

 // entete octave
 fprintf(file, "# name: %s\n", name_var) ;
 fprintf(file, "# type: matrix\n") ;
 fprintf(file, "# rows: %d\n",dim_Ligne) ;
 fprintf(file, "# columns: 1\n") ;

 // donnees
 for (i=0;i<dim_Ligne;i++){
   fprintf(file," %.16e", input[i]) ;fprintf(file, "\n") ;  //" %f",input[i]); fprintf(file,"\n");
 }  fclose(file);
}


// vecteur entier
void saveVecI2oct(int *input, int dim_Ligne, char* name_var, char* name_file, int erase_file)
{
 FILE *file ;
 int i ;

 if (erase_file == 1)
   {
     file=fopen(name_file,"w");
     fprintf(file, "# Created by Physics and Images Processing Group \n") ;
   }
 else
   file=fopen(name_file,"a");

 // entete octave
 fprintf(file, "# name: %s\n", name_var) ;
 fprintf(file, "# type: matrix\n") ;
 fprintf(file, "# rows: %d\n",dim_Ligne) ;
 fprintf(file, "# columns: 1\n") ;

 // donnees
 for (i=0;i<dim_Ligne;i++){
   fprintf(file," %d", input[i]) ;fprintf(file, "\n") ;  //" %f",input[i]); fprintf(file,"\n");
 }  fclose(file);
}






// Matrice real
void saveMat2oct(double **input, int dim_Ligne, int dim_Col, char* name_var, char* name_file, int erase_file)
{
 FILE *file ;
 int j,i ;

 if (erase_file == 1)
   {
     file=fopen(name_file,"w");
     fprintf(file, "# Created by Physics and Images Processing Group \n") ;
   }
 else
   file=fopen(name_file,"a");

 // entete octave
 fprintf(file, "# name: %s\n", name_var) ;
 fprintf(file, "# type: matrix\n") ;
 fprintf(file, "# rows: %d\n",dim_Ligne) ;
 fprintf(file, "# columns: %d\n", dim_Col) ;

 // donnees
 for (i=0;i<dim_Ligne;i++) {
   for(j=0;j<dim_Col;j++){
       fprintf(file, " %.16e", input[i][j]) ;
   }
   fprintf(file, "\n") ;
 }

 fclose(file);
}



void stat_image( double **image, double *variance, double *moy, int haut, int larg)
{
  int i,j;
  double sum=0., var=0.;
  
  for (i=0;i<haut;i++)
    for (j=0;j<larg;j++)
      {
	sum=sum+image[i][j];
	var=var+image[i][j]*image[i][j];
      }
  *moy=sum/haut/larg;
  *variance=var/haut/larg - (*moy)*(*moy);
}//fin stat_image


double d_round( double reel)
{
  return floor(reel + 0.5);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void d_annul_mat(int idim, int jdim,double **mat)
{
int i,j;
  for (i=0;i<idim;i++)
    for (j=0;j<jdim;j++)
	{
	  mat[i][j]=0.;
	}
}

double d_max( double a, double b)
{
double ret=0.;

if (a>b){
ret=a; }
else{ 
ret=b; }
return ret;
}

double d_max_1d( double *a, int dim)
{
double ret;
int k=0;

ret=a[0];

for(k=1;k<dim;k++)
{
ret= d_max( a[k], ret);
}
return ret;
}



void d_addmat(int idim, int jdim,double **sortie, double **in1, double **in2)
{
int i,j;
  for (i=0;i<idim;i++)
    for (j=0;j<jdim;j++)
	{
	  sortie[i][j]=in1[i][j]+in2[i][j];
	}
}



int compare_doubles(const void *a, const void *b)
     {
        double *da = ( double *) a;
        double *db = ( double *) b;
     
       return (*da > *db) - (*da < *db);
     }

