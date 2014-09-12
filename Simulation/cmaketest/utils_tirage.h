/***************************************************************
Nom du fichier: /.../qtk_genere.h
Descriptif    :	 Bibliothèque pour bruitage d'une image bruitée par bruit sub-poissonnien .
Auteur        : J.FADE .
Date          : 10/03/07
Nom du coordinateur informatique: F Galland .
****************************************************************/

// fonction gammln
// GAMMLN   Natural log of the complete Gamma function.
//  GAMMLN(X) returns the log of the gamma function
//   for every element of X.

//   Useful in formulas involving, e.g., gamma(x)/gamma(y)
//     for large x and y.

// Peter R. Shaw, Woods Hole Oceanographic Institution
// Woods Hole, MA 02543
// (508) 457-2000 ext. 2473  pshaw@aqua.whoi.edu
// Converted from the Fortran subroutine "GAMMLN" in:
// Numerical Recipes, Press et al., Cambridge, 1986.
extern float gammln(float x);

//==========================================================
// 14 : GENERATEUR DE POISSON PAR LA METHODE DE REJECTION (format DOUBLE)
// Descriptif            : cette fonction bruite une image de reference
//			   par un bruit de Poisson.
// L'image a bruiter est transmise d'un seul bloc
//==========================================================
//ENTREES :
//			alpha		: tel que alpha*image[i] = esperance de la loi
//                      larg, haut	: dimensions de l'image de sortie
//			im_orig		: image a bruiter
//                      init		: initialisation ou non de la graine
//
//SORTIES :		im_bruit	: image bruitee
//==========================================================
//Auteur                :M.Guillaume
//Date                  :19/11/98
//==========================================================
//ATTENTION: LA GRAINE ALEATOIRE N"EST JAMAIS INITIALISEE ICI
//==========================================================
extern void poisson_rejection_d(double alpha,double **im_orig,int larg,int haut, double **im_bruit);


//==========================================================
// 15 : GENERATEUR DE POISSON PAR LA METHODE DIRECTE
// Descriptif            : cette fonction bruite une image de reference
//			   par un bruit de Poisson.
// L'image a bruiter est transmise d'un seul bloc
//==========================================================
//ENTREES :
//                      alpha           : tel que alpha*image[i] = esperance de la loi
//                      larg, haut      : dimensions de l'image de sortie
//                      im_orig         : image a bruiter
//                      init            : initialisation ou non de la graine
//
//SORTIES :             im_bruit        : image bruitee
//==========================================================
//Auteur                :M.Guillaume
//Date                  :19/11/98
//==========================================================
//ATTENTION 1: LA GRAINE ALEATOIRE N"EST JAMAIS INITIALISEE ICI
//ATTENTION 2: Cette methode est plus rapide mais peu fiable
//		pour les valeurs faibles de alpha
//=======================================================
extern void poisson_direct_d( double alpha, double **im_orig,int larg,int haut, double **im_bruit);


/*----------------------------------------------------------------
11 : nom de la fonction :bruit_exponentiel_d(m,idim,jdim,init,sortie)
Descriptif    	: Genere un tableau de variables aleatoires blanches
			et de pdf exponentielle (images en format double) .
entree          : m		moyenne de la distribution
		  idim,jdim	dimensions du tableau
retour          : sortie.
remarques       : ne pas confondre cette pdf avec la loi de Laplace !!
Auteur       	: O. Germain
-------------------------------------------------------------------*/
extern void bruit_exponentiel_d(double m,int idim,int jdim, double **sortie);




/*----------------------------------------------------------------
12 : nom de la fonction :bruit_gamma_d(m,ordre,idim,jdim,init,sortie)
Descriptif    	: Genere un tableau de variables aleatoires blanches 
			et de pdf exponentielle. (images en format double)
			entree          : m		moyenne de la distribution
	        : ordre		ordre de la loi Gamma
		idim,jdim	dimensions du tableau
		retour          : sortie.
		remarques       : ATTENTION, ceci ne vaut que pour des ordres entiers....
Auteur       	: V.PAGE
-------------------------------------------------------------------*/
extern void bruit_gamma_d(double m,int ordre,int idim,int jdim,double **tampon, double **sortie);



/*=========================================================================*/
/*                                                                         */
/* 01 :               TIRAGE POISSON                                       */
/*                                                                         */
/*=========================================================================*/
/* descriptif : Tire une  variable aleatoire poisson                       */
/*              de moyenne m                                               */
/*                                                                         */  
/* L'algo utilisé est donné par les numerical recipes p296                 */
/*=========================================================================*/
/* entree :  m     : parametres de la loi et proba de réussite          */
/* sortie :  tirage : realisation de la v.a. binomia                       */
/*=========================================================================*/
/* auteur : F.Julien                                                       */
/* date   : 03/03/08                                                       */
/* Coordinateur : F. G.                                                     */
/*=========================================================================*/
extern double  tirage_poisson(double m) ;



/*=========================================================================*/
/*                                                                         */
/* 01 :               BRUITAGE IMAGE PAR BRUIT POISSON                     */
/*                                                                         */
/*=========================================================================*/
/* descriptif : Bruite une image avec un bruit      poissonien             */
/*   selon la valeur de l'intensité du pixel : methode rejection/directe   */
/*=========================================================================*/
/* entree : im_bruit: image a bruiter                                      */
/*          larg,haut: dimensions de l'image                               */
/*	    alpha    : facteur multiplicatif des intensités                */
/* sortie : im_bruit : image bruitée                                       */
/*=========================================================================*/
/* auteur : F.Julien                                                       */
/* date   : 08/03/08                                                       */
/* Coordinateur : F. G.                                                     */
/*=========================================================================*/
extern void bruite_poisson( double **im_orig, int larg, int haut, double alpha, double **im_bruit);


/*=========================================================================*/
/*                                                                         */
/* 01 :               TIRAGE BINOMIAL                                      */
/*                                                                         */
/*=========================================================================*/
/* descriptif : Tire une  variable aleatoire binomiale                     */
/*              de parametres N et pp                                      */
/* P(x)=C^N_x pp^x (1-pp)^(1-x)                                            */  
/* La moy est n*pp                                                         */
/* L'algo utilisé est donné par les numerical recipes p296                 */
/*=========================================================================*/
/* entree : n ,pp     : parametres de la loi et proba de réussite          */
/* sortie :  tirage : realisation de la v.a. binomia                       */
/*=========================================================================*/
/* auteur : F.Julien                                                       */
/* date   : 03/06/07                                                       */
/* Coordinateur : F. G.                                                     */
/*=========================================================================*/
extern double  tirage_binomial(double pp, int n);






/*=========================================================================*/
/*                                                                         */
/* 01 :               BRUITAGE IMAGE PAR SOUS POISSON BINOMIAL             */
/*                                                                         */
/*=========================================================================*/
/* descriptif : Bruite une image avec un bruit sous-poissonien             */
/*  modèle binomial de facteur de Fano min F au plus intense de l'image    */
/*=========================================================================*/
/* entree : im_bruit: image a bruiter                                      */
/*          larg,haut: dimensions de l'image                               */
/*	    alpha    : facteur multiplicatif des intensités                */
/*          fano_factor : Facteur de Fano au plus intense de l'image       */
/* sortie : im_bruit : image bruitée                                       */
/*=========================================================================*/
/* auteur : F.Julien                                                       */
/* date   : 08/03/08                                                       */
/* Coordinateur : F. G.                                                     */
/*=========================================================================*/

extern void bruite_sub_poisson_binomial( double **im_orig, int larg, int haut, double alpha, double fano_factor, double **im_bruit);



/*=========================================================================*/
/*                                                                         */
/* 01 :               TIRAGE BRUIT CAUCHY                                  */
/*                                                                         */
/*=========================================================================*/
/* descriptif : Tire une  variable aleatoire de Cauchy                     */
/*              de médiane med et de  parametre scale                      */
/*=========================================================================*/
/* entree : med ,scale     : mediane et scale parameter  gaussien          */
/* sortie :  tirage : realisation de la v.a.  cauchy                       */
/*=========================================================================*/
/* auteur : F.Julien                                                       */
/* date   : 03/06/07                                                       */
/* Coordinateur : F. G.                                                     */
/*=========================================================================*/
extern float  tirage_cauchy(float med,float scale);







/*----------------------------------------------------------------
nom de la fonction :tirage_loi_K(m,sortie)
Descriptif    	: Genere un tirage de loi K .
entree          : m		moyenne de la distribution
                  ord_M         ordres de la distribution  
		  ord_L

                  Loi K   Multiplication de 2 tirages Gamma 
		  une de moyenne unitaire et ordre ord_L
		  autre de moyenne m et ordre ord_M
retour          : sortie.

Auteur       	: J.Fade
-------------------------------------------------------------------*/

extern float tirage_loi_K(float m, int ord_M, int ord_L);




/*----------------------------------------------------------------
nom de la fonction :tirage_weibull(m,sortie)
Descriptif    	: Genere un tirage de pdf weibull.
entree          : m		moyenne de la distribution
                              1/m*exp(-t/m)
retour          : sortie.
remarques       : ne pas confondre cette pdf avec la loi de Laplace !!
Auteur       	: J.Fade
-------------------------------------------------------------------*/

extern float tirage_weibull(float lambda, float ordre);






/*----------------------------------------------------------------
nom de la fonction :tirage_bernoulli(eta)
Descriptif    	: Genere un tirage de variable bernoulli.
entree          : eta		moyenne de la distribution
                      P(1)=eta ; P(0)=1-eta
retour          : sortie.
Auteur       	: J.Fade
-------------------------------------------------------------------*/

extern float tirage_bernoulli(float eta);







/*----------------------------------------------------------------
nom de la fonction :tirage_exponentiel(m,sortie)
Descriptif    	: Genere un tirage de pdf exponentielle.
entree          : m		moyenne de la distribution
		  
retour          : sortie.
remarques       : ne pas confondre cette pdf avec la loi de Laplace !!
Auteur       	: J.Fade
-------------------------------------------------------------------*/

extern float tirage_exponentiel(float m);






/******----------------------------------------------------------------
       tirage_gamma(m,ordre,sortie)
       Descriptif      : Genere un tirage de pdf Gamma
       entree          : m		moyenne de la distribution
       : ordre		ordre de la loi Gamma
       
       retour          : tirage
       remarques       : ATTENTION, ceci ne vaut que pour des ordres entiers....
       Auteur          : J.Fade
- ----- -------------------------------------------------------------------*/

extern float tirage_gamma(float m, int ordre);

/******----------------------------------------------------------------
       tirage_gauss(m,sig,sortie)
       Descriptif      : Genere un tirage gaussien 
       entree          : m		moyenne de la distribution
       : sig 		ecart-type de la loi Gamma
       
       retour          : tirage
       Auteur          : J.Fade
- ----- -------------------------------------------------------------------*/

float tirage_Gauss(float m, float sig);



/******----------------------------------------------------------------
       tirage_gamma_gen(m,ordre,sortie)
       Descriptif      : Genere un tirage de pdf Gamma ordre non entier
       entree          : m		moyenne de la distribution
       : ordre		ordre de la loi Gamma
       
       retour          : tirage
       remarques       : ATTENTION, ceci ne vaut que pour des ordres entiers....
       Auteur          : J.Fade
- ----- -------------------------------------------------------------------*/

extern float tirage_gamma_gen(float m, float ordre);



/*----------------------------------------------------------------
nom de la fonction : generation_sub_poisson(alpha,im_orig,larg,haut,init,sqz_factor,integ_time,im_bruit)
GENERATEUR DE BRUIT SUB-POISSONNIEN
 Descriptif            : cette fonction bruite une image de reference
			   par un bruit de sub-Poisson.
 L'image a bruiter est transmise d'un seul bloc
-------------------------------------------------------------------
ENTREES :
			alpha		: tel que alpha*image[i] = esperance de la loi
			im_orig		: image a bruiter
                         larg, haut  	: dimensions de l'image de sortie
			init		: initialisation ou non de la graine
                         sqz_factor        : facteur de squeezing
                         integ_time        : temps d'integration

SORTIES :		im_bruit	: image bruitee
--------------------------------------------------------------------
Auteur                :J.Fade	
Date                  :10/03/2007
-------------------------------------------------------------------*/
extern int generation_sub_poisson_pseudo_photons(float alpha,float **im_orig,int larg,int haut,
			    float sqz_factor, float integ_time, float **im_bruit);




/*----------------------------------------------------------------
nom de la fonction : generation_sub_poisson_gamma(alpha,im_orig,larg,haut,init,sqz_factor,integ_time,im_bruit)
GENERATEUR DE BRUIT SUB-POISSONNIEN

UTILISEE POUR TESTER  LA STATISTIQUE DES TEMPS D'ATTENTES ENTRES ARRIVEES PHOTONIQUES

 Descriptif            : cette fonction bruite une image de reference
			   par un bruit de sub-Poisson.
 L'image a bruiter est transmise d'un seul bloc
-------------------------------------------------------------------
ENTREES :
			alpha	       : tel que alpha*image[i] = esperance de la loi
			im_orig	       : image a bruiter
                         larg, haut        : dimensions de l'image de sortie
			init	       : initialisation ou non de la graine
                         sqz_factor        : facteur de squeezing
                         integ_time        : temps d'integration

SORTIES :		im_bruit	: image bruitee
--------------------------------------------------------------------
Auteur                :J.Fade	
Date                  :10/03/2007
-------------------------------------------------------------------*/
extern void generation_sub_poisson_gamma(float alpha,float **im_orig,int larg,int haut,
					 float sqz_factor, float integ_time, float **im_bruit);







/*----------------------------------------------------------------
nom de la fonction : generation_sub_poisson_loiK(alpha,im_orig,larg,haut,init,sqz_factor,integ_time,im_bruit)
GENERATEUR DE BRUIT SUB-POISSONNIEN
-*/

extern void generation_sub_poisson_loiK(float alpha,float **im_orig,int larg,int haut,
				 float sqz_factor, float integ_time, float **im_bruit);





/*----------------------------------------------------------------
nom de la fonction : generation_sub_poisson_weibull(alpha,im_orig,larg,haut,init,sqz_factor,integ_time,im_bruit)
GENERATEUR DE BRUIT SUB-POISSONNIEN
*/

extern void generation_sub_poisson_weibull(float alpha,float **im_orig,int larg,int haut,
				    float sqz_factor, float integ_time, float **im_bruit);





/*----------------------------------------------------------------
nom de la fonction : generation_sub_poisson_lognorm(alpha,im_orig,larg,haut,init,sqz_factor,integ_time,im_bruit)
GENERATEUR DE BRUIT SUB-POISSONNIEN
*/

extern void generation_sub_poisson_lognorm(float alpha,float **im_orig,int larg,int haut,
				    float sqz_factor, float integ_time, float **im_bruit);




/*---------------------------------------------------------------- nom
de la fonction : generation_gaussien(alpha,im_orig,larg,haut,snr,
im_bruit, itot , waist) GENERATEUR DE BRUIT SUB-POISSONNIEN
*/
extern void generation_sub_poisson_gamma_non_entier(float alpha,float **im_orig,int larg,int haut,
						    float sqz_factor, float integ_time, float **im_bruit);



/*----------------------------------------------------------------
nom de la fonction : generation_sub_poisson_gauss_cauchy(alpha,im_orig,larg,haut,init,sqz_factor,integ_time,im_bruit)
GENERATEUR DE BRUIT SUB-POISSONNIEN
*/
extern void generation_sub_poisson_cauchy(float alpha,float **im_orig,int larg,int haut,
				    float sqz_factor, float integ_time, float **im_bruit);



