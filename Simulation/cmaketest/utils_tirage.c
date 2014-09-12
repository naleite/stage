/***************************************************************
Nom du fichier: /.../qtk_genere.c
Descriptif    :	 Bibliothèque pour bruitage d'une image bruitée par bruit sub-poissonnien .
Auteur        : J.FADE .
Date          : 10/03/07
Nom du coordinateur informatique: F Galland .
****************************************************************/


# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include "utils_fade.h"

# define PI2     6.28318531
# define PI     3.1415926535



//#include <gsl/gsl_rng.h> 
//#include <gsl/gsl_randist.h>
/* #include <gsl/gsl_cdf.h> */
/* #include <gsl/gsl_sf_psi.h> */
/* #include <gsl/gsl_roots.h> */



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

float gammln(float xx)
{
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
  int j;
  
  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp+log(2.50662827465*ser);
}

/* //========================================================== */
/* // 14 : GENERATEUR DE POISSON PAR LA METHODE DE REJECTION (format DOUBLE) */
/* // Descriptif            : cette fonction bruite une image de reference */
/* //			   par un bruit de Poisson. */
/* // L'image a bruiter est transmise d'un seul bloc */
/* //========================================================== */
/* //ENTREES : */
/* //			alpha		: tel que alpha*image[i] = esperance de la loi */
/* //                      larg, haut	: dimensions de l'image de sortie */
/* //			im_orig		: image a bruiter */
/* //                      init		: initialisation ou non de la graine */
/* // */
/* //SORTIES :		im_bruit	: image bruitee */
/* //========================================================== */
/* //Auteur                :M.Guillaume */
/* //Date                  :19/11/98 */
/* //========================================================== */
/* //ATTENTION: LA GRAINE ALEATOIRE N"EST JAMAIS INITIALISEE ICI */
/* //========================================================== */
void poisson_rejection_d(double alpha,double **im_orig,int larg,int haut, double **im_bruit)
{

double   xm;
static double    sq,alxm,oldm=(-1.0);
static double   g;
double           em,y;
double          t;
int             i,j;




// Boucle portant sur les pixels de l'image a bruiter
//==========================================================
for(i=0;i<haut;i++)
  for(j=0;j<larg;j++)
    {

    xm=alpha*im_orig[i][j];

    /* Si xm est le meme qu'au precedent tirage, on ne */
    /* recalcule pas le critere g                      */

    if (xm != oldm)
      {
      oldm=xm;
      sq=sqrt(2.0*xm);

      if (xm!=0) alxm=log(xm);
      else alxm=-1000;

      g=xm*alxm-gammln(xm+1.0);

      }

    do
      {
      do
        {
        y=tan(M_PI*drand48());
        em=sq*y+xm;
        } while(em<0.0);

      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      } while (drand48()>t);

    /* Valeur placee au pixel i,j de l'image bruitee */
    im_bruit[i][j] = em;
    // fprintf(stdout,"%f ",em);
    } /* Fin de la boucle de parcours des pixels */

} /* fin poisson_rejection */






/* //========================================================== */
/* // 15 : GENERATEUR DE POISSON PAR LA METHODE DIRECTE */
/* // Descriptif            : cette fonction bruite une image de reference */
/* //			   par un bruit de Poisson. */
/* // L'image a bruiter est transmise d'un seul bloc */
/* //========================================================== */
/* //ENTREES : */
/* //                      alpha           : tel que alpha*image[i] = esperance de la loi */
/* //                      larg, haut      : dimensions de l'image de sortie */
/* //                      im_orig         : image a bruiter */
/* //                      init            : initialisation ou non de la graine */
/* // */
/* //SORTIES :             im_bruit        : image bruitee */
/* //========================================================== */
/* //Auteur                :M.Guillaume */
/* //Date                  :19/11/98 */
/* //========================================================== */
/* //ATTENTION 1: LA GRAINE ALEATOIRE N"EST JAMAIS INITIALISEE ICI */
/* //ATTENTION 2: Cette methode est plus rapide mais peu fiable */
/* //		pour les valeurs faibles de alpha */
/* //======================================================= */
void poisson_direct_d( double alpha, double **im_orig,int larg,int haut, double **im_bruit)
{
static double    oldm=(-1.0);
static double   g;
double           em,xm;
int             i,j;
double          t;


// Boucle portant sur les pixels de l'image a bruiter
//==========================================================
for(i=0;i<haut;i++)
  for(j=0;j<larg;j++)
    {

    xm=alpha*im_orig[i][j];
   
    /* Si xm est le meme qu'au precedent tirage, on ne */
    /* recalcule pas le critere g                      */

    if (xm != oldm)
      {
	oldm=xm;
	g=exp(-xm);
	//if ( isnormal(g)  ){exit(0);}
	//	fprintf(stdout,"%f ",log10(g));
      }

    /* On tire plusieurs nombres t de facon uniforme.     */
    /* On multiplie ces nombres successivement jusqu'a ce */
    /* que leur produit depasse g. Le nombre de tirages   */
    /* effectues-1 (em) est issu d'une loi de Poisson.    */

    em=-1.0;
    t=1.0;

    do
      {
      ++em;
      t *= drand48();
      } while (t>g);

    /* Valeur placee au pixel i,j de l'image bruitee */
    im_bruit[i][j] = em;
  

    } /* Fin de la boucle de parcours des pixels */

} /* fin poisson_direct */







/*----------------------------------------------------------------
11 : nom de la fonction :bruit_exponentiel_d(m,idim,jdim,init,sortie)
Descriptif    	: Genere un tableau de variables aleatoires blanches
			et de pdf exponentielle (images en format double) .
entree          : m		moyenne de la distribution
		  idim,jdim	dimensions du tableau
		  init		initialisation (1) ou non (0) de la graine.
retour          : sortie.
remarques       : ne pas confondre cette pdf avec la loi de Laplace !!
Auteur       	: O. Germain
-------------------------------------------------------------------*/
void bruit_exponentiel_d(double m,int idim,int jdim, double **sortie)
{
int	i,j;
//struct timeb temps;
double z;


/*----------------------------------------------------*/
/* Initialisation du generateur de nombres aleatoires */
/*----------------------------------------------------*/

/*if(init==1)
{
  time(&(temps.time));
  srand48(temps.time);
}*/

/*-------------------*/
/* Tirage du tableau */
/*-------------------*/

for (i = 0; i < idim; i++)
{
  for (j = 0; j < jdim; j++)
  {
    z = drand48();
    sortie [i][j] =  -1.0 * m * log(1.0 - z);
  }  /* next j */
}  /* next i */

}  /* end of sub bruit_exponentiel */




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
void bruit_gamma_d(double m,int ordre,int idim,int jdim,double **tampon, double **sortie)
{
    int i;
    d_annul_mat(idim,jdim,sortie);
    for (i=1;i<=ordre;i++){
        bruit_exponentiel_d(m/ordre,idim,jdim,tampon);
        d_addmat(idim,jdim,sortie,sortie,tampon);
    }
} /* fin bruit_gamma */




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
double  tirage_poisson(double m)      // ! Pompée intégralement sur les numerical recipes !
{
  double  tirage,y,t;
  static double sq,alxm,g,oldm=(-1.0);
    
  
  if (m < 12.0)              // Methode directe pour les n faibles
    {
      if (m != oldm) {       // Calcule l'exponentielle si m nouveau
	oldm=m;
	g=exp(-m);
      }
      tirage = -1.0;
      t=1.0;
      do {
	++tirage;
	t *= drand48();
      } while (t> g);
    }
  else                 // Dans le reste des cas, il faut implanter une methode de rejection
    {
      if (m != oldm) {
	oldm=m;
	sq=sqrt(2.0*m);
	alxm=log(m);
	g=m*alxm-gammln(m+1.0);
      }
      do {
	do {
	  y=tan( M_PI*drand48() );
	  tirage=sq*y+m;             // Pour une fonction lorentzienne f=c/(1+(x-x0)/a), on tire une réalisation de
	                         // la v.a. avec U uniforme [0,1] : X= a * tan( Pi*U ) +x0;
	} while (tirage < 0.0 ); // Elimine cas ne conduisant pas à une réalisation possible
	tirage=floor(tirage);              // On  rend la réalisation entière
	t=0.9*(1.0+y*y)*exp(tirage*alxm-gammln(tirage+1.0)-g);
	} while ( drand48() > t); // Critère de rejet (calculé : cf Numerical Recipes)
      }
      

  return tirage;
}



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

void bruite_poisson( double **im_orig, int larg, int haut, double alpha, double **im_bruit)
{
  int i,j;
  double tmp; 
  
  for(i=0;i<haut;i++)
    for(j=0;j<larg;j++) 
      {
	tmp=im_orig[i][j]*alpha;
	im_bruit[i][j]=tirage_poisson(tmp); 
      }
} /* end of bruite_poisson */



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
double  tirage_binomial(double pp, int n)      // ! Pompée intégralement sur les numerical recipes !
{
  double  p,moy,tirage,sq,angle,em,y,t,g;
  int j;
  static int nold=(-1);
  static double pold=(-1.0),pc,plog,pclog,en,oldg;
  
  
  p=( pp<= 0.5 ? pp : 1.0-pp);   // invariance de la pdf binomial si p->1-p, alors, n -> N-n; On se limite à p<0.5
  moy=n*p;
  
  if (n<25)              // Methode ultra classique pour les n faibles
    {
      tirage=0,0;
      for (j=1;j<=n;j++)
	{
	  if ( drand48() < p) ++tirage;
	}
    }
  else if (moy < 1.0)   // Pour les n + forts mais moins d'un événement probable tous les 25 essais,
    {                   // on approxime bien la distrib par une realisation poisson;
      g=exp(-moy);
      t=1.0;
      for (j=0;j<=n;j++)
       	{
	  t *= drand48();
	  if (t < g) break;
	}
      tirage=( j <= n ? j : n);
    }
  else                  // Dans le reste des cas, il faut implanter une methode de rejection
    {
      if (n != nold) {
	en=n;
	oldg=gammln(en+1.0);
	nold=n;
      } if (p !=pold ) {
	pc=1.0-p;
	plog=log(p);
	pclog=log(pc);
	pold=p;
      }
      /*     en=n;             // La fonction de comparaison/domination est une fonction Lorentzienne, de mode moy, et largeur
      oldg=gammln(en+1.0);   // sq=sqrt(2.0*moy*(1-p));
      pc=1.0-p;
      plog=log(p);
      pclog=log(pc);*/
      
      sq=sqrt(2.0*moy*pc);
      do {
	do {
	  angle=M_PI*drand48();
	  y=tan(angle);
	  em=sq*y+moy;             // Pour une fonction lorentzienne f=c/(1+(x-x0)/a), on tire une réalisation de
	                           // la v.a. avec U uniforme [0,1] : X= a * tan( Pi*U ) +x0;
	} while (em < 0.0 || em >= (en + 1.0)); // Elimine cas ne conduisant pas à une réalisation possible binomial

	em=floor(em);              // On  rend la réalisation entière

	t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)- gammln(en-em+1.0)+em*plog+(en-em)*pclog);
      } while ( drand48() > t); // Critère de rejet (calculé : cf Numerical Recipes)
      tirage=em;
    }
  
  if ( p!= pp) tirage=n-tirage;      // Gestion de la symétrie de la loi binomiale si p>0.5
  return tirage;
}


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

void bruite_sub_poisson_binomial( double **im_orig, int larg, int haut, double alpha, double fano_factor, double **im_bruit)
{
  int i,j;
  int N_zero=0.0;
  double maj=0.0;
  
  
  for(i=0;i<haut;i++)
    for(j=0;j<larg;j++) 
      {
	maj= ( maj <= im_orig[i][j] ? im_orig[i][j] : maj );
      }
  
  N_zero = d_round( alpha*maj / (1-fano_factor) );  
  fprintf(stderr,"Valeur max de l'illumination image d'entree %g \n",alpha*maj);
  fprintf(stderr,"Valeur de N_zero retenue %d \n", N_zero);
  fprintf(stderr,"Valeur générée max de l'illumination image d'entree %g \n", (1-fano_factor)*N_zero);
  fprintf(stderr,"Attention: Facteur de correction de l'intensité: %g \n", (1-fano_factor)*N_zero/(alpha*maj));

  for(i=0;i<haut;i++)
    for(j=0;j<larg;j++) 
      {
	im_bruit[i][j]=tirage_binomial((1-fano_factor)*im_orig[i][j]/maj,N_zero); 
      }
} /* end of bruite_sub_poisson_binomial */




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

float  tirage_cauchy(float med,float scale) 
{ 
  
  double z; 
  float   tirage=0.;
  z=drand48();

  tirage=(float) med + scale*tan(M_PI*z);

  return tirage;
}

/*----------------------------------------------------------------
nom de la fonction :tirage_exponentiel(m,sortie)
Descriptif    	: Genere un tirage de pdf exponentielle.
entree          : m		moyenne de la distribution
                              1/m*exp(-t/m)
retour          : sortie.
remarques       : ne pas confondre cette pdf avec la loi de Laplace !!
Auteur       	: J.Fade
-------------------------------------------------------------------*/

float tirage_exponentiel(float m)
{
  float tirage=0;
  double z;
  
  z = drand48();
  tirage= (float) -1.0 * m * log(1.0 - z);
  return(tirage);    
}  /* end of tirage_exponentiel */






/******----------------------------------------------------------------
       tirage_gamma(m,ordre,sortie)
       Descriptif      : Genere un tirage de pdf Gamma
       entree          : m		moyenne de la distribution
       : ordre		ordre de la loi Gamma
       
       retour          : tirage
       remarques       : ATTENTION, ceci ne vaut que pour des ordres entiers....
       Auteur          : J.Fade
- ----- -------------------------------------------------------------------*/

float tirage_gamma(float m, int ordre)
{
  float tirage=0;
  float tampon=0;
  int i;
  
  for (i=1;i<=ordre;i++)
    {
      tampon=tirage_exponentiel((float) m/ordre);
      tirage += tampon;

    }
  return(tirage);
} 
/* fin tirage_gamma */


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

float tirage_loi_K(float m, int ord_M, int ord_L)
{
  float tirage=0;
   
  tirage= tirage_gamma( (float) 1.0,ord_M)*tirage_gamma( (float) m,ord_L);
   
  return(tirage);    
}  /* end of tirage_loi_K */



/*----------------------------------------------------------------
nom de la fonction :tirage_weibull(m,sortie)
Descriptif    	: Genere un tirage de pdf weibull.
entree          : m		moyenne de la distribution
                  ordre         ordre de la distribution  
            Necessite le calcul d'un parametre L=m/Gamma(1+1/ordre)

                  Loi Weibull   ordre/L*(t/L)^(ordre-1)*exp(-(t/L)^ordre)
retour          : sortie.

Auteur       	: J.Fade
-------------------------------------------------------------------*/

float tirage_weibull(float m, float ordre)
{
  float tirage=0;
  double z,lgam;

  lgam = lgamma(1+1/ordre);  //Calcul du log de la fonciton Gamma necessaire 
                             /* pour réajuster la moyenne */
  // Ce calcul pourra être mutualisé dans l'avenir puisque l'ordre est identique pour 1 image donnée

  z = drand48();
  //  tirage= (float) m/signgam*exp(-lgam)*pow((-log(1.0 - z)),1/ordre);
  tirage= (float) m/signgam*exp(-lgam)*pow((-log(z)),1/ordre);
  return(tirage);    
}  /* end of tirage_weibull */




/*----------------------------------------------------------------
nom de la fonction :tirage_bernoulli(eta)
Descriptif    	: Genere un tirage de variable bernoulli.
entree          : eta		moyenne de la distribution
                      P(1)=eta ; P(0)=1-eta
retour          : sortie.
Auteur       	: J.Fade
-------------------------------------------------------------------*/

float tirage_bernoulli(float eta)
{
  double z;
  
  z = drand48();
  if ( z <= eta ) z=1.0;
  else z=0.;
  return((float) z);    
}  /* end of tirage_bernoulli */







/******----------------------------------------------------------------
       tirage_gauss(m,sig,sortie)
       Descriptif      : Genere un tirage gaussien 
       entree          : m		moyenne de la distribution
       : sig 		ecart-type de la loi Gamma
       
       retour          : tirage
       Auteur          : J.Fade
- ----- -------------------------------------------------------------------*/

float tirage_Gauss(float m, float sig)
{
  float tirage=0;
  double z1,dmoy;

  dmoy= (double)m;
   
  z1 = sqrt(-2.0*log( drand48() )); 
   tirage= (float) (dmoy + sig * z1 * cos(2*M_PI *drand48()));

  return(tirage);
} 

/* fin tirage_gamma */



/*!
 * @brief  Genere un tirage de variable aléatoire gaussienne 
 * @param[in] m : moyenne de la distribution
 * @param[in] sig: ecart-type de la distribution
 * @return tirage gaussien 
 * @author  J.Fade
 * @date    07/08/2007
 */
float tirage_gauss(float m, float sig)
{
  float tirage=0;
  double z1,dmoy;

  dmoy= (double) m;
   
  z1 = sqrt(-2.0*log( drand48() )); 
  tirage= (float) (dmoy + sig * z1 * cos(2*M_PI *drand48()));

  return(tirage);
} /* fin tirage_gauss */

/*!
 * @brief  Genere un tirage de variable aléatoire lognormale 
 * @detail La moyenne d'une telle v.a. vaut m*exp(sig^2/2) \n La variance vaut quant à elle m^2*exp(sig^2)*(exp(sig^2)-1). \n L'implantation et les valeurs des parametres sont tires de Evans, 'statistical distributions', p.129.
 * @param[in] m : médiane de la distribution
 * @param[in] sig: paramètre de forme la distribution
 * @return tirage lognormal 
 * @author J. Fade
 * @date   07/08/2007
 */
float tirage_lognormal(float m, float sig)
{
  float tirage=0;
 
  tirage= (float) m*exp(  sig*tirage_gauss(0,1));

  return(tirage);
} /* fin tirage_lognormal */


/******----------------------------------------------------------------
       tirage_gamma_gen(m,ordre,sortie)
       Descriptif      : Genere un tirage de pdf Gamma ordre non entier
       entree          : m		moyenne de la distribution
       : ordre		ordre de la loi Gamma
       
       retour          : tirage
       remarques       : ATTENTION, ceci ne vaut que pour des ordres entiers....
       Auteur          : J.Fade
- ----- -------------------------------------------------------------------*/

float tirage_gamma_gen(float m, float ordre)
{
  float tirage=0;
  float am,e,s,v1,v2,y;
  
  
  do {
    do {
      do {
	v1=drand48();
	v2=2.0*drand48()-1.0;
      } while (v1*v1+v2*v2 > 1.0);
      y=v2/v1;
      am=ordre-1;
      s=sqrt(2.0*am+1.0);
      tirage=s*y+am;
    } while (tirage <= 0.0);
    e=(1.0+y*y)*exp(am*log(tirage/am)-s*y);
  } while (drand48() > e);
  
  return(tirage*m/ordre);
} 

/* fin tirage_gamma_gen */





/*----------------------------------------------------------------
nom de la fonction : generation_sub_poisson_pseudo_photons(alpha,im_orig,larg,haut,sqz_factor,integ_time,im_bruit)
GENERATEUR DE BRUIT SUB-POISSONNIEN
 Descriptif            : cette fonction bruite une image de reference
			   par un bruit de sub-Poisson en utilisant un modèle
de sous photons tires avec une loi exponentielle. 
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


int generation_sub_poisson_pseudo_photons(float alpha,float **im_orig,int larg,int haut,
			    float sqz_factor, float integ_time, float **im_bruit)
     
{
  float 	       t,xm;
  int             i=0,j=0,k,m,hyp,cnt=0;
  


  for(i=0;i<haut;i++) 
    { 
      for(j=0;j<larg;j++)                        /* double boucle sur les pixels de l'image */
	{
	  k=0;t=0;

	 xm=alpha*im_orig[i][j]+1e-6;                /* xm = flux du process (sub) poisson */
  
	  if (sqz_factor >= 2 )                  /* Distinction pour q>1 des différents régimes de démarrage de la simul*/
	    {
	      hyp= (int)( (double)sqz_factor*drand48()); /* tirage de 'hyp': l'une des q hypothèses equiprobables    */
	                                                 /* i.e. entre 0 et q-1 pseudo-photons avant le 1er comptage */

		  for (m=0;m<=hyp;m++) 
		    { t += tirage_exponentiel((float) 1/xm/sqz_factor );    /* Tirage des 'hyp' pseudo photons avant */
		    }                        
                               /* premier comptage                      */
		  k=1;                                                      /* incré du tps tau=tps du 1er vrai phot */
	    
	    }

	  while (t< integ_time)       /* Comptage des arrivées de photons durant la fenetre integration T- tau*/
	    {
	      t += tirage_gamma(1/xm,(int) sqz_factor);    /* incrémentation du tps */
	      k++;                                   /* Incrémentation du process de comptage */ 
	    }

	  /*nb de phot détectés dans la plage T  en i,j */
	  im_bruit[i][j]=k-1;

	  	    cnt+=k-1;  //permet dans un but de vérification de compter le nombre total de photons générés.
	}
  
    
    }
 

  return(cnt);  
} /*  Fin generation_sub_poisson_pseudo_photons */






/*----------------------------------------------------------------
nom de la fonction : generation_sub_poisson_gamma(alpha,im_orig,larg,haut,sqz_factor,integ_time,im_bruit)
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

void generation_sub_poisson_gamma(float alpha,float **im_orig,int larg,int haut,
			    float sqz_factor, float integ_time, float **im_bruit)
     /* PREMIER ESSAI : RETROUVER LES COURBES OBTENUES EN GAMMA PRECEDENTES */
{
  float 	 t,xm,tir=0,gd_fenetre;
  float           dates[50000];
  int             i,j,k,cnt,cntd,kmax;
  


  for(i=0;i<haut;i++)      // ADAPTATION DE LA GRANDE FENETRE AU FLUX/* début correction avec dates */
    {
      for(j=0;j<larg;j++)            // Boucle sur chaque pixel
	{
	  
	  k=0;t=0;kmax=0;cnt=0;cntd=0;
	  xm=alpha*im_orig[i][j]+1e-6;    // xm= flux de photons à génerer
	  
	  gd_fenetre=2*integ_time;                     // Choix d'une fenêtre plus grande pour eviter synchronisations
	  if (5/xm > gd_fenetre) { gd_fenetre= 5/xm;
	  }

	  while (t< gd_fenetre) /* Comptage des arrivées de photons durant la fenetre integration k*T=integ_time */
	    {

	      tir=tirage_gamma(1/xm,sqz_factor);
	      t += tir;  k++;   /* incrémentation du tps */
	      dates[k]=t; /* enregistrement de la date dans tableau dates */
	    }
          kmax=k;
	  
	  tir= gd_fenetre/4+((double)(gd_fenetre*3/4-integ_time)*drand48()); //Selection d'une fenetre tempo au hasard
	  
          for(k=0;k<kmax;k++)              // Selection des dates d'arrivée convenables
	    {
	      if (dates[k]> tir && dates[k] < tir + integ_time ) {          // FENETRE HASARD
		cnt++;
	      }
	    }
	  im_bruit[i][j]=cnt;
	  
	} // next j
    } //next i
  
} /*  Fin generation_sub_poisson_gamma */


/*----------------------------------------------------------------
nom de la fonction : generation_sub_poisson_loiK(alpha,im_orig,larg,haut,sqz_factor,integ_time,im_bruit)
GENERATEUR DE BRUIT SUB-POISSONNIEN

UTILISEE POUR TESTER  LA STATISTIQUE DES TEMPS D'ATTENTES ENTRES ARRIVEES PHOTONIQUES

 Descriptif            : cette fonction bruite une image de reference
			   par un bruit de sub-Poisson temps d'attente Loi K.
 L'image a bruiter est transmise d'un seul bloc
! ! Les paramètres de la loi K ne sont tabulés que pour des squeezing d'ordre 2 3 4 5 10 ! ! !
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
void generation_sub_poisson_loiK(float alpha,float **im_orig,int larg,int haut,
			    float sqz_factor, float integ_time, float **im_bruit)

{
  float 	 t,xm,tir=0,gd_fenetre;
  float           dates[50000];
  int ord_L=0,ord_M=0;
  int             i,j,k,cnt,kmax;
  


  /* //////////////////////////// */
  /* // Cas special et temporaire des loi K */
  /* //  */
    if (sqz_factor > 1.99999 && sqz_factor < 2.000001) {ord_L=8 ; ord_M=3;}
    else  if (sqz_factor > 2.99999 && sqz_factor < 3.000001) {ord_L=15 ; ord_M=4;}
    else if  (sqz_factor > 3.99999 && sqz_factor < 4.000001) {ord_L=24 ; ord_M=5;}
    else if (sqz_factor > 4.99999 && sqz_factor < 5.000001) {ord_L=35 ; ord_M=6;}
    else if (sqz_factor > 9.99999 && sqz_factor < 10.000001) {ord_L=120 ; ord_M=11;}
    else {fprintf(stderr,"\n Parametres non tabulés pour un tel sqz_factor !"); exit(0);}
  /* //////////////////////////// */


  for(i=0;i<haut;i++)      // ADAPTATION DE LA GRANDE FENETRE AU FLUX/* début correction avec dates */
    {
      for(j=0;j<larg;j++)            // Boucle sur chaque pixel
	{
	  
	  k=0;t=0;kmax=0;cnt=0;
	  xm=alpha*im_orig[i][j]+1e-6;    // xm= flux de photons à génerer
	  
	  gd_fenetre=2*integ_time;                     // Choix d'une fenêtre plus grande pour eviter synchronisations
	  if (5/xm > gd_fenetre) { gd_fenetre= 5/xm;
	  }
	  
	  while (t< gd_fenetre) /* Comptage des arrivées de photons durant la fenetre integration k*T=integ_time */
	    {

	      tir=tirage_loi_K(1/xm,ord_M,ord_L);
	      t += tir;  k++;   /* incrémentation du tps */
	      dates[k]=t; /* enregistrement de la date dans tableau dates */
	    }
          kmax=k;
	  
	  tir= gd_fenetre/4+((double)(gd_fenetre*3/4-integ_time)*drand48()); //Selection d'une fenetre tempo au hasard
	  
          for(k=0;k<kmax;k++)              // Selection des dates d'arrivée convenables
	    {
	      if (dates[k]> tir && dates[k] < tir + integ_time ) {          // FENETRE HASARD
		cnt++;
	      }
	    }
	  im_bruit[i][j]=cnt;
	  
	} // next j
    } //next i
  
} /*  Fin generation_sub_poisson_loiK  */



/*----------------------------------------------------------------
nom de la fonction : generation_sub_poisson_weibull(alpha,im_orig,larg,haut,sqz_factor,integ_time,im_bruit)
GENERATEUR DE BRUIT SUB-POISSONNIEN
UTILISEE POUR TESTER  LA STATISTIQUE DES TEMPS D'ATTENTES ENTRES ARRIVEES PHOTONIQUES

 Descriptif            : cette fonction bruite une image de reference
			   par un bruit de sub-Poisson temps d'attente Weibull.
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

void generation_sub_poisson_weibull(float alpha,float **im_orig,int larg,int haut,
			    float sqz_factor, float integ_time, float **im_bruit)
     /* BRUITAGE AVEC DES TEMPS D'ATTENTE WEIBULL  */
{
  float 	 t,xm,tir=0,gd_fenetre;
  float           dates[50000];
  int             i,j,k,cnt,cntd,kmax;
  

  for(i=0;i<haut;i++)      // ADAPTATION DE LA GRANDE FENETRE AU FLUX/* début correction avec dates */
    {
      for(j=0;j<larg;j++)            // Boucle sur chaque pixel
	{
	  
	  k=0;t=0;kmax=0;cnt=0;cntd=0;
	  xm=alpha*im_orig[i][j]+1e-6;    // xm= flux de photons à génerer
	  
	  gd_fenetre=2*integ_time;                     // Choix d'une fenêtre plus grande pour eviter synchronisations
	  if (5/xm > gd_fenetre) { gd_fenetre= 5/xm;
	  }
	  
	  while (t< gd_fenetre) /* Comptage des arrivées de photons durant la fenetre integration k*T=integ_time */
	    {

	      tir=tirage_weibull(1/xm,sqz_factor);
	      t += tir;  k++;   /* incrémentation du tps */
	      dates[k]=t; /* enregistrement de la date dans tableau dates */
	    }
          kmax=k;
	  
	  tir= gd_fenetre/4+((double)(gd_fenetre*3/4-integ_time)*drand48()); //Selection d'une fenetre tempo au hasard
	  
          for(k=0;k<kmax;k++)              // Selection des dates d'arrivée convenables
	    {
	      if (dates[k]> tir && dates[k] < tir + integ_time ) {          // FENETRE HASARD
		cnt++;
	      }
	    }
	  im_bruit[i][j]=cnt;
	  
	} // next j
    } //next i
  
} /*  Fin generation_sub_poisson_weibull */





/*----------------------------------------------------------------
nom de la fonction : generation_sub_poisson_lognorm(alpha,im_orig,larg,haut,sqz_factor,integ_time,im_bruit)
GENERATEUR DE BRUIT SUB-POISSONNIEN
UTILISEE POUR TESTER  LA STATISTIQUE DES TEMPS D'ATTENTES ENTRES ARRIVEES PHOTONIQUES

 Descriptif            : cette fonction bruite une image de reference
			   par un bruit de sub-Poisson temps d'attente Lognormal.
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

void generation_sub_poisson_lognorm(float alpha,float **im_orig,int larg,int haut,
			    float sqz_factor, float integ_time, float **im_bruit)
     /* BRUITAGE AVEC DES TEMPS D'ATTENTE lognormal*/
{
  float 	  t,xm,tir=0,gd_fenetre;
  float           dates[50000];
  float           median=0,sig=0;
  int             i,j,k,cnt,kmax;
  

  for(i=0;i<haut;i++)      // ADAPTATION DE LA GRANDE FENETRE AU FLUX/* début correction avec dates */
    {
      for(j=0;j<larg;j++)            // Boucle sur chaque pixel
	{
	  
	  k=0;t=0;kmax=0;cnt=0;
	  xm=alpha*im_orig[i][j]+1e-6;    // xm= flux de photons à génerer

	  sig= sqrt(log(1+ 1/sqz_factor)); // Calcul des parametres de la loi lognormale
	  median = 1/xm*exp(-sig*sig/2);
	  
	  gd_fenetre=2*integ_time;                     // Choix d'une fenêtre plus grande pour eviter synchronisations
	  if (5/xm > gd_fenetre) { gd_fenetre= 5/xm;
	  }
	  
	  while (t< gd_fenetre) /* Comptage des arrivées de photons durant la fenetre integration k*T=integ_time */
	    {
	      tir= tirage_lognormal(median,sig);//(float) exp( tirage_Gauss(mu,sig)); 
	      t += tir;  k++;   /* incrémentation du tps */
	      dates[k]=t; /* enregistrement de la date dans tableau dates */
	    }
          kmax=k;
	  
	  tir= gd_fenetre/4+((double)(gd_fenetre*3/4-integ_time)*drand48()); //Selection d'une fenetre tempo au hasard
	  
          for(k=0;k<kmax;k++)              // Selection des dates d'arrivée convenables
	    {
	      if (dates[k]> tir && dates[k] < tir + integ_time ) {          // FENETRE HASARD
		cnt++;
	      }
	    }
	  im_bruit[i][j]=cnt;
	  
	} // next j
    } //next i
  
} /*  Fin generation_sub_poisson_lognorm */



/*----------------------------------------------------------------
nom de la fonction : generation_sub_poisson_cauchym(alpha,im_orig,larg,haut,sqz_factor,integ_time,im_bruit)

ATTENTION : N'A JAMAIS VRAIMENT MARCHE CAR NECESSITE ENORMEMENT DE TIRAGE POUR FAIRE LA DESYNCHRO

 Descriptif            : cette fonction bruite une image de reference
			   par un bruit de sub-Poisson temps d'attente Cauchy.
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
void generation_sub_poisson_cauchy(float alpha,float **im_orig,int larg,int haut,
			    float sqz_factor, float integ_time, float **im_bruit)
     /* BRUITAGE AVEC DES TEMPS D'ATTENTE cauchy*/
{
  float 	 t,xm,tir=0,gd_fenetre;
  float           dates[50000];
  int             i,j,k,cnt,cntd,kmax;

    FILE *tirage;
  
   tirage=fopen("waittime.dat","w");

  for(i=0;i<haut;i++)      // ADAPTATION DE LA GRANDE FENETRE AU FLUX/* début correction avec dates */
    {
      for(j=0;j<larg;j++)            // Boucle sur chaque pixel
	{
	  
	  k=0;t=0;kmax=0;cnt=0;cntd=0;
	  xm=alpha*im_orig[i][j]+1e-6;    // xm= flux de photons à génerer

	  
	  gd_fenetre=2*integ_time;                     // Choix d'une fenêtre plus grande pour eviter synchronisations
	  if (5/xm > gd_fenetre) { gd_fenetre= 5/xm;
	  }

	  while (t< gd_fenetre) /* Comptage des arrivées de photons durant la fenetre integration k*T=integ_time */
	    {
	      tir=tirage_cauchy(1/xm,2);
	      fprintf(tirage,"%f \n",tir);
	      
	      if (tir> 0){
		t += tir;  k++;   /* incrémentation du tps */
		//	      	fprintf(stderr,"%d ",k);
		dates[k]=t; /* enregistrement de la date dans tableau dates */
	      }

	    }
	  
          kmax=k;
	  
	  tir= gd_fenetre/4+((double)(gd_fenetre*3/4-integ_time)*drand48()); //Selection d'une fenetre tempo au hasard
	  
          for(k=0;k<kmax;k++)              // Selection des dates d'arrivée convenables
	    {
	      if (dates[k]> tir && dates[k] < tir + integ_time ) {          // FENETRE HASARD
		cnt++;
	      }
	    }
	  im_bruit[i][j]=cnt;
	  
	} // next j
    } //next i
  fprintf(stderr,"BREAK");
    fclose(tirage);
} /*  Fin generation_sub_poisson_cauchy*/







/*----------------------------------------------------------------
nom de la fonction : generation_sub_poisson_gamma_non_entier(alpha,im_orig,larg,haut,sqz_factor,integ_time,im_bruit)
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

void generation_sub_poisson_gamma_non_entier(float alpha,float **im_orig,int larg,int haut,
			    float sqz_factor, float integ_time, float **im_bruit)
     /* PREMIER ESSAI : RETROUVER LES COURBES OBTENUES EN GAMMA PRECEDENTES */
{
  float 	 t,xm,tir=0,gd_fenetre;
  float           dates[50000];
  int             i,j,k,cnt,cntd,kmax;
  


  for(i=0;i<haut;i++)      // ADAPTATION DE LA GRANDE FENETRE AU FLUX/* début correction avec dates */
    {
      for(j=0;j<larg;j++)            // Boucle sur chaque pixel
	{
	  
	  k=0;t=0;kmax=0;cnt=0;cntd=0;
	  xm=alpha*im_orig[i][j]+1e-6;    // xm= flux de photons à génerer
	  
	  gd_fenetre=2*integ_time;                     // Choix d'une fenêtre plus grande pour eviter synchronisations
	  if (5/xm > gd_fenetre) { gd_fenetre= 5/xm;
	  }

	  while (t< gd_fenetre) /* Comptage des arrivées de photons durant la fenetre integration k*T=integ_time */
	    {

	      tir=tirage_gamma_gen(1/xm,sqz_factor);
	      t += tir;  k++;   /* incrémentation du tps */
	      dates[k]=t; /* enregistrement de la date dans tableau dates */
	    }
          kmax=k;
	  
	  tir= gd_fenetre/4+((double)(gd_fenetre*3/4-integ_time)*drand48()); //Selection d'une fenetre tempo au hasard
	  
          for(k=0;k<kmax;k++)              // Selection des dates d'arrivée convenables
	    {
	      if (dates[k]> tir && dates[k] < tir + integ_time ) {          // FENETRE HASARD
		cnt++;
	      }
	    }
	  im_bruit[i][j]=cnt;
	  
	} // next j
    } //next i
  
} /*  Fin generation_sub_poisson_gamma_non_entier */



/* /\*----------------------------------------------------------------*\/ */
/* /\* nom de la fonction : generation_sub_poisson_fisher(alpha,im_orig,larg,haut,init,sqz_factor,integ_time,im_bruit,r) *\/ */
/* /\* GENERATEUR DE BRUIT SUB-POISSONNIEN *\/ */
/* /\* *\\/ *\/ */

/* void generation_sub_poisson_fisher(float alpha,float **im_orig,int larg,int haut, */
/* 			    float sqz_factor, float integ_time, float **im_bruit, gsl_rng *r) */
/*      /\* BRUITAGE AVEC DES TEMPS D'ATTENTE cauchy*\/ */
/* { */
/*   float 	 t,xm,tir=0,gd_fenetre; */
/*   float           dates[50000]; */
/*   int             i,j,k,cnt,cntd,kmax; */
/*   double          test; */
  
/*   FILE *tirage; */
  
/*   tirage=fopen("waittime.dat","w"); */
  
/*   for(i=0;i<haut;i++)      // ADAPTATION DE LA GRANDE FENETRE AU FLUX/\* début correction avec dates *\/ */
/*     { */
/*       for(j=0;j<larg;j++)            // Boucle sur chaque pixel */
/* 	{ */
	  
/* 	  k=0;t=0;kmax=0;cnt=0;cntd=0;  */
	  
/* 	  xm=alpha*im_orig[i][j];    // xm= flux de photons à génerer */
	  
	  
/* 	  gd_fenetre=2*integ_time;                     // Choix d'une fenêtre plus grande pour eviter synchronisations */
/* 	  if (5/xm > gd_fenetre) { gd_fenetre= 5/xm; */
/* 	  } */
	  
/* 	  while (t< gd_fenetre) /\* Comptage des arrivées de photons durant la fenetre integration k*T=integ_time *\/ */
/* 	    {   */
/* 	      tir =(float)  1.0/( (double) xm )  *  gsl_ran_fdist ( r , 2.0, 4.0); */
/* 	      fprintf(stderr,"TEST %f \n",tir); */
/* 	      fprintf(tirage,"%f \n",tir); */
	      
/* 	      t += tir;  k++;   /\* incrémentation du tps *\/ */
	      
/* 	      dates[k]=t; /\* enregistrement de la date dans tableau dates *\/ */
	      
	      
/* 	    } */
	  
/*           kmax=k; */
	  
/* 	  tir= gd_fenetre/4+((double)(gd_fenetre*3/4-integ_time)*drand48()); //Selection d'une fenetre tempo au hasard */
	  
/*           for(k=0;k<kmax;k++)              // Selection des dates d'arrivée convenables */
/* 	    { */
/* 	      if (dates[k]> tir && dates[k] < tir + integ_time ) {          // FENETRE HASARD */
/* 		cnt++; */
/* 	      } */
/* 	    } */
/* 	  im_bruit[i][j]=cnt; */
/* 	  fprintf(stderr," k= %d ",cnt); */
	  
/* 	} // next j */
/*     } //next i */
/*   fprintf(stderr,"BREAK"); */
/*   fclose(tirage); */
/* } /\*  Fin generation_sub_poisson_fisher*\/ */














