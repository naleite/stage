#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utils_fade.h"
#include "utils_tirage.h"
#include <unistd.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

#define MAXIT 1000
//Fonctions utils rajoutees.
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}
//fin fonctions utils
//==========================================================
//==========================================================

/* Generation d'image // et |_ (image_X et image_Y)
 Entrees:	Nw,Nw_: nombre de pixels que w ou w_ a,
 Imoy_w:	Intensite moyenne de w
 Imoy_w_:---------------------w_,
 DOP_w,	DOP_w_: DOP de w et w_,
 reals:	Le nombre de realisations
 ordre:	Ordre
 Sorties: deux tableau R*(Nw+Nw_) qui contient les donnees generees pour I// et I|_ (image_X et image_Y.)
 */
void generer_image_XY(int Nw, double Imoy_w, double DOP_w, int Nw_, double Imoy_w_,double DOP_w_, int ordre, int reals, double **image_X, double **image_Y)

{
	printf("Generer images X et Y en cours: Nw= %d, Nw_=%d,Imoy_w et w_=%f %f,DOP_w et w_= %f %f, ordre=%d, nb reals=%d ...\n",Nw,Nw_,Imoy_w,Imoy_w_,DOP_w,DOP_w_,ordre,reals);
	//boucle pour generation d'image:
	for(int i=0;i<reals;i++){
		//Boucle pour w:
		for (int j=0;j<Nw;j++){
			image_X[i][j]=tirage_gamma((1+DOP_w)*Imoy_w/2, ordre); //I//
			image_Y[i][j]=tirage_gamma((1-DOP_w)*Imoy_w/2, ordre); //I|_
		}
		
		//Boucle pour w_
		for (int j=0;j<Nw_;j++){
			image_X[i][j+Nw]=tirage_gamma((1+DOP_w_)*Imoy_w_/2, ordre); //I//
			image_Y[i][j+Nw]=tirage_gamma((1-DOP_w_)*Imoy_w_/2, ordre); //I|_
		}
		
	}
	printf("Generation des images X et Y reussit!\n\n");
}

/*
 Generation image Atot, A_ICEO, A_delta
 Entree: image_X,image_Y: Les images de ?? et |_
 Sortie: 3 images: A_tot, A_ICEO,A_delta
 */
void generer_3image(int taille, int reals,double **image_X,double **image_Y, double **image_A_tot, double **image_A_ICEO, double **image_A_delta){
	printf("Generer images en cours: la taille=%d, nb de reals=%d... \n",taille,reals);
	for(int i=0;i<reals;i++){
		for (int j=0; j<taille; j++) {
			
			double x=image_X[i][j];
			double y=image_Y[i][j];
			image_A_tot[i][j]=x+y;
			image_A_ICEO[i][j]=(x-y)/(x+y);
			image_A_delta[i][j]=atanh(image_A_ICEO[i][j]);
			
			//printf("Image_tot[%d][%d]:\t %f \n",i,j,image_A_tot[i][j]);
			//printf("Image_ICEO[%d][%d]:\t %f \n",i,j,image_A_ICEO[i][j]);
			//printf("Image_delta[%d][%d]:\t %f \n",i,j,image_A_delta[i][j]);
		}
		
	}
	printf("Fin de generation des 3 images\n\n");
}//fin de generation 3 images

//==================FONCTIONS DE TEST DE DETECTION DE TOT ET DELTA==================================

/*test_tot: Pour les tests de detections R_Itot
 Entrees:	tab_tot: tableau de Image generer totales de R* (Nw+Nw_),
 reals: nombre de realisations
 Nw:		Nombre de pixels de w
 Nw_		Nombre de pixels de w_
 Sorties: tableau test_tot de R
 */

double *test_tot(int reals, int nbw, int nbw_, double **tab_tot){
	printf("Test de detection: Intensite totale\n");
	double *tab_res=d_alloue_1d(reals);
	double *temp_w=d_alloue_1d(nbw);
	double *temp_w_=d_alloue_1d(nbw_);
	
	for (int i=0;i<reals;i++){
		//Separation de w et w_
		for(int j=0;j<nbw;j++){
			temp_w[j]=tab_tot[i][j];
		}
		for(int j=0;j<nbw_;j++){
			temp_w_[j]=tab_tot[i][j+nbw];
		}
		//Calculer le moy de w et w_
		double moy_w=d_moyenne(nbw, temp_w);
		double moy_w_=d_moyenne(nbw_, temp_w_);
		
		//calculer le test detection de I_tot
		tab_res[i]=(moy_w-moy_w_)*(moy_w-moy_w_);
		//printf("Test_tot[%d]:  %f \n",i,tab_res[i]);
	}
	printf("Test de detection1 se termine.\n\n");
	return	tab_res;
}

/*test_delta: Pour les tests de detections R_delta
 Entrees:	tab_tot: tableau de Image generer totales de R* (Nw+Nw_),
 reals:		nombre de realisations
 Nw:		Nombre de pixels de w
 Nw_		Nombre de pixels de w_
 Sorties: tableau test_reals de taille R
 */
double *test_delta(int reals, int nbw, int nbw_, double **tab_delta){
	printf("Test de detection: Delta...\n");
	//double *tab_res;//=d_alloue_1d(reals);
	double *tab_res=d_alloue_1d(reals);
	double A=1;//nbw*nbw_/(nbw_+nbw);
	for (int i=0; i<reals; i++) {
		double imoy_w=moyenne(0, nbw, tab_delta[i]);
		double imoy_w_=moyenne(nbw_, nbw+nbw_, tab_delta[i]);
		tab_res[i]=A*(imoy_w-imoy_w_)*(imoy_w-imoy_w_);
		//printf("Test_delta[%d]:  %f \n",i,tab_res[i]);
	}
	printf("Test de detection2 se termine.\n\n");
	return	tab_res;
}

//=============TEST DE DETECTION LRT=============================================
//===============================================================================

/*log_Proba_p(ro,p)----Fonction utile qui calcule le ln_proba de p(i,j)
	Entrees:	ro: les donnees simulees
	p:  fixee
	Sorites:	double res: les log_probas de donnees ro
 */
double log_Proba_p(double ro,double p){
	//printf("args ro=%f, p=%f",ro,p);
	double proba=((1-p*p)/2)*(1/(1-p*ro)/(1-p*ro));
	//printf("proba: %f\n",proba);
	double res=log(proba);
	//printf("res: %f\n",res);
	return res;
}

/*test_LRT, cette fonction calcule le LRT sur p, ainsi que s
	Entrees:	tab_p: tableau de Image genere totales de R* (Nw+Nw_),
	reals:		nombre de realisations
	Nw:		Nombre de pixels de w
	Nw_:	Nombre de pixels de w_
	Sorties: tableau test_LRT de R
 */
double *test_LRT(double DOP_F, double DOP_w_, double DOP_w,int reals, int nbw, int nbw_, double **tab){
	printf("Test de detection: LRT_p...\n");
	double *tab_res=d_alloue_1d(reals);
	//double *temp_h1_w=d_alloue_1d(nbw);
	//double *temp_h0_w_=d_alloue_1d(nbw_);
	
	for (int i=0; i<reals; i++) {
		
		double sum_pr_h1_w=0.;
		double sum_pr_h0_w_=0.;
		double sum_pr_h0_F=0.;//
		
		for(int j=0;j<nbw;j++){
			//temp_h1_w[j]=tab1[i][j];
			sum_pr_h1_w+=log_Proba_p(tab[i][j], DOP_w);
			//printf("%f\n",sum_pr_h1_w);
		}
		for(int j=0;j<nbw_;j++){
			//temp_h0_w_[j]=tab0[i][j+nbw];
			sum_pr_h0_w_+=log_Proba_p(tab[i][j+nbw], DOP_w_);
			//printf("%f\n",sum_pr_h0_w_);
		}
		
		for (int j=0; j<nbw+nbw_; j++) {
			sum_pr_h0_F+=log_Proba_p(tab[i][j], DOP_F);
		}
		tab_res[i]=sum_pr_h1_w+sum_pr_h0_w_-sum_pr_h0_F;
		//printf("LRT[%d]: %f \n",i,tab_res[i]);
	}
	printf("Test de detection3 se termine.\n\n");
	//free(temp_h0_w_);
	//free(temp_h1_w);
	return tab_res;
}

/*test_GLRT, cette fonction calcule le GLRT sur p
 Entrees:	tab_p: tableau de Image genere totales de R* (Nw+Nw_),
 reals:		nombre de realisations
 Nw:		Nombre de pixels de w
 Nw_:	Nombre de pixels de w_
 Sorties: tableau test_LRT de R
 */
double *test_GLRT(int reals, int nbw, int nbw_, double **tab){
	printf("Test de detection: GLRT_p...\n");
	double *tab_res=d_alloue_1d(reals);
	double tmp_DOP=0.;
	//double *temp_h1_w=d_alloue_1d(nbw);
	//double *temp_h0_w_=d_alloue_1d(nbw_);
	
	for (int i=0; i<reals; i++) {
		
		double sum_pr_h1_w=0.;
		double sum_pr_h0_w_=0.;
		double sum_pr_h0_F=0.;//
		
		tmp_DOP = moyenne(0,nbw,tab[i]);
		
		for(int j=0;j<nbw;j++){
			//temp_h1_w[j]=tab1[i][j];
			sum_pr_h1_w+=log_Proba_p(tab[i][j],tmp_DOP);
			//printf("%f\n",sum_pr_h1_w);
		}
		
		tmp_DOP = moyenne(nbw,nbw+nbw_,tab[i]);

		for(int j=0;j<nbw_;j++){
			//temp_h0_w_[j]=tab0[i][j+nbw];
			sum_pr_h0_w_+=log_Proba_p(tab[i][j+nbw], tmp_DOP);
			//printf("%f\n",sum_pr_h0_w_);
		}
		
		tmp_DOP = d_moyenne(nbw+nbw_,tab[i]);
		
		for (int j=0; j<nbw+nbw_; j++) {
			sum_pr_h0_F+=log_Proba_p(tab[i][j], tmp_DOP);
		}
		tab_res[i]=sum_pr_h1_w+sum_pr_h0_w_-sum_pr_h0_F;
		//printf("LRT[%d]: %f \n",i,tab_res[i]);
	}
	
	//free(temp_h0_w_);
	//free(temp_h1_w);
	printf("Test de detection4 se termine.\n\n");
	return tab_res;
}

//======================Test de detection de LRTs==========================
/*log_Proba_up
	Entrees:	s: s(i,j), la donnees du tableau_tot
				u:	Ui, moyenne d'intensite
				DOP: degres du polar
	Sorties:	le resultat de log-vraisemblemse.
 */
double log_Proba_up(double s,double u,double DOP){
	//printf("s=%f, u=%f, DOP=%f \n\n",s,u,DOP);
	double proba;
	if (DOP==0.) {
		proba=(2/u)*(2/u)*s*exp(-2*s/u);
		//printf("DOP==0 \n");
	}
	
	else {
		proba=(1/DOP/u)*(exp(-2*s/(1+DOP)/u)-exp(-2*s/(1-DOP)/u));
	}
	//printf("proba_up: %f\n",proba);
	double res=log(proba);
	//printf("log_up_res: %f\n",res);
	return res;

}

/*test_LRT_s,	Cette fonction calcule le LRT_s
	Entrees:	moy_h1_w, Moyenne d'intensite de w de H1
				moy_h0, Moyenne d'intensite de w_ de H0
				DOP_h1_w, DOP de w de H1
				DOP_h0, DOP de h0
				nbw,nbw_: nombre de w et w_
				reals, le nombre de realisation
				**tab, Les donnees simulees.
	Sorties:	double *tab_res, les resultat de test de detection LRTs
 
 */
double *test_LRT_s(double moy_h1_w,double moy_h0,double DOP_h1_w,double DOP_h0,int reals,int nbw,int nbw_,double **tab){
	printf("Test de detection: LRT_s...\n");
	double *tab_res=d_alloue_1d(reals);
	
	for (int i=0; i<reals; i++) {
		
		double sum_pr_h1_w=0.;
		double sum_pr_h0_w_=0.;
		double sum_pr_h0_F=0.;//
		
		for(int j=0;j<nbw;j++){
			//temp_h1_w[j]=tab1[i][j];
			sum_pr_h1_w+=log_Proba_up(tab[i][j],moy_h1_w, DOP_h1_w);
			//printf("%f\n",sum_pr_h1_w);
		}
		for(int j=0;j<nbw_;j++){
			//temp_h0_w_[j]=tab0[i][j+nbw];
			sum_pr_h0_w_+=log_Proba_up(tab[i][j+nbw], moy_h0, DOP_h0);
			//printf("%f\n",sum_pr_h0_w_);
		}
		
		for (int j=0; j<nbw+nbw_; j++) {
			sum_pr_h0_F+=log_Proba_up(tab[i][j],moy_h0, DOP_h0);
		}
		tab_res[i]=sum_pr_h1_w+sum_pr_h0_w_-sum_pr_h0_F;
		//printf("LRT[%d]: %f \n",i,tab_res[i]);
	}

	printf("Test de detection5 se termine.\n\n");
	return tab_res;
}
//==============================FIN de test_LRT_s===================================================

//=============================Fonctions de GLRT_S=================================
/*calculer beta2 estime*/

double estime_p(double estime_u, double var_tot){
	double beta2=(2*var_tot/estime_u/estime_u)-1;
	if (beta2<0.) {
		beta2=0.;
	}
	else if(beta2>1){
		beta2=1;
	}
	double res=sqrt(beta2);
	
	return res;
}

/*test_GLRT, cette fonction calcule le GLRT sur p, ainsi que s
 Entrees:	tab_p: tableau de Image genere totales de R* (Nw+Nw_),
 reals:		nombre de realisations
 Nw:		Nombre de pixels de w
 Nw_:	Nombre de pixels de w_
 tab_tot: tableau tot
 tab_p : tableau ICEO
 Sorties: tableau test_LRT de R
 */
double *test_GLRT_s(int reals, int nbw, int nbw_, double **tab_tot){
	printf("Test de detection: GLRT_s...\n");
	double *tab_res=d_alloue_1d(reals);
	double tmp_DOP=0.;
	double tmp_u=0.;
	double tmp_var=0.;
	
	for (int i=0; i<reals; i++) {
		
		double sum_pr_h1_w=0.;
		double sum_pr_h0_w_=0.;
		double sum_pr_h0_F=0.;//
		
		
		tmp_u=moyenne(0,nbw,tab_tot[i]);
		tmp_var = variance(0,nbw,tab_tot[i],tmp_u);
		tmp_DOP = estime_p(tmp_u,tmp_var);
		
		for(int j=0;j<nbw;j++){
			//temp_h1_w[j]=tab1[i][j];
			sum_pr_h1_w+=log_Proba_up(tab_tot[i][j],tmp_u,tmp_DOP);
			//printf("%f\n",sum_pr_h1_w);
		}
		

		tmp_u=moyenne(nbw,nbw+nbw_,tab_tot[i]);
		tmp_var = variance(nbw,nbw+nbw_,tab_tot[i],tmp_u);
		tmp_DOP = estime_p(tmp_u,tmp_var);
		
		for(int j=0;j<nbw_;j++){
			//temp_h0_w_[j]=tab0[i][j+nbw];
			sum_pr_h0_w_+=log_Proba_up(tab_tot[i][j+nbw],tmp_u, tmp_DOP);
			//printf("%f\n",sum_pr_h0_w_);
		}
		
		
		tmp_u=d_moyenne(nbw+nbw_,tab_tot[i]);
		tmp_var = d_variance(nbw+nbw_,tab_tot[i],tmp_u);
		tmp_DOP = estime_p(tmp_u,tmp_var);
		
		for (int j=0; j<nbw+nbw_; j++) {
			sum_pr_h0_F+=log_Proba_up(tab_tot[i][j],tmp_u, tmp_DOP);
		}
		
		
		tab_res[i]=sum_pr_h1_w+sum_pr_h0_w_-sum_pr_h0_F;
		//printf("LRT[%d]: %f \n",i,tab_res[i]);
	}
	
	//free(temp_h0_w_);
	//free(temp_h1_w);
	printf("Test de detection5 se termine.\n\n");
	return tab_res;
}

//============================FIN du test GLRTs=====================================================


//test_var_gauss
double *test_var_gausse(int reals, int nbw, int nbw_, double **tab_tot){
	double *tab_res=d_alloue_1d(reals);
	
	for (int i=0; i<reals; i++) {
		double var_w=variance(0, nbw, tab_tot[i], moyenne(0, nbw, tab_tot[i]));
		double var_w_=variance(nbw, nbw+nbw_, tab_tot[i], moyenne(nbw, nbw+nbw_, tab_tot[i]));
		double var_F=d_variance(nbw+nbw_, tab_tot[i], d_moyenne(nbw+nbw_, tab_tot[i]));
		
		tab_res[i]=-nbw/2*log(var_w)-nbw_/2*log(var_w_)+(nbw+nbw_)/2*log(var_F);
	}
	
	return tab_res;
}

//test detection sur intensite
double *test_intensite(int reals, int nbw, int nbw_, double **tab_tot){
	double *tab_res=d_alloue_1d(reals);
	double A=nbw*nbw_/(nbw_+nbw);
	for (int i=0; i<reals; i++) {
		double imoy_w=moyenne(0, nbw, tab_tot[i]);
		double imoy_w_=moyenne(nbw_, nbw+nbw_, tab_tot[i]);
		tab_res[i]=A*(imoy_w-imoy_w_)*(imoy_w-imoy_w_);
	}
	return tab_res;
}
//fin test varairance gaussian

//test_gausse_multi
double *test_gausse_multi(int reals, int nbw,int nbw_, double **tab_test){
	double *tab_res=d_alloue_1d(reals);
	double tmp_u_w;
	double tmp_v_w;
	double tmp_u_w_;
	double tmp_v_w_;
	double tmp_u_F;
	double tmp_v_F;
	for (int i=0; i<reals; i++) {
		tmp_u_w=moyenne(0, nbw, tab_test[i]);
		tmp_v_w=variance(0, nbw, tab_test[i], tmp_u_w);//calculer moyenne et variance de w
		
		tmp_u_w_=moyenne(nbw, nbw+nbw_, tab_test[i]);
		tmp_v_w_=variance(nbw, nbw+nbw_, tab_test[i], tmp_u_w_);//calculer moyenne et variance de w_
		
		tmp_u_F=d_moyenne(nbw+nbw_, tab_test[i]);
		tmp_v_F=d_variance(nbw+nbw_, tab_test[i],tmp_u_F);//calculer moyenne et variance de F
		
		double tmp1=tmp_v_w/2/tmp_u_w/tmp_u_w-nbw*log(tmp_u_w);
		double tmp2=tmp_v_w_/2/tmp_u_w_/tmp_u_w_-nbw_*log(tmp_u_w_);
		double tmp3=tmp_v_F/2/tmp_u_F/tmp_u_F-(nbw+nbw_)*log(tmp_u_F);//caculer chaque
		
		tab_res[i]=tmp1+tmp2-tmp3;//resultat de la ligne.
	}
	
	return tab_res;
}//fin test gausse multi

//log_moyenne: sert la fonction test_diff_logmoment et test_quot_logmoment
double log_moyenne(int debut, int fin, double *tab){
	double res;
	double sum=0;
	int card=fin-debut;
	for (int i=debut;i<fin;i++){
		sum+=log(tab[i]);
	}
	res=sum/card;
	return res;
}

//test_diff_logmoment
double *test_diff_logmoment(int reals, int nbw,int nbw_, double **tab_test){
	double *tab_res=d_alloue_1d(reals);
	for (int i=0; i<reals; i++) {
		double lm_w=log_moyenne(0, nbw, tab_test[i]);
		double lm_w_=log_moyenne(nbw, nbw+nbw_, tab_test[i]);
		//printf("lm_w=%f,\nlm_w_=%f.\n",lm_w,lm_w_);
		tab_res[i]=(lm_w-lm_w_)*(lm_w-lm_w_);
	}
	
	return tab_res;
}//fin test_diff_logmoment

//test_quot_logmoment
double *test_quot_logmoment(int reals, int nbw,int nbw_, double **tab_test){
	double *tab_res=d_alloue_1d(reals);
	for (int i=0; i<reals; i++) {
		double lm_w=log_moyenne(0, nbw, tab_test[i]);
		double lm_w_=log_moyenne(nbw, nbw+nbw_, tab_test[i]);
		tab_res[i]=(lm_w/lm_w_)*(lm_w/lm_w_);
	}
	
	return tab_res;
}//fin test_quot_logmoment



//==========================Calculer AUC============================================================
double calculer_auc(double *vect_x,double *vect_y,int nb){
	double aire=0.;
	for (int i=0;i<nb-1; i++) {
		aire+=((vect_y[i]+vect_y[i+1])*(vect_x[i]-vect_x[i+1])/2);
		//printf("aire=%f\n",aire);
	}
	aire+=(vect_x[nb-1]*vect_y[nb-1])/2;
	return aire;
}


//========TEST DE DETECTION---LRT_GAMMA=============================
double d_sum_index(int start,int end,double *tab){
	double sum=0.;
	for (int i=start; i<end; i++) {
		sum+=tab[i];
	}
	
	return sum;
}
double d_lnsum_index(int start,int end,double *tab){
	double sum=0.;
	for (int i=start; i<end; i++) {
		sum+=log(tab[i]);
	}
	
	return sum;
}

//param l: Ordre de la simulation.
double *test_LRT_gamma(int reals,int l, double h1_moy_w, double h0_moy_w_ ,int nbw,int nbw_, double **tab_test){
	double *res;
	res=d_alloue_1d(reals);
	for (int i=0; i<reals; i++) {
		
		double item_w=0.;
		double item_w_=0.;
		double item_F=0.;
		
		item_w = -d_sum_index(0, nbw, tab_test[i])*l/h1_moy_w + (l-1)*d_lnsum_index(0, nbw, tab_test[i]) - nbw*l*log(h1_moy_w/l) - nbw*gsl_sf_lngamma(l);
		item_w_ = -d_sum_index(nbw, nbw+nbw_, tab_test[i])*l/h0_moy_w_ + (l-1)*d_lnsum_index(nbw, nbw+nbw_, tab_test[i]) - nbw_*l*log(h0_moy_w_/l) - nbw_*gsl_sf_lngamma(l);
		item_F = -d_sum_index(0, nbw+nbw_, tab_test[i])*l/h0_moy_w_ + (l-1)*d_lnsum_index(0, nbw+nbw_, tab_test[i]) - (nbw+nbw_)*l*log(h0_moy_w_/l) - (nbw+nbw_)*gsl_sf_lngamma(l);
		
		res[i]=item_w+item_w_-item_F;
	}
	return res;
		
		
	}
//======fin de la fonction LRT_GAUSSIAN
//======Routines which returns the function value and the first derivative of the function
double estime_u_et_y(int start,int end,double *tab){
	double y=0.0;
	double es_u=0.0;
	double sum=0.0;
	int nb=end-start;
	//Estimer u
	for (int i=start;i<end;i++){
		sum+=tab[i];
	}
	es_u=sum/nb;
	//Calculer 1/N*sigma(ln Si)
	double log_sum=0.0;
	for (int i=start;i<end;i++){
		log_sum+=log(tab[i]);
	}
	y=log(es_u)-log_sum/nb;
	return y;
}



void estime_l(double l, double y,double *res,double *deriv){
	double digramme=gsl_sf_psi(l);
	double trigramme=gsl_sf_psi_1(l);
	double r=log(l)-digramme-y;
	double d=1/l-trigramme;
	*res=r;
	*deriv=d;
	//printf("res=%f, deriv=%f\n",r,d);
}
//========NEWTON METHOD=================

/*Using a combination of Newton-Raphson and bisection, find the root of a function bracketed between x1 and x2. The root, returned as the function value rtsafe, will be refined until its accuracy is known within ±xacc. funcd is a user-supplied routine that returns both the function value and the first derivative of the function.*/
//double Newton(void (*funcd)(double, double, double *, double *), double x1, double x2, double xacc)
double Newton(double l0,double y,double x1, double x2, double xacc)
{
	void nrerror(char error_text[]);
	int j;
	double df=0,dx=0,dxold=0,f=0,fh=0,fl=0;
	double temp,xh,xl,rts;
	//(*funcd)(x1,y,&fl,&df);
	estime_l(x1, y, &fl, &df);
	estime_l(x2, y, &fh, &df);
	//printf("fl=%f,df=%f,fh=%f\n",fl,df,fh);
	//(*funcd)(x2,y,&fh,&df);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
	{
		printf("fl=%f,fh=%f\n\n",fl,fh);
		nrerror("Root must be bracketed in rtsafe");}
	if (fl == 0.0)
		return x1;
	if (fh == 0.0)
		return x2;
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	//rts=0.5*(x1+x2);	/*Initialize the guess for root, the “stepsize before last,” and the last step.*/
	//rts=(3-y+sqrt((y-3)*(y-3)+24*y))/12/y;
	rts=l0;
	dxold=fabs(x2-x1);
	dx=dxold;
	estime_l(rts,y,&f,&df);
	for (j=1;j<=MAXIT;j++) {	//Loop over allowed iterations.

		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)|| (fabs(2.0*f) > fabs(dxold*df))){ // or not decreasing fast enough.
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts)
				return rts;
		}
		else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts)
				return rts;
			//Change in root is negligible. Newton step acceptable. Take it.
		}
		if (fabs(dx) < xacc)
			return rts;
		estime_l(rts,y,&f,&df);
		if (f < 0.0)// Maintain the bracket on the root.
			xl=rts;
		else
			xh=rts;
	}
	nrerror("Maximum number of iterations exceeded in rtsafe");
	
	return INFINITY;// Never get here.return inf. ATTENTION!!!!!!ICI!!!!!!! RETURN INF? OU 0.0?????
}

//=======FIN NEWTON METHOD==============

//=======GLRT-gamma====================
double *test_GLRT_gamma(int reals, int nbw, int nbw_, double **tab_tot,double x1,double x2,double xacc){
	double *tab_res=d_alloue_1d(reals);
	double es_uw=0.0;
	double es_uw_=0.0;
	double es_uF=0.0;
	double Lw;
	double Lw_;
	double LF;
	for (int i=0; i<reals; i++) {
		es_uw=d_sum_index(0, nbw, tab_tot[i])/nbw;
		es_uw_=d_sum_index(nbw, nbw+nbw_, tab_tot[i])/nbw_;
		es_uF=d_sum_index(0, nbw+nbw_, tab_tot[i])/(nbw+nbw_);
		double yw=estime_u_et_y(0, nbw, tab_tot[i]);
		double yw_=estime_u_et_y(nbw, nbw+nbw_, tab_tot[i]);
		double yF=estime_u_et_y(0, nbw+nbw_, tab_tot[i]);
		double lw0=(3-yw+sqrt((yw-3)*(yw-3)+24*yw))/12/yw;
		double lw_0=(3-yw+sqrt((yw_-3)*(yw_-3)+24*yw_))/12/yw_;
		double lF0=(3-yF+sqrt((yF-3)*(yF-3)+24*yF))/12/yF;
		//double resw,resw_,resF,derivw,derivw_,derivF;
		Lw=Newton(lw0,yw,x1, x2, xacc);
		Lw_=Newton(lw_0,yw_,x1, x2, xacc);
		LF=Newton(lF0,yF,x1, x2, xacc);
		//printf("\nLw=%f,Lw-=%f,LF=%f\n",Lw,Lw_,LF);
		double item_w=0.;
		double item_w_=0.;
		double item_F=0.;
		
		item_w = -nbw*Lw+ (Lw-1)*d_lnsum_index(0, nbw, tab_tot[i]) - nbw*Lw*log(es_uw/Lw) - nbw*gsl_sf_lngamma(Lw);
		
		item_w_ = -nbw_*Lw_ + (Lw_-1)*d_lnsum_index(nbw, nbw+nbw_, tab_tot[i]) - nbw*Lw_*log(es_uw_/Lw_) - nbw_*gsl_sf_lngamma(Lw_);
		
		item_F =-(nbw+nbw_)*LF + (LF-1)*d_lnsum_index(0, nbw+nbw_, tab_tot[i]) - (nbw+nbw_)*LF*log(es_uF/LF) - (nbw+nbw_)*gsl_sf_lngamma(LF);
		
		tab_res[i]=item_w+item_w_-item_F;


	}
	
	
	return tab_res;
}

//=======================FONCTION GENERATION DE COURBE COR==========================================
//==================================================================================================

//Declaration de fonction d_min et d_min_1d qui calcule le min de tableau.
double d_min( double a, double b)
{
	double ret=0.;
	
	if (a<b){
		ret=a; }
	else{
		ret=b; }
	return ret;
}

double d_min_1d( double *a, int dim)
{
	double ret;
	int k=0;
	
	ret=a[0];
	
	for(k=1;k<dim;k++)
	{
		ret= d_min( a[k], ret);
	}
	return ret;
}

//Calcule histogrmme en fonction de donnees.
void calcul_histogram(double *vecteur, int size, double min, double max, int pt_nb_bin, double *histo)
{
	int k,j;
	double delta=0;
	
	j=0;
	
	delta = (double) (max -min) / pt_nb_bin;
	// fprintf(stderr,"\n  %g ",delta);
	// printf("\n %d ",pt_nb_bin);printf("\n %g ",min);printf("\n %g ",max);
	
	for (k=0;k< size;k++)
    {
		vecteur[k] =( vecteur[k] - min)/delta;
		     //printf(stderr,"\n %g",vecteur[k]);
    }
	
	for (k=0;k< size;k++)
    {
		j= (int) vecteur[k];
		histo[j]++;
    }
	
} /* end of calcul_histogram */

//Nouvelle fonction de Calculer_pd_pfa
void  calculer_Pd_Pfa(int nb_niv, double *tab_P, double *tab_Q, double *hist_P, double *hist_Q, int nb_pts)
{
	int k=0;
	
	/* Normalisation des vecteurs pour avoir Pfa et Pd max de 1 */
	for(k=0;k<nb_niv;k++)
	{
		tab_P[k]  =  hist_P[k]/nb_pts;
		tab_Q[k]  =  hist_Q[k]/nb_pts;
		//fprintf(stderr,"(%f) %f", tab_P[k], tab_Q[k]);
	}
	/* Construction des courbes Pf Pfa */
	
	for(k=2;k<nb_niv+1;k++)
	{
		tab_P[nb_niv - k]  +=  tab_P[nb_niv+1-k];
		tab_Q[nb_niv - k]  +=  tab_Q[nb_niv+1-k];
		
	}
	
}

/*Generer_COR
 Entrees:	h0_test et h1_test: le tableau qui contient les resultats de test de detections.
 reals:				le nombre de realisations
 nom_fichier[]: String de nom du fichier
 Sortie: le courbe COR
 */
void generer_COR(double *h0_test,double *h1_test,int reals,char nom_fichier[]){
	printf("Generation courbe COR en cours...\n");
	//Calculer pt_nb_bin
	int pt_nb_bin = (int) reals/10;
	double Vmax = d_max( d_max_1d(h0_test, reals), d_max_1d(h1_test, reals));
	double Vmin = d_min( d_min_1d(h0_test, reals), d_min_1d(h1_test, reals));
	double taille_bin;
	
	taille_bin = (Vmax - Vmin)/ (double) pt_nb_bin;
	
	//printf("nb bin:  %d \n",pt_nb_bin);
	//printf("Val min:  %f \n",Vmin);
	//printf("Val max:  %f \n",Vmax);
	//printf("Taille bin:  %f \n",taille_bin);
	
	//Allouer les memoires
	double *histo_h0=d_alloue_1d(pt_nb_bin);
	double *histo_h1=d_alloue_1d(pt_nb_bin);
	double *vect_Pfa=d_alloue_1d(pt_nb_bin);
	double *vect_Pd=d_alloue_1d(pt_nb_bin);
	double *abscisse=d_alloue_1d(pt_nb_bin+1);
	
	//Creation abscisse histo
	for (int k=0; k<=pt_nb_bin; k++)
    {
		abscisse[k]=Vmin + (double) k * taille_bin;
    }
	
	
	//Calculer de l'histo de H0
	calcul_histogram(h0_test, reals, Vmin, Vmax, pt_nb_bin, histo_h0);
	
	//Calculer de l'histo de H1
	calcul_histogram(h1_test, reals, Vmin, Vmax, pt_nb_bin, histo_h1);
	
	//calculer des courbes Pd/Pfa
	calculer_Pd_Pfa(pt_nb_bin, vect_Pfa, vect_Pd, histo_h0, histo_h1, reals);
	
	//calculer AUC
	double auc[1];
	auc[0]=calculer_auc(vect_Pfa, vect_Pd, pt_nb_bin);
	printf("AUC de COR = %f .\n",auc[0]);
	//sauvegarder
	saveVec2oct(abscisse,pt_nb_bin,"Abs",nom_fichier,1);
	saveVec2oct(histo_h0,pt_nb_bin,"Hist_H0",nom_fichier,0);
	saveVec2oct(histo_h1,pt_nb_bin,"Hist_H1",nom_fichier,0);
	saveVec2oct(vect_Pfa,pt_nb_bin,"Vect_Pfa",nom_fichier,0);
	saveVec2oct(vect_Pd,pt_nb_bin,"Vect_Pd",nom_fichier,0);
	saveVec2oct(auc,1,"AUC",nom_fichier,0);
//	for (int i=0;i<pt_nb_bin;i++){
//		printf("Pfa[%d]:  %f \n",i,vect_Pfa[i]);
//	}
//	
//	for (int i=0;i<pt_nb_bin;i++){
//		printf("Pd[%d]:  %f \n",i,vect_Pd[i]);
//	}
	
	free(histo_h1);
	free(histo_h0);
	free(vect_Pfa);
	free(vect_Pd);
	free(abscisse);
	
//	histo_h0=NULL;
//	histo_h1=NULL;
//	vect_Pd=NULL;
//	vect_Pfa=NULL;
//	abscisse=NULL;
	printf("Un courbe COR cree!\n\n");
}