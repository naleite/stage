#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utils_fade.h"
#include "utils_tirage.h"
#include <unistd.h>


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
			
			///////////////////////////////////////PROBLEM ICI
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
	double A=nbw*nbw_/(nbw_+nbw);
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
	
	histo_h0=NULL;
	histo_h1=NULL;
	vect_Pd=NULL;
	vect_Pfa=NULL;
	abscisse=NULL;
	printf("Un courbe COR cree!\n\n");
}
