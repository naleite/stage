/*
 Le programme de la simulation
 Commencer le 15 mai 2014 a IPR
 */

# include <stdio.h>
# include <time.h>
# include <math.h>
# include <stdlib.h>
# include <getopt.h>
# include <string.h>
# include "librairie.h"
# include <unistd.h>

int verifier(int N, int n, double I,double i,double P,double p,int O,int R, char *flag){
	int nb_error=0;
	if (N<=0||n<=0) {
		printf("ERROR: Nombre de pixels doit etre positive. E.x N=100\n");
		nb_error++;
		//return 1;
	}
	if(I<=0||i<=0){
		printf("ERROR: Intensite moyenn doit etre positive. E.x I=100\n");
		nb_error++;
		//return 2;
	}
	
	if (P>=1||p>=1||P<=0||p<=0) {
		printf("ERROR: DOP doit etre (0,1),(0<DOP<1). E.x p=0.8\n");
		nb_error++;
		//return 3;
	}
	if (O<=0) {
		printf("ERROR: Ordre doit etre positive. E.x O=1\n");
		nb_error++;
		//return 4;
	}
	if (R<=0){
		printf("ERROR: Nombre de realisation doit etre positive. E.x R=1000\n");
		nb_error++;
//		return 5;
	}
	if (strlen(flag)!=10) {
		printf("ERROR: Trop de Tests(%lu) choisis. \n",strlen(flag));
		nb_error++;
//		return 6;
	}
	return nb_error;
}


void affiche_params(int N, int n, double I,double i,double P,double p,int O,int R,char *flag,int aff){
	if(aff==1){
		//system("clear");
		printf("\n----------------TABLEAU DE PARAMETRES DE LA SIMULATION ------------------\n");
		printf("Nombre de realisations:		%d \n",R);
		printf("Ordre de la simulation:		%d \n",O);
		printf("Nombre de pixels cible:		%d \n",N);
		printf("Nombre de pixels fond:		%d \n",n);
		printf("Intensite moyenne de cible :	%f \n",I);
		printf("Intensite moyenne de fond  :	%f \n",i);
		printf("DOP de cible :			%f \n",P);
		printf("DOP de fond :			%f \n",P);
		printf("Lest tests choisis:\n");
		int nbtest=10;
		char *test_string[nbtest]; //
		test_string[0]="Test de delta.";
		test_string[1]="Test de LRTp.";
		test_string[2]="Test de GLRTp.";
		test_string[3]="Test de LRTs.";
		test_string[4]="Test de GLRTs.";
		test_string[5]="Test de variance Gaussian.";
		test_string[6]="Test de Intensite.";
		test_string[7]="Test de Gausse Multi.";
		test_string[8]="Test de Diff Logmoment.";
		test_string[9]="Test de Quot Logmoment.";
		for (int k=0; k<nbtest; k++) {
			if (flag[k]=='y'||flag[k]=='Y') {
				printf("\t%s\n",test_string[k]);
			}
		}
		printf("-------------------------------------------------------------------------\n");
		printf("\nVerifiez les parametres et tapez 'y' ou 'entree' pour continuer la simulation...\n");
		
		char ch;
		scanf("%c",&ch);
		if (ch!='y'&&ch!='Y'&&ch!=10) {
			printf("La simulatioin se termine.\n");
			exit(0);
		}
	}
}
void usage(FILE *fp_out, char **argv)
{
	fprintf(fp_out,"\n USAGE :%s",argv[0]);
	fprintf(fp_out,"\n           -N 'Taille de cible, Nw' : Le nombre de pixels que le cible contient.");
	fprintf(fp_out,"\n           -n 'Taille de fond, Nw_' : Le nombre de pixels que le fond contient.");
	fprintf(fp_out,"\n           -R 'Reals' : Le nombre de realisations");
	fprintf(fp_out,"\n           -O 'Ordre' :L'ordre.");
	fprintf(fp_out,"\n           -I 'moy cible' : valeur moyenne du intensite de cible");
	fprintf(fp_out,"\n           -i 'moy fond' : valeur moyenne du intensite de fond");
	fprintf(fp_out,"\n           -P 'DOP de cible' : Degres de polarisation de cible ");
	fprintf(fp_out,"\n           -p 'DOP de fond' : Degres de polarisation de fond  ");
	fprintf(fp_out,"\n           -t Choisir les tests(saisir 10 'y' ou 'n').Ex: yyyyynnyyy\n");
	int nbtest=10;
	char *test_string[nbtest]; //
	test_string[0]="Test de delta.";
	test_string[1]="Test de LRTp.";
	test_string[2]="Test de GLRTp.";
	test_string[3]="Test de LRTs.";
	test_string[4]="Test de GLRTs.";
	test_string[5]="Test de variance Gaussian.";
	test_string[6]="Test de Intensite.";
	test_string[7]="Test de Gausse Multi.";
	test_string[8]="Test de Diff Logmoment.";
	test_string[9]="Test de Quot Logmoment.";
	for (int k=0; k<nbtest; k++) {
		
			printf("\t\t%d) %s\n",k+1,test_string[k]);
		
	}

	fprintf(fp_out,"\n           -h -? fournit cette page de manuel");
	fprintf(fp_out,"\n           -v Voir les arguments pour cette simulation.");
	fprintf(fp_out,"\n CrÈe une image N*R contenant R realisations de tirages sous poissoniens");
	fprintf(fp_out,"\n d'intensitÈ comportant C pixels cible et N-C pixels de fond");
	fprintf(fp_out,"\n Par defaut N=R=256, Fano = 1 (poisson) et m egale a 10.");
	fprintf(fp_out,"\n Par dÈfaut, la cible fait N/2 ");
	fprintf(fp_out,"\n Image de sortie au format srf ou bin");
	fprintf(fp_out,"\n     La sortie se fait toujours sur la sortie standard \n");
	fprintf(fp_out,"\n Exemple : -v -N 300 -n 300 -I 200 -i 200 -P 0.9 -p 0.2 -O 1 -R 500 -t yyyyyyyynn \n ");
	exit(0);
	
}


/*************P.Principal**************/
//==================================================================================================
//==================================================================================================

int main(int argc, char * argv[])
{
	int c;
	extern char *optarg;
	char *chaine="t:N:n:R:O:I:i:P:p:h?vdf:z";
	//int nb_tot_test=10;
	char *test_flag="yyyyyyyyyy";
	
	int Nw=400; //Nb de pixels des images w et w_ d'hypo H0 et H1.
	int Nw_=400;
	
	int ordre=1;
	int reals=1000;
	char *prefix="";
	int zip=0;
	//arguments de H1
	double h1_Imoyco_w=200;
	double h1_Imoyco_w_=200;
	double h1_DOP_w=0.99;
	double h1_DOP_w_=0.25;
	
	int affiche=0;
	int nb_test=0;
	if (argc==1) {
		system("clear");
		printf("Vous ne saisissez aucun parametre.\n");
		printf("\tTapez -h ou -? pour le manual.\n");
		printf("\tTapez -d pour lancer la simulation par defaut.\n\n");
		exit(1);
	}
	while ((c=getopt(argc, argv, chaine))!=EOF) {
		switch (c) {
			case 'N':
				Nw=atoi(optarg);break;
			case 'n':
				Nw_=atoi(optarg);break;
			case 'R':
				reals=atoi(optarg);break;
			case 'O':
				ordre=atoi(optarg);break;
			case 'I':
				h1_Imoyco_w=atof(optarg);break;
			case 'i':
				h1_Imoyco_w_=atof(optarg);break;
			case 'P':
				h1_DOP_w=atof(optarg);break;
			case 'p':
				h1_DOP_w_=atof(optarg);break;
			case 'h':
			case '?':
				usage(stdout, argv);break;
			case 'v':
				//affiche_params(Nw, Nw_, h1_DOP_w, h1_DOP_w_, h1_DOP_w, h1_DOP_w_, ordre, reals);
				affiche=1;
				break;
			case 't':
				test_flag=optarg;
				//printf("%s",optarg);
				break;
			case 'd':
				printf("La simulation utilise les parametres par defaut.\n");
				affiche_params(Nw, Nw_, h1_Imoyco_w, h1_Imoyco_w_, h1_DOP_w, h1_DOP_w_, ordre, reals,test_flag,1);
				break;
			case 'f':
				prefix=optarg;
				break;
			case 'z':
				zip=1;break;
			default:
				break;
		}
	}
	system("clear");
	int bon_param=verifier(Nw, Nw_, h1_Imoyco_w, h1_Imoyco_w_, h1_DOP_w, h1_DOP_w_, ordre, reals,test_flag);
	if (bon_param>0) {
		printf("\nLe program trouve %d error(s) parametres.\nMerci de verifier les parametres et de relancer la simulation.\n",bon_param);
		exit(1);
	}
	affiche_params(Nw, Nw_, h1_Imoyco_w, h1_Imoyco_w_, h1_DOP_w, h1_DOP_w_, ordre, reals,test_flag,affiche);
	
	
	//arguments de H0
	double h0_Imoyco_w_=h1_Imoyco_w_;
	double h0_Imoyco_w =h0_Imoyco_w_;
	double h0_DOP_w_=h1_DOP_w_;
	double h0_DOP_w =h0_DOP_w_;
	
	//================GENERATIONS ET TESTS POUR H0==================//
	double **h0_tab_X=d_alloue_2d(reals, Nw+Nw_);
	double **h0_tab_Y=d_alloue_2d(reals, Nw+Nw_);
	
	double **h0_tab_tot=d_alloue_2d(reals, Nw+Nw_);
    double **h0_tab_ICEO=d_alloue_2d(reals, Nw+Nw_);
	double **h0_tab_delta=d_alloue_2d(reals, Nw+Nw_);
	
	//Generer des images // et |_ (X ET Y) pour w et w_ de H0;
	generer_image_XY(Nw, h0_Imoyco_w, h0_DOP_w, Nw_, h0_Imoyco_w_, h0_DOP_w_, ordre, reals, h0_tab_X, h0_tab_Y);
	
	//generer Itot, p et delta de image en fonction de Image_X et Image_Y
	generer_3image(Nw+Nw_, reals, h0_tab_X, h0_tab_Y, h0_tab_tot, h0_tab_ICEO, h0_tab_delta);
	
	//================GENERATIONS ET TESTS POUR H1==================//
	
	double **h1_tab_X=d_alloue_2d(reals, Nw+Nw_);//Generation de image//
	double **h1_tab_Y=d_alloue_2d(reals, Nw+Nw_);
	
	double **h1_tab_tot=d_alloue_2d(reals, Nw+Nw_);//tot
    double **h1_tab_ICEO=d_alloue_2d(reals, Nw+Nw_);//ICEO
	double **h1_tab_delta=d_alloue_2d(reals, Nw+Nw_);//Delta
	
	//Generer des images // et |_ (X ET Y) pour w et w_ de H1;
	generer_image_XY(Nw, h1_Imoyco_w, h1_DOP_w, Nw_, h1_Imoyco_w_, h1_DOP_w_, ordre, reals, h1_tab_X, h1_tab_Y);
	
	//generer Itot, ICEO et Delta de image en fonction de Image_X et Image_Y
	generer_3image(Nw+Nw_, reals, h1_tab_X, h1_tab_Y, h1_tab_tot, h1_tab_ICEO, h1_tab_delta);
	
	
/*===========================TEST DE DETECTIONS ET GENERATION COURBE COR===============================*/
	if(test_flag[0]=='y'||test_flag[0]=='Y'){
	//Pour test delat
		double *h0_test_delta;//Test de dectection de delta pour H0
		double *h1_test_delta;//Test de dectection de delta pour H1
		h0_test_delta=test_delta(reals, Nw, Nw_, h0_tab_delta); //test2 de detection pour Delta
		h1_test_delta=test_delta(reals, Nw, Nw_, h1_tab_delta);
		printf("COR_delta:\n");
		char nom_fiche[100];
		strcat(nom_fiche, prefix);
		strcat(nom_fiche, "_Cor_delta.dat");
		//printf("%s\n",nom_fiche);
		//generer_COR(h0_test_delta, h1_test_delta, reals,"Cor_delta.dat");
		generer_COR(h0_test_delta, h1_test_delta, reals,nom_fiche);
		nb_test++;
	}//test1
	
	if(test_flag[1]=='y'||test_flag[1]=='Y'){
    //Test de LRT p
		double *h0_test_LRTp;
		double *h1_test_LRTp;
		h0_test_LRTp=test_LRT(h0_DOP_w_, h1_DOP_w_, h1_DOP_w, reals, Nw, Nw_, h0_tab_ICEO);
		h1_test_LRTp=test_LRT(h0_DOP_w_, h1_DOP_w_, h1_DOP_w, reals, Nw, Nw_, h1_tab_ICEO);
		printf("COR_LRTp:\n");
		char nom_fiche[100];
		strcat(nom_fiche, prefix);
		strcat(nom_fiche, "_Cor_LRTp.dat");
		generer_COR(h0_test_LRTp, h1_test_LRTp, reals,nom_fiche);
		nb_test++;
	}//test2
	
	if(test_flag[2]=='y'||test_flag[2]=='Y'){
	//Test de GLRT p
		double *h0_test_GLRTp;
		double *h1_test_GLRTp;
		h0_test_GLRTp=test_GLRT(reals, Nw, Nw_, h0_tab_ICEO);
		h1_test_GLRTp=test_GLRT(reals, Nw, Nw_, h1_tab_ICEO);
		printf("COR_GLRTp:\n");
		char nom_fiche[100];
		strcat(nom_fiche, prefix);
		strcat(nom_fiche, "_Cor_GLRTp.dat");
		generer_COR(h0_test_GLRTp, h1_test_GLRTp, reals,nom_fiche);
		nb_test++;
	}//test3
	
	if(test_flag[3]=='y'||test_flag[3]=='Y'){
	//Test de LRT_s
		double *h0_test_LRT_s;
		double *h1_test_LRT_s;
		h0_test_LRT_s=test_LRT_s(h1_Imoyco_w, h0_Imoyco_w_, h1_DOP_w, h0_DOP_w_, reals, Nw,Nw_,h0_tab_tot);
		h1_test_LRT_s=test_LRT_s(h1_Imoyco_w, h0_Imoyco_w_, h1_DOP_w, h0_DOP_w_, reals, Nw,Nw_,h1_tab_tot);
		printf("COR_LRTs:\n");
		char nom_fiche[100];
		strcat(nom_fiche, prefix);
		strcat(nom_fiche, "_Cor_LRTs.dat");
		generer_COR(h0_test_LRT_s, h1_test_LRT_s, reals,nom_fiche);
		nb_test++;
	}//test4
	
	if(test_flag[4]=='y'||test_flag[4]=='Y'){
	//Test de GLRT_S
		double *h0_test_GLRT_s;
		double *h1_test_GLRT_s;
		h0_test_GLRT_s=test_GLRT_s(reals, Nw, Nw_, h0_tab_tot);
		h1_test_GLRT_s=test_GLRT_s(reals, Nw, Nw_, h1_tab_tot);
		printf("COR_GLRTs:\n");
		char nom_fiche[100];
		strcat(nom_fiche, prefix);
		strcat(nom_fiche, "_Cor_GLRTs.dat");
		generer_COR(h0_test_GLRT_s, h1_test_GLRT_s, reals,nom_fiche);
		nb_test++;
	}//test5
	
	if(test_flag[5]=='y'||test_flag[5]=='Y'){
	//Test de variance gaussian
		double *h0_test_var_gausse;
		double *h1_test_var_gausse;
		h0_test_var_gausse=test_var_gausse(reals, Nw, Nw_, h0_tab_tot);
		h1_test_var_gausse=test_var_gausse(reals, Nw, Nw_, h1_tab_tot);
		printf("COR_Var_Gausse:\n");
		char nom_fiche[100];
		strcat(nom_fiche, prefix);
		strcat(nom_fiche, "_Cor_Var_Gausse.dat");
		generer_COR(h0_test_var_gausse, h1_test_var_gausse, reals, nom_fiche);
		nb_test++;
	}//test6
	
	if(test_flag[6]=='y'||test_flag[6]=='Y'){
	//Test de l'intensite
		double *h0_test_intensite;
		double *h1_test_intensite;
		h0_test_intensite=test_intensite(reals, Nw, Nw_, h0_tab_tot);
		h1_test_intensite=test_intensite(reals, Nw, Nw_, h1_tab_tot);
		printf("COR_Intensite:\n");
		char nom_fiche[100];
		strcat(nom_fiche, prefix);
		strcat(nom_fiche, "_Cor_intensite.dat");
		generer_COR(h0_test_intensite, h1_test_intensite, reals,nom_fiche);
		nb_test++;
	}//test7
	
	if(test_flag[7]=='y'||test_flag[7]=='Y'){
		//Test de gausse multi
		double *h0_test_gausse_multi;
		double *h1_test_gausse_multi;
		h0_test_gausse_multi=test_gausse_multi(reals, Nw, Nw_, h0_tab_tot);
		h1_test_gausse_multi=test_gausse_multi(reals, Nw, Nw_, h1_tab_tot);
		printf("COR_Gausse_Multi:\n");
		char nom_fiche[100];
		strcat(nom_fiche, prefix);
		strcat(nom_fiche, "_Cor_Gausse_Multi.dat");
		generer_COR(h0_test_gausse_multi, h1_test_gausse_multi, reals, nom_fiche);
		nb_test++;
	}//test8
	
	if(test_flag[8]=='y'||test_flag[8]=='Y'){
		//Test de diff logmoment
		double *h0_test_diff_lm=test_diff_logmoment(reals, Nw, Nw_, h0_tab_tot);
		double *h1_test_diff_lm=test_diff_logmoment(reals, Nw, Nw_, h1_tab_tot);
		printf("COR_Diff_logmoment:\n");
		char nom_fiche[100];
		strcat(nom_fiche, prefix);
		strcat(nom_fiche, "_Cor_Diff_Logmoment.dat");
		generer_COR(h0_test_diff_lm, h1_test_diff_lm, reals, nom_fiche);
		nb_test++;
	}//test9
	
	if(test_flag[9]=='y'||test_flag[9]=='Y'){
		//Test de Quot logmoment
		double *h0_test_quot_lm=test_quot_logmoment(reals, Nw, Nw_, h0_tab_tot);
		double *h1_test_quot_lm=test_quot_logmoment(reals, Nw, Nw_, h1_tab_tot);
		printf("COR_Quot_logmoment:\n");
		char nom_fiche[100];
		strcat(nom_fiche, prefix);
		strcat(nom_fiche, "_Cor_Quot_Logmoment.dat");
		generer_COR(h0_test_quot_lm, h1_test_quot_lm, reals,nom_fiche);
		nb_test++;
	}//test10
	
	
	char buffer[100];
	getcwd(buffer, sizeof(buffer));
	printf("Generation COR reussit.\nLes donnees sont sauvegardees sous \n%s avec nom du fichier comme\"%s_Cor_test.dat\" \n",buffer,prefix);
	
	
	
//	for (int i=0; i<reals; i++) {
//		printf("h0_GLRTs[%d]: %f \n",i,h0_test_GLRT_s[i]);
//	}
//	
//	for (int i=0; i<reals; i++) {
//		printf("h1_GLRTs[%d]: %f \n",i,h1_test_GLRT_s[i]);
//	}
	
	
	
	//printf("test log 10=%f\n",log(10));
	
	//Afficher des resultat de LRT p
//	printf("exp2=%f \n",exp(2));
//	printf("exp2=%f \n",exp2(2));
//	printf("log2=%f \n",log(-2));
	
	//============Des nettoyages=========================
	d_efface_2d(reals, Nw+Nw_, h0_tab_X);
	d_efface_2d(reals, Nw+Nw_, h0_tab_Y);
	d_efface_2d(reals, Nw+Nw_, h0_tab_tot);
	d_efface_2d(reals, Nw+Nw_, h0_tab_ICEO);
	d_efface_2d(reals, Nw+Nw_, h0_tab_delta);
	
	d_efface_2d(reals, Nw+Nw_, h1_tab_X);
	d_efface_2d(reals, Nw+Nw_, h1_tab_Y);
	d_efface_2d(reals, Nw+Nw_, h1_tab_tot);
	d_efface_2d(reals, Nw+Nw_, h1_tab_ICEO);
	d_efface_2d(reals, Nw+Nw_, h1_tab_delta);

//	free(h0_test_delta);
//	free(h0_test_itot);
//	free(h0_test_LRTp);
//	free(h0_test_GLRTp);
//	free(h0_test_GLRT_s);
//	free(h0_test_LRT_s);
//
//	free(h1_test_delta);
//	free(h1_test_itot);
//	free(h1_test_LRTp);
//	free(h1_test_GLRTp);
//	free(h1_test_GLRT_s);
//	free(h1_test_LRT_s);
	
	
	if (zip==1) {
		char command[100];
		char *cmd="tar -cf ";
		char *zip=prefix;
		strcat(command, cmd);
		strcat(command, zip);
		strcat(command, ".tar *.dat");
		printf("Les fichiers .dat sont archives. %s.tar\n",prefix);
		system(command);

	}
	
	return 0;
}

