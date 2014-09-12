//
//  main.c
//  UnitTest
//
//  Created by NALEITE on 14-6-19.
//  Copyright (c) 2014年 NALEITE. All rights reserved.
//

#include <stdio.h>
#include "librairie.h"
#include <stdlib.h>
int main(int argc, const char * argv[])
{

	int Nw=250; //Nb de pixels des images w et w_ d'hypo H0 et H1.
	int Nw_=250;
	
	int ordre=1;
	int reals=1000;
	//arguments de H1
	double h1_Imoyco_w=1000;
	double h1_Imoyco_w_=1000;
	double h1_DOP_w=0.5;
	double h1_DOP_w_=0.01;
	
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


	double *h0_test_LRT_gamma=test_LRT_gamma(reals, ordre, h1_Imoyco_w, h0_Imoyco_w_, Nw, Nw_, h0_tab_tot);
	double *h1_test_LRT_gamma=test_LRT_gamma(reals, ordre, h1_Imoyco_w, h0_Imoyco_w_, Nw, Nw_, h1_tab_tot);
	printf("COR_LRT_gamma:\n");
	generer_COR(h0_test_LRT_gamma, h1_test_LRT_gamma, reals, "Cor_LRT_gamma.dat");
	free(h0_test_LRT_gamma);
	free(h1_test_LRT_gamma);
	
	double *h0_test_GLRT_gamma=test_GLRT_gamma(reals, Nw, Nw_, h0_tab_tot, 0.01, 5, 0.01);
	double *h1_test_GLRT_gamma=test_GLRT_gamma(reals, Nw, Nw_, h1_tab_tot, 0.01, 5, 0.01);
	printf("COR_GLRT_gamma:\n");
	generer_COR(h0_test_GLRT_gamma, h1_test_GLRT_gamma, reals, "Cor_GLRT_gamma.dat");
	free(h0_test_GLRT_gamma);
	free(h1_test_GLRT_gamma);
    return 0;
}

