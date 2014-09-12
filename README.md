stage
=====
Le sujet du stage est l’évaluation des performances de détection en imagerie polarimétrique par contraste de speckle. 

L’imagerie polarimétrique consiste à mesurer et analyser l’état de polarisation de la lumière rétrodiffusée par un objet en tout point de l’image. Ce type d’imagerie trouve des applications de plus en plus nombreuses biomédical, industrie, télédétection. Les méthodes usuelles nécessitent d’acquérir et traiter plusieurs images. Une méthode “computationnelle” fondée sur l’acquisition d’une unique image a été étudiée récemment à l’IPR. Elle est basée sur l’analyse des statistiques d’intensité lumineuse de la figure de speckle créée par la diffusion de l’éclairement actif cohérent utilisé pour illuminer l’objet imagé. Le stage vise à évaluer rigoureusement les performances d’une telle approche “computationnelle” pour des problèmes de détection/discrimination de matériaux aux propriétés polarimétriques distinctes.

Sous une forme simplifiée, ce problème consiste à distinguer les 2 hypothèses suivantes: les deux images ont mêmes statistiques (hypothèse H0) ou non (hypothèse H1). Pour faire cela, il faut comparer deux approches polarimétriques: une est basée sur l’image OSC, dont les statistiques d’intensité ont des distributions exponentielle; l’autre approche est basée sur une image d’intensité de speckle dont la statistique d’intensité dépend de degré de polarisation (DOP). Ensuite, il faut construire différents tests de détection. Un test de détection est une fonction des données mesurées que l’on compare avec un seuil de décision — λ. La qualité globale du test se mesure avec AUC (Area Under Curve), qui est l’aire sous la courbe COR. Une courbe COR est une courbe paramétrique qui décrit l’évolution de la probabilité de détection en fonction de la probabilité de fausse alarme2. 
Donc, le but du stage est de simuler numériquement des signaux optiques réalistes afin de comparer différents tests statistiques de détection en générant les courbes COR et  calculant l’AUC. La partie la plus importante du travail a été implémentations des tests de détection. Au cours du stages, 12 tests de détection ont été implémentés:

Delta(image OSC)
LRTp (OSC)
GLRTp (OSC)	
LRTs (speckle)
GLRTs (image speckle)
Var Gausse (speckle)
Intensité (speckle)
Gausse Multi (speckle)
La différence du logmoment (speckle)
Le Quotient du logmoment (specle)	
LRT gamma (speckle)	GLRT gamma (speckle)

Les détails mathématiques de ces différents tests ne seront pas détaillés dans ce rapport faute de place. Les 3 premiers sont basés sur l’images d’OSC et les 9 autres sur l’image d’intensité de speckle.

