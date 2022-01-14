#  Optimisation combinatoire pour le *Production Routing Problem*

Dans le cadre du projet en optimisation combinatoire, le projet consiste à résoudre le PRP par des méthodes heuristiques et exactes.  

La documentation sur l'utilisation des fonctions sont dans les fichiers .jl. Mais, en particulier le fichier main.jl contient les principaux fonctions qui permettent d'obtenir les solutions qui sont placées dans le fichier Results ou Example_results selon les fonctions que vous utilisez.  

Le rapport du projet est dans le fichier Rapport_MAOA_PRP_FU_LIU.pdf  
  
Pour pouvoir lancer les fonctions, veuillez installer Julia et Cplex.

Pour obtenir tous les résultats des instances de type A pour 6 véhicules:  
• Lancez Julia sur un terminal ayant le path sur le production_routing_problem  
• Entrez sur le terminal : include("main.jl")  
• Puis entrez : Provide_full_results("A", 6)  
Les résultats (les logs et les graphes de chemins pour les VRP) seront placés dans le dossier Example_results.  

Pour obtenir tous les résultats par classe des instances de type A pour 6 véhicules:  
• Lancez Julia sur un terminal ayant le path sur le production_routing_problem  
• Entrez sur le terminal : include("main.jl")  
• Puis entrez : Provide_class_results("A", 6)  
Les résultats (les logs et les graphes de chemins pour les VRP) seront placés dans le dossier Results.  

Pour tout question sur le fonctionnement du repository, veuillez envoyer un mail à :  
Vincent Fu : vincent.fu@etu.sorbonne-universite.fr  
  
![Exemple de solution pour le VRP](./Results/exact/PDI_exact_A_050_ABS12_50_5/VRP_exact_A_050_ABS12_50_5_p5.pdf)
<embled src="./Results/exact/PDI_exact_A_050_ABS12_50_5/VRP_exact_A_050_ABS12_50_5_p5.pdf" type="application/pdf" width="100%">
</object>
