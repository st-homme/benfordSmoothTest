# require("R/probabilite_benford.R")
# require("R/quantile_benford.R")
# library("polynom")
# library("gtools")

#' Tk du test de Benford
#' Une fonction pour appliquer le smooth test à des données
#'
#' @param data un vecteur numerique
#' @param K un nombre entre 1 et 7
#' @param digits le digit a testé
#' @param optimal_t_k boolean
#' @param alpha niveau de confiance
#' @param nb_repetition Monte carlo repetition
#' @return la fonction renvoie un vecteur contenant les tk
#' @import gtools
#' @export
BenfordSmooth.test <- function(data, K, digits= 1, alpha=0.05,
                               nb_repetition= 5000,optimal_t_k=TRUE ){
  support_vector <- c(1:9)
  probabilite_theorique <- generer_probabilite_theorique(digits)
  result_empirique <- calcul_tk(data, K,support_vector, probabilite_theorique)

  taille <- length(data)
  probabilite_empirique <- generer_probabilite_empirique(data, support_vector)
  result_empirique <- c(result_empirique,calcul_tk_widehat(result_empirique, K, taille , base=10))
  p_value_vector <-c()
  if(taille <= 100){
    quantile_vector <- calcul_quantile_monte_carlo(alpha, taille=1000, nb_repetition,
                                                    digits, K, support_vector,
                                                    base_smooth=10, optimal_t_k)
    p_value_vector <- calcul_p_value(taille ,nb_repetition, quantile_vector,
                                     K,probabilite_theorique,
                                     probabilite_empirique,base_smooth=10,
                                     optimal_t_k,support_vector)
  }else{
    p_value_vector < c()
    for(i in 1:K){
      p_value_vector <- c(p_value_vector, pchisq(result_empirique[i], df=i, lower.tail=TRUE) )
    }
    p_value_vector <- c(p_value_vector, pchisq(result_empirique[K+1], df=1,lower.tail=TRUE) )
    # p_value_vector <- pchisq(result_empirique, df=K, lower.tail=FALSE)
  }
  result <- rbind(result_empirique,p_value_vector )
  row.names(result) <- c("Tk", "P-value")
  colnames(result) <- c(1:K, "K data driven")
  return(result)
}
###############################################################################################################

#' Polynome Hj
#' Une fonction pour appliquer le smooth test à des données
#'
#' @param j numerique entre 1 et 7
#' @return la fonction renvoie le polynome
#' @import polynom
#' @export
get_polynome_h_for_j <- function(j){
  if(j==1){
    return(polynomial(c(-1.3979030247488394121040493613182696526739584777368 ,
                        0.40633916736200686478542965064543693283491024599939)))
  }else if(j==2){
    return(polynomial(c( 2.2835458942827613994342486254403905701270610581437,
                         -1.6127530382510774426950563380168084894310786771510 ,
                         0.18247002110725297036924261140018110452702111183442)))
  }else if(j==3){
    return(polynomial(c(-4.0815121188328689936236326398218860591304380838753 ,
                        4.5719372742522112044947303129821412996094568232480 ,
                        -1.2052506737490148470311505196645954099462916875905,
                        0.086173299060539497073071021550861682063512004101864)))
  }else if(j==4){
    return(polynomial(c( 8.0795029241386908450486312960577659502034469468917,
                         -12.094575500356753355982353234285254134889081144404 ,
                         5.1951355207630563189562812915032734479956498201195,
                         -0.82485600201281199967265529487511547787181262658796 ,
                         0.043126384085651298652336797627506539669328872876724 )))
  }else if(j==5){
    return(polynomial(c( -18.108494394772161751612380555693712085844319263352 ,
                         33.142284947991963143915043695725705210025788289119 ,
                         - 19.722974020860465175272441917939990432050441230630,
                         5.0173393268200010622998332406879860011024286123698 ,
                         - 0.56656631604890066847536348651878387350781426664251 ,
                         0.023336822888880408751909327047967310412185600689676 )))
  }else if(j==6){
    return(polynomial(c( 47.539025123688956939839276197426162757171984295486,
                         -100.97034588448909037802964553469158020757299454309 ,
                         75.855344772076234632745120190525279742211301165629,
                         -26.682825198088079673933396924974441588864816295226 ,
                         4.7568124723412390834620445451992011968266182408605 ,
                         -0.41562928281524381714954870899838795772527009926521 ,
                         0.014116816834808702805742210360195029761782308937628 )))
  }else if(j==7){
    return(polynomial(c(-155.23869158405472059098634159455986865585205117244 ,
                        370.17852496798326169953247150561067308461284772334,
                        -331.11599730173645571307956788170061765775873248699,
                        147.72553389848838714703346742044011676782185532865,
                        -36.176614461752994684988881308483105755882446054477,
                        4.9340195510159607707278807129264094913710664275380,
                        - 0.35100285248652007809095369751269381430321388096641,
                        0.010139198243335377835234192471425857563997234141646)))
  }else{
    return(NA)

  }

}

#' Calcul des Uj
#' Une fonction pour calculer les Uj
#'
#' @param j numerique entre 1 et 7
#' @param data vecteur numérique
#' @param taille numerique
#' @param polynome.h.j venant la fonction h_for_j
#' @return la fonction renvoie le polynome
#' @importFrom stats predict
#' @export
calcul_U_for_j <- function(j, data, taille, polynome.h.j){
  return((1/sqrt(taille)) * sum(predict(polynome.h.j, data)))
}

#' Calcul des tk
#' Une fonction pour calculer les tk
#'
#' @param K numerique entre 1 et 7
#' @param data vecteur numérique
#' @param support.vector numerique
#' @param proba.vector venant la fonction h_for_j
#' @return la fonction renvoie le tk
#' @export
calcul_tk<- function(data, K, support.vector=c(1:9), proba.vector=rep(1/9, 9)){
  #print(polynome.h)
  resultat <- vector()
  taille <- length(data)
  for(j in 1:K ){
    if (j == 1 ){
      resultat[j] = (calcul_U_for_j(j, data, taille, get_polynome_h_for_j(j))^2)
    }else{
      resultat[j] = resultat[j-1] +  (calcul_U_for_j(j, data, taille,get_polynome_h_for_j(j) )^2)
    }
  }
  return(resultat)
}

##################### Gestion optimum Calcul de  $ T_{\widehat{k}$####################################
#' Calcul du tk chapeau
#' Une fonction pour calculer le tk chapeau
#'
#' @param K numerique entre 1 et 7
#' @param Tki vecteur de tk
#' @param taille taille de l'echantillon
#' @param base base logarithmique
#' @return la fonction renvoie le meilleur tk
#' @export
calcul_tk_widehat <- function(Tki, K ,taille, base){
  return(Tki[which((Tki - (1:K) * log(taille, base))== max(Tki - (1:K) * log(taille, base)) )])
}

######################### Gestion Permutation ########################"
#' Calcul du tk chapeau
#' Une fonction pour calculer le tk chapeau
#'
#' @param K numerique entre 1 et 7
#' @param data vecteur numerique
#' @return la fonction renvoie le polynome

#' @export
calcul_permutations_tk_matrix <- function(data, K){
  base <-  10
  vector.permutations <-  permutations(n=3, r=3, v=c(1:K))
  #vector.permutations <- matrix(c(1, 2, 3, 2, 1, 3, 2 , 3 ,1) , nrow = 3, ncol = 3, byrow = TRUE)
  vector.uj <- calcul_permutations_uj(data, K)

  res <- vector()
  n  <- length(data)
  #print(vector.permutations)
  #Ti.matrix <- matrix(NA , factorial(3), K)
  for(i in 1: nrow(vector.permutations)){
    Ti <- calcul_permutations_tk( c(vector.uj[c(vector.permutations[i,],4:K)]), K)
    res[i] <- calcul_tk_widehat(Ti, K, n , base)
  }

  return(res)
}

# ################################################
#' Calcul du tk chapeau
#' Une fonction pour calculer le tk chapeau
#'
#' @param K numerique entre 1 et 7
#' @param uj.vector vecteur de uj
#' @return la fonction renvoie le polynome
#' @export
calcul_permutations_tk<- function(uj.vector,K){
  return(cumsum(uj.vector^2))
}

# ##### Calcul de  $T_{permutations}$
#' Calcul du tk chapeau
#' Une fonction pour calculer le tk chapeau
#'
#' @param K numerique entre 1 et 7
#' @param data vecteur numerique
#' @return la fonction renvoie le polynome
#' @export
calcul_permutations_uj <- function(data, K){
  resultat <- vector()
  taille <- length(data)
  for(j in 1:K ){
    resultat[j] <- (calcul_U_for_j(j, data, taille, get_polynome_h_for_j(j)))
  }
  return(resultat)
}

