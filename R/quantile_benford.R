#' Calcul de la table de quantile Smooth test par Monte carlo
#' Une fonction pour calculer la table de quantile à la Benford
#'
#' @param alpha niveau de confiance
#' @param taille taille de l echantillon
#' @param nb_repetition nombre de repetititons
#' @param digits FSD
#' @param K nombre  de smooth test
#' @param support_vector vecteur numerique
#' @param base_smooth base
#' @param optimal_t_k boolean
#' @return la fonction renvoie une matrice
#' @import stats
#' @export
calcul_quantile_monte_carlo <- function(alpha, taille, nb_repetition, digits=1,
                            K=4, support_vector=c(1:9),
                            base_smooth=10, optimal_t_k=TRUE){
  nb_row <- K
  name_row <- paste0("T",1:K)
  if (optimal_t_k){
    nb_row <- nb_row + 1
    name_row <- c(name_row, "T_k")
  }
  tmp_test_simuler <- matrix(NA, nb_row, nb_repetition)

  row.names(tmp_test_simuler) <- name_row

  proba_theorique  <- generer_probabilite_theorique(digits)

  # Si freq.cumule.digits.data
  # P freq.digits.data
  for(compteur in 1:nb_repetition){
    data <-  generer_data_benford(n = taille, support_vector_to_use = support_vector,
                                  digits = digits, base = base_smooth)
    proba_emprique <- generer_probabilite_empirique(data, support_vector)
    result <- calcul_tk(data, K,support_vector, proba_theorique)
    if(optimal_t_k){
      result <- c(result,calcul_tk_widehat(result, K, taille , base_smooth))
    }
    tmp_test_simuler[, compteur] <- result
  }
  return(apply(tmp_test_simuler, 1, function(x){quantile(x, (1-alpha))}))
}



#' Calcul de la p-value
#' Une fonction pour calculer la p-value
#'
#' @param probabilite_theorique vecteur de probabilite theorique
#' @param probabilite_empirique vecteur de probabilite empirique
#' @param taille taille de l echantillon
#' @param nb_repetition nombre de repetititons
#' @param K nombre  de smooth test
#' @param support_vector vecteur numerique
#' @param base_smooth Base
#' @param optimal_t_k boolean
#' @param quantile_vector vecteur quantile
#' @return la fonction renvoie un vecteur de p-value
#' @export
calcul_p_value <- function(taille ,nb_repetition, quantile_vector,
                                   K=4,probabilite_theorique,
                                   probabilite_empirique,base_smooth,
                                   optimal_t_k=FALSE,support_vector=c(1:9)){
  nb_row <- K

  name_row <- paste0("S", 1:K)

  if (optimal_t_k){
    nb_row <- nb_row + 1
    name_row <- c(name_row, "S_k")
  }

  puissance_vect<- rep(0, nb_row)
  puissance_matrice <- matrix(0, nb_row, nb_repetition)
  #progress.bar <- txtProgressBar(min = 0, max = length(taille.vector), style = 3)

  for (a in 1: nb_repetition){
      data <- sample(support_vector, size = taille,
                     prob = probabilite_empirique,
                     replace = TRUE )

      result <- calcul_tk(data, K, support_vector, probabilite_theorique)
      if(optimal_t_k){
        result <- c(result,calcul_tk_widehat(result, K, taille , base_smooth))
      }

      puissance_matrice[1:nb_row, a] <- as.numeric(result > quantile_vector)

    }
  p_value <- apply(puissance_matrice, 1, sum)
  p_value <- p_value / nb_repetition
  return(p_value)
}

