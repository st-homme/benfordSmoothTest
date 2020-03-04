
.calcul_p_value <- function(data_empirique, taille, nb_repetition, digits=1,
                            K=5, support_vector=c(1:9),
                            base_smooth=10){
  nb_row <- K+1
  name_row <- paste0("T",1:K)
  name_row <- c(name_row, "T_k")

  tmp_test_simuler <- matrix(NA, nb_row, nb_repetition)

  row.names(tmp_test_simuler) <- name_row

  proba_theorique  <- .generer_probabilite_theorique(digits)

  # Si freq.cumule.digits.data
  # P freq.digits.data
  puissance_matrice <- matrix(0, nb_row, nb_repetition)
  result_empirique <- .calcul_tk(data_empirique, K)
  result_empirique <- c(result_empirique,.calcul_tk_widehat(result_empirique, K, taille , base_smooth))
  pb <- txtProgressBar(min = 1, max = nb_repetition, style = 3)
  for(compteur in 1:nb_repetition){

    data <-  .generer_data_benford(n = taille, support_vector_to_use = support_vector,
                                  digits = digits, base = base_smooth)
    result_benford <- .calcul_tk(data, K)
    result_benford <- c(result_benford,.calcul_tk_widehat(result_benford, K, taille , base_smooth))

    puissance_matrice[1:nb_row, compteur] <- as.numeric(result_empirique > result_benford)

    # tmp_test_simuler[, compteur] <- result
    setTxtProgressBar(pb, compteur)
  }
  close(pb)
  p_value <- apply(puissance_matrice, 1, sum)

  p_value <- p_value / nb_repetition
  return(p_value)
}



