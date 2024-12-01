# rm(list=ls())
# 
# n=3000
# n_X=20
# R=2
# taus=0.5
# factor1=1
# factor2=0.5 
# perc_conf=0.25
# 
# taus = rep(taus, R)

sim_data <- function(n, n_X, R, taus, factor1, factor2, perc_conf) {
  
  # 1. 공변량 생성 (절반은 정규 분포, 절반은 베르누이 분포)
  n_X_actual <- n_X - R  # 실제 공변량 개수 계산
  X_cont <- matrix(rnorm(n * n_X_actual / 2), n, n_X_actual / 2)
  X_bin <- matrix(rbinom(n * n_X_actual / 2, 1, 0.3), n, n_X_actual / 2)
  X <- cbind(X_cont, X_bin)
  colnames(X) <- paste0("X", 1:n_X_actual)
  
  # 2. 하위 그룹 변수 생성
  S <- matrix(rbinom(n * R, 1, 0.25), n, R)
  colnames(S) <- paste0("S", 1:R)
  
  # 3. 치료 지표 생성
  alpha_r <- -2
  alpha_s <- rep(1, R)
  alpha_x <- rep(0, n_X_actual)  # 실제 공변량 개수 사용
  
  # 연속형 변수와 이진형 변수에 대해 각각 회귀 계수 설정
  n_conf_per_type <-  floor(perc_conf * n_X / 2)
  alpha_x[1:n_conf_per_type] <- seq(0.25 * factor1, 0.5 * factor1, length.out = n_conf_per_type)
  alpha_x[(n_X_actual/2 + 1):(n_X_actual/2 + n_conf_per_type)] <- seq(0.25 * factor1, 0.5 * factor1, length.out = n_conf_per_type)
  
  #alpha_xs <- -alpha_x * factor2
  alpha_xs <- rep(-alpha_x * factor2, times = R)  # alpha_xs를 n_X_actual * R 길이로 확장
  
  # 상호 작용 항 생성
  interaction_terms <- model.matrix( ~ X * S)[, -c(1:(n_X + 1))] 
  colnames(interaction_terms) <- paste0(rep(colnames(S), each = ncol(X)), ":", rep(colnames(X), times = ncol(S)))
  
  lin_pred <- alpha_r + S %*% alpha_s + X %*% alpha_x + interaction_terms %*% alpha_xs
  pZ <- plogis(lin_pred)
  Z <- rbinom(n, 1, pZ)
  
  # 4. 결과 변수 생성
  Y <- 0 + X %*% alpha_x + S %*% rep(0.8, R) - Z + (S * Z) %*% taus + rnorm(n)
  
  # true_X와 just_X 구분
  true_X_index <- which(alpha_x != 0)
  true_X <- X[, true_X_index, drop = FALSE]
  just_X <- X[, -true_X_index, drop = FALSE]
  
  # prop_X 생성
  prop_X <- cbind(X, S, interaction_terms)
  
  # 데이터프레임 생성 및 반환
  list(
    Y = Y, 
    Z = Z, 
    pZ = pZ, 
    true_X = true_X, 
    just_X = just_X, 
    S = S, 
    prop_X = prop_X
  )
}
