setwd("C:/Users/apple/OneDrive/바탕 화면/4-1/tjb/fpj")

install.packages(spdep)
install.packages("glmnet")

library(tidyverse)
library(rsample)
library(pROC)
library(MASS)
library(spdep)
library(glmnet)


df <- read_csv("wildfire(ver251204).csv")
df <- df %>%
  mutate(
    fire_occ = as.integer(CNT > 0),
    logBA1   = log1p(BA)
  )
theme_set(theme_bw(base_size = 12))

ggplot(df, aes(x = CNT)) +
  geom_histogram(bins = 60) +
  scale_x_continuous(limits = c(0, quantile(df$CNT, 0.99, na.rm = TRUE))) +
  labs(title = "CNT 분포(상위 1% 절단)", x = "산불 발생 빈도(CNT)", y = "관측치 수")

df %>%
  mutate(is_zero = if_else(CNT == 0, "CNT=0", "CNT>0")) %>%
  count(is_zero) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = is_zero, y = prop)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "CNT의 0 비율", x = NULL, y = "비율")

cnt_mv <- df %>%
  group_by(month) %>%
  summarise(mean_cnt = mean(CNT), var_cnt = var(CNT), .groups = "drop")

ggplot(cnt_mv, aes(x = mean_cnt, y = var_cnt)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(title = "월별 CNT 평균-분산 비교", x = "평균", y = "분산")

ggplot(df, aes(x = BA)) +
  geom_histogram(bins = 60) +
  scale_x_continuous(limits = c(0, quantile(df$BA, 0.99, na.rm = TRUE))) +
  labs(title = "BA 분포(상위 1% 절단)", x = "산불 규모(BA)", y = "관측치 수")

df %>%
  dplyr::select(BA, logBA1) %>%
  tidyr::pivot_longer(cols = c(BA, logBA1),
                      names_to = "type", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 60) +
  facet_wrap(~type, scales = "free_x") +
  labs(title = "산불 규모 분포: BA vs log(1+BA)",
       x = NULL, y = "관측치 수")

month_occ <- df %>%
  group_by(month) %>%
  summarise(fire_rate = mean(fire_occ), .groups = "drop")

ggplot(month_occ, aes(x = month, y = fire_rate)) +
  geom_line() + geom_point() +
  labs(title = "월별 산불 발생 확률", x = "month", y = "P(fire_occ=1)")

month_cnt <- df %>%
  group_by(month) %>%
  summarise(mean_cnt = mean(CNT), .groups = "drop")

ggplot(month_cnt, aes(x = month, y = mean_cnt)) +
  geom_line() + geom_point() +
  labs(title = "월별 평균 산불 발생 빈도(CNT)", x = "month", y = "mean(CNT)")

month_ba <- df %>%
  filter(BA > 0) %>%
  group_by(month) %>%
  summarise(mean_ba = mean(BA), median_ba = median(BA), .groups = "drop")

ggplot(month_ba, aes(x = month, y = mean_ba)) +
  geom_line() + geom_point() +
  labs(title = "월별 평균 산불 규모(BA, 발생시)", x = "month", y = "mean(BA | BA>0)")

df <- df %>%
  mutate(
    lon_bin = ntile(lon, 4) - 1,
    lat_bin = ntile(lat, 3) - 1,
    region  = lon_bin * 3 + lat_bin
  )

region_sum <- df %>%
  group_by(region) %>%
  summarise(
    fire_rate   = mean(fire_occ),
    mean_cnt    = mean(CNT),
    mean_ba_pos = mean(BA[BA > 0], na.rm = TRUE),
    .groups = "drop"
  )
ggplot(region_sum, aes(x = factor(region), y = fire_rate)) +
  geom_col() +
  labs(title = "region별 산불 발생 확률", x = "region", y = "fire_rate")

ggplot(df, aes(lon, lat, z = fire_occ)) +
  stat_summary_2d(fun = mean, bins = 40) +
  coord_equal() +
  scale_fill_viridis_c() +
  labs(title = "공간 격자별 산불 발생 확률(2D bin 평균)", fill = "fire_rate")

df %>%
  sample_n(min(nrow(df), 5000)) %>%
  ggplot(aes(lon, lat, color = factor(fire_occ))) +
  geom_point(alpha = 0.3, size = 0.7) +
  coord_equal() +
  labs(title = "산불 발생 여부의 공간 분포(표본)", color = "fire_occ")


set.seed(1)
split <- initial_split(df, prop = 0.8, strata = fire_occ)
train <- training(split)
test  <- testing(split)

m_logit <- glm(
  fire_occ ~ poly(lon,2) + poly(lat,2) + factor(month) +
    lc_crop + lc_tree + lc_flooded + lc_urban + lc_bare + lc_water +
    clim_u + clim_v + temp + evap + prec,
  data = train,
  family = binomial()
)


prob <- predict(m_logit, newdata = test, type = "response")
roc_obj <- roc(test$fire_occ, prob)

plot(roc_obj, main = paste0("ROC Curve (AUC = ", round(auc(roc_obj), 3), ")"))
auc(roc_obj)

test %>%
  mutate(prob = prob) %>%
  ggplot(aes(x = prob, fill = factor(fire_occ))) +
  geom_histogram(bins = 40, alpha = 0.5, position = "identity") +
  labs(title = "발생/비발생별 예측확률 분포", x = "predicted probability", fill = "fire_occ")


m_pois <- glm(
  CNT ~ poly(lon,2) + poly(lat,2) + factor(month) +
    lc_crop + lc_tree + lc_flooded + lc_urban + lc_bare + lc_water +
    clim_u + clim_v + temp + evap + prec,
  data = train,
  family = poisson()
)

m_nb <- glm.nb(
  CNT ~ poly(lon,2) + poly(lat,2) + factor(month) +
    lc_crop + lc_tree + lc_flooded + lc_urban + lc_bare + lc_water +
    clim_u + clim_v + temp + evap + prec,
  data = train
)

disp_pois <- deviance(m_pois) / df.residual(m_pois)
disp_pois

AIC(m_pois, m_nb)

pred_pois <- predict(m_pois, newdata = test, type = "response")
pred_nb   <- predict(m_nb,   newdata = test, type = "response")

perf <- tibble(
  model = c("Poisson", "NegBin"),
  MAE   = c(mean(abs(test$CNT - pred_pois)),
            mean(abs(test$CNT - pred_nb))),
  RMSE  = c(sqrt(mean((test$CNT - pred_pois)^2)),
            sqrt(mean((test$CNT - pred_nb)^2)))
)
perf

tibble(obs = test$CNT,
       Poisson = pred_pois,
       NegBin  = pred_nb) %>%
  pivot_longer(cols = c(Poisson, NegBin), names_to = "model", values_to = "pred") %>%
  ggplot(aes(x = pred, y = obs)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~model, scales = "free") +
  labs(title = "관측 CNT vs 예측 CNT", x = "predicted", y = "observed")

# Poisson 과산포 지표
disp_dev_pois <- deviance(m_pois) / df.residual(m_pois)

pear_pois <- residuals(m_pois, type = "pearson")
disp_pear_pois <- sum(pear_pois^2) / df.residual(m_pois)

c(disp_dev_pois = disp_dev_pois,
  disp_pear_pois = disp_pear_pois)


par(mfrow = c(1,2))

plot(fitted(m_pois), residuals(m_pois, type="pearson"),
     main="Poisson: Residuals vs Fitted",
     xlab="Fitted values", ylab="Pearson residuals",
     pch=16, cex=0.6)
abline(h=0, lty=2)

plot(fitted(m_nb), residuals(m_nb, type="pearson"),
     main="NegBin: Residuals vs Fitted",
     xlab="Fitted values", ylab="Pearson residuals",
     pch=16, cex=0.6)
abline(h=0, lty=2)

par(mfrow = c(1,1))


obs_zero <- mean(train$CNT == 0)

mu_pois <- predict(m_pois, type = "response")
exp_zero_pois <- mean(exp(-mu_pois))

mu_nb <- predict(m_nb, type = "response")
theta <- m_nb$theta
exp_zero_nb <- mean((theta / (theta + mu_nb))^theta)

c(obs_zero = obs_zero,
  exp_zero_pois = exp_zero_pois,
  exp_zero_nb   = exp_zero_nb)


set.seed(1)
sim_pois <- replicate(200, simulate(m_pois)[[1]])
sim_nb   <- replicate(200, simulate(m_nb)[[1]])

sim_zero_pois <- colMeans(sim_pois == 0)
sim_zero_nb   <- colMeans(sim_nb == 0)

par(mfrow=c(1,2))
hist(sim_zero_pois, breaks=30, main="Poisson: simulated zero rate",
     xlab="simulated zero proportion")
abline(v=obs_zero, lty=2)

hist(sim_zero_nb, breaks=30, main="NegBin: simulated zero rate",
     xlab="simulated zero proportion")
abline(v=obs_zero, lty=2)
par(mfrow=c(1,1))


coords <- cbind(train$lon, train$lat)
knn <- knearneigh(coords, k = 8)
nb  <- knn2nb(knn)
lw  <- nb2listw(nb, style = "W")

res_nb <- residuals(m_nb, type = "pearson")
moran.test(res_nb, lw)


train$res_nb <- residuals(m_nb, type = "pearson")


ggplot(train, aes(lon, lat, color = res_nb)) +
  geom_point(alpha = 0.5, size = 0.8) +
  coord_equal() +
  scale_color_gradient2() +
  labs(title = "NB 모형 Pearson 잔차의 공간 분포", color = "residual")


set.seed(1)

m <- min(nrow(train), 3000)
idx <- sample(seq_len(nrow(train)), m)

coords <- as.matrix(train[idx, c("lon", "lat")])
e <- train$res_nb[idx]

k <- 8
D <- as.matrix(dist(coords))
diag(D) <- Inf

nn_idx <- apply(D, 1, function(row) order(row)[1:k])

e_lag <- sapply(seq_len(m), function(i) mean(e[nn_idx[, i]]))

plot(e, e_lag,
     xlab = "Pearson residual (e_i)",
     ylab = "Spatial lag (mean of kNN residuals)",
     main = "Spatial dependence check via kNN lag")
abline(lm(e_lag ~ e), lty = 2)

cor(e, e_lag)


df_pos <- df %>% filter(BA > 0)

set.seed(1)
split_pos <- initial_split(df_pos, prop = 0.8)
train_pos <- training(split_pos)
test_pos  <- testing(split_pos)


m_lm <- lm(
  logBA1 ~ poly(lon,2) + poly(lat,2) + factor(month) +
    lc_crop + lc_tree + lc_flooded + lc_urban + lc_bare + lc_water +
    clim_u + clim_v + temp + evap + prec,
  data = train_pos
)

summary(m_lm)


par(mfrow=c(1,2))
plot(fitted(m_lm), resid(m_lm), main="Residuals vs Fitted", xlab="Fitted", ylab="Residuals")
abline(h=0, lty=2)
qqnorm(resid(m_lm)); qqline(resid(m_lm))
par(mfrow=c(1,1))


df_pos <- df %>% filter(BA > 0)

set.seed(1)
split_pos <- initial_split(df_pos, prop = 0.8)
train_pos <- training(split_pos)
test_pos  <- testing(split_pos)


m_lm <- lm(
  logBA1 ~ poly(lon,2) + poly(lat,2) + factor(month) +
    lc_crop + lc_tree + lc_flooded + lc_urban + lc_bare + lc_water +
    clim_u + clim_v + temp + evap + prec,
  data = train_pos
)

summary(m_lm)


par(mfrow=c(1,2))
plot(fitted(m_lm), resid(m_lm),
     main="Residuals vs Fitted",
     xlab="Fitted", ylab="Residuals")
abline(h=0, lty=2)

qqnorm(resid(m_lm)); qqline(resid(m_lm))
par(mfrow=c(1,1))


pred_lm <- predict(m_lm, newdata = test_pos)
rmse_lm <- sqrt(mean((test_pos$logBA1 - pred_lm)^2))
rmse_lm


num_vars <- train_pos %>%
  dplyr::select(lc_crop, lc_tree, lc_flooded, lc_urban, lc_bare, lc_water,
                clim_u, clim_v, temp, evap, prec) %>%
  drop_na()

C <- cor(num_vars)

as.data.frame(as.table(C)) %>%
  ggplot(aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  coord_equal() +
  scale_fill_gradient2() +
  labs(title = "설명변수 상관구조(히트맵)", x = NULL, y = NULL, fill = "cor")


X_train <- model.matrix(
  logBA1 ~ poly(lon,2) + poly(lat,2) + factor(month) +
    lc_crop + lc_tree + lc_flooded + lc_urban + lc_bare + lc_water +
    clim_u + clim_v + temp + evap + prec,
  data = train_pos
)[, -1]
y_train <- train_pos$logBA1

X_test <- model.matrix(
  logBA1 ~ poly(lon,2) + poly(lat,2) + factor(month) +
    lc_crop + lc_tree + lc_flooded + lc_urban + lc_bare + lc_water +
    clim_u + clim_v + temp + evap + prec,
  data = test_pos
)[, -1]
y_test <- test_pos$logBA1


X_train <- model.matrix(
  logBA1 ~ poly(lon,2) + poly(lat,2) + factor(month) +
    lc_crop + lc_tree + lc_flooded + lc_urban + lc_bare + lc_water +
    clim_u + clim_v + temp + evap + prec,
  data = train_pos
)[, -1]
y_train <- train_pos$logBA1

X_test <- model.matrix(
  logBA1 ~ poly(lon,2) + poly(lat,2) + factor(month) +
    lc_crop + lc_tree + lc_flooded + lc_urban + lc_bare + lc_water +
    clim_u + clim_v + temp + evap + prec,
  data = test_pos
)[, -1]
y_test <- test_pos$logBA1


set.seed(1)
cv_ridge <- cv.glmnet(X_train, y_train, alpha = 0)
plot(cv_ridge)


set.seed(1)
cv_lasso <- cv.glmnet(X_train, y_train, alpha = 1)
plot(cv_lasso)


pred_ridge <- predict(cv_ridge, s = "lambda.1se", newx = X_test)
pred_lasso <- predict(cv_lasso, s = "lambda.1se", newx = X_test)

rmse_ridge <- sqrt(mean((y_test - pred_ridge)^2))
rmse_lasso <- sqrt(mean((y_test - pred_lasso)^2))

tibble(
  model = c("OLS(lm)", "Ridge", "Lasso"),
  RMSE  = c(rmse_lm, rmse_ridge, rmse_lasso)
)


coef_lasso <- coef(cv_lasso, s = "lambda.1se")
coef_lasso[coef_lasso[,1] != 0, , drop = FALSE]