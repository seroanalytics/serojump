devtools::load_all()

# Case Study 1

models_all <- readRDS(here::here("hpc", "sim", "models_all.RDS"))

pp_model_case1 <- models_all[[1]]$model

p1 <- plotSero(pp_model_case1)
p2 <- plotPriorPredictive(pp_model_case1)
p3 <- plotPriorInfection(pp_model_case1)

p1
ggsave(here::here("outputs", "figs", "supp", "pp_case1_sero.png"))
p2
ggsave(here::here("outputs", "figs", "supp", "pp_case1_ab.png"))
p3
ggsave(here::here("outputs", "figs", "supp", "pp_case1_infs.png"))


# Case Study 2

models_all <- readRDS(here::here("hpc", "transvir", "models_all.RDS"))

pp_model_case1 <- models_all[[1]]$model

p1 <- plotSero(pp_model_case1)
p2 <- plotPriorPredictive(pp_model_case1)
p3 <- plotPriorInfection(pp_model_case1)


p1
ggsave(here::here("outputs", "figs", "supp", "pp_case2_sero.png"))
p2
ggsave(here::here("outputs", "figs", "supp", "pp_case2_ab.png"))
p3
ggsave(here::here("outputs", "figs", "supp", "pp_case2_infs.png"))

