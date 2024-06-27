get_infection_dates <- function() {
        # Get information on KNOWN infection dates
    gambia_inf_dates <- read.csv(file = here::here("data", "transvir", "David_inf_data.csv") )

    gambia_inf_dates_time <- gambia_inf_dates %>% 
        dplyr::mutate(first_inf_time = dplyr::case_when(
            First_positive_date < dmy("01/07/2021") ~ "Pre-Delta",
            First_positive_date >= dmy("01/07/2021") & First_positive_date < dmy("05/12/2021")~ "Delta",
            First_positive_date >= dmy("05/12/2021") & First_positive_date < dmy("01/06/2022")~ "Omicron"
            )
        ) %>%
        dplyr::mutate(second_inf_time = dplyr::case_when(
            second_positive_date < dmy("01/07/2021") ~ "Pre-Delta",
            second_positive_date >= dmy("01/07/2021") & second_positive_date < dmy("05/12/2021")~ "Delta",
            second_positive_date >= dmy("05/12/2021") & second_positive_date < dmy("01/06/2022")~ "Omicron"
        ) 
        )%>% 
        dplyr::mutate(third_inf_time = dplyr::case_when(
            third_positive_date < dmy("01/07/2021") ~ "Pre-Delta",
            third_positive_date >= dmy("01/07/2021") & third_positive_date < dmy("05/12/2021")~ "Delta",
            third_positive_date >= dmy("05/12/2021") & third_positive_date < dmy("01/06/2022")~ "Omicron"
        )
    ) 

    gambia_known_inf_temp <- gambia_inf_dates_time %>% mutate(
        known_pre_delta = 
            case_when(
                second_inf_time == "Pre-Delta"~ second_positive_date,
                first_inf_time == "Pre-Delta"~ First_positive_date
                ),
    known_delta = 
            case_when(
                third_inf_time == "Delta"~ second_positive_date,
                second_inf_time == "Delta"~ second_positive_date,
                first_inf_time == "Delta"~ First_positive_date
                ),
    known_omicron = 
            case_when(
                third_inf_time == "Omicron"~ third_positive_date,
                second_inf_time == "Omicron"~ second_positive_date,
                first_inf_time == "Omicron"~ First_positive_date
                )
        ) %>% dplyr::select(Participant_ID, known_pre_delta, known_delta, known_omicron)
    gambia_known_inf_temp
}

get_data_titre_model_wave1 <- function() {
    odd_people <- c("17-197B", "07-077B", "41-481E", "43-509K", "34-399H", "41-483J", "43-512H")
    gambia_pvn_raw <- read.csv(file = here::here("data", "transvir", "Pseudovirus_data_V1_V2_V3.csv") ) %>% filter(!Participant_ID %in% odd_people)

    start_date <- dmy("17-03-2020")
    start_day <- yday(start_date)

    gambia_pvnt_b1 <- gambia_pvn_raw %>% dplyr::select(Participant_ID, V1_date, contains("ptna_Ancestral_V1")) %>% mutate(
        V1_date = as.numeric(dmy(V1_date) - start_date)) %>% rename(pid = Participant_ID, time = V1_date, titre = ptna_Ancestral_V1) %>% 
        filter(!is.na(titre))
    gambia_pvnt_b0 <- gambia_pvnt_b1 %>% dplyr::select(pid) %>% mutate(time = 0, titre = 0)
    ids <- gambia_pvnt_b1 %>% pull(pid) %>% unique
    gambia_pvnt <- bind_rows(gambia_pvnt_b0, gambia_pvnt_b1) %>% filter(pid %in% ids) %>% mutate(id = as.numeric(factor(pid))) %>%
        arrange(id, time) %>% 
        mutate(titre = case_when(
            titre < 1 ~ 1,
            TRUE ~ titre
        )) %>% mutate(titre = log10(titre))
    log10(20)
    gambia_pvnt <- gambia_pvnt %>% mutate(pid = factor(id), id = as.numeric(factor(pid)))
    gambia_pvnt

}

get_exposures_wave1 <- function() {
    # Get known infection dates
    gambia_inf_dates <- get_infection_dates()
    gambia_inf_dates_pd <- gambia_inf_dates %>% filter(!is.na(known_pre_delta)) %>% mutate(inf_days = as.numeric(ymd(known_pre_delta) - dmy("17-03-2020")) + 1) %>%
        select(pid = Participant_ID, inf_days)

    # Get vaccine dates
    gambia_meta <- read.csv(file = here::here("data", "transvir","full_demog_data.csv") )
    gambia_meta_short <- gambia_meta %>% dplyr::select(pid = Participant_ID, vax_date) %>% mutate(vax_days = as.numeric(ymd(vax_date) - dmy("17-03-2020")) + 1) 

    data_titre <- get_data_titre_model_wave1()
    known_exposures <- data_titre %>% left_join(gambia_inf_dates_pd) %>% left_join(gambia_meta_short) %>% group_by(id) %>% mutate(t = row_number()) %>% ungroup %>% pivot_wider(names_from = "t", values_from = c("titre", "time")) %>% 
        mutate(inf_time = case_when(inf_days > time_1 & inf_days <= time_2 ~ inf_days, TRUE~-1)) %>% as.data.frame %>% 
        mutate(vax_time = case_when(vax_days > time_1 & vax_days <= time_2 ~ vax_days, TRUE~-1)) %>%
        select(pid, id, time_1, time_2, titre_1, titre_2, inf_time, vax_time)
        
    known_exposures
}

get_exp_prior_wave1 <- function() {
    start_date <- ymd("2020-03-17")

    gambia_daily_cases <- read.csv(file = here::here("data", "transvir", "gambia_covid_daily.csv") )
    gambia_daily_cases_ed <- gambia_daily_cases %>% filter(ymd(date) < dmy("01/07/2021")) %>% mutate(
        date_time = as.numeric(ymd(date) - start_date)) %>% mutate(week = date_time %/% 14) %>% 
        summarise(total = sum(new_cases, na.rm = TRUE), .by = week) %>% 
        mutate(prop = total / sum(total)) %>% dplyr::select(week, prop)
    exp_prior <- data.frame(
        day = 1:(nrow(gambia_daily_cases_ed)*14),
        week = gambia_daily_cases_ed$week %>% map(~rep(.x, 14)) %>% unlist,
        prob = gambia_daily_cases_ed$prop %>% map(~rep(.x / 14, 14)) %>% unlist
    )
    exp_prior
}

get_data_titre_model_wave2_pvnt <- function() {
    odd_people <- c("17-197B", "07-077B", "41-481E", "43-509K", "34-399H", "41-483J", "43-512H")
    gambia_pvn_raw <- read.csv(file = here::here("data", "transvir", "Pseudovirus_data_V1_V2_V3.csv") ) %>% filter(!Participant_ID %in% odd_people)

    gambia_iga_raw <- read.csv(file = here::here("data", "transvir", "IgA_data_V1_V2_V3.csv") ) %>% filter(!Participant_ID %in% odd_people)


    start_date <- gambia_pvn_raw %>% dplyr::select(Participant_ID, V1_date, contains("ptna_B.1.617.2_V1")) %>% pull(V1_date) %>% dmy %>% min  #"2021-03-02"
    start_day <- yday(start_date)
    gambia_pvnt_b0 <- gambia_pvn_raw %>% dplyr::select(Participant_ID, V1_date, contains("ptna_B.1.617.2_V1")) %>% mutate(
        V1_date = as.numeric(dmy(V1_date) - start_date + 1)) %>% rename(pid = Participant_ID, time = V1_date, titre = ptna_B.1.617.2_V1) %>% 
        filter(!is.na(titre))
    gambia_pvnt_b1 <- gambia_pvn_raw %>% dplyr::select(Participant_ID, V2_date, contains("ptna_B.1.617.2_V2")) %>% mutate(
        V2_date = as.numeric(dmy(V2_date) - start_date + 1)) %>% rename(pid = Participant_ID, time = V2_date, titre = ptna_B.1.617.2_V2) %>% 
        filter(!is.na(titre))
    ids <- gambia_pvnt_b1 %>% pull(pid)
    gambia_pvnt <- bind_rows(gambia_pvnt_b0, gambia_pvnt_b1) %>% filter(pid %in% ids) %>% mutate(id = as.numeric(factor(pid))) %>%
        arrange(id, time) %>% 
        mutate(titre = case_when(titre <= 40 ~ 6.3245, TRUE~titre)) %>%
        mutate(titre = case_when(
            titre < 1 ~ 1,
            TRUE ~ titre
        )) %>% mutate(titre = pmax(log10(titre), 0))

    gambia_pvnt
}

get_data_titre_model_wave2_pvnt_iga <- function() {
    odd_people <- c("17-197B", "07-077B", "41-481E", "43-509K", "34-399H", "41-483J", "43-512H")
    gambia_pvn_raw <- read.csv(file = here::here("data", "transvir", "Pseudovirus_data_V1_V2_V3.csv") ) %>% filter(!Participant_ID %in% odd_people)


    start_date <- gambia_pvn_raw %>% dplyr::select(Participant_ID, V1_date, contains("ptna_B.1.617.2_V1")) %>% pull(V1_date) %>% dmy %>% min  #"2021-03-02"
    start_day <- yday(start_date)
    gambia_pvnt_b0 <- gambia_pvn_raw %>% dplyr::select(Participant_ID, V1_date, contains("ptna_B.1.617.2_V1")) %>% mutate(
        V1_date = as.numeric(dmy(V1_date) - start_date + 1)) %>% rename(pid = Participant_ID, time = V1_date, titre = ptna_B.1.617.2_V1) %>% 
        filter(!is.na(titre))
    gambia_pvnt_b1 <- gambia_pvn_raw %>% dplyr::select(Participant_ID, V2_date, contains("ptna_B.1.617.2_V2")) %>% mutate(
        V2_date = as.numeric(dmy(V2_date) - start_date + 1)) %>% rename(pid = Participant_ID, time = V2_date, titre = ptna_B.1.617.2_V2) %>% 
        filter(!is.na(titre))
    ids <- gambia_pvnt_b1 %>% pull(pid)
    gambia_pvnt <- bind_rows(gambia_pvnt_b0, gambia_pvnt_b1) %>% filter(pid %in% ids) %>% mutate(id = as.numeric(factor(pid))) %>%
        arrange(id, time) %>% 
        mutate(titre = case_when(titre <= 40 ~ 6.3245, TRUE~titre)) %>%
        mutate(titre = case_when(
            titre < 1 ~ 1,
            TRUE ~ titre
        )) %>% mutate(titre = pmax(log10(titre), 0)) %>% rename(sVNT = titre)

    gambia_iga_raw <- read.csv(file = here::here("data", "transvir", "IgA_data_V1_V2_V3.csv") ) %>% filter(!Participant_ID %in% odd_people)

    gambia_iga_b0 <- gambia_iga_raw %>% dplyr::select(Participant_ID, V1_date, contains("IgA_B.1.617.2_V1")) %>% mutate(
        V1_date = as.numeric(dmy(V1_date) - start_date + 1)) %>% rename(pid = Participant_ID, time = V1_date, titre = IgA_B.1.617.2_V1) %>% 
        filter(!is.na(titre))

    gambia_iga_b1 <- gambia_iga_raw %>% dplyr::select(Participant_ID, V2_date, contains("IgA_B.1.617.2_V2")) %>% mutate(
        V2_date = as.numeric(dmy(V2_date) - start_date + 1)) %>% rename(pid = Participant_ID, time = V2_date, titre = IgA_B.1.617.2_V2) %>% 
        filter(!is.na(titre))

    ids <- gambia_pvnt_b1 %>% pull(pid)

    gambia_iga <- bind_rows(gambia_iga_b0, gambia_iga_b1) %>% filter(pid %in% ids) %>% mutate(id = as.numeric(factor(pid))) %>%
        arrange(id, time)  %>% mutate(titre = pmax(log10(titre), 0))  %>% rename(IgA = titre)   

    gambia_pvnt_both <- gambia_pvnt %>% left_join(gambia_iga)
    full_ids <- gambia_pvnt_both$pid[!complete.cases(gambia_pvnt_both)]
    gambia_pvnt_both_clean <- gambia_pvnt_both %>% filter(!pid %in% full_ids) %>% mutate(id = as.numeric(factor(pid, levels = unique(pid))))
    gambia_pvnt_both_clean
}

get_exposures_wave2 <- function() {
    odd_people <- c("17-197B", "07-077B", "41-481E", "43-509K", "34-399H", "41-483J", "43-512H")
    gambia_pvn_raw <- read.csv(file = here::here("data", "transvir", "Pseudovirus_data_V1_V2_V3.csv") ) %>% filter(!Participant_ID %in% odd_people)

    start_date <- gambia_pvn_raw %>% dplyr::select(Participant_ID, V1_date, contains("ptna_Ancestral_V1")) %>% pull(V1_date) %>% dmy %>% min  #"2021-03-02"

    gambia_inf_dates <- get_infection_dates()
    gambia_inf_dates_pd <- gambia_inf_dates %>% filter(!is.na(known_pre_delta)) %>% mutate(inf_days = as.numeric(ymd(known_pre_delta) - start_date) + 1) %>%
        dplyr::select(pid = Participant_ID, time = inf_days) %>% mutate(exposure_type = "predelta")
    gambia_inf_dates_d <- gambia_inf_dates %>% filter(!is.na(known_delta)) %>% mutate(inf_days = as.numeric(ymd(known_delta) - start_date) + 1) %>%
        dplyr::select(pid = Participant_ID, time = inf_days) %>% mutate(exposure_type = "delta")

    # Get vaccine dates
    gambia_meta <- read.csv(file = here::here("data", "transvir","full_demog_data.csv") )
    gambia_meta_short <- gambia_meta %>% dplyr::select(pid = Participant_ID, vax_date) %>% mutate(vax_days = as.numeric(ymd(vax_date) - ymd("2021-03-02")) + 1) %>% 
        select(pid, time = vax_days) %>% mutate(exposure_type = "vax")

    data_titre <- get_data_titre_model_wave2_pvnt_iga()

    known_exposures <- data_titre %>% select(pid, id) %>% unique %>% left_join(bind_rows(gambia_inf_dates_d, gambia_meta_short, gambia_inf_dates_pd) %>% filter(!is.na(time)) ) %>% 
        filter(!is.na(exposure_type))

    #known_exposures <- data_titre %>% left_join(gambia_inf_dates_pd) %>% left_join(gambia_inf_dates_d) %>% left_join(gambia_meta_short) %>%
    #    group_by(id) %>% mutate(t = row_number()) %>% ungroup %>% 
     #   pivot_wider(names_from = "t", values_from = c("titre", "time")) %>% as.data.frame %>% filter(pid != "34-399H") %>%
     #   mutate(inf_pd_time = case_when(inf_days_pd > time_1 & inf_days_pd <= time_2 ~ inf_days_pd, TRUE~-1)) %>% as.data.frame %>% 
     #   mutate(inf_d_time = case_when(inf_days_d > time_1 & inf_days_d <= time_2 ~ inf_days_d, TRUE~-1)) %>%
     #   mutate(vax_time = case_when(vax_days > time_1 & vax_days <= time_2 ~ vax_days, TRUE~-1)) %>%
     #   dplyr::select(pid, id, time_1, time_2, titre_1, titre_2, vax_days, inf_pd_time, inf_d_time, vax_time) 
        
    known_exposures
}


get_exp_prior_wave2 <- function() {

    odd_people <- c("17-197B", "07-077B", "41-481E", "43-509K", "34-399H", "41-483J", "43-512H")
    gambia_pvn_raw <- read.csv(file = here::here("data", "transvir", "Pseudovirus_data_V1_V2_V3.csv") ) %>% filter(!Participant_ID %in% odd_people)

    start_date <- gambia_pvn_raw %>% dplyr::select(Participant_ID, V1_date, contains("ptna_Ancestral_V1")) %>% pull(V1_date) %>% dmy %>% min  #"2021-03-02"
    end_date <- gambia_pvn_raw %>% dplyr::select(Participant_ID, V2_date, contains("ptna_Ancestral_V2")) %>% pull(V2_date) %>% dmy %>% max(na.rm = TRUE)  #"2021-03-02"

    gambia_daily_cases <- read.csv(file = here::here("data", "transvir", "gambia_covid_daily.csv") )
    gambia_daily_cases_ed <- gambia_daily_cases %>% filter(date >= start_date &  date < end_date) %>% mutate(
        date_time = as.numeric(ymd(date) - start_date)) %>% mutate(week = date_time %/% 14) %>% 
        summarise(total = sum(new_cases, na.rm = TRUE), .by = week) %>% 
        mutate(prop = total / sum(total)) %>% dplyr::select(week, prop)
    # Truncate the exposure prior such thae probab of infection is 0 before fornight 7 and zero after day 19
    gambia_daily_cases_ed <- gambia_daily_cases_ed %>% mutate(prop = case_when(week < 7 ~ 0, week > 19 ~ 0, TRUE ~ prop))
    exp_prior <- data.frame(
        day = 1:(nrow(gambia_daily_cases_ed)*14),
        week = gambia_daily_cases_ed$week %>% map(~rep(.x, 14)) %>% unlist,
        prob = gambia_daily_cases_ed$prop %>% map(~rep(.x / 14, 14)) %>% unlist
    )

    exp_prior
}