 #       date < dmy("01/07/2021") ~ "Pre-Delta",
 #       date >= dmy("01/07/2021") &  date < dmy("05/12/2021")~ "Delta",
 #       date >= dmy("05/12/2021") & date < dmy("01/06/2022")~ "Omicron"


get_data_titre_nih_prior2022 <- function(year, subtype_i) {
    require(lubridate)
    #year_i <- 2020
    #subtype_i <- "H1"

    nih_hai_raw <- read.csv(file = here::here("data", "data 2024", "serology.csv") )
    nih_hai_raw_20 <- nih_hai_raw %>% filter(year == year_i, subtype == subtype_i) %>% 
        select(pid, year, day, vax_inf, subtype, virus, titre) %>%
         mutate(titre = log2(titre / 5)) %>%
        pivot_wider(names_from = "virus", values_from = "titre") 

    biomarkers <- names(nih_hai_raw_20)[c(6, 7)]

    nih_hai_date_20 <- nih_hai_raw_20 %>% left_join(add_dates, by = c("pid", "year", "day", "vax_inf")) %>%
    select(c(pid, day, vax_inf, date, biomarkers)) %>% filter(!is.na(date))

    nih_bleed_raw <- read.csv(file = here::here("data", "nih_2024_XX", "bleed-dates.csv") )

    add_dates <- bind_rows(
        nih_bleed_raw %>% mutate(vax_inf = "V", date = dmy(date)),
    )

    start_date <- nih_hai_date_20$date %>% min
    end_date <- nih_hai_date_20$date %>% max

    data_titre <- nih_hai_date_20 %>% mutate(time = as.numeric(date - start_date) + 1) %>% 
        mutate(pid = factor(pid), id = as.numeric(factor(pid))) %>% 
        select(1, 8, 7, 5, 6) %>% arrange(id, time) %>% filter_all(all_vars(!is.na(.))) %>%
        summarise(across(all_of(biomarkers), mean), .by = c(pid, id, time)) 


    time_diff_too_small <- data.frame(
        id = 1:(data_titre$id %>% max),
        start_time = data_titre %>% group_by(id) %>% filter(time == min(time)) %>% unique %>% pull(time),
        end_time = data_titre %>% group_by(id) %>% filter(time == max(time)) %>% unique %>% pull(time)
    ) %>% mutate(time_diff = end_time - start_time) %>% filter(time_diff <= 7) %>% pull(id)

    data_titre <- data_titre %>% filter(!id %in% time_diff_too_small) %>% mutate(pid = factor(pid), id = as.numeric(factor(pid))) 

    pid_key <- data_titre %>% select(pid, id) %>% distinct
    # Get exposures
    # 2023 flu summary: https://www.health.gov.au/sites/default/files/2023-12/aisr-2023-national-influenza-season-summary.pdf
    vax_exp <- data_titre %>% group_by(id) %>% filter(time == min(time)) %>% ungroup %>% select(pid, id, time) %>% mutate(exposure_type = "vax")

    exposure_data <- vax_exp

    require(readxl)
    flu_incidence <- read_excel(path = here::here("data", "nih_2024_XX", "aus_flu_records_flunet.xlsx") )
    week_cases <- flu_incidence %>% select(ISO_SDATE, INF_ALL, INF_NEGATIVE) %>% 
        mutate(date = ymd(substr(ISO_SDATE, 1, 10))) %>%
        mutate(INF_ALL = as.numeric(INF_ALL)) %>%
        mutate(INF_NEGATIVE = as.numeric(INF_NEGATIVE)) %>%
        mutate(INF_PROB = INF_ALL / INF_NEGATIVE) %>%
        filter(date > start_date - 14 & date < end_date ) %>%
        select(date, INF_PROB) %>% mutate(INF_PROB = ifelse(is.na(INF_PROB), 0, INF_PROB))
        
    exp_prior <- data.frame(
        day = 1:(nrow(week_cases)*7),
        week = week_cases$date %>% map(~rep(.x, 7)) %>% unlist,
        prob = week_cases$INF_PROB %>% map(~rep(.x / 7, 7)) %>% unlist
    )

     list(data_titre = data_titre, exposure_data  = exposure_data, exp_prior = exp_prior)
   
}

get_data_titre_nih_2023_h3 <- function() {
    # Read in the titre data and exposure data
    nih_hai_raw <- read.csv(file = here::here("data", "data 2024", "serology.csv") )
    nih_hai_raw_24 <- nih_hai_raw %>% filter(year == 2023, subtype == "H3") %>% 
        select(pid, year, day, vax_inf, subtype, virus, titre) %>%
         mutate(titre = log2(titre / 5)) %>%
        pivot_wider(names_from = "virus", values_from = "titre") 

    biomarkers <- names(nih_hai_raw_24)[c(6, 7)]

    # Find known infections
    nih_bleed_raw <- read.csv(file = here::here("data", "nih_2024_XX", "bleed-dates.csv") )
    nih_inf_dates_raw <- read.csv(file = here::here("data", "nih_2024_XX", "postinf-bleed-dates.csv") )

    add_dates <- bind_rows(
        nih_bleed_raw %>% mutate(vax_inf = "V", date = dmy(date)),
        nih_inf_dates_raw %>% select(pid, year, day, date = bleed_date) %>% mutate(vax_inf = "I", date = ymd(date))
    )

    nih_hai_date_24 <- nih_hai_raw_24 %>% left_join(add_dates, by = c("pid", "year", "day", "vax_inf")) %>%
        select(1, 3, 4, 8, 6, 7) %>% filter(!is.na(date))

    start_date <- nih_hai_date_24$date %>% min
    end_date <- nih_hai_date_24$date %>% max

    data_titre <- nih_hai_date_24 %>% mutate(time = as.numeric(date - start_date) + 1) %>% 
        mutate(pid = factor(pid), id = as.numeric(factor(pid))) %>% 
        select(1, 8, 7, 5, 6) %>% arrange(id, time) %>% filter_all(all_vars(!is.na(.))) %>%
        summarise(across(all_of(biomarkers), mean), .by = c(pid, id, time)) 
    
    data_titre %>% group_by(id) %>% filter(time == min(time)) %>% unique %>% as.data.frame %>% group_by(id) %>% mutate(r = row_number()) %>% as.data.frame %>% filter(r == 2)

    time_diff_too_small <- data.frame(
        id = 1:(data_titre$id %>% max),
        start_time = data_titre %>% group_by(id) %>% filter(time == min(time)) %>% unique %>% pull(time),
        end_time = data_titre %>% group_by(id) %>% filter(time == max(time)) %>% unique %>% pull(time)
    ) %>% mutate(time_diff = end_time - start_time) %>% filter(time_diff <= 7) %>% pull(id)

    data_titre <- data_titre %>% filter(!id %in% time_diff_too_small) %>% mutate(pid = factor(pid), id = as.numeric(factor(pid))) 


    pid_key <- data_titre %>% select(pid, id) %>% distinct
    # Get exposures
    # 2023 flu summary: https://www.health.gov.au/sites/default/files/2023-12/aisr-2023-national-influenza-season-summary.pdf
    vax_exp <- data_titre %>% group_by(id) %>% filter(time == min(time)) %>% ungroup %>% select(pid, id, time) %>% mutate(exposure_type = "vax")

    require(readxl)

    
    nih_inf_raw <- read_excel(path = here::here("data", "nih_2024_XX", "2022_2023_Flu_Swabs.xlsx") )
    h3_exp <- nih_inf_raw %>% filter(year == 2023) %>% filter(grepl("H3", swab_virus)) %>% select(pid, samp_date) %>% 
        mutate(exposure_type = "h3_2023") %>% mutate(time = as.numeric(ymd(samp_date) - start_date) + 1) %>% 
        select(pid, time, exposure_type) %>% left_join(pid_key)

    # get post-infection serology
    nih_inf_dates_raw <- read.csv(file = here::here("data", "nih_2024_XX", "postinf-bleed-dates.csv") )

    exposure_data <- bind_rows(vax_exp, h3_exp) 

    # Read in the exposure prior data
    flu_incidence <- read_excel(path = here::here("data", "nih_2024_XX", "aus_flu_records_flunet.xlsx") )
    week_cases <- flu_incidence %>% select(ISO_SDATE, INF_ALL, INF_NEGATIVE) %>% 
        mutate(date = ymd(substr(ISO_SDATE, 1, 10))) %>%
        mutate(INF_ALL = as.numeric(INF_ALL)) %>%
        mutate(INF_NEGATIVE = as.numeric(INF_NEGATIVE)) %>%
        mutate(INF_PROB = INF_ALL / INF_NEGATIVE) %>%
        filter(date > start_date - 14 & date < end_date ) %>%
        select(date, INF_PROB)
        
    exp_prior <- data.frame(
        day = 1:(nrow(week_cases)*7),
        week = week_cases$date %>% map(~rep(.x, 7)) %>% unlist,
        prob = week_cases$INF_PROB %>% map(~rep(.x / 7, 7)) %>% unlist
    )
    list(data_titre = data_titre, exposure_data  = exposure_data, exp_prior = exp_prior)
}


get_data_titre_nih_2023_h1 <- function() {
    # Read in the titre data and exposure data
    nih_hai_raw <- read.csv(file = here::here("data", "data 2024", "serology.csv") )
    nih_hai_raw_24_h1 <- nih_hai_raw %>% filter(year == 2023, subtype == "H1") %>% 
        select(pid, year, day, vax_inf, subtype, virus, titre) %>%
         mutate(titre = log2(titre / 5)) %>%
        pivot_wider(names_from = "virus", values_from = "titre") 

    biomarkers_h1 <- names(nih_hai_raw_24_h1)[c(6, 7)]

    nih_hai_raw_24_h1 %>% filter(vax_inf == "I")
    nih_bleed_raw <- read.csv(file = here::here("data", "nih_2024_XX", "bleed-dates.csv") )
    nih_inf_dates_raw <- read.csv(file = here::here("data", "nih_2024_XX", "postinf-bleed-dates.csv") )

    add_dates <- bind_rows(
        nih_bleed_raw %>% mutate(vax_inf = "V", date = dmy(date)),
        nih_inf_dates_raw %>% select(pid, year, day, date = bleed_date) %>% mutate(vax_inf = "I", date = ymd(date))
    )

    nih_hai_date_24_h1 <- nih_hai_raw_24_h1 %>% left_join(add_dates, by = c("pid", "year", "day", "vax_inf")) %>%
        select(1, 3, 4, 8, 6, 7) %>% filter(!is.na(date))

    start_date <- nih_hai_date_24_h1$date %>% min
    end_date <- nih_hai_date_24_h1$date %>% max

    data_titre_h1 <- nih_hai_date_24_h1 %>% mutate(time = as.numeric(date - start_date) + 1) %>% 
        mutate(pid = factor(pid), id = as.numeric(factor(pid))) %>% 
        select(1, 8, 7, 5, 6) %>% arrange(id, time) %>% filter_all(all_vars(!is.na(.))) %>%
        summarise(across(all_of(biomarkers_h1), mean), .by = c(pid, id, time)) 
    
    data_titre_h1 %>% group_by(id) %>% filter(time == min(time)) %>% unique %>% as.data.frame %>% group_by(id) %>% mutate(r = row_number()) %>% as.data.frame %>% filter(r == 2)

    time_diff_too_small <- data.frame(
        id = 1:(data_titre_h1$id %>% max),
        start_time = data_titre_h1 %>% group_by(id) %>% filter(time == min(time)) %>% unique %>% pull(time),
        end_time = data_titre_h1 %>% group_by(id) %>% filter(time == max(time)) %>% unique %>% pull(time)
    ) %>% mutate(time_diff = end_time - start_time) %>% filter(time_diff <= 7) %>% pull(id)

    data_titre_h1 <- data_titre_h1 %>% filter(!id %in% time_diff_too_small) %>% mutate(pid = factor(pid), id = as.numeric(factor(pid))) 


    pid_key_h1 <- data_titre_h1 %>% select(pid, id) %>% distinct
    # Get exposures
    # 2023 flu summary: https://www.health.gov.au/sites/default/files/2023-12/aisr-2023-national-influenza-season-summary.pdf
    vax_exp_h1 <- data_titre_h1 %>% group_by(id) %>% filter(time == min(time)) %>% ungroup %>% select(pid, id, time) %>% mutate(exposure_type = "vax")

    require(readxl)

    nih_inf_raw <- read_excel(path = here::here("data", "nih_2024_XX", "2022_2023_Flu_Swabs.xlsx") )
    h1_exp <- nih_inf_raw %>% filter(year == 2023) %>% filter(grepl("H1", swab_virus)) %>% select(pid, samp_date) %>% 
        mutate(exposure_type = "h1_2023") %>% mutate(time = as.numeric(ymd(samp_date) - start_date) + 1) %>% 
        select(pid, time, exposure_type) %>% left_join(pid_key_h1)

    # get post-infection serology
    nih_inf_dates_raw <- read.csv(file = here::here("data", "nih_2024_XX", "postinf-bleed-dates.csv") )

    exposure_data_h1 <- bind_rows(vax_exp_h1, h1_exp) 

    # Read in the exposure prior data
    flu_incidence <- read_excel(path = here::here("data", "nih_2024_XX", "aus_flu_records_flunet.xlsx") )
    week_cases <- flu_incidence %>% select(ISO_SDATE, INF_ALL, INF_NEGATIVE) %>% 
        mutate(date = ymd(substr(ISO_SDATE, 1, 10))) %>%
        mutate(INF_ALL = as.numeric(INF_ALL)) %>%
        mutate(INF_NEGATIVE = as.numeric(INF_NEGATIVE)) %>%
        mutate(INF_PROB = INF_ALL / INF_NEGATIVE) %>%
        filter(date > start_date - 14 & date < end_date ) %>%
        select(date, INF_PROB)
        
    exp_prior <- data.frame(
        day = 1:(nrow(week_cases)*7),
        week = week_cases$date %>% map(~rep(.x, 7)) %>% unlist,
        prob = week_cases$INF_PROB %>% map(~rep(.x / 7, 7)) %>% unlist
    )
    list(data_titre = data_titre_h1, exposure_data  = exposure_data_h1, exp_prior = exp_prior)
}
