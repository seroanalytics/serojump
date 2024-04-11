 #       date < dmy("01/07/2021") ~ "Pre-Delta",
 #       date >= dmy("01/07/2021") &  date < dmy("05/12/2021")~ "Delta",
 #       date >= dmy("05/12/2021") & date < dmy("01/06/2022")~ "Omicron"


get_data_titre_nih_2023 <- function() {
    # Read in the titre data and exposure data
    nih_hai_raw <- read.csv(file = here::here("data", "nih_2024", "serology.csv") )
    nih_hai_raw_24 <- nih_hai_raw %>% filter(year == 2023, subtype == "H3", virus == "A/Darwin/06/2021" ) 

    nih_hai_raw_24 %>% filter(vax_inf == "I")
    nih_bleed_raw <- read.csv(file = here::here("data", "nih_2024", "bleed-dates.csv") )
    nih_inf_dates_raw <- read.csv(file = here::here("data", "nih_2024", "postinf-bleed-dates.csv") )

    add_dates <- bind_rows(
        nih_bleed_raw %>% mutate(vax_inf = "V", date = dmy(date)),
        nih_inf_dates_raw %>% select(pid, year, day, date = bleed_date) %>% mutate(vax_inf = "I", date = ymd(date))
    )

    nih_hai_date_24 <- nih_hai_raw_24 %>% left_join(add_dates, by = c("pid", "year", "day", "vax_inf")) %>%
        select(pid, day, virus, titre, vax_inf, date) %>% filter(!is.na(date))

    start_date <- nih_hai_date_24$date %>% min
    end_date <- nih_hai_date_24$date %>% max

    data_titre <- nih_hai_date_24 %>% mutate(time = as.numeric(date - start_date) + 1) %>% 
        mutate(titre = log2(titre / 5)) %>% 
        mutate(pid = factor(pid), id = as.numeric(factor(pid))) %>% 
        select(pid, id, time, titre, biomarker = virus) %>% arrange(id, time) %>% filter(!is.na(titre))
    

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

    
    nih_inf_raw <- read_excel(path = here::here("data", "nih_2024", "2022_2023_Flu_Swabs.xlsx") )
    h3_exp <- nih_inf_raw %>% filter(year == 2023) %>% filter(grepl("H3", swab_virus)) %>% select(pid, samp_date) %>% 
        mutate(exposure_type = "h3_2023") %>% mutate(time = as.numeric(ymd(samp_date) - start_date) + 1) %>% 
        select(pid, time, exposure_type) %>% left_join(pid_key)

    # get post-infection serology
    nih_inf_dates_raw <- read.csv(file = here::here("data", "nih_2024", "postinf-bleed-dates.csv") )

    exposure_data <- bind_rows(vax_exp, h3_exp) 

    # Read in the exposure prior data
    flu_incidence <- read_excel(path = here::here("data", "nih_2024", "aus_flu_records_flunet.xlsx") )
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
