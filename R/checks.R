check_titre <- function(data_titre) {
    check_time(data_titre)
}

check_time <- function(data_titre) {
    time_diff_too_small <- data.frame(
        id = 1:(data_titre$id %>% max),
        start_time = data_titre %>% group_by(id) %>% filter(time == min(time)) %>% unique %>% pull(time),
        end_time = data_titre %>% group_by(id) %>% filter(time == max(time)) %>% unique %>% pull(time)
    ) %>% mutate(time_diff = end_time - start_time) %>% filter(time_diff <= 7)

    if(nrow(time_diff_too_small) > 0) {
        stop(paste("Error: Time difference between first and last titre measurement is less than 7 days for the following id(s):", paste(time_diff_too_small$id, collapse = ", "), ". Please remove these values and recode the data inputs."))
    }
}