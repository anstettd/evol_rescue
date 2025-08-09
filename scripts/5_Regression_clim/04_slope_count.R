##################################################################################
## Get per population counts for parts of the S distribution
## Author Daniel Anstett
## 
## 
## Last Modified Oct 31, 2023
###################################################################################

obs_env_unique <- read_csv("data/snp_change_2/slope_obs_all_unique.csv") %>% 
  filter(SE<5) %>% mutate(abs_slope = abs(Slope)) %>% filter(Site!=12)



slope_counts <- obs_env_unique %>%
  mutate(
    slope_bin = case_when(
      Slope >= 0 & Slope <= 0.2 ~ "0–0.2",
      Slope > 0.2 & Slope <= 0.5 ~ "0.2–0.5",
      Slope > 0.5 ~ ">0.5",
      TRUE ~ "Missing/Out of range"
    )
  ) %>%
  group_by(Site, slope_bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = slope_bin,
    values_from = count,
    values_fill = 0
  )

slope_counts

min(slope_counts$"0–0.2")
max(slope_counts$"0–0.2")

min(slope_counts$"0.2–0.5")
max(slope_counts$"0.2–0.5")

min(slope_counts$"0.5")
max(slope_counts$"0.5")



