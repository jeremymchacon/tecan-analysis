

rm(list = ls())
source('tecan_functions.R')
OD_path = "25oct16_OD.csv"
CFP_path = "25oct16_cfp.csv"
RFP_path = "25oct16_rfp.csv"
YFP_path = "25oct16_yfp.csv"
map_path = "25oct16_platemap_short.csv"

# note -- had to remove the extra rows
OD = read_tecan(OD_path, "OD")
CFP = read_tecan(CFP_path, "CFP")[, c("cycle","well","CFP")]
RFP = read_tecan(RFP_path, "RFP")[, c("cycle","well","RFP")]
YFP = read_tecan(YFP_path, "YFP")[, c("cycle","well","YFP")]

# easy enough to merge together
OD = merge(merge(merge(OD, CFP), YFP), RFP)

OD = apply_well_map(OD, map_path, "treatment")
OD <- filter(OD, !is.na(treatment))
OD <- filter(OD, treatment != "")
OD <- get_all_growth_rates_per_hr(OD, "OD", 12, 0.98,gr_name = "window_gr", cycle_name = "window_cycle")

get_baranyi_rate(OD, "B2", "OD")
plot_well(OD, "B2")

OD <- get_all_baranyi_rates(OD, "OD", "baranyi_gr","baranyi_lag", "baranyi_ymax", "baranyi_y0",30)
OD$hours = OD$time_s / 3600

OD <- generate_baranyi_prediction(results,  "baranyi_gr","baranyi_lag", "baranyi_ymax", "baranyi_y0")


OD %>%
  filter(!(treatment %in% c("a", "s"))) %>%
  ggplot(aes(x = time_s / 3600, y = OD, group = well))+
  geom_point(shape = 1, color = "#AAAAAA")+
  geom_line(aes( y= baranyi_pred))+
  facet_wrap(well~treatment)+
  theme_bw()

OD %>%
  filter(!(treatment %in% c("a", "s"))) %>%
  ggplot(aes(x = time_s / 3600, y = OD, group = well))+
  geom_point(size = 0.2, color = "#AAAAAA")+
  geom_line(aes( y= baranyi_pred))+
  facet_wrap(~treatment)+
  theme_bw()

OD %>%
  filter(!(treatment %in% c("a", "s"))) %>%
  ggplot(aes(x = time_s / 3600, y = OD, group = well))+
  geom_point(size = 0.2, color = "#AAAAAA")+
  geom_line(aes( y= baranyi_pred))+
  facet_wrap(~treatment)+
  scale_y_log10()+
  theme_bw()

a = OD %>%
  filter(!(treatment %in% c("a", "s"))) %>%
  group_by(well, treatment) %>%
  summarize(gr = first(window_gr)) %>%
  ggplot(aes(x = treatment, y= gr))+
  geom_jitter(width = 0.1)+
  theme_bw()+
  labs(y = "linear fit in window gr")

b = OD %>%
  filter(!(treatment %in% c("a", "s"))) %>%
  group_by(well, treatment) %>%
  summarize(gr = first(baranyi_gr)) %>%
  ggplot(aes(x = treatment, y= gr))+
  geom_jitter(width = 0.1)+
  theme_bw()+
  labs(y = "baryani gr")

gridExtra::grid.arrange(a,b, nrow = 1)


a = OD %>%
  filter(!(treatment %in% c("a", "s"))) %>%
  group_by(well, treatment) %>%
  summarize(gr = first(window_gr)) %>%
  group_by(treatment) %>%
  summarize(r = mean(gr), r_sd = sd(gr), n = length(gr)) %>%
  ggplot(aes(x = treatment, y= r,
             ymin = r - r_sd / sqrt(n-1),
             ymax = r + r_sd / sqrt(n-1)))+
  geom_point()+
  geom_errorbar()+
  theme_bw()+
  labs(y = "linear fit in window gr")

b = OD %>%
  filter(!(treatment %in% c("a", "s"))) %>%
  group_by(well, treatment) %>%
  summarize(gr = first(baranyi_gr)) %>%
  group_by(treatment) %>%
  summarize(r = mean(gr), r_sd = sd(gr), n = length(gr)) %>%
  ggplot(aes(x = treatment, y= r,
             ymin = r - r_sd / sqrt(n-1),
             ymax = r + r_sd / sqrt(n-1)))+
  geom_point()+
  geom_errorbar()+
  theme_bw()+
  labs(y = "baryani gr")

gridExtra::grid.arrange(a,b, nrow = 1)