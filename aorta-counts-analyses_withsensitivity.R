#### PACKAGES ####
library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)
library(tidyverse)
library(readxl)

#### READ FILE ####
setwd("~/Desktop/TERVA/TERVA2DATA/VALIDATIONS_2023/Miran tulokset 25022024/")
cell_counts_adipose <- read_excel(
  "~/Desktop/TERVA/TERVA2DATA/VALIDATIONS_2023/Miran tulokset 25022024/results-reannotated/cell_counts_adipose.xlsx"
)

summary(cell_counts_adipose)

#### NEW VARIABLES ####
cell_counts_adipose <- cell_counts_adipose %>%
  mutate(
    Status = case_when(
      startsWith(`ROI ID`, "LD") ~ "Obese",
      startsWith(`ROI ID`, "PL") ~ "Non-obese",
      TRUE ~ NA_character_
    )
  )

adipose <- cell_counts_adipose %>%
  mutate(
    # main analysis variables you already used
    postototal = `positive counts` / total_count,         # fraction PI16+ of all detected nuclei/cells
    postoarea  = (`positive counts` / area_mm2) / 100,    # your existing scaled density
    
    # sensitivity analysis (more interpretable units)
    pi16_per_1000 = 1000 * `positive counts` / total_count,
    pi16_per_mm2  = `positive counts` / area_mm2
  ) %>%
  filter(!is.na(Status))

#### POS OF TOTAL (fraction PI16+/total) ####
postotal.beta <- glmmTMB(
  postototal ~ Status,
  data = adipose,
  family = beta_family(link = "logit")
)

sim.base <- simulateResiduals(fittedModel = postotal.beta, plot = TRUE)
testCategorical(sim.base, catPred = adipose$Status)

Anova(postotal.beta, type = "III")
emmeans(postotal.beta, "Status", type = "response")

#### POS OF AREA ####
postoarea.normal <- glmmTMB(
  postoarea ~ Status,
  data = adipose,
  family = gaussian(link = "identity")
)
sim.base <- simulateResiduals(fittedModel = postoarea.normal, plot = TRUE)
testCategorical(sim.base, catPred = adipose$Status)

postoarea.lognormal <- glmmTMB(
  postoarea ~ Status,
  data = adipose,
  family = lognormal(link = "log")
)
sim.base <- simulateResiduals(fittedModel = postoarea.lognormal, plot = TRUE)
testCategorical(sim.base, catPred = adipose$Status)

postoarea.beta <- glmmTMB(
  postoarea / 100 ~ Status,
  data = adipose,
  family = beta_family(link = "logit")
)
sim.base <- simulateResiduals(fittedModel = postoarea.beta, plot = TRUE)
testCategorical(sim.base, catPred = adipose$Status)

Anova(postoarea.beta, type = "III")
emm_area <- emmeans(postoarea.beta, "Status", type = "response")
emm_area
p1 <- plot(emm_area, cex.lab = 5, cex.axis = 5, comparisons = TRUE)

ggsave(plot = p1, filename = "PI16groupmeanresponse.tiff",
       device = "tiff", dpi = 600, width = 15, height = 4, units = "cm")
ggsave(plot = p1, filename = "PI16groupmeanresponse.png",
       device = "png", dpi = 600, width = 15, height = 4, units = "cm")

#### SENSITIVITY ANALYSIS: PI16+ per 1,000 nuclei/cells ####
m_per1000 <- glmmTMB(
  pi16_per_1000 ~ Status,
  data = adipose,
  family = gaussian()
)
Anova(m_per1000, type = "III")
emmeans(m_per1000, "Status", type = "response")

#### SUPPLEMENTARY FIGURE (area-normalized + nucleus-normalized) ####
p_area <- ggplot(adipose, aes(Status, pi16_per_mm2)) +
  geom_point(position = position_jitter(width = 0.1), size = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, , color = "blue", shape = 18) +
  labs(y = "PI16+ cells per mm²", x = NULL) +
  theme_classic()

p_nuc <- ggplot(adipose, aes(Status, pi16_per_1000)) +
  geom_point(position = position_jitter(width = 0.1), size = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "blue", shape = 18) +
  labs(y = "PI16+ cells per 1,000 nuclei", x = NULL) +
  theme_classic()

# side-by-side:
library(patchwork)
p_sens <- p_area + p_nuc
ggsave("PI16_sensitivity_area_vs_nuclei.png", p_sens, dpi = 600, width = 16, height = 7, units = "cm")
