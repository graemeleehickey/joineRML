library(ggplot2)
library(joineRML)
library(reshape2)
library(hexSticker)

data(pbc2)
pbc2 <- subset(pbc2, id != "228")
pbc2 <- subset(pbc2, select = c(id, year, serBilir, albumin))
pbc2 <- melt(pbc2, id.vars = c("id", "year"))
#pbc2$variable <- relevel(pbc2$variable, ref = "albumin")

p <- ggplot(aes(x = year, y = log(value)), data = pbc2) +
  geom_line(aes(group = id), alpha = 0.5) +
  geom_line(data = subset(pbc2, id == "11"), colour = "red") +
  facet_wrap( ~ variable, scales = "free_y") +
  labs(x = "", y = "") +
  theme_void() +
  theme(strip.background = element_blank(),
        strip.text = element_blank())

sticker(p, package = "joineRML",
        p_size = 8,
        s_x = 1, s_y = .8, s_width = 1.5, s_height = 0.9,
        #p_color = "#FFFFFFDD",
        h_color = "#C76730",
        h_fill = "#3090C7")
