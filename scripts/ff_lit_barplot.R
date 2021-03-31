

library("viridis")

dta <- read.csv2("raw-data/lit_rev_data.csv")

p_increase <- table(dta$p_increase)
QoI <- table(dta$QoI)
QoI_type <- table(dta$QoI_type)
QoI_plot <- table(dta$QoI_plot)
only_direction <- table(dta$only_direction)

barplot(p_increase)

plot <- data.frame("direct_only" = only_direction["1"],
                   "p_increase" = p_increase["1"],
                   "QoI" = QoI["1"])


plot <- as.matrix(plot)

height_EV <- QoI_type["EV"]
height_FD <- QoI_type["FD"]
height_PV <- QoI_type["PV"]


bgcol <- "grey95"
barcol <- "grey50"

QoI_cols <- viridis(4)


pdf("figures/lit_plot.pdf",
    width = 15,
    height = 5)

par(mar = c(5, 6, 1, 1))
barplot(rev(plot),
        horiz = T,
        names.arg = c("Quantity of\n Interest",
                      "% increase",
                      "direction\n only"),
        las = 1,
        col = barcol,
        border = barcol,
        xlim = c(0,32),
        axes = F)
axis(1, 
     at = seq(0, 25, 5))
polygon(x = c(0, 0, height_EV, height_EV),
        y = c(0.35+0.15, 0.35-0.1, 0.35-0.1, 0.35+0.15),
        col = QoI_cols[1],
        border = QoI_cols[1])
polygon(x = c(0, 0, height_FD,  height_FD),
        y = c(0.7+0.125, 0.7-0.125, 0.7-0.125, 0.7+0.125),
        col = QoI_cols[2],
        border = QoI_cols[2])
polygon(x = c(0, 0, height_PV, height_PV),
        y = c(1.05+0.1, 1.05-0.15, 1.05-0.15, 1.05+0.1),
        col = QoI_cols[3],
        border = QoI_cols[3])


text(x = rep(26.5, 3),
     y = c(0.35, 0.7, 1.05),
     labels = c("Expected Values", "First Difference", "Predicted Values"),
     pos = 4,
     cex = 1.3,
     col = QoI_cols)

dev.off()
