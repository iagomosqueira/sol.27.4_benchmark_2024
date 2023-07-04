# data.R - DESC
# /home/mosqu003/NEW/Backlog/sol274_BENCHMARK+ices/sol.27.4_benchmark/data.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(icesTAF)
mkdir("data")

library(FLSAM)

load("boot/data/data.rda")

type(indices$BTS) <- "number"
type(indices$SNS) <- "number"

save(stock, indices, refpts, file="data/data.rda", compress="xz")

# CONVERT to stockassessment

data <- FLSAM2SAM(FLStocks(residual=stock), indices)

save(data, file="data/sam.rda", compress="xz")


