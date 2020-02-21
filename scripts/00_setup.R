# Simulating Quantities of Interest when the Dependent Variable is logged

##### Setup

## ----setup----
## Save package names as a vector of strings
pkgs <-
  c("MASS",
    "scales")

## Install uninstalled packages
lapply(pkgs[!(pkgs %in% installed.packages())], install.packages)

## Load all packages to library and adjust options
print(sapply(pkgs, require, character.only = TRUE))

rm(list = ls())