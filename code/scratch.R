### TBC Test with endothelial cells for now



```{r}
oldColnames <- pb$endothelial %>% colnames
newColnames <- data.frame(sample_ID = oldColnames) %>%
  dplyr::left_join(map_df, by = "sample_ID") %>%
  dplyr::pull(sample_treatment)

colnames(pb$endothelial) <- newColnames

endo_veh <- pb$endothelial %>%
  as.data.frame %>%
  dplyr::select(contains("Vehicle"))

endo_ppf <- pb$endothelial%>%
  as.data.frame %>%
  dplyr::select(contains("Pl-47615"))
```

## Z-scores

```{r}
endo_veh_z <- apply(endo_veh, 2, cmapR::robust_zscore)
endo_ppf_z <- apply(endo_ppf, 2, cmapR::robust_zscore)

removeLargeZscores <- function(x){

  x[x > 10] <- 10
  x[x < -10] <- -10

  return(x)

}

endo_veh_z <- apply(endo_veh_z, 2, removeLargeZscores)
endo_ppf_z <- apply(endo_ppf_z, 2, removeLargeZscores)
```

## Distil (perform weighted average of columns)

```{r}
endo_veh_z_dis <- cmapR::distil(endo_veh_z,
                                dimension = "col")
endo_ppf_z_dis <- cmapR::distil(endo_ppf_z,
                                dimension = "col")

```

### Initial comparison between PPF and Vehicle

```{r}
diff_ppf_and_veh <- endo_ppf_z_dis$values -
  endo_veh_z_dis$values

diff_ppf_and_veh %>% hist(breaks=100)
summary(diff_ppf_and_veh)

endo_veh_z_dis$values %>% length
endo_veh_z_dis$weights %>%
  hist(breaks=50)

endo_ppf_z_dis$weights %>%
  hist(breaks=50)
```


## Each cell type


