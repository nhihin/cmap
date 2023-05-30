data(drugs10)

drugs10

drugs_orig <- unique(drugs)
drugs <- unique(tolower(drugs10))
# load dtlink_db_clue_sti.db stored in AnnotationHub
conn <- signatureSearch:::load_sqlite("EH3228")


dtlink_db <- dbGetQuery(conn, 'SELECT * FROM dtlink_db')
dtlink_clue <- dbGetQuery(conn, 'SELECT * FROM dtlink_clue')
dtlink_sti <- dbGetQuery(conn, 'SELECT * FROM dtlink_sti')
dtlink <- dbGetQuery(conn, 'SELECT * FROM dtlink')
dbDisconnect(conn)

dtl <- dtlink
dt <- unique(do.call(rbind, dtl))
dtlist <- split(dt$t_gn_sym, dt$drug_name)

## Extract Colnames

```{r}
all_drugbank_tables %>% lapply(colnames)

all_drugbank_tables$products %>% View

indication_tables <- list()
for(table in names(all_drugbank_tables)){
  if("indication_id" %in% colnames(all_drugbank_tables[[table]])){
    indication_tables[[table]] <- all_drugbank_tables[[table]]
  }
}

indication_tables %>% lapply(colnames)

indications <- indication_tables$indication_drugs %>%
  full_join(indication_tables$indication_conditions, by = "indication_id") %>%
  full_join(indication_tables$indication_categories, by = "indication_id")

View(indications)
dim(indications)
```

## Join indications to conditions

```{r}
indications2 <- indications %>%
  full_join(all_drugbank_tables$conditions, by = c("condition_id"="id"))

dim(indications2)
```

## Join indications onto drug info

```{r}
indications3 <- indications2 %>%
  full_join(all_drugbank_tables$drugs, by = c("drug_id"="id"))
dim(indications3)
```

## Join indications

```{r}

```




## Get condition tables

```{r}
condition_tables <- list()
for(table in names(all_drugbank_tables)){
  if("condition_id" %in% colnames(all_drugbank_tables[[table]])){
    condition_tables[[table]] <- all_drugbank_tables[[table]]
  }
}
condition_tables %>% lapply(colnames)
condition_tables$indication_conditions %>% View
```

## Tables with ID

```{r}
anyid_tables <- list()
for(table in names(all_drugbank_tables)){
  if("id" %in% colnames(all_drugbank_tables[[table]])){
    anyid_tables[[table]] <- all_drugbank_tables[[table]]
  }
}
anyid_tables %>% lapply(colnames)
anyid_tables$structured_indications%>% View
```


```{r}
indications <- indication_tables %>% purrr::reduce(dplyr::left_join, by = "indication_id")
View(indications)
```



```{r}
conditions <- file.path(drugbank_dir, "conditions.csv") %>%
  readr::read_csv()

View(conditions)

drugs <- file.path(drugbank_dir, "drugs.csv") %>%
  readr::read_csv()

View(drugs)

drug_synonyms <- file.path(drugbank_dir, "drug_synonyms.csv") %>%
  readr::read_csv()

View(drug_synonyms)
```

## Join relevant tables

```{r}
drugs2Conditions <- dplyr::inner_join(drugs, conditions, by = "id")
```

