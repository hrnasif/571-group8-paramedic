library(dplyr)
absolute <- readr::read_csv("data/qpcr_measurements.csv")
relative <- readr::read_csv("data/raw_counts.csv")

#### Checking that the data sets are wholly aligned before subsetting
# Particpant and Hours in Study are the only features that appear in both tables

sum(setdiff(absolute$Participant, relative$Participant)) == 0
sum(setdiff(absolute$Hours_In_Study, relative$Hours_In_Study)) == 0


######### Obtaining absolute abundance data #
absolute <- absolute[,-1] %>% rename(subject_id = Participant) %>% 
  group_by(subject_id) %>% mutate(sample_id = 1:n())

drop_cols = c("Hours_In_Study", "Total_16S")
absolute <- absolute %>% relocate(subject_id, sample_id, Gardnerella_vaginalis,
                                         Lactobacillus_crispatus, Lactobacillus_iners,
                                         Lactobacillus_jensenii, Megasphaera, BVAB2,
                                         Atopobium_vaginae) %>% select(-one_of(drop_cols))

qpcr_ind <- 3:ncol(absolute)
colnames(absolute)[qpcr_ind] <- colnames(example_qPCR_data[2:ncol(example_qPCR_data)])

# Scaling as done in paramedic data
absolute[,qpcr_ind] <- round(absolute[,qpcr_ind]/1000)


########### Obtaining relative abundance data #
colnames(relative) <- gsub("\\[|\\]|\\(|\\)", "", colnames(relative))
colnames(relative) <- gsub("\\s|\\/", "_", colnames(relative))

relative <- relative[,-1] %>% rename(subject_id = Participant, 
                                Lactobacillus_crispatus = Lactobacillus_crispatus_helveticus,
                                Megasphaera.genomosp..type_1 =  Megasphaera_genomosp._type_1,
                                BVAB2..species = BVAB2_species) %>% 
  group_by(subject_id) %>% mutate(sample_id = 1:n()) %>% 
  relocate(subject_id, sample_id, Gardnerella_vaginalis,
                  Lactobacillus_crispatus, Lactobacillus_iners,
                  Lactobacillus_jensenii, Megasphaera.genomosp..type_1, BVAB2..species,
                  Atopobium_vaginae) %>% 
  select(-one_of("Hours_In_Study"))

colnames(relative)[qpcr_ind] = colnames(absolute)[qpcr_ind]


### Removing observations where M < 1000
M_min = 1000
M <- apply(relative[,-c(1,2)], 1, sum)

absolute <- absolute[M >= M_min, ]
relative <- relative[M >= M_min, ]

### Subsetting
set.seed(123)
subjs <- sample(unique(absolute$subject_id), 3, replace = F)

absolute_dt <- absolute %>% filter(subject_id %in% subjs & sample_id <= 10)
relative_dt <- relative %>% filter(subject_id %in% subjs & sample_id <= 10)

saveRDS(absolute_dt, file = "data/group8_qpcr_data.rda")
saveRDS(relative_dt, file = "data/group8_16S_data.rda")

# reference
load("data/example_16S_data.rda")
load("data/example_qPCR_data.rda")
