#**********************************************
# Read-in CSV data
#**********************************************

# Haematocrit
haem <- read.csv("haem_raw.csv")
haem <- haem[c("id", "haematocrit", "yearse", "gender", "weight", "age",
               "failure", "followyear")]
names(haem) <- c("id", "haematocrit", "years", "gender", "weight", "age",
                 "failure", "fuyears")
haem$haematocrit <- 100 * haem$haematocrit

# Proteinuria
prot <- read.csv("prot_raw.csv")
prot$proteinuria <- as.numeric(prot$proteinuria >= 1)
prot <- prot[c("PATID", "proteinuria", "yearse",
               "failure", "followyear")]
names(prot) <- c("id", "proteinuria", "years",
                 "failure", "fuyears")

# eGFR
gfr <- read.csv("gfr_raw.csv")
names(gfr) <- c("id", "failure", "failure10", "years", "gfr", "fuyears")
gfr <- gfr[c("id", "gfr", "years", "failure", "fuyears")]

#------------------------------------------------------------------

#**********************************************
# Re-order and restrict to common patients
#**********************************************

ind.gfr <- gfr$id %in% unique(haem$id)
gfr <- gfr[ind.gfr, ]
gfr <- gfr[order(gfr$id, gfr$years), ]
ind.prot <- prot$id %in% unique(haem$id)
prot <- prot[ind.prot, ]
prot <- prot[order(prot$id, prot$years), ]
haem <- haem[order(haem$id, haem$years), ]

#------------------------------------------------------------------

#**********************************************
# Covariate data
#**********************************************

# Weight
sp <- sapply(split(haem$weight, haem$id), "[", 1)
prot$weight <- rep(sp, tapply(prot$id, prot$id, length))
gfr$weight <- rep(sp, tapply(gfr$id, gfr$id, length))

# Age
sp <- sapply(split(haem$age, haem$id), "[", 1)
prot$age <- rep(sp, tapply(prot$id, prot$id, length))
gfr$age <- rep(sp, tapply(gfr$id, gfr$id, length))

# Gender
sp <- sapply(split(haem$gender, haem$id), "[", 1)
prot$gender <- rep(sp, tapply(prot$id, prot$id, length))
gfr$gender <- rep(sp, tapply(gfr$id, gfr$id, length))

# Event time (for heamoglobin data)
sp <- sapply(split(gfr$fuyears, gfr$id), "[", 1)
haem$fuyears <- rep(sp, tapply(haem$id, haem$id, length))
haem <- haem[, c(1:3, 7:8, 4:6)]

#------------------------------------------------------------------

#**********************************************
# Finalise datasets
#**********************************************

# Convert IDs to factors
gfr$id <- factor(gfr$id)
haem$id <- factor(haem$id)
prot$id <- factor(prot$id)

# Remove missing values
gfr <- gfr[!is.na(gfr$gfr), ]
haem <- haem[!is.na(haem$haematocrit), ]
prot <- prot[!is.na(prot$proteinuria), ]

# Remove biomarkers measured at or after event time
sp <- split(haem, haem$id)
lsp <- lapply(sp, function (x) x[x$years < x$fuyears, ])
haem <- do.call("rbind", lsp)
sp <- split(gfr, gfr$id)
lsp <- lapply(sp, function (x) x[x$years < x$fuyears, ])
gfr <- do.call("rbind", lsp)
sp <- split(prot, prot$id)
lsp <- lapply(sp, function (x) x[x$years < x$fuyears, ])
prot <- do.call("rbind", lsp)

# Clean the rownames
row.names(prot) <- seq_len(nrow(prot))
row.names(haem) <- seq_len(nrow(haem))
row.names(gfr) <- seq_len(nrow(gfr))

#------------------------------------------------------------------

#**********************************************
# Survival data
#**********************************************

surv <- data.frame(
    id = names(tapply(gfr$id, gfr$id, "[", 1)),
    fuyears = tapply(gfr$fuyears, gfr$id, "[", 1),
    failure = tapply(gfr$failure, gfr$id, "[", 1),
    weight = tapply(gfr$weight, gfr$id, "[", 1),
    age = tapply(gfr$age, gfr$id, "[", 1),
    gender = factor(tapply(gfr$gender, gfr$id, "[", 1), levels = 1:2,
                    labels = c("female", "male"))
)
row.names(surv) <- seq_len(nrow(surv))

renal <- list("prot" = prot, "haem" = haem, "gfr" = gfr, "surv" = surv)
save(renal, file = "renal.Rdata")
