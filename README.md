# geneticAnalysisR 1.2.0

An R wrapper for **BLUPF90+**: pedigree preparation, Gaussian EM imputation,
REML estimation, full-pedigree breeding values, prediction SEs, accuracy,
heritability, genetic correlations and inbreeding. This release removes the
sommer implementation and dependency. Model fitting uses the external BLUPF90+
program; the R package prepares inputs, runs isolated jobs and validates outputs.

## Install

```r
install.packages(c("Matrix", "processx", "curl", "openssl"))
install.packages("geneticAnalysisR_1.2.0.tar.gz", repos=NULL, type="source")
library(geneticAnalysisR)

# Explicit download; fitting never downloads or installs programs.
info <- install_blupf90()
# For a custom location, retain the path for later sessions:
Sys.setenv(BLUPF90_PLUS=info$executable)
```

The package source and repository do not contain BLUPF90+ or Intel binaries.
Their licensing and runtime requirements apply separately. Publisher name/email
in DESCRIPTION remain placeholders for the repository owner to fill before release.

## Operating systems

| Platform | R wrapper | External executable/setup |
|---|---|---|
| Linux x86-64 | Supported and tested here | Official ELF binary; optional private pinned Intel runtime installation through Python/pip. |
| Windows x86-64 | Portable process runner and paths; CI configured | Official `.exe`; install compatible Intel Fortran/MKL runtime DLLs separately. |
| Intel macOS | Portable runner and paths; CI configured | Official Intel Mach-O binary, if available; compatible Intel runtime required. |
| Apple Silicon macOS | R code is portable | Intel binary requires Rosetta and compatible runtimes. Explicitly opt in with `install_blupf90(allow_emulation=TRUE)`. A locally supplied compatible binary can also be used. |
| Other architectures | R preparation may work | Automatic downloads are rejected; supply a compatible executable. |

Automatic installation selects platform-specific official URLs, verifies HTTPS,
checks the executable format/architecture, records SHA256, and optionally tests
startup. Windows/macOS runtime installation is manual. No administrator access,
Gatekeeper changes or automatic Rosetta installation is performed. Upstream
availability and compatibility may vary: [official software notes](https://nce.ads.uga.edu/software/).
Only Linux native execution was validated for this release on this device;
Windows/macOS CI checks R code without downloading the external solver.

Use a preinstalled trusted binary on any supported host:

```r
Sys.setenv(BLUPF90_PLUS="C:/tools/blupf90/blupf90+.exe") # Windows example
# Sys.setenv(BLUPF90_PLUS="/opt/blupf90/blupf90+")        # macOS/Linux example
```

## Quick start and scalable traits

```r
pedigree <- DF[, c("ID", "Sire", "Dam")]
# Optional US-format birth dates:
# pedigree$Birth_Date <- as.Date(as.character(DF$Born), "%m/%d/%Y")
phenotypes <- DF[, c("ID", "Sex", "Age_6Month", "BW6", "BL6", "CG6")]
# Fixed covariates are not imputed. Keep the full pedigree.
phenotypes <- phenotypes[complete.cases(phenotypes[c("Sex", "Age_6Month")]), ]

one <- run_genetic_analysis(
  pedigree, phenotypes,
  traits=c("BW6", "BL6", "CG6"),
  fixed=~Age_6Month+Sex
)
one$ebv
one$heritability
one$genetic_correlations
one$inbreeding
```

`traits` accepts one, two, three or more numeric columns. All trait pairs are
returned. The default resource guard is 32 traits, adjustable through
`blupf90_control=list(max_traits=...)`; model size and sample information still
limit sensible analyses. Older calls with `backend="blupf90"` continue to work;
other engines are rejected. Use `blupf90_control`, not the retired `model_control`.

## Choose EM predictor traits

Leaving `predictor_traits=NULL` preserves the existing joint multivariate Gaussian
EM: each missing outcome uses the other observed selected outcomes in that row.
Rows missing **all selected outcome traits** are excluded first, even if auxiliary
predictors are observed. They remain in the pedigree for prediction.

For more than two outcomes, or for a single outcome, select auxiliary numerical
traits explicitly. Predictor columns must be present in the phenotype frame;
they do not become genetic-analysis outcomes unless also listed in `traits`.

```r
# Shared predictor pool, excluding the target itself:
one <- run_genetic_analysis(
  pedigree, phenotypes,
  traits=c("BW6", "BL6", "CG6"), fixed=~Age_6Month+Sex,
  predictor_traits=c("BW6", "BL6")
)

# Or define predictors separately for every outcome:
one <- run_genetic_analysis(
  pedigree, phenotypes,
  traits=c("BW6", "BL6", "CG6"), fixed=~Age_6Month+Sex,
  predictor_traits=list(BW6="BL6", BL6="BW6", CG6=c("BW6", "BL6"))
)
```

A named list must cover every outcome. `character()` selects an intercept-only
model for that target. Each selected-predictor target has its own observed-data
Gaussian EM distribution over that target and its predictors. These separate
models do **not** imply a single shared joint covariance across all outcomes.
Observed predictor values are used; filled values are not recycled as observations.
If all chosen predictors are missing in a row, its target prediction is the fitted
marginal mean. Targets and predictors need enough nonconstant observations and
pairwise overlap (at least three); singular or unidentified models fail explicitly.

Each reported cycle performs two Gaussian E/M updates. The E-step retains
conditional missing-value covariance, with no random error injected. The defaults
require at least ten cycles, correlation change <=1e-6 and relative parameter
change <=1e-8; `em_control` can change `max_iter`, `tol`, `min_iter`, `cor_tol`.
Complete data skip iteration. The legacy `method="sequential"` remains available
only without explicit predictor selection.

```r
one$imputation$diagnostics
# Default joint mode:
one$imputation$correlation_history
# Selected-predictor mode: inspect each fitted target model separately.
one$imputation$target_models$CG6$correlation_history
one$imputation$target_models$CG6$correlation_pairs
```

Imputation residual variance is phenotypic conditional variance, in squared trait
units; it is distinct from the animal-model residual variance. Imputed values are
conditional means. The downstream fit treats them as single imputations, so its
SEs and accuracy exclude imputation uncertainty and estimates may be biased.

## Maternal and other random effects

The direct animal genetic effect is always included with the pedigree A matrix.
Add terms through a named `random` list:

```r
one <- run_genetic_analysis(
  pedigree, phenotypes, traits=c("BW6", "BL6"), fixed=~Age_6Month+Sex,
  random=list(
    maternal=pedigree_effect("Dam", correlate_with="animal"),
    maternal_environment=iid_effect("Dam")
  ),
  blupf90_control=list(keep_files=TRUE, timeout=900)
)

# Other independent groups, when Herd is a phenotype column:
# random=list(herd=iid_effect("Herd"), litter=iid_effect("Litter"))

one$random_solutions$maternal
one$random_solutions$maternal_environment
one$random_covariances
one$variance_ratios
```

`pedigree_effect(column)` fits a pedigree genetic effect at the IDs in that
column, using the full A inverse. `Sire` and `Dam` are looked up from the matched
pedigree; other ID columns come from phenotypes. Maternal EBVs and prediction SEs,
PEVs, reliability and accuracy are returned for the full pedigree.
`iid_effect(column)` fits independent group levels with a separate covariance
across traits. Group solutions contain `ID`, `trait`, `EBV` (the generic solution
column), `SE` and `effect`; an environment solution is not a genetic breeding value.

Set `correlate_with="animal"` on the maternal `Dam` term to jointly estimate
direct and maternal genetics. Leaving it NULL retains an independent maternal
model, allowing explicit comparison of the two models. Other random terms remain
independent of this joint genetic group and each other.

For k traits the joint covariance is a 2k by 2k matrix, ordered
`direct:trait1,...,direct:traitk,maternal:trait1,...,maternal:traitk`:
`[G_direct, C_dm; t(C_dm), G_maternal]`. All k-squared entries of C_dm are
estimated, including direct trait i with maternal trait j; C_dm itself need not
be symmetric. Positive and negative covariances are allowed, subject to the
admissibility of the full joint matrix. One maternal term can correlate with
direct genetics. Correlated terms are automatically made consecutive in native
effect order and their ID maps retain the full pedigree.

```r
one$direct_maternal_covariance   # Rows: direct traits; columns: maternal traits
one$joint_genetic_covariance    # Full 2k by 2k genetic matrix
one$direct_maternal_parameters  # Covariance, covariance_SE, correlation, correlation_SE
```

The covariance and correlation SEs use BLUPF90+'s native `se_covar_function`
approximation. The final BLUP pass carries the estimated cross-effect covariance
as well as all variance components, so both direct and maternal EBVs, SEs and
accuracy reflect the correlated model. Genetic correlations in
`genetic_correlations` still describe the direct genetic effect across traits;
direct-maternal correlations are in `direct_maternal_parameters`.

Arbitrary correlated groups, random regressions, genomic models and repeated
phenotype rows per animal are not exposed in this release. Maternal permanent
environment is supported through offspring sharing a dam; individual repeated
records still require a future repeated-record interface.

Missing dam/group IDs in fitting rows are errors: exclude those phenotype rows
or supply scientifically valid known groups. Unknown mothers are never pooled
into a shared fictional dam. IID terms need at least two groups and replicated
levels. Maternal genetic and maternal environment can be difficult to separate;
use informative family structure and multiple generations. A richer model is not
a guarantee of identifiability or convergence.

Direct heritability uses the reference phenotypic variance
`P[t] = G_direct[t,t] + G_maternal[t,t] + C_dm[t,t] + sum(G_other[t,t]) + R[t,t]`
when the correlated maternal model is selected. Thus `h2_direct=G_direct[t,t]/P[t]`
and `h2_maternal=G_maternal[t,t]/P[t]`. `variance_ratios` reports the maternal
heritability and other component proportions with their native approximate SEs.
For independent terms C_dm is zero; without a maternal term its variance is absent.
The ratios need not sum to one when covariance contributes to P.

This follows the reference non-inbred offspring/dam relationship A_id,dam=0.5:
the covariance contribution to an offspring's variance is 2*0.5*C_dm=C_dm.
This is a conventional reference-population heritability, not an individual
record's realized variance in an inbred pedigree. The prediction model itself
uses the full pedigree relationships. Native SE expressions use the identical
denominator; estimates are not silently clipped to [0,1]. See the
[BLUPF90+ maternal variance examples](https://nce.ads.uga.edu/wiki/doku.php?id=readme.blupf90plus)
and [cross-trait/cross-effect matrix order](https://nce.ads.uga.edu/wiki/doku.php?id=faq.parameter).

## Outputs, accuracy and audit

`ebv` contains full-pedigree direct EBV, SE, PEV, Reliability, Accuracy and
Accuracy_percent for each trait. `PEV=SE^2`,
`Reliability=1-PEV/((1+F)*G_animal[t,t])`, and
`Accuracy=sqrt(Reliability)`. F is a fraction; `F_percent` is also returned.
Maternal accuracy uses its own genetic covariance. Covariance-parameter and
imputation uncertainty are not included in prediction SEs/accuracy.

The pipeline also returns `completed_data`, `pedigree`, `imputation`, covariance
matrices, parameter tables and `diagnostics`. Use `keep_model=TRUE` and
`blupf90_control=list(keep_files=TRUE)` to retain job metadata and native files.
Every job has a separate directory, original-to-numeric ID maps for each effect,
parameter map, sorted pedigree, sparse A inverse, numeric phenotypes and logs.
The runner performs REML variance estimation then final BLUP prediction with all
converged component covariances. Missing output or failed convergence is an error.

Native random groups follow the [BLUPF90 parameter specification](https://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=blupf90_all.pdf)
and [user-file inverse convention](https://nce.ads.uga.edu/wiki/doku.php?id=user_defined_files_for_covariances_of_random_effects).

## Data guide: pedigree, phenotype and fixed effects

The analysis takes **two ordinary R data frames**, not a database connection,
file path, or raw BLUPF90 parameter file:

```r
result <- run_genetic_analysis(
  pedigree=pedigree, data=phenotypes,
  traits=c("BW", "WW", "YW"), fixed=~Sex+Herd,
  backend="blupf90"
)
```

`pedigree` supplies relationships for all animals to predict. `phenotypes`
supplies measured traits and fixed predictors for the animals included in fitting.
Their row orders can differ: the package matches by `ID`, not row position.

### Pedigree data dictionary

| Column | Required? | R type | Meaning and rules |
|---|---|---|---|
| `ID` | Yes | Character | Permanent, unique animal identifier. Never missing, empty or `"0"`. |
| `Sire` | Yes | Character or NA | Biological sire's ID; added as a founder if its row is absent; NA if unknown. |
| `Dam` | Yes | Character or NA | Biological dam's ID; added as a founder if its row is absent; NA if unknown. |
| `Birth_Date` | No | Date or ISO character | Animal's own birth date, `YYYY-MM-DD`. Missing dates are allowed. |
| Other columns | No | Any appropriate type | Metadata such as breed or sex; not automatically included in the model. |

Keep IDs as text when importing. `"001"` and `"1"` are different identifiers;
converting numeric-looking IDs to numbers can destroy leading zeros. Use stable,
globally unique IDs when combining farms. Keep spelling/case consistent and
remove surrounding whitespace in both tables before matching. Pedigree sorting
trims pedigree IDs, but phenotype IDs are not automatically trimmed.

You can supply just `ID`, `Sire` and `Dam` rows for the recorded animals.
The package automatically adds one row for each named sire/dam absent from `ID`,
with missing Sire, Dam, Birth_Date and other metadata. Added parents sort above
their offspring and receive EBVs/SEs and F values. Each added parent is treated
as an unrelated, non-inbred founder; its unrecorded ancestors are not recovered.
If you have the parent's real ancestry, include those rows to retain that
information. Leave a parent unknown only when its identity is unknown.
Keeping known parent IDs preserves half-sib/full-sib relationships. Ancestor-only
animals belong in the pedigree and still receive EBVs/SEs and F values.

Empty strings and `"0"` are accepted as unknown **parent** codes and normalized to
NA. In a database, use SQL `NULL` instead. Do not create an animal with ID `"0"`.
An explicit parent birth date must precede the child's. Birth dates control sorting;
they do not automatically become fixed effects. The package checks self-parentage,
identical sire/dam and cycles, but does not cross-check parent sex metadata.

### Automatic conversion of an ID/Sire/Dam table 

For example, only the offspring need appear in your supplied table:

```r
pedigree <- data.frame(
  ID=c("C001", "C002"),
  Sire=c("S001", "S001"),
  Dam=c("D001", "D002")
)
sort_pedigree(pedigree)
#     ID Sire  Dam Birth_Date
# 1 D001 <NA> <NA>       <NA>
# 2 D002 <NA> <NA>       <NA>
# 3 S001 <NA> <NA>       <NA>
# 4 C001 S001 D001       <NA>
# 5 C002 S001 D002       <NA>
```

There is no need to run `sort_pedigree()` separately: the pipeline performs this
conversion before BLUPF90+ fitting. Missing birth dates are allowed. If offspring dates
are supplied, missing-date founders use their oldest offspring's date for priority;
with no dates, ID breaks ties while parents still precede children.

Your BW6 call can therefore remain:

```r
one <- run_genetic_analysis(pedigree, phenotypes, traits="BW6",
  fixed=~Age_6Month+Sex, backend="blupf90")
one$pedigree                          # completed, sorted pedigree
one$diagnostics$added_founder_count   # number of inserted parent rows
one$diagnostics$added_founder_ids     # original IDs of inserted parents
```

`BW6` and `Age_6Month` must be numeric columns in your phenotype table, and Sex
must be a complete categorical predictor. Automatically added parents need no
phenotype rows. Duplicate IDs, self-parentage, identical sire/dam and cycles are
still errors. Resource limits apply to the **expanded** pedigree, including added
founders. After updating the package, restart an R session that already loaded an
older version before retrying this call.

### Phenotype and covariate data dictionary

| Column | Required? | R type | Meaning and rules |
|---|---|---|---|
| `ID` | Yes | Character | One unique row per animal; must occur in the full pedigree. |
| Trait columns | At least one | Numeric | Finite measurements or NA, with exactly the names supplied in `traits`. |
| `Sex`, `Herd`, `Batch`, etc. | If used in `fixed` | Factor or character | Categorical effects; convert numeric codes to factors when they represent groups. |
| `Age`, `Dam_Age`, etc. | If used in `fixed` | Numeric | Continuous predictors, complete and finite. |
| Other columns | No | Appropriate type | Retained metadata; not fitted unless explicitly referenced. |

The labels BW, WW and YW below are examples, not reserved names:

| Example trait | Meaning | Example unit | Model role |
|---|---|---|---|
| `BW` | Birth weight | kg | Missing values use other observed traits in bidirectional EM. |
| `WW` | Weaning weight at the defined recording/adjustment age | kg | May use observed BW/YW to predict missing WW. |
| `YW` | Yearling weight at the defined recording/adjustment age | kg | May use observed BW/WW to predict missing YW. |

Define each trait consistently across animals and sources: unit, age/window,
measurement method, and whether values are raw or pre-adjusted. Convert units
before analysis; do not combine kg and lb in one column. The package does not
perform age standardization or unit conversion. EBVs and EBV SEs use the trait's
input unit; covariance entries use the product of the corresponding units.
Heritability, genetic correlation and F are dimensionless; `F_percent` is percent.

The vector `traits` defines model/output order. Default bidirectional EM uses
other observed traits in each row in both directions; observation counts set
processing priority. For a one-way biological chain, use
`em_control=list(method="sequential", trait_order="supplied")`. For example:

```r
traits <- "BW"                    # one trait
traits <- c("BW", "WW")           # two traits
traits <- c("BW", "WW", "YW")     # three traits
# More traits: append exact numeric column names; joint EM uses all observed traits.
```

Use actual R `NA` for missing observations. Zero is a real response value, not a
missing code. Convert external sentinel codes such as -999 to NA only when the
source data dictionary defines them as missing; do this in trait columns only.
Text such as `"185 kg"` is not a numeric measurement. Each selected trait needs
observed variation and enough observed records for its EM regression. A trait
with all missing values cannot be recovered by this implementation.

Fixed predictors must be complete: EM imputes response traits, not Sex/Herd/Age.
Numeric Herd codes would be treated as a continuous regression unless converted
to factors. One-level categorical predictors and confounded/saturated designs can
be unestimable; omit an unvarying effect rather than inventing extra categories.
All selected traits share the supplied per-animal covariate values, with separate
coefficients per trait. If an animal changes herd between recording stages, a
single `Herd` column must have a clearly defined interpretation; this version does
not accept different per-trait fixed-effect formulas or repeated-record models.

### Example CSV inputs and loading in R

The following small dataset demonstrates formatting and relationships. It is
**not** evidence that a stable multi-trait REML analysis is possible with only six
phenotyped animals. Use the larger simulated example for a runnable model fit.

Save this as `pedigree.csv`:

```csv
ID,Sire,Dam,Birth_Date
S001,,,
S002,,,
D001,,,
D002,,,
A001,S001,D001,2022-01-10
A002,S001,D002,2022-01-12
A003,S002,D001,2022-01-15
A004,S002,D002,2022-01-18
A005,S001,D001,2022-01-20
A006,S002,D002,2022-01-22
```

Save this as `phenotypes.csv`:

```csv
ID,Sex,Herd,BW,WW,YW
A001,M,H1,32,185,310
A002,F,H1,30,,295
A003,M,H2,34,210,350
A004,F,H2,29,175,
A005,F,H1,31,190,320
A006,M,H2,33,205,340
```

```r
pedigree <- read.csv("pedigree.csv", colClasses="character",
  na.strings=c("", "NA"), check.names=FALSE)
pedigree$Birth_Date <- as.Date(pedigree$Birth_Date, format="%Y-%m-%d")

phenotypes <- read.csv("phenotypes.csv",
  colClasses=c("character", "character", "character", "numeric", "numeric", "numeric"),
  na.strings=c("", "NA"), check.names=FALSE)
phenotypes$ID <- trimws(phenotypes$ID)
phenotypes$Sex <- factor(phenotypes$Sex)
phenotypes$Herd <- factor(phenotypes$Herd)
traits <- c("BW", "WW", "YW")

# Quick checks and preparation; no external solver is needed for these steps.
stopifnot(!anyDuplicated(pedigree$ID), !anyDuplicated(phenotypes$ID))
stopifnot(all(phenotypes$ID %in% pedigree$ID))
sorted_pedigree <- sort_pedigree(pedigree)
completed <- sequential_em(phenotypes, traits)
job <- prepare_blupf90(sorted_pedigree, completed$data, traits, fixed=~Sex+Herd)
```

Validate original date text before conversion: invalid strings can become NA in
`as.Date()`. You can also pass the original ISO character column directly to the
package, which rejects malformed nonmissing dates. CSV column names are
case-sensitive. Spreadsheet files can be imported with an appropriate reader,
then converted to these same data frames; the core functions do not read Excel
files or database files themselves.

### Repeated measurements and cohort selection

A source database may hold many measurements per animal. Before fitting, define
one reproducible record per **animal and selected trait**: for example a specified
age window with a documented rule to choose among measurements. Do not use
`distinct(ID)` or an arbitrary first row to discard observations. Do not silently
average biologically different stages. Reshape the chosen records into one wide
row per ID. Individual repeated-record modeling requires an extension to the current
interface. Maternal permanent environment through shared dams is supported.

The pipeline now automatically excludes phenotype rows with no observed value
in any selected trait, while keeping those animals in the pedigree. The equivalent
manual filter is shown below; it is no longer required before calling the pipeline:

```r
has_observation <- rowSums(!is.na(as.matrix(phenotypes[traits]))) > 0L
phenotypes <- phenotypes[has_observation, , drop=FALSE]
```

Partially missing response rows can remain for bidirectional EM. With `impute=FALSE`,
fully missing rows are still excluded, and selected responses in retained rows
must already be complete; that flag is not an
observed-data/missing-response REML mode. Database selection should retain the
known ancestors of the selected animals. Do not filter the pedigree to just IDs
with phenotype rows, or reject founder rows merely because their birth dates are
unknown. Sorting can add referenced missing parents as founders, but cannot
restore their own ancestry if those parent records were removed by an SQL query.

## Database integration and schema

### Which data sources can be used?

The program is database-independent **after extraction into R data frames**.
Database connectors are optional application-layer dependencies, not built into
the analysis package. Changing the database does not change the animal model.

| Data source | R connection/import route | What to supply to the program |
|---|---|---|
| SQLite `.sqlite` / `.db` | `DBI` + [RSQLite](https://rsqlite.r-dbi.org/reference/SQLite.html) | Query results converted to the two input data frames. |
| PostgreSQL | `DBI` + [RPostgres](https://rpostgres.r-dbi.org/) | The same input data frames. |
| MySQL or MariaDB | `DBI` + [RMariaDB](https://rmariadb.r-dbi.org/reference/dbConnect-MariaDBDriver-method.html) | The same input data frames. |
| ODBC-accessible database | `DBI` + [odbc](https://odbc.r-dbi.org/) and the database's driver | The same input data frames; SQL syntax depends on the server. |
| CSV / TSV exports | `read.csv()` / `read.delim()` | Typed pedigree and wide phenotype data. |
| Other stores, APIs or JSON | Your existing connector and flattening/pivot step | Ordinary data frames with the documented columns. |

A database file is not the same as a BLUPF90 job file. The package generates its
own native numeric data, sparse inverse and ID mappings after input validation.

### Flexible relational schema for adding traits

A wide phenotype table with BW/WW/YW columns is sufficient for a fixed trait set.
For a changing trait set, the following **example SQLite schema** stores each
trait as a row. Adding a trait then adds a catalog/measurement row instead of
requiring another SQL table column. Adapt types/constraints to your database;
the SQL below is a proposed application schema, not a migration executed by the
package.

```sql
PRAGMA foreign_keys = ON;

CREATE TABLE animals (
  animal_id TEXT PRIMARY KEY NOT NULL,
  sire_id TEXT REFERENCES animals(animal_id),
  dam_id TEXT REFERENCES animals(animal_id),
  birth_date TEXT,                         -- ISO YYYY-MM-DD or NULL
  sex TEXT CHECK (sex IN ('F', 'M')),       -- optional metadata
  CHECK (length(trim(animal_id)) > 0 AND animal_id <> '0'),
  CHECK (sire_id IS NULL OR sire_id <> animal_id),
  CHECK (dam_id IS NULL OR dam_id <> animal_id),
  CHECK (sire_id IS NULL OR dam_id IS NULL OR sire_id <> dam_id)
);

CREATE TABLE phenotype_context (
  analysis_id TEXT NOT NULL,
  animal_id TEXT NOT NULL REFERENCES animals(animal_id),
  herd TEXT NOT NULL,
  PRIMARY KEY (analysis_id, animal_id)
);

CREATE TABLE trait_catalog (
  trait_code TEXT PRIMARY KEY NOT NULL,
  trait_name TEXT NOT NULL,
  unit TEXT NOT NULL,
  sequence_order INTEGER NOT NULL UNIQUE CHECK (sequence_order > 0)
);

CREATE TABLE trait_values (
  analysis_id TEXT NOT NULL,
  animal_id TEXT NOT NULL,
  trait_code TEXT NOT NULL REFERENCES trait_catalog(trait_code),
  measured_on TEXT,                        -- ISO date or NULL
  trait_value REAL,                        -- NULL means missing
  PRIMARY KEY (analysis_id, animal_id, trait_code),
  FOREIGN KEY (analysis_id, animal_id)
    REFERENCES phenotype_context(analysis_id, animal_id)
);

CREATE INDEX animals_sire_idx ON animals(sire_id);
CREATE INDEX animals_dam_idx ON animals(dam_id);
CREATE INDEX trait_values_analysis_trait_idx
  ON trait_values(analysis_id, trait_code);
```

`animals` is the pedigree/master list. `phenotype_context` defines membership and
a shared fixed-effect context for an analysis cohort. `trait_catalog` records
meaning, units and intended ordering. `trait_values` holds one selected response
per animal/trait/cohort; its primary key prevents duplicate selected records.
Raw repeated measurements can be retained in a separate source table before this
selection step. Unknown parent IDs and missing responses are SQL `NULL`, not the
literal text `"NA"` and not a response of zero.

Insert founder/parent rows before children when using these immediate foreign-key
constraints. Automatic founder completion occurs in R and does not insert rows
into your source database. This example SQL schema therefore still needs parent
rows before child inserts; an offspring-only R table can instead come from a CSV
or an existing staging table without those parent foreign keys. If storing EBVs
with a foreign key to the source animal table, account for any added founder IDs
before inserting those results. For existing databases, either query compatible views or rename the
returned columns; you do not need to restructure the original database. These SQL
constraints do not detect multi-generation pedigree cycles or validate biological
measurement ranges; package validation and source-data checks still apply.
SQLite foreign keys must be enabled on each connection; see the
[SQLite foreign-key documentation](https://www.sqlite.org/foreignkeys.html).

### Read an existing SQLite database and pivot traits into program inputs

This example assumes the preceding schema has already been created and populated
in `genetics.sqlite`. It only reads data; it does not create a schema or insert
records. Optional dependencies: `install.packages(c("DBI", "RSQLite"))`.

```r
read_analysis_inputs <- function(db_path, analysis_id, traits) {
  stopifnot(is.character(traits), length(traits) > 0L,
    !anyNA(traits), !anyDuplicated(traits))
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname=db_path,
    flags=RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(con), add=TRUE)
  DBI::dbBegin(con)  # Keep these SQLite reads in one consistent snapshot.
  on.exit(DBI::dbRollback(con), add=TRUE, after=FALSE)

  # Retain the full master pedigree, including animals without measurements.
  pedigree <- DBI::dbGetQuery(con, '
    SELECT animal_id AS "ID", sire_id AS "Sire", dam_id AS "Dam",
           birth_date AS "Birth_Date"
    FROM animals ORDER BY animal_id')
  pedigree$ID <- as.character(pedigree$ID)
  # Keep ISO text here so the package can validate malformed dates.

  context <- DBI::dbGetQuery(con, '
    SELECT c.animal_id AS "ID", a.sex AS "Sex", c.herd AS "Herd"
    FROM phenotype_context c JOIN animals a ON a.animal_id = c.animal_id
    WHERE c.analysis_id = ? ORDER BY c.animal_id', params=list(analysis_id))

  catalog <- DBI::dbGetQuery(con,
    'SELECT trait_code, unit, sequence_order FROM trait_catalog')
  stopifnot(all(traits %in% catalog$trait_code))
  # Caller order governs model output; EM uses observed counts by default.

  placeholders <- paste(rep("?", length(traits)), collapse=",")
  sql <- paste0('
    SELECT animal_id AS "ID", trait_code, trait_value FROM trait_values
    WHERE analysis_id = ? AND trait_code IN (', placeholders, ')')
  values <- DBI::dbGetQuery(con, sql,
    params=c(list(analysis_id), as.list(traits)))
  stopifnot(!anyDuplicated(context$ID),
    !anyDuplicated(values[c("ID", "trait_code")]),
    all(values$ID %in% context$ID), is.numeric(values$trait_value))

  phenotypes <- context
  for (trait in traits) {
    chosen <- values[values$trait_code == trait, , drop=FALSE]
    phenotypes[[trait]] <- chosen$trait_value[match(phenotypes$ID, chosen$ID)]
  }
  observed <- rowSums(!is.na(as.matrix(phenotypes[traits]))) > 0L
  phenotypes <- phenotypes[observed, , drop=FALSE]
  phenotypes$ID <- as.character(phenotypes$ID)
  phenotypes$Sex <- factor(phenotypes$Sex)
  phenotypes$Herd <- factor(phenotypes$Herd)
  list(pedigree=pedigree, phenotypes=phenotypes, traits=traits,
    trait_dictionary=catalog[catalog$trait_code %in% traits, , drop=FALSE])
}

inputs <- read_analysis_inputs("genetics.sqlite", analysis_id="cohort_2022",
  traits=c("BW", "WW", "YW"))

# Run on a scientifically suitable, sufficiently informative cohort.
result <- run_genetic_analysis(inputs$pedigree, inputs$phenotypes,
  traits=inputs$traits, fixed=~Sex+Herd, backend="blupf90")
```

A missing `trait_values` row becomes NA after matching, just like a stored NULL.
If none of the animals has observed records for a requested trait, select another
trait set or collect more data; EM cannot invent that trait's distribution.
The example deliberately fetches the full pedigree for clarity. For a large
master database, extract the target animals **and their complete known ancestor
closure** into a consistent snapshot, then respect the backend's resource limits.

DBI passes values through `params`; do not concatenate analysis IDs or animal IDs
into SQL. Placeholder syntax varies by driver: the example uses SQLite `?`,
whereas RPostgres uses positional `$1`, `$2`, etc. Change the SQL placeholders when
changing drivers; `dbGetQuery()` does not translate SQL dialects. See the
[DBI parameter-binding reference](https://dbi.r-dbi.org/reference/dbBind.html).
Keep SQL identifiers and fixed formulas server-defined. Connection credentials
belong in your application configuration, outside the README and source code.

### Connecting other databases

For PostgreSQL, an application can establish its connection using configuration
provided outside the source code:

```r
con <- DBI::dbConnect(RPostgres::Postgres(),
  host=Sys.getenv("PGHOST"), dbname=Sys.getenv("PGDATABASE"),
  user=Sys.getenv("PGUSER"), password=Sys.getenv("PGPASSWORD"))
# Use PostgreSQL SQL/placeholders and extract the documented data frames.
# Disconnect when the application has finished using this connection:
DBI::dbDisconnect(con)
```

The SQLite extraction function above is not a universal database adapter. For
PostgreSQL replace `?` with `$1`, `$2`, etc.; preserve case using quoted aliases
such as `AS "ID"`. For other servers adapt their quoting, date types and transaction
isolation as appropriate. Preserve text IDs and SQL NULL values during extraction.
A shared application can manage connections/pools separately from the analysis
worker; the statistical functions only receive the extracted data frames.

### Which additional livestock data are used?

Breed, birth year, management group and dam age can be included as complete fixed
predictors when scientifically appropriate and estimable. They are not included
automatically. Measurement dates, laboratory methods and recording ages can be
stored for selection/audit even if they are not fitted. SNP/genotype tables are
not used by this pedigree-only module: passing genotype columns does not activate
GBLUP or single-step genomic analysis. Maternal genetic effects, permanent
environment effects and repeated-record models likewise need an explicit model
extension rather than extra database columns alone.

### Before fitting and when storing results

Check that both extracts describe the same snapshot/cohort; every phenotype ID
matches the expanded pedigree; added founders reflect the intended unknown
ancestry assumption; trait units and selection
rules are consistent; one phenotype row remains per ID; and fixed predictors are
complete and estimable. At least three phenotype rows are required by the package,
but that minimum is **not** a sufficient sample size for reliable genetic parameter
estimation. Connectivity, families, missingness and the number of fitted parameters
also matter. Do not combine species or incompatible breeding populations merely
because their databases share column names.

If saving results back to a database, use a unique analysis run ID rather than
overwriting the raw phenotypes. Useful result tables are:

| Result table | Suggested key | Values from the returned object |
|---|---|---|
| Analysis run | `run_id` | Backend, ordered trait list, fixed formula, package/native versions, cohort selection, unit dictionary, convergence, imputation flag and SE method. |
| EBVs | `run_id, ID, trait` | `result$ebv`: EBV, SE, PEV, Reliability, Accuracy and Accuracy_percent. |
| Heritabilities | `run_id, trait` | `result$heritability`: h2 and SE. |
| Genetic correlations | `run_id, trait1, trait2` | `result$genetic_correlations`: rg and SE. |
| Inbreeding | `run_id, ID` | `result$inbreeding`: F and F_percent. |

Store the ordered trait list and units with each run, and retain source-measurement
identifiers/selection rules for reproducibility. Preserve the original observed
traits separately from `result$completed_data` so imputed values are not mistaken
for measurements. The core package does not automatically read or write any
application database. Use your own DBI transactions and access controls for result
storage.


## Public API

- `run_genetic_analysis()`: complete pipeline.
- `sort_pedigree()`, `pedigree_precision()`, `pedigree_inbreeding()`, `pedigree_relationship()`.
- `sequential_em()`: default joint EM or selected-predictor target models.
- `pedigree_effect()`, `iid_effect()`: extra random-effect definitions.
- `prepare_blupf90()`, `write_blupf90_job()`, `fit_blupf90_model()`, `read_blupf90_results()`.
- `install_blupf90()`, `add_ebv_accuracy()`, `reshape_traits()`.

For Shiny/Plumber use trusted server-controlled executable paths and formulas,
isolated job directories, resource limits and background jobs. Never let a web
client select arbitrary executables. Data frames are passed explicitly; no global
working-directory changes or implicit network installation occurs during fitting.

## Development and release

```r
install.packages(c("testthat", "roxygen2", "covr"))
roxygen2::roxygenise()
testthat::test_local()
# Set BLUPF90_PLUS to include real solver integration tests.
covr::package_coverage()
```

Run `R CMD build .` and `R CMD check --no-manual geneticAnalysisR_1.2.0.tar.gz`.
CI checks R code on Linux, Windows and macOS; a separate Linux job installs the
external engine and exercises native tests. See VALIDATION.md for local evidence.

Upload the repository ZIP contents to your own GitHub repository. Fill maintainer
metadata in DESCRIPTION and LICENSE first. Do not upload private animal datasets,
installed R libraries, native executables/runtimes or retained analysis jobs.
