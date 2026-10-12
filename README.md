# geneticAnalysisR 1.1.2

A modular R animal-model package with selectable **BLUPF90+** and **sommer**
backends. BLUPF90+ is the default and does not require sommer to be installed.
Supply any ordered vector of numeric trait names: one, two, three, or more.
The same result schema provides full-pedigree EBVs with prediction-error SEs,
accuracy and reliability; heritability and SE for each trait; genetic correlation
and SE for every pair; and exact inbreeding coefficients as fractions and percentages.

## Contents

- [Release features](#what-is-included-in-release-112)
- [BW6 quick example](#quick-example-six-month-weight)
- [R package installation](#install-the-r-package)
- [Linux BLUPF90+ installation](#install-blupf90-explicitly)
- [Backend and trait selection](#select-the-backend-and-trait-count)
- [Pedigree and phenotype data guide](#data-guide-pedigree-phenotype-and-fixed-effects)
- [Database integration](#database-integration-and-schema)
- [Missing-data selection and EM ordering](#missing-measurement-selection-and-adaptive-em-112)
- [Backend comparison](#how-the-backends-differ)
- [EBV accuracy and reliability](#ebv-accuracy-and-reliability-included-in-112)
- [API and audit files](#public-api-and-audit-files)
- [Web-backend use and scientific limits](#web-backend-use-and-scientific-limits)
- [Development and CI](#development-and-ci)
- [Upload to GitHub](#upload-this-release-to-github)

## What is included in release 1.1.2?

- Selectable BLUPF90+ (default) and optional sommer backends.
- One, two, three or more traits through an ordered `traits` vector.
- Bidirectional Gaussian EM for missing responses, with an optional legacy chain.
- Automatic pedigree completion: a three-column ID/Sire/Dam table is sufficient;
  absent named parents become unknown-ancestry founders above their offspring.
- Exact sparse pedigree A-inverse and inbreeding for the native backend.
- EBV, SE, PEV, Reliability, Accuracy and Accuracy_percent for each animal/trait.
- Trait heritabilities and pairwise genetic correlations with SEs.
- Linux installation helper, isolated native jobs, diagnostics, unit/integration
  tests, generated R help pages and GitHub Actions validation.
- CSV/data dictionaries, relational database schema and DBI extraction examples.

This 1.1.2 release consolidates all the improvements developed in this project.
The R package name is `geneticAnalysisR`; the suggested GitHub repository name is
`genetic-analysis-r`. Those names can differ: `library()` uses the R package name.

## Quick example: six-month weight

```r
library(geneticAnalysisR)
# Install BLUPF90+ once if it is not already available:
# install_blupf90()

# These are columns in your existing data frame DF.
pedigree <- DF[c("ID", "Sire", "Dam")]
phenotypes <- DF[c("ID", "Age_6Month", "Sex", "BW6")]
phenotypes$Sex <- factor(phenotypes$Sex)

one <- run_genetic_analysis(pedigree, phenotypes, traits="BW6",
  fixed=~Age_6Month+Sex, backend="blupf90")
head(one$ebv[c("ID", "trait", "EBV", "SE", "Accuracy", "Accuracy_percent", "Reliability")])
one$heritability
one$inbreeding
one$diagnostics$added_founder_ids
```

DF is your own dataset, not a bundled object. ID/Sire/Dam identify animals;
BW6 and Age_6Month must be numeric, and Sex must be complete and estimable.
If your date column is called `Born`, you may include it as `Birth_Date` explicitly:
`pedigree$Birth_Date <- DF$Born`, provided its values are Date or ISO YYYY-MM-DD.
Otherwise dates can be omitted. See the detailed data guide below for valid types,
missing values, units and database inputs. For a self-contained simulated fit,
use the bundled `example_script.R` after installing the package and backend.

## Install the R package

```r
install.packages(c("Matrix", "processx", "curl"))
install.packages("geneticAnalysisR_1.1.2.tar.gz", repos=NULL, type="source")
# Optional, only for the sommer backend:
install.packages("sommer")
```

After uploading this repository to GitHub, use
`remotes::install_github("YOUR_ACCOUNT/genetic-analysis-r")`.
R >= 4.1 is required. Linux x86-64 is the target for automatic BLUPF90+ installation.
Replace the placeholder maintainer name/email in DESCRIPTION before publishing.

## Install BLUPF90+ explicitly

```r
library(geneticAnalysisR)
installed <- install_blupf90()
# The default installation is found automatically on later runs.
# Custom locations also work:
installed <- install_blupf90(directory="/your/writable/tools/blupf90/bin")
Sys.setenv(BLUPF90_PLUS=installed$executable)
```

The helper downloads the official Linux ELF executable over verified HTTPS,
records its SHA256, and checks `--help`. When Intel MKL/OpenMP libraries are
missing, Python3/pip installs pinned 2025.3.0 runtimes into a private subdirectory.
A launcher sets library paths for its child process. No administrator privileges,
system Python changes, global library changes or shell-profile edits are needed.
Installation can use roughly 1 GB. Use `install_runtimes=FALSE` if your machine
already provides suitable runtimes. Existing executables require explicit
`overwrite=TRUE`. Downloads never occur during genetic analysis.

The upstream server currently omits intermediate TLS certificates. The helper
can supplement them with bundled public certificates, first verified against the
machine's existing CA roots. It never disables TLS verification. An administrator
can alternatively provide `ca_bundle` with a valid trust chain.

BLUPF90+ is external software: its own usage/licensing terms apply. Its executable
and Intel libraries are not distributed inside this repository or R archive.
See [official distribution](https://nce.ads.uga.edu/wiki/doku.php?id=distribution)
and [Intel MKL installation](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html?operatingsystem=linux).

## Select the backend and trait count

```r
# Single trait
one <- run_genetic_analysis(pedigree, phenotypes, traits="Birth",
  fixed=~Sex+Herd, backend="blupf90")

# Two traits
two <- run_genetic_analysis(pedigree, phenotypes,
  traits=c("Birth","Weaning"), fixed=~Sex+Herd, backend="blupf90")

# Three or any larger trait set; EM defaults to most-observed-first ordering
many <- run_genetic_analysis(pedigree, phenotypes,
  traits=c("Birth","Weaning","Yearling","Mature","Carcass"),
  fixed=~Sex+Herd, backend="blupf90",
  blupf90_control=list(threads=4L, n_iters=1000L, timeout=600))

# Switch the same inputs to sommer
alternative <- run_genetic_analysis(pedigree, phenotypes,
  traits=c("Birth","Weaning","Yearling"), fixed=~Sex+Herd,
  backend="sommer", model_control=list(n_iters=200L))

many$ebv
many$heritability
many$genetic_correlations
many$inbreeding
many$diagnostics
```

These blocks use your own `pedigree` and `phenotypes`. For a runnable simulated
example, source `system.file("example_script.R", package="geneticAnalysisR")`.
A `k`-trait analysis returns `k` heritabilities, `k*(k-1)/2` correlations, and
`n_pedigree*k` EBV rows. A single trait returns an empty correlation table.

## Input requirements and analysis workflow

Pedigree: required `ID`, `Sire`, `Dam`; optional `Birth_Date` as R Date or ISO text.
Unknown parents use NA, empty strings or "0". Named parents missing from ID
are automatically added as founders with unknown ancestry.
Duplicate IDs, self-parentage, cycles and impossible explicit dates are errors.
Founders are unrelated and non-inbred. Founders sort first; undated founders use
the date of their oldest dated direct offspring. Parent precedence overrides
chronology; birth date and then ID prioritize eligible offspring. The new
priority-queue sorter costs O(n log n), avoiding repeated full-pedigree scans.

Phenotypes: one row per ID, numeric finite traits or NA, complete fixed predictors.
Fixed formulas are one-sided; model.matrix handles factors, transformations and
interactions, with separate coefficients for each trait. Repeated records and
permanent environmental effects are outside this module's current model.

Bidirectional joint Gaussian EM is now the default. Completely unmeasured rows
are excluded first. Each iteration performs a forward and reverse EM update,
using conditional means and missing-value covariances. For two traits, missing
X uses observed Y and missing Y uses observed X. Observed data are never replaced.
EM runs at least 10 cycles before checking correlation change; parameter stability
is also required. Set `em_control=list(method="sequential")` for the older chaining
method; most-observed-first/supplied ordering remains available for that method.
Model trait/output order remains the caller's supplied order.

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

### Automatic conversion of an ID/Sire/Dam table (included in 1.1.2)

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
conversion for both backends. Missing birth dates are allowed. If offspring dates
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
row per ID. Repeated-measurement/permanent-environment modeling requires an
extension to the current model.

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

## Missing-measurement selection and adaptive EM (1.1.2)

| X | Y | Default action |
|---|---|---|
| Observed | Observed | Keep both observations unchanged. |
| Observed | NA | Predict missing Y from observed X. |
| NA | Observed | Predict missing X from observed Y. |
| NA | NA | Exclude phenotype row before EM/fitting; keep its pedigree animal. |

For more traits, a row is excluded only when all selected traits are missing.
For one trait, missing rows are excluded and there is no response left to impute.
Your original DF/phenotypes table is not overwritten. The returned completed_data
contains retained rows; excluded IDs/indices/rows are available for audit.

### Forward/reverse cycles and stopping criteria

Default `method="bidirectional"` estimates a joint multivariate normal phenotype
distribution by EM. Each cycle contains two E/M updates with forward and reversed
missingness-pattern processing. Missing X is conditional on observed Y, and missing
Y on observed X. With three or more traits, missing values are conditional on the
other **observed** traits in their row. Provisional fills are never treated as
measured predictors. Conditional missing-value covariance is included in the
second moments, avoiding hard-imputation variance collapse during EM estimation.
This is phenotype-distribution EM, not joint pedigree animal-model REML.

Every incomplete-data run performs at least `min_iter=10` full cycles before
checking convergence. The cycle criterion is the maximum absolute change across
pairwise correlations, considering both estimated-distribution and completed-data
correlations. Default `cor_tol=1e-6`; a relative parameter-change guard (`tol=1e-8`)
must also pass, so stable correlation alone cannot hide drifting means/variances.
`max_iter=500` is the default limit; nonconvergence raises an error. Data already
complete after selection skip EM and report zero iterations.

Initial nonmissing counts determine processing priority; ties keep supplied order.
This does not impose one-way dependence in joint EM. Model/output trait order stays
as supplied. Scaling/centering is internal to EM for numerical stability; final
fills and distribution parameters are restored to the original measurement units.
Initialization may shrink pairwise correlations toward identity to obtain a
positive-definite starting covariance. Estimated updates are not silently repaired;
singular/redundant trait models raise errors.

```r
one <- run_genetic_analysis(
  pedigree, phenotypes, traits=c("BW6", "BL1Y"),
  fixed=~Age_6Month+Sex, backend="blupf90",
  em_control=list(method="bidirectional", min_iter=10,
    cor_tol=1e-6, tol=1e-8, max_iter=500)
)
one$imputation$diagnostics
one$imputation$correlation_history
one$imputation$correlation_pairs
one$imputation$estimated_correlation
one$imputation$completed_correlation
one$imputation$distribution_mean
one$imputation$distribution_covariance
one$imputation$excluded_ids
```

Diagnostics iterations are **global forward/reverse cycles**, including for a trait
whose values were already complete. `corr_diff` is the cycle correlation change;
`delta` is relative parameter change. For incomplete-data runs, residual_variance
is the fitted full conditional variance of that trait given all the other traits.
The actual row-specific conditional variance may be larger when fewer other
traits are observed. Correlation history includes iteration, corr_diff,
parameter_delta, observed-data log_likelihood and estimated phenotype correlations.
Columns r_T1_T2, etc. are mapped to names by correlation_pairs. These are phenotype
correlations used for EM convergence, not the animal model's genetic correlations.
The estimated distribution correlation and correlation of conditional-mean-filled
data are different quantities and may differ.

The pairwise covariance must be informed by overlapping measured records. The
implementation requires at least three observed values per trait and three
overlapping observed records per pair when imputation is needed; this is a
computational guard, not a guarantee of adequate scientific sample size. Fully
unobserved selected traits, fewer than three retained rows, and unidentifiable
covariances raise errors. Missing fixed effects (Age/Sex) still require correction
or explicit exclusion. Single-imputation SE/accuracy limitations still apply.

### Optional earlier sequential method

```r
em_control=list(method="sequential", trait_order="most_observed",
  predictors="previous")
# Or use trait_order="supplied", predictors="all" for a specified chain.
```

Only this optional method uses an intercept-only first-trait imputation. Its
conditional regression estimates can already be an EM fixed point, so one
iteration can be legitimate. `min_iter` and `cor_tol` apply to the bidirectional
method, not this legacy method. The default now fills the first trait from the
other observed trait rather than a mean-only rule.

## How the backends differ

| Property | BLUPF90+ | sommer |
|---|---|---|
| Required installation | External Linux binary/runtimes | Optional R package |
| Pedigree representation | Sparse exact Henderson A-inverse | Dense A matrix |
| Model data layout | Native wide multi-trait files | Trait-major long data |
| Estimation | AI-REML after configurable EM warmup | mmer REML |
| EBV SE | Native square root of PEV | Exact marginal-equation square root of PEV |
| h2/rg SE | Native Monte-Carlo `se_covar_function` approximation | First-order delta method |
| Default pedigree guard | 100000 | 2000 |

Both fit `Cov(a)=G ⊗ A` and within-animal `Cov(e)=E ⊗ I`, with unstructured trait
covariance matrices. BLUPF90+ model files are generated from the trait vector;
there are no fixed two/four-trait blocks. Its sparse inverse is written using
`RANDOM_TYPE user_file`, preserving exact inbreeding corrections without rounded
pedigree codes. Pedigree IDs are mapped to consecutive numeric codes and restored
in the returned EBV table.

A native VCE pass estimates G and E and computes parameter SEs. A final BLUP pass
uses those converged covariance matrices for EBVs/SEs, avoiding the one-iteration
lag in VCE solutions. Results require successful execution, verified REML
convergence, admissible covariance matrices, complete trait/ID solution coverage
and compatible output headers. Unknown output formats fail explicitly.
The parameter SE method is recorded in diagnostics; SEs from the two methods
need not be identical. See the
[author's multi-trait tutorial](https://masuday.github.io/blupf90_tutorial/quicktour_mt.html)
and [REML/SE tutorial](https://masuday.github.io/blupf90_tutorial/vc_aireml.html).

BLUPF90+'s default `max_traits=32` is a configurable resource guard, not a model
assumption; raise it via `blupf90_control`. With k traits there are `k*(k+1)`
variance/covariance parameters; information-matrix memory grows rapidly.
Pedigree inbreeding uses an iterative bounded ancestral-pair cache, not dense A.
Worst-case inbreeding cost remains quadratic for highly connected pedigrees;
`max_pairs` bounds cache size. Larger counts are not automatically faster.
Benchmark your actual connected pedigree and traits using `benchmark_backends.R`.

The tested sommer version is 4.4.87. Its mmer interface is deprecated upstream;
output compatibility is checked, and no future-version support is assumed.

## EBV accuracy and reliability (included in 1.1.2)

Every fitted backend now includes these columns in `result$ebv`:

| Column | Meaning | Scale/unit |
|---|---|---|
| `EBV` | Predicted additive genetic merit. | Trait's input unit. |
| `SE` | Prediction-error standard deviation. | Trait's input unit. |
| `PEV` | Prediction-error variance, `SE^2`. | Squared trait unit. |
| `Reliability` | Model-based squared correlation between predicted and true breeding value. | 0 to 1. |
| `Accuracy` | Square root of Reliability. | 0 to 1. |
| `Accuracy_percent` | `100 * Accuracy`, for display. | 0 to 100 percent. |

For animal i and trait t, the fitted prior additive genetic variance is
`(1 + F_i) * G_tt`. The common calculation for both backends is:

```text
PEV_it         = SE_it^2
Reliability_it = 1 - PEV_it / ((1 + F_i) * G_tt)
Accuracy_it    = sqrt(Reliability_it)
```

Here F is the inbreeding **fraction**, not F_percent, and G_tt is the estimated
additive genetic variance for that trait, not its heritability or phenotypic
variance. The package's SE already has the prediction-error SD scale; no extra
residual-variance multiplication is needed. This follows the reliability/accuracy
definition in the [official BLUPF90+ documentation](https://nce.ads.uga.edu/wiki/doku.php?id=readme.blupf90plus).

For example, genetic variance 25, F=0 and SE=3 give PEV=9, reliability=0.64 and
accuracy=0.80 (80%). With F=0.25, the same SE and genetic variance use prior variance
31.25, giving reliability=0.712 and accuracy approximately 0.844 (84.4%). Accuracy
is a prediction-quality metric, not the EBV magnitude or a confidence interval;
an 80% accuracy does not mean an 80% probability that an individual EBV is correct.

```r
# New analyses: accuracy is calculated automatically.
one <- run_genetic_analysis(pedigree, phenotypes, traits="BW6",
  fixed=~Age_6Month+Sex, backend="blupf90")
head(one$ebv[c("ID", "trait", "EBV", "SE", "Accuracy", "Accuracy_percent", "Reliability")])

# Existing result from an earlier version: update it without another model fit.
one$ebv <- add_ebv_accuracy(one$ebv, one$genetic_covariance, one$inbreeding)
```

Genetic variance of zero gives undefined accuracy/reliability (NA). Tiny negative
reliability caused by rounding is set to zero using a default dimensionless
tolerance of 1e-4. A larger negative value means PEV exceeds the supplied prior
variance: the helper warns and returns NA for those rows rather than presenting
an unsupported accuracy. Trait and animal mappings are validated explicitly.
Model-based accuracy conditions on estimated covariance components and completed
phenotypes. It does not propagate variance-component estimation uncertainty or
single-imputation uncertainty; imputation can make apparent accuracy optimistic.

## Public API and audit files

- `run_genetic_analysis(...)`: backend-neutral high-level pipeline.
- `add_ebv_accuracy(ebv, genetic_covariance, inbreeding)`: add reliability/accuracy
  to an existing result without refitting.
- `fit_blupf90_model(...)`: native sparse backend, no sommer calls.
- `prepare_blupf90(...)`: generate data, sparse inverse, ID map and parameter text.
- `write_blupf90_job(job, directory)`: write a prepared job to a fresh directory.
- `read_blupf90_results(job, log_file, solutions_file)`: validate/import native results.
- `install_blupf90(...)`: explicit Linux installation helper.
- `pedigree_precision(...)`: sparse A-inverse and exact inbreeding.
- `pedigree_inbreeding(...)`: exact F without dense A.
- `sort_pedigree`, `pedigree_relationship`, `sequential_em`,
  `reshape_traits`, and `fit_animal_model` remain exported.

`blupf90_control` accepts executable, job_root, keep_files, timeout, threads,
n_iters, tol, em_rounds, max_pairs and max_traits. `model_control` remains exclusive
to sommer. For retained job files use `keep_files=TRUE` and an existing `job_root`;
the actual job is a unique subdirectory. Audit files include `phenotypes.txt`,
`ainverse.txt`, `pedigree.txt`, `analysis.par`, `prediction.par`, `id_map.csv`,
`parameter_map.csv`, `inbreeding.csv`, native solutions and process logs.
The pedigree file is for audit; the native model uses `ainverse.txt`.
Set `keep_model=TRUE` to retain backend metadata in the returned pipeline object.

## Web-backend use and scientific limits

Analysis jobs do not change R's working directory. External processes receive
isolated job directories, a timeout, and explicit thread counts. Deleted job files
are the default. Keep executable paths, formula objects and filesystem paths
server-controlled: do not accept arbitrary ones from an HTTP client. Do not install
software on incoming requests. Limit input size, workers and retained-file access.
No hosting service, authentication layer or endpoint is deployed by this package.

Single imputation can bias variances, h2 and rg, and all model SEs exclude
imputation uncertainty. The pipeline warns and marks imputed fits. Gaussian EM
assumes ignorable missingness and multivariate-normal conditional relationships;
it is not an animal-model EM or multiple-imputation procedure. With complete data use
`impute=FALSE`. Unidentified/nonconverged/indefinite covariance estimates are
rejected. Validate statistical assumptions on your own livestock data.

## Development and CI

```r
testthat::test_local()
covr::package_coverage()
# To include real external-program tests:
Sys.setenv(BLUPF90_PLUS="/path/to/blupf90+")
testthat::test_local()
```

The CI job explicitly installs BLUPF90+ before testing, and tests the native
one-, two-, three- and five-trait paths. Additional parser/preparation tests cover
six and ten traits. Known-variance native EBVs and SEs are checked against exact
R calculations on an inbred pedigree. Tests can run without the executable;
only external-program integration tests are then skipped and must not be reported
as native validation. See VALIDATION.md for this delivery's actual results.

## Upload this release to GitHub

1. Extract `genetic-analysis-r-1.1.2.zip`. Upload the contents of its
   `genetic-analysis-r` folder as the repository root: DESCRIPTION, NAMESPACE,
   README.md, R/, tests/, inst/, man/ and the hidden .github/ directory.
2. Edit `Authors@R` in DESCRIPTION with your preferred author/maintainer name and
   email. The current `Package Maintainer <maintainer@example.org>` is a placeholder.
   After changing author metadata, regenerate documentation with
   `roxygen2::roxygenise()` so the package help page also identifies you.
3. Replace `YOUR_ACCOUNT` in this README with your GitHub username and adjust the
   repository name if you choose another one. The install command will then be
   `remotes::install_github("YOUR_ACCOUNT/genetic-analysis-r")`.
4. Include LICENSE and LICENSE.md, and review the copyright holder text before
   publishing. The R package uses MIT; external BLUPF90+/Intel software retains
   its own licensing terms and is downloaded separately by users.
5. Keep livestock records, credentials, installed R libraries, downloaded binaries,
   native job folders and runtime libraries outside the public repository. They
   are not included in the release ZIP.
6. After uploading, review the GitHub Actions result. When ready, create a release
   tag `v1.1.2` and attach `geneticAnalysisR_1.1.2.tar.gz` as its installable R source
   archive. Local validation is recorded in VALIDATION.md; it does not substitute
   for a successful remote CI run on your repository.

If you use Git locally, run these commands inside the extracted package folder
against a newly created, empty GitHub repository. Replace YOUR_ACCOUNT first:

```sh
git init -b main
git add .
git commit -m "Release geneticAnalysisR 1.1.2"
git remote add origin https://github.com/YOUR_ACCOUNT/genetic-analysis-r.git
git push -u origin main
git tag v1.1.2
git push origin v1.1.2
```

The `.tar.gz` archive is installed with `install.packages(..., repos=NULL,
type="source")`. It is not the GitHub repository ZIP. The ZIP contains the complete
repository, including its README and CI configuration; the tarball follows the
package's R-build exclusions. This delivery prepares files for your own upload.
