dbboost='''PRAGMA cache_size = 20000;'''
dbcreate='''
CREATE TABLE project (
    "name" TEXT NOT NULL,
    "description" TEXT,
    "kind" TEXT,
    "tmp" TEXT,
    "creation" TEXT,
    "last" TEXT,
    "genome" TEXT,
    "phenome" TEXT,
    "pangenome" INTEGER   DEFAULT (0),
    "kegg" REAL
);
CREATE TABLE organism (
    "org_id" TEXT NOT NULL,
    "name" TEXT,
    "description" TEXT,
    "file" TEXT,
    "mutant" INTEGER,
    "reference" TEXT,
    "genome" TEXT   DEFAULT ('none'),
    "phenome" TEXT   DEFAULT ('none'),
    "mkind" TEXT,
	"color" TEXT
);
CREATE TABLE "protein" (
    "prot_id" TEXT NOT NULL,
    "org_id" TEXT NOT NULL,
    "description" TEXT,
    "sequence" TEXT
);
CREATE TABLE "mapko" (
    "prot_id" TEXT,
    "ko_id" TEXT,
    "indirect" INTEGER DEFAULT (0)
);
CREATE TABLE "ortholog" (
    "group_id" TEXT,
    "prot_id" TEXT
);
CREATE TABLE "ko" (
    "ko_id" TEXT NOT NULL,
    "name" TEXT,
    "description" TEXT,
    "analyzed" INTEGER   DEFAULT (0)
);
CREATE TABLE "reaction" (
    "re_id" TEXT NOT NULL,
    "name" TEXT,
    "description" TEXT,
    "enzyme" TEXT
);
CREATE TABLE "rpair" (
    "rp_id" TEXT NOT NULL,
    "co1" TEXT,
    "co2" TEXT,
    "kind" TEXT 
);
CREATE TABLE "compound" (
    "co_id" TEXT NOT NULL,
    "name" TEXT,
    "description" TEXT
);
CREATE TABLE "pathway" (
    "path_id" TEXT NOT NULL,
    "name" TEXT,
    "description" TEXT,
    "html" TEXT
);
CREATE TABLE "ko_react" (
    "ko_id" TEXT NOT NULL,
    "re_id" TEXT NOT NULL
);
CREATE TABLE "react_comp" (
    "re_id" TEXT NOT NULL,
    "co_id" TEXT NOT NULL
);
CREATE TABLE "react_path" (
    "re_id" TEXT NOT NULL,
    "path_id" TEXT NOT NULL
);
CREATE TABLE "comp_path" (
    "co_id" TEXT NOT NULL,
    "path_id" TEXT NOT NULL
);
CREATE TABLE "rpair_react" (
    "rp_id" TEXT NOT NULL,
    "re_id" TEXT NOT NULL
);
CREATE TABLE "pathmap" (
    "path_id" TEXT NOT NULL,
    "png" BLOB,
    "html" TEXT
);
CREATE TABLE biolog (
    "plate_id" TEXT NOT NULL,
    "well_id" TEXT NOT NULL,
    "concentration" INTEGER,
    "zero_well_id" TEXT,
    "chemical" TEXT,
    "cas_id" TEXT,
    "co_id" TEXT,
    "moa" TEXT
, "category" INTEGER
);
CREATE TABLE biolog_exp (
    "plate_id" TEXT NOT NULL,
    "well_id" TEXT NOT NULL,
    "org_id" TEXT NOT NULL,
    "replica" INTEGER NOT NULL,
    "activity" INTEGER,
    "zero" INTEGER,
    "min" REAL,
    "max" REAL,
    "height" REAL,
    "plateau" REAL,
    "slope" REAL,
    "lag" REAL,
    "area" REAL,
    "v" REAL,
    "y0" REAL,
    "model" TEXT,
    "source" TEXT
);
CREATE TABLE biolog_exp_det (
    "plate_id" TEXT NOT NULL,
    "well_id" TEXT NOT NULL,
    "org_id" TEXT NOT NULL,
    "replica" INTEGER NOT NULL,
    "times" TEXT,
    "signals" TEXT
);
CREATE TABLE biolog_purged_exp (
    "plate_id" TEXT NOT NULL,
    "well_id" TEXT NOT NULL,
    "org_id" TEXT NOT NULL,
    "replica" INTEGER NOT NULL,
    "activity" INTEGER,
    "zero" INTEGER,
    "min" REAL,
    "max" REAL,
    "height" REAL,
    "plateau" REAL,
    "slope" REAL,
    "lag" REAL,
    "area" REAL,
    "v" REAL,
    "y0" REAL,
    "model" TEXT,
    "source" TEXT
);
CREATE TABLE biolog_purged_exp_det (
    "plate_id" TEXT NOT NULL,
    "well_id" TEXT NOT NULL,
    "org_id" TEXT NOT NULL,
    "replica" INTEGER NOT NULL,
    "times" TEXT,
    "signals" TEXT
);
CREATE UNIQUE INDEX project_id ON project(name ASC);
CREATE UNIQUE INDEX "organism_id" on organism (org_id ASC);
CREATE UNIQUE INDEX "protein_id" on protein (prot_id ASC, org_id ASC);
CREATE UNIQUE INDEX "mapko_id" on mapko (prot_id ASC, ko_id ASC);
CREATE UNIQUE INDEX "ortholog_id" on ortholog (group_id ASC, prot_id ASC);
CREATE UNIQUE INDEX "ko_id" on ko (ko_id ASC);
CREATE UNIQUE INDEX "reaction_id" on reaction (re_id ASC);
CREATE UNIQUE INDEX "rpair_id" on rpair (rp_id ASC);
CREATE UNIQUE INDEX "compound_id" on compound (co_id ASC);
CREATE UNIQUE INDEX "pathway_id" on pathway (path_id ASC);
CREATE UNIQUE INDEX "koreact_id" on ko_react (ko_id ASC, re_id ASC);
CREATE UNIQUE INDEX "reactcomp_id" on react_comp (re_id ASC, co_id ASC);
CREATE UNIQUE INDEX "reactpath_id" on react_path (re_id ASC, path_id ASC);
CREATE UNIQUE INDEX "comppath_id" on comp_path (co_id ASC, path_id ASC);
CREATE UNIQUE INDEX "rpairreact_id" on rpair_react (rp_id ASC, re_id ASC);
CREATE UNIQUE INDEX "pathmap_id" on pathmap (path_id ASC);
CREATE UNIQUE INDEX "biolog_id" on biolog (plate_id ASC, well_id ASC);
CREATE UNIQUE INDEX "biologexp_id" on biolog_exp (plate_id ASC, well_id ASC, org_id ASC, replica ASC);
CREATE UNIQUE INDEX "biologexpdet_id" on biolog_exp_det (plate_id ASC, well_id ASC, org_id ASC, replica ASC);
CREATE UNIQUE INDEX "biologpurgedexp_id" on biolog_purged_exp (plate_id ASC, well_id ASC, org_id ASC, replica ASC);
CREATE UNIQUE INDEX "biologpurgedexpdet_id" on biolog_purged_exp_det (plate_id ASC, well_id ASC, org_id ASC, replica ASC);
'''
