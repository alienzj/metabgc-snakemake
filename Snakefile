#!/usr/bin/env snakemake

shell.executable("bash")

configfile: "config.yaml"

include: "rules/build.smk"
include: "rules/identify.smk"
include: "rules/quantify.smk"
include: "rules/cluster.smk"
