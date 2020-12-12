#!/bin/bash
# run from home dirsctory, not from python

mkdir toy_pipeline
cp -r python toy_pipeline
rm -rf toy_pipeline/python/__pycache__ toy_pipeline/python/region_summary.xlsx
mkdir toy_pipeline/report
cp report/report.tex  report/report.bib report/plos2009.bst toy_pipeline/report
cp -r report/figures toy_pipeline/report
cp -r fastqc toy_pipeline
find toy_pipeline/fastqc -name "*.zip" | xargs rm -f
mkdir toy_pipeline/task_dna
mv toy_pipeline/python/README.md toy_pipeline/
zip -r toy_pipeline.zip toy_pipeline/
