#!/bin/bash
# One command to rebuild the report after ANY change to scoring/results or to the text (2026-10-13).
#   bash report/build/rebuild_report.sh            build; stop if the text and the results disagree
#   bash report/build/rebuild_report.sh --fix      first rewrite the section 6 tables of the text from the results
# Order: appendix tables -> section 6 numbers -> consistency check -> Word file -> check of the Word file itself
# -> copy into the reproducibility package, new zip, copy to Google Drive.
# New models (e.g. CryoZeta) must be scored first: score_benchmark.py -> usalign_fixed_build_audit.py -> apply_fixed_audit.py.
set -e
R=/Users/anony/Downloads/MICA_GAIVI_run; PY=$R/scoring/.venv/bin/python; cd $R
SOFFICE=(soffice); BACKUP_DIR=""                      # LibreOffice command for the PDF render; machine-specific overrides below
[ -f report/build/local_env.sh ] && . report/build/local_env.sh
$PY scoring/make_tables.py > report/build/last_make_tables.log 2>&1
$PY report/build/make_section6.py > report/build/last_section6_facts.log 2>&1
$PY report/build/make_figures.py                           # Figures 1 and 2 are redrawn from the current results at every build
[ "$1" = "--fix" ] && $PY report/build/sync_report_tables.py --fix
$PY report/build/sync_report_tables.py                     # exits 1 on any inconsistency
( cd report/build && node build_report.js ) | tail -n 2
# the Word file itself: the per-case TM table of section 6.1 and appendix Table 1 must be identical cell for cell
$PY - <<'PY'
import re, sys, zipfile
x = zipfile.ZipFile("/Users/anony/Downloads/MICA_GAIVI_run/report/MICA_benchmark_report.docx").read("word/document.xml").decode()
tabs = []   # per-case tables with 5 method columns, in document order: the first is section 6.1, the second appendix Table 1 (TM-score)
for tb in re.findall(r"<w:tbl>.*?</w:tbl>", x, flags=re.S):
    rows = {}
    for tr in re.findall(r"<w:tr[ >].*?</w:tr>", tb, flags=re.S):
        c = ["".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", tc)).strip() for tc in re.findall(r"<w:tc>.*?</w:tc>", tr, flags=re.S)]
        if len(c) == 9 and re.fullmatch(r"[1-6]", c[1]) and all(re.fullmatch(r"0\.\d{3}", v) for v in c[4:]): rows[c[0]] = c[4:]
    if len(rows) == 24: tabs.append(rows)
bad = {k: (tabs[0][k], tabs[1].get(k)) for k in tabs[0] if tabs[0][k] != tabs[1].get(k)} if len(tabs) >= 2 else "tables not found"
if bad: print("WORD FILE CHECK FAILED:", bad); sys.exit(1)
print("Word file: section 6.1 and appendix Table 1 agree for all 24 cases x 5 methods")
PY
P=reproducibility_package/report
cp report/MICA_benchmark_report.docx report/REPORT_DRAFT.md report/RESULTS_TABLES.md report/raw_results_24cases.csv $P/
cp report/build/build_report.js report/build/make_figures.py report/build/make_section6.py report/build/sync_report_tables.py report/build/sync_package_cryozeta.py report/build/validate_report.py report/build/rebuild_report.sh $P/build/
cp report/figures/image1.png report/figures/image2.png $P/figures/
# result files keep the package's names: CryoZeta.json -> CryoZeta_interim.json, CryoZeta_partial.json -> CryoZeta_partial_models.json
for f in scoring/results/*.json; do
  n=$(basename $f); [ $n = CryoZeta.json ] && n=CryoZeta_interim.json; [ $n = CryoZeta_partial.json ] && n=CryoZeta_partial_models.json
  cp $f reproducibility_package/scoring/results/$n
done
cp scoring/usalign_fixed_build_audit.json scoring/make_tables.py reproducibility_package/scoring/
$PY report/build/sync_package_cryozeta.py           # CryoZeta models, run records and manifest rows; exits 1 if a file is not the scored one
rm -f MICA_benchmark_reproducibility_package.zip && zip -qr -X MICA_benchmark_reproducibility_package.zip reproducibility_package -x '*.DS_Store'
D="${BACKUP_DIR:-}"                                    # optional second copy; set in report/build/local_env.sh
if [ -n "$D" ] && [ -d "$D" ]; then cp report/MICA_benchmark_report.docx report/REPORT_DRAFT.md report/RESULTS_TABLES.md report/raw_results_24cases.csv MICA_benchmark_reproducibility_package.zip "$D"/ && echo "backup copy updated"; else echo "no backup folder set"; fi
ls -la report/MICA_benchmark_report.docx MICA_benchmark_reproducibility_package.zip | awk '{print $5, $6, $7, $8, $9}'
# render the Word file for the page checks, then run every check and write report/qa/QA_RECORD.md (exit 1 on any FAIL)
rm -rf report/qa/render && mkdir -p report/qa/render && cp report/MICA_benchmark_report.docx report/qa/render/
( cd report/qa/render && "${SOFFICE[@]}" --headless --convert-to pdf MICA_benchmark_report.docx > /dev/null 2>&1; rm -f MICA_benchmark_report.docx ) || true
$PY report/build/validate_report.py
[ -n "$D" ] && [ -d "$D" ] && cp report/qa/QA_RECORD.md "$D"/
echo REBUILD-OK
