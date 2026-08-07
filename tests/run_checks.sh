#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

python_cmd="${FDRREG_PYTHON:-python}"
"$python_cmd" -m unittest discover -s tests -p 'test_*.py' -v

while IFS= read -r -d '' file; do "$python_cmd" -m py_compile "$file"; done < <(find scripts -name '*.py' -print0)

for file in scripts/run_pipeline.sh scripts/check_inputs.sh scripts/imputation/*.sh scripts/ldsc/*.sh scripts/real/v7/*.sh scripts/real/v8/*.sh; do
  [[ -f "$file" ]] && bash -n "$file"
done

bash scripts/run_pipeline.sh --help >/dev/null
[[ "$(bash scripts/run_pipeline.sh --stage list | wc -l)" -eq 16 ]]
bash scripts/imputation/run_dist_imputation.sh --help >/dev/null

if command -v Rscript >/dev/null 2>&1; then
  while IFS= read -r -d '' file; do Rscript -e "parse(file='$file')" >/dev/null; done < <(find scripts -path '*/legacy' -prune -o \( -name '*.R' -o -name '*.r' \) -print0)
  Rscript tests/test_fdreg.R
  Rscript tests/test_imputation.R
else
  echo "SKIP: Rscript is not installed; R parse checks were not run." >&2
fi

echo "All available checks passed."
