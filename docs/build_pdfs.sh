#!/usr/bin/env bash
# Build the PDF user guide and theory guide from the Markdown sources (pandoc + xelatex).
set -euo pipefail
export LANG=C.UTF-8 LC_ALL=C.UTF-8
cd "$(dirname "$0")"
mkdir -p pdf
VERSION=$(python3 -c "import re,pathlib;print(re.search(r'__version__\s*=\s*[\"\x27]([^\"\x27]+)',pathlib.Path('../rvestudio/__init__.py').read_text()).group(1))")
build () {  # $1 source, $2 output, $3 title, $4 subtitle
  python3 - "$1" > /tmp/rvestudio_doc_$$.md <<'PY'
import re, sys
s = open(sys.argv[1], encoding='utf-8').read()
s = re.sub(r'\A# .*\n', '', s, count=1)                                  # title comes from the metadata
s = re.sub(r'\n## Contents\n.*?\n---\n', '\n', s, count=1, flags=re.S)   # pandoc writes the table of contents
s = re.sub(r'\]\((USER_GUIDE|THEORY|API_REFERENCE|FE_VERIFICATION|DEVELOPER_GUIDE)\.md(#[^)]*)?\)', r'](\1.md)', s)
sys.stdout.write(s)
PY
  pandoc /tmp/rvestudio_doc_$$.md -o "$2" --pdf-engine=xelatex --toc --toc-depth=2 --shift-heading-level-by=-1 \
    --resource-path=.:.. -V geometry:margin=2.2cm -V mainfont="DejaVu Serif" -V sansfont="DejaVu Sans" \
    -V monofont="DejaVu Sans Mono" -V fontsize=10pt -V colorlinks=true -V linkcolor=blue -V urlcolor=blue \
    -V toccolor=black -V documentclass=article -V papersize=a4 \
    --metadata-file=<(printf 'title: "%s"\nsubtitle: "%s"\nauthor: "RVE Studio %s"\ndate: "%s"\n' "$3" "$4" "$VERSION" "$(date +%Y-%m-%d)") \
    -H <(printf '%s\n' '\usepackage{float}' '\let\origfigure\figure' '\let\endorigfigure\endfigure' \
         '\renewenvironment{figure}[1][2]{\origfigure[H]}{\endorigfigure}' '\usepackage{fvextra}' \
         '\DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\{\}}' \
         '\setlength{\emergencystretch}{3em}' '\let\oldtoc\tableofcontents' '\renewcommand{\tableofcontents}{\oldtoc\newpage}')
  rm -f /tmp/rvestudio_doc_$$.md
}
build USER_GUIDE.md pdf/RVE_Studio_User_Guide.pdf "RVE Studio — User Guide" "Composite cellular core design: homogenisation, FE design and validation, design databases"
build THEORY.md     pdf/RVE_Studio_Theory_Guide.pdf "RVE Studio — Theory and Modelling Guidance" "Strain-energy homogenisation, wall models, verified limits and mesh convergence"
mkdir -p ../rvestudio/docs && cp pdf/RVE_Studio_User_Guide.pdf pdf/RVE_Studio_Theory_Guide.pdf ../rvestudio/docs/   # opened from Help
echo "Written docs/pdf/RVE_Studio_User_Guide.pdf and docs/pdf/RVE_Studio_Theory_Guide.pdf"
