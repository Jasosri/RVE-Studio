#!/usr/bin/env bash
# Publish RVE Studio as a public GitHub repository and trigger the first release (and, with Zenodo linked, a DOI).
#
#   tools/publish_github.sh <github-owner> [repo-name]
#
# Before running: (1) `gh auth login` as the owner; (2) on https://zenodo.org/account/settings/github/ sign in with
# GitHub and switch the repository ON as soon as it exists (step 3 below pauses for this). Zenodo archives every
# GitHub *release* published after the switch is on.
set -euo pipefail
OWNER=${1:?usage: tools/publish_github.sh <github-owner> [repo-name]}
REPO=${2:-RVE-Studio}
cd "$(dirname "$0")/.."
VERSION=$(sed -n "s/^__version__ *= *['\"]\([^'\"]*\)['\"].*/\1/p" rvestudio/__init__.py)

# 1. fill the repository placeholders
grep -rl "OWNER/RVE-Studio" README.md CITATION.cff pyproject.toml installer docs 2>/dev/null | xargs -r sed -i "s#OWNER/RVE-Studio#${OWNER}/${REPO}#g"

# 2. commit (the repository is already initialised with the release commit)
git add -A
git diff --cached --quiet || git commit -m "Set repository owner to ${OWNER}/${REPO}"

# 3. create the public repository and push (a clone of the .gitbundle has the bundle file as "origin": drop it)
git remote remove origin 2>/dev/null || true
gh repo create "${OWNER}/${REPO}" --public --source . --remote origin \
   --description "RVE Studio: composite cellular core design - homogenisation, FE validation and design of laminated-composite cellular cores" \
   --homepage "https://doi.org/10.1016/j.engstruct.2023.117106"
git branch -M main
git push -u origin main
read -r -p "Switch the repository ON in Zenodo (https://zenodo.org/account/settings/github/), then press Enter... "

# 4. tag and push: the Release workflow tests, builds wheel/sdist/PDFs and publishes the GitHub release
git tag -f -a "v${VERSION}" -m "RVE Studio ${VERSION}"   # move the tag to the commit with the owner filled in
git push origin "v${VERSION}"
echo "Release v${VERSION} pushed. When Zenodo has archived it, put the concept DOI into README.md (ZENODO_DOI) and"
echo "CITATION.cff (doi:), commit and push."
