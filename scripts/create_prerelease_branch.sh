#!/usr/bin/env bash
set -euo pipefail

# Create and publish a pre-release branch with current local updates.
# It will:
#   1) create/switch to the target branch
#   2) stage all non-ignored changes
#   3) commit (if needed)
#   4) push to origin and set upstream
#
# Usage:
#   bash scripts/create_prerelease_branch.sh [branch_name] [commit_message]
#
# Example:
#   bash scripts/create_prerelease_branch.sh prerelease/vb-clustering \
#     "Prepare pre-release branch for VB clustering integration"

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "Error: run this script inside a git repository."
  exit 1
fi

BRANCH_NAME="${1:-prerelease/vb-clustering-$(date +%Y%m%d-%H%M)}"
COMMIT_MSG="${2:-Prepare pre-release branch for VB clustering integration}"

echo "Target branch: ${BRANCH_NAME}"

# Create/switch branch from current HEAD while preserving local changes.
if git show-ref --verify --quiet "refs/heads/${BRANCH_NAME}"; then
  git switch "${BRANCH_NAME}"
else
  git switch -c "${BRANCH_NAME}"
fi

# Stage everything except ignored files (e.g. notebooks/tests if ignored).
git add -A

if git diff --cached --quiet; then
  echo "No staged changes to commit."
else
  git commit -m "${COMMIT_MSG}"
fi

git push -u origin "${BRANCH_NAME}"

echo
echo "Done. Branch pushed: ${BRANCH_NAME}"

# Helpful compare URL if upstream exists.
if git remote get-url upstream >/dev/null 2>&1; then
  UPSTREAM_URL="$(git remote get-url upstream)"
  REPO_PATH="$(echo "${UPSTREAM_URL}" | sed -E 's#(https://github.com/|git@github.com:)([^.]+)(\\.git)?#\\2#')"
  echo "Compare URL: https://github.com/${REPO_PATH}/compare/${BRANCH_NAME}?expand=1"
fi
