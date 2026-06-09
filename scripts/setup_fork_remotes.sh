#!/usr/bin/env bash
set -euo pipefail

# Configure remotes for fork-based development:
# - upstream -> original repository URL (current origin URL)
# - origin   -> your GitHub fork URL
#
# Usage:
#   # Option 1 (no gh required): pass your fork URL directly
#   bash scripts/setup_fork_remotes.sh "https://github.com/<you>/XClone.git"
#
#   # Option 2 (with gh authenticated): auto-create/find fork
#   bash scripts/setup_fork_remotes.sh
#
# Optional second arg:
#   bash scripts/setup_fork_remotes.sh "<fork_url>" "<upstream_url>"

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "Error: run this script inside a git repository."
  exit 1
fi

UPSTREAM_URL="${2:-$(git remote get-url origin)}"
FORK_URL="${1:-}"

if [[ -z "${FORK_URL}" ]]; then
  if command -v gh >/dev/null 2>&1 && gh auth status >/dev/null 2>&1; then
    REPO_SLUG="$(gh repo view --json nameWithOwner --jq .nameWithOwner)"
    FORK_OWNER="$(gh api user --jq .login)"
    FORK_SLUG="${FORK_OWNER}/$(gh repo view --json name --jq .name)"
    FORK_URL="https://github.com/${FORK_SLUG}.git"

    echo "Current repository: ${REPO_SLUG}"
    echo "Upstream URL:       ${UPSTREAM_URL}"
    echo "Fork target:        ${FORK_SLUG}"

    # Create fork if missing (safe to re-run).
    if ! gh repo view "${FORK_SLUG}" >/dev/null 2>&1; then
      echo "Creating fork ${FORK_SLUG} ..."
      gh repo fork "${REPO_SLUG}" --clone=false
    else
      echo "Fork already exists: ${FORK_SLUG}"
    fi
  else
    cat <<'EOF'
Error: fork URL not provided, and gh CLI is unavailable/not authenticated.

Either:
  1) Install/authenticate gh, then rerun:
     bash scripts/setup_fork_remotes.sh
  2) Provide your fork URL directly:
     bash scripts/setup_fork_remotes.sh "https://github.com/<you>/XClone.git"
EOF
    exit 1
  fi
fi

# Ensure upstream remote points to the original repository.
if git remote get-url upstream >/dev/null 2>&1; then
  git remote set-url upstream "${UPSTREAM_URL}"
else
  git remote add upstream "${UPSTREAM_URL}"
fi

# Set origin to your fork.
git remote set-url origin "${FORK_URL}"

echo
echo "Remotes configured:"
git remote -v
echo
echo "Done. You can now push branches to your fork (origin) and open PRs to upstream."
