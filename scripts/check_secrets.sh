#!/bin/bash
# Lightweight secret/PII regression check over the tracked worktree.
#
# Scope limit, stated plainly: this is pattern-matching against a short list
# of known credential shapes, not entropy-based detection — a secret in an
# unlisted shape (a generic API key, a JWT, a DB connection string) will sail
# straight through. This is a cheap, no-dependency backstop against the known
# shapes, not a replacement for the manual audits already done on this repo,
# and not comparable to a dedicated tool like gitleaks/trufflehog.
#
# Run from the repo root: bash scripts/check_secrets.sh
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

# This script's own pattern definitions don't currently self-match any of the
# patterns below (checked by hand: each pattern's literal regex syntax in the
# array below breaks its own match, e.g. the literal "(" right after "AKIA["
# in the pattern text isn't itself a char in [0-9A-Z]). So this exclude isn't
# actually load-bearing today — it's kept anyway as cheap insurance, since a
# less carefully-written pattern added later easily could self-match without it.
EXCLUDE=':(exclude)scripts/check_secrets.sh'

PATTERNS=(
  '-----BEGIN (RSA|OPENSSH|EC|DSA|PGP) PRIVATE KEY-----'
  'AKIA[0-9A-Z]{16}'
  'gh[pousr]_[A-Za-z0-9]{36,}'
  'github_pat_[A-Za-z0-9_]{70,}'
  'xox[baprs]-[A-Za-z0-9-]{10,}'
  '/Users/[A-Za-z0-9_.-]+'
  '/home/[A-Za-z0-9_.-]+'
)

found=0
for pattern in "${PATTERNS[@]}"; do
  set +e
  # -e explicitly marks $pattern as the pattern, not an option flag — needed
  # because a pattern starting with "-" (e.g. "-----BEGIN ...") would
  # otherwise be misparsed as a git-grep flag and fail with a usage error.
  output=$(git grep -InE -e "$pattern" -- . "$EXCLUDE" 2>&1)
  status=$?
  set -e

  if [ "$status" -eq 0 ]; then
    echo "::error::Potential secret/PII pattern matched: $pattern"
    echo "$output"
    found=1
  elif [ "$status" -gt 1 ]; then
    # status 1 just means "no match" (the common case); anything higher is a
    # real tooling error (e.g. a bad pattern) and must not be swallowed as if
    # it were a clean "nothing found".
    echo "::error::git grep itself failed (exit $status) on pattern: $pattern"
    echo "$output"
    found=1
  fi
done

if [ "$found" -eq 1 ]; then
  echo
  echo "One or more suspicious patterns (or a tooling error) found above. Review before merging."
  exit 1
fi

echo "No secret/PII patterns found."
