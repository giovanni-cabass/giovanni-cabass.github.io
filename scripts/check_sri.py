#!/usr/bin/env python3
"""Verify every pinned CDN resource's SRI hash still matches what's declared
in the source templates. Run from the repo root: python3 scripts/check_sri.py
"""
import base64
import concurrent.futures
import glob
import hashlib
import re
import sys
import time
import urllib.request

# Auto-discovered, not hardcoded: any new _includes/_layouts file that ever
# declares its own integrity="..." tag is picked up automatically, instead of
# silently going unchecked until someone remembers to add it to a list here.
FILES = sorted(glob.glob("_includes/*.html") + glob.glob("_layouts/*.html"))

TAG_RE = re.compile(r"<(?:script|link)\b[^>]*>", re.IGNORECASE)
URL_RE = re.compile(r'''(?:src|href)=["'](https://[^"']+)["']''', re.IGNORECASE)
INTEGRITY_RE = re.compile(
    r'''integrity=["']sha(256|384|512)-([^"']+)["']''', re.IGNORECASE
)

HASH_FUNCS = {"256": hashlib.sha256, "384": hashlib.sha384, "512": hashlib.sha512}

MAX_ATTEMPTS = 5
INITIAL_RETRY_DELAY_SECONDS = 2
REQUEST_TIMEOUT_SECONDS = 10  # these are small CSS/JS files; 10s is generous
MAX_CONCURRENT_FETCHES = 16  # fetches are independent and I/O-bound, so run
# them in parallel threads instead of one after another: the GIL is released
# during network I/O, so this bounds total wall time by the slowest single
# resource's retry chain instead of the sum over every resource.


def gha_warning(message, title):
    """Emit a GitHub Actions warning annotation (visible in the PR/run UI
    even on a passing job), in addition to the plain log line."""
    escaped = message.replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")
    print(f"::warning title={title}::{escaped}")


def fetch(url):
    """Fetch url, retrying with exponential backoff on transient network
    errors. Raises on final failure."""
    last_error = None
    delay = INITIAL_RETRY_DELAY_SECONDS
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(url, timeout=REQUEST_TIMEOUT_SECONDS) as resp:
                return resp.read()
        except Exception as e:
            last_error = e
            if attempt < MAX_ATTEMPTS:
                print(f"  [{url}] attempt {attempt}/{MAX_ATTEMPTS} failed ({e}), retrying in {delay}s...")
                time.sleep(delay)
                delay *= 2
    raise last_error


def main():
    # failures = security-relevant: a hash we *did* check no longer matches,
    # or a local file we expected to read couldn't be read. These must fail CI.
    # warnings = we couldn't even reach the CDN to check it (its outage, not
    # ours) — surfaced loudly but must not block an unrelated merge/deploy.
    failures = []
    warnings = []
    tasks = []  # (path, url, algo, expected_b64)

    for path in FILES:
        try:
            with open(path, encoding="utf-8") as f:
                content = f.read()
        except OSError as e:
            failures.append(f"{path}: could not read file ({e})")
            continue

        for tag in TAG_RE.findall(content):
            url_m = URL_RE.search(tag)
            integ_m = INTEGRITY_RE.search(tag)
            if not (url_m and integ_m):
                continue

            url = url_m.group(1)
            algo, expected_b64 = integ_m.group(1), integ_m.group(2)
            tasks.append((path, url, algo, expected_b64))
            print(f"Checking {url} ({path}) ...")

    def verify_one(task):
        _path, url, algo, expected_b64 = task
        try:
            data = fetch(url)
        except Exception as e:
            return "warning", f"{url}: fetch failed after {MAX_ATTEMPTS} attempts ({e})"

        digest = HASH_FUNCS[algo](data).digest()
        actual_b64 = base64.b64encode(digest).decode()
        if actual_b64 != expected_b64:
            return "failure", (
                f"{url}: SRI MISMATCH — declared sha{algo}-{expected_b64}, "
                f"actual sha{algo}-{actual_b64}"
            )
        return "ok", None

    if tasks:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=min(MAX_CONCURRENT_FETCHES, len(tasks))
        ) as pool:
            for kind, message in pool.map(verify_one, tasks):
                if kind == "warning":
                    warnings.append(message)
                elif kind == "failure":
                    failures.append(message)

    checked = len(tasks)
    print(f"\nChecked {checked} pinned resources across {len(FILES)} files.")
    if warnings:
        print("WARNINGS (CDN unreachable, not verified this run):")
        for w in warnings:
            print(" -", w)
            gha_warning(w, title="SRI-check-CDN-unreachable")

        if len(warnings) == checked:
            # Every single resource was unreachable: nothing was actually
            # verified this run, so this isn't "one flaky CDN" anymore — it's
            # the whole check being a no-op. Don't let that pass silently.
            failures.append(
                f"all {checked} pinned resources were unreachable this run — "
                "nothing could be verified, refusing to pass"
            )

    if failures:
        print("FAILURES:")
        for f in failures:
            print(" -", f)
        sys.exit(1)

    if warnings:
        print("No SRI mismatches, but some resources could not be verified (see warnings above).")
    else:
        print("All SRI hashes verified OK.")


if __name__ == "__main__":
    main()
