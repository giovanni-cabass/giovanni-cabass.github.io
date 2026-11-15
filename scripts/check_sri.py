#!/usr/bin/env python3
"""Verify every pinned CDN resource's SRI hash still matches what's declared
in the source templates. Run from the repo root: python3 scripts/check_sri.py
"""
import base64
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

MAX_ATTEMPTS = 3
RETRY_DELAY_SECONDS = 2


def fetch(url):
    """Fetch url, retrying on transient network errors. Raises on final failure."""
    last_error = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                return resp.read()
        except Exception as e:
            last_error = e
            if attempt < MAX_ATTEMPTS:
                print(f"  attempt {attempt}/{MAX_ATTEMPTS} failed ({e}), retrying...")
                time.sleep(RETRY_DELAY_SECONDS)
    raise last_error


def main():
    failures = []
    checked = 0

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
            checked += 1
            print(f"Checking {url} ({path}) ...")

            try:
                data = fetch(url)
            except Exception as e:
                failures.append(f"{url}: fetch failed after {MAX_ATTEMPTS} attempts ({e})")
                continue

            digest = HASH_FUNCS[algo](data).digest()
            actual_b64 = base64.b64encode(digest).decode()
            if actual_b64 != expected_b64:
                failures.append(
                    f"{url}: SRI MISMATCH — declared sha{algo}-{expected_b64}, "
                    f"actual sha{algo}-{actual_b64}"
                )

    print(f"\nChecked {checked} pinned resources across {len(FILES)} files.")
    if failures:
        print("FAILURES:")
        for f in failures:
            print(" -", f)
        sys.exit(1)

    print("All SRI hashes verified OK.")


if __name__ == "__main__":
    main()
