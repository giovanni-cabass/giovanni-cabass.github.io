#!/usr/bin/env python3
"""Verify the CSP's hardcoded inline-script sha256 hashes still match what's
actually shipped. Scans the *built* _site/ output (not the source templates):
neither inline script currently contains any Liquid, so source and built
output hash identically today, but nothing enforces that staying true, and
scanning the built output has no such blind spot.

No network calls — unlike check_sri.py, everything needed is already on disk
once the site is built. Run from the repo root, after `bundle exec jekyll
build`: python3 scripts/check_csp_hashes.py
"""
import base64
import glob
import hashlib
import re
import sys
from html.parser import HTMLParser

PAGE_GLOB = "_site/**/*.html"


class PageScanner(HTMLParser):
    """Collects every inline <script> block's raw text and the CSP meta
    tag's content="..." attribute from a single built HTML page.

    convert_charrefs is left at its default (True): verified directly that
    this doesn't affect <script> content either way — html.parser already
    treats <script>/<style> as raw text and never decodes character
    references inside them, matching real browser behavior, regardless of
    this flag. (Confirmed with a one-off test before writing this.)
    """

    def __init__(self):
        super().__init__()
        self.inline_scripts = []
        self.csp_content = None
        self._capturing = False
        self._buf = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "script" and "src" not in attrs:
            # Only scripts with no src= have inline content subject to a
            # hash-source check; a <script src="..."> tag's text content is
            # empty, and hashing an empty string collides with the
            # well-known SHA-256-of-empty-string — meaningless noise that
            # must never enter the comparison.
            self._capturing = True
            self._buf = []
        elif tag == "meta" and (attrs.get("http-equiv") or "").lower() == "content-security-policy":
            self.csp_content = attrs.get("content")

    def handle_data(self, data):
        if self._capturing:
            self._buf.append(data)

    def handle_endtag(self, tag):
        if tag == "script" and self._capturing:
            self.inline_scripts.append("".join(self._buf))
            self._capturing = False


def script_src_hashes(csp_content):
    """Extract the sha256-... values from just the script-src directive."""
    for directive in csp_content.split(";"):
        parts = directive.strip().split()
        if parts and parts[0] == "script-src":
            return {f"sha256-{h}" for h in re.findall(r"'sha256-([A-Za-z0-9+/=]+)'", directive)}
    return set()


def hash_of(text):
    digest = hashlib.sha256(text.encode("utf-8")).digest()
    return "sha256-" + base64.b64encode(digest).decode()


def main():
    pages = sorted(glob.glob(PAGE_GLOB, recursive=True))
    if not pages:
        print(f"No pages matched {PAGE_GLOB} — did you run `bundle exec jekyll build` first?")
        sys.exit(1)

    failures = []
    warnings = []
    all_declared = {}  # page -> frozenset of declared hashes
    found_anywhere = set()
    mismatched_pages = {}  # hash -> list of pages where it wasn't declared

    for page in pages:
        with open(page, encoding="utf-8") as f:
            content = f.read()

        scanner = PageScanner()
        scanner.feed(content)

        if scanner.csp_content is None:
            # A page with no CSP meta tag at all (e.g. a future layout that
            # bypasses head.html) is itself the failure, not a "0 declared
            # hashes" case to silently compare against — fail loudly and
            # name the page, rather than letting a later .split() crash on
            # None stand in for a real error message.
            failures.append(f"{page}: no CSP meta tag found on this page at all")
            continue

        declared = script_src_hashes(scanner.csp_content)
        all_declared[page] = frozenset(declared)

        for script_text in scanner.inline_scripts:
            h = hash_of(script_text)
            found_anywhere.add(h)
            if h not in declared:
                mismatched_pages.setdefault(h, []).append(page)

    # Report mismatches grouped by hash, not by page: one edit typically
    # breaks the same hash on every page at once (the theme-flash script
    # renders on all of them), so reporting it once per page would just be
    # noise repeating the same root cause.
    for h, pages_hit in sorted(mismatched_pages.items()):
        failures.append(
            f"Inline <script> hash {h} (found on {len(pages_hit)}/{len(pages)} pages, "
            f"e.g. {pages_hit[0]}) is NOT in the CSP script-src allow-list on those pages. "
            f"If this script's content changed intentionally, update the hash in "
            f"_includes/head.html's CSP <meta> tag to this value."
        )

    # Sanity check: every page shares one static CSP today (all pages chain
    # through _layouts/default.html -> _includes/head.html), so their
    # declared sets should be identical. This is a tripwire tied to today's
    # architecture, not a permanent invariant — if the site ever
    # deliberately varies CSP per page, relax or remove this check then,
    # rather than being surprised by it failing on a legitimate change.
    distinct_sets = set(all_declared.values())
    if len(distinct_sets) > 1:
        failures.append(
            "Pages declare different CSP script-src hash sets, but every page is "
            "expected to share the same one (all routed through head.html): "
            + ", ".join(f"{p}: {sorted(h)}" for p, h in all_declared.items())
        )

    # Orphaned hashes: checked globally, not per-page. Every page's declared
    # set includes the MathJax hash even on pages that never render MathJax
    # (same shared-CSP reason as above) — so "declared but not found on
    # *this* page" would misfire on every non-post page, permanently, for a
    # hash that's legitimately just used elsewhere on the site.
    all_declared_hashes = set().union(*all_declared.values()) if all_declared else set()
    for h in sorted(all_declared_hashes - found_anywhere):
        warnings.append(f"CSP declares {h} but no inline <script> anywhere on the site hashes to it (stale/orphaned entry?)")

    print(f"Scanned {len(pages)} built pages.")

    if warnings:
        print("\nWARNINGS (non-fatal):")
        for w in warnings:
            print(" -", w)

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(" -", f)
        sys.exit(1)

    print("All inline <script> hashes verified OK.")


if __name__ == "__main__":
    main()
