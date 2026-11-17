# giovanni-cabass.github.io

## Development note

This site was built with help from Claude Code: page structure, styling, and content were developed collaboratively with it as a coding assistant.

(Earlier commits were missing the git trailer so GitHub never recognized Claude Code as co-author.)

Design inspired by [al-folio](https://github.com/alshedivat/al-folio) (also credited in the footer) — built independently from scratch (vanilla HTML/CSS/JS, no theme gem or `remote_theme` dependency), not a fork of its code.

## Notes

- The CV embed on the Contact page doesn't use an iframe `sandbox` attribute: Chrome's built-in PDF viewer needs script execution to render inline, and `sandbox="allow-same-origin"` (no `allow-scripts`) made Chrome block the embed outright. The PDF is verified to contain no embedded JavaScript, `/OpenAction`, or `/Launch` actions — it does carry 3 ordinary `/URI` hyperlink actions (email, personal site, InspireHEP), no different in kind from a normal `<a href>` link. If that ever changes, re-check before relying on this. See: [Chromium issue 413851](https://groups.google.com/a/chromium.org/g/chromium-bugs/c/z5hA0H6LFws) and [this write-up of the same failure mode](https://lightrun.com/answers/cburgmer-rasterizehtml-js-documents-frame-is-sandboxed-and-the-allow-scripts-permission-is-not-set) for why `allow-scripts` is needed; [WebKit's bug discussion](https://bugs.webkit.org/show_bug.cgi?id=218086) for why `allow-scripts`+`allow-same-origin` together isn't the right fix either.
- The CI "Verify pinned CDN resources still match their SRI hashes" step tolerates a CDN being transiently unreachable: it retries with backoff, then logs a GitHub Actions `::warning::` annotation and still passes, rather than failing the whole build over someone else's outage. It only fails if *every* pinned resource is unreachable (nothing could be verified) or if a fetched resource's hash actually mismatches what's declared — that's a real tampered/changed resource, not a network blip. See `scripts/check_sri.py`.
