# giovanni-cabass.github.io

## Development note

This site was built with help from Claude Code: page structure, styling, and content were developed collaboratively with it as a coding assistant.

(Earlier commits were missing the git trailer so GitHub never recognized Claude Code as co-author.)

## Notes

- The CV embed on the Contact page doesn't use an iframe `sandbox` attribute: Chrome's built-in PDF viewer needs script execution to render inline, and `sandbox="allow-same-origin"` (no `allow-scripts`) made Chrome block the embed outright. The PDF is verified to contain no embedded JavaScript, `/OpenAction`, or launch/URI actions — if that ever changes, re-check before relying on this.
