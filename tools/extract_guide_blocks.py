"""Extract the Python blocks of the user guide, with the script each belongs to.

This is the parser for the documentation verification procedure. It reads
docs/USER_GUIDE.md, finds every fenced block opened with ```python, and writes
them in guide order to a JSON file as a list of objects with three keys:
``label`` (the last line of prose above the block), ``script`` (the file name
the label assigns the block to, or null) and ``body`` (the code itself).

It relies on the label contract of the guide. The paragraph of prose directly
above a code block says what the block is, in one of these forms:

    Start a new file named `x.py`.
    Continues `x.py`. Add these lines at the end of the file.
    The whole of `x.py`:
    Settings rather than Python, save as `path`:

The first three assign a Python block to a script; a label may wrap over more
than one line. The fourth marks a settings block, which is not fenced as
python and so is not extracted. A block with no label above it is the printed
output of the block before it, and is fenced as plain text, so it is not
extracted either. Rewording a label drops its block out of the check, which
shows up here as a block whose script is null.

A prose only edit to the guide is verified by running this before and after
the edit and confirming that the two JSON files are identical.

Usage:

    python tools/extract_guide_blocks.py docs/USER_GUIDE.md out.json
"""

import json
import re
import sys
from pathlib import Path

LABELS = [
    re.compile(r"Start a new file named `([^`]+)`\."),
    re.compile(r"Continues `([^`]+)`\. Add these lines at the end of the file\."),
    re.compile(r"The whole of `([^`]+)`:"),
]


def extract_blocks(lines):
    """Return the Python blocks in ``lines`` as dicts of label, script and body."""
    blocks = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("```python"):
            j = i - 1
            while j >= 0 and not lines[j].strip():
                j -= 1
            k = j
            while k > 0 and lines[k - 1].strip() and not lines[k - 1].startswith("```"):
                k -= 1
            paragraph = " ".join(line.strip() for line in lines[k : j + 1])
            script = None
            for pattern in LABELS:
                match = pattern.search(paragraph)
                if match:
                    script = match.group(1)
            end = i + 1
            while not lines[end].startswith("```"):
                end += 1
            blocks.append(
                {
                    "label": lines[j].strip() if j >= 0 else "",
                    "script": script,
                    "body": "\n".join(lines[i + 1 : end]),
                }
            )
            i = end
        i += 1
    return blocks


def main(argv):
    if len(argv) != 3:
        sys.exit("usage: python tools/extract_guide_blocks.py GUIDE.md OUT.json")
    guide, out = Path(argv[1]), Path(argv[2])
    blocks = extract_blocks(guide.read_text(encoding="utf-8").splitlines())
    out.write_text(json.dumps(blocks, indent=1) + "\n", encoding="utf-8")
    scripts = {block["script"] for block in blocks if block["script"]}
    unlabelled = sum(1 for block in blocks if not block["script"])
    print(
        f"{len(blocks)} Python blocks, {len(scripts)} scripts, {unlabelled} unlabelled"
    )


if __name__ == "__main__":
    main(sys.argv)
