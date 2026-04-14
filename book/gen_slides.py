"""
gen_slides.py — generate revealjs slide .qmd files from book chapters.

For each ch*.qmd in the book directory:
  - strip its front matter (title:, format:, execute:, etc.)
  - write a slide-ready copy into _slides/ that declares only title
    (the _slides/_quarto.yml supplies the revealjs format for the whole project)

Also creates symlinks in _slides/ so notebooks can resolve relative paths
(%run _common.py, from python.L3 import *, image references) at execution time.

Run from the book directory:  python3 gen_slides.py
"""
import re, os, glob

BOOK_DIR   = os.path.dirname(os.path.abspath(__file__))
SLIDES_DIR = os.path.join(BOOK_DIR, "_slides")
os.makedirs(SLIDES_DIR, exist_ok=True)

# ── Symlinks for execution-time resources ─────────────────────────────────────
# Notebooks use %run _common.py and `from python.X import *`, so _slides/ must
# appear to contain these resources when Quarto runs notebooks from there.
SYMLINK_TARGETS = [
    "_common.py",   # shared setup run by every chapter
    "python",       # lecture module directory
    "images",       # images directory (if present)
]

for name in SYMLINK_TARGETS:
    src = os.path.join(BOOK_DIR, name)
    link = os.path.join(SLIDES_DIR, name)
    if os.path.exists(src) and not os.path.exists(link):
        os.symlink(os.path.join("..", name), link)
        print(f"  symlinked: _slides/{name} -> ../{name}")

# Symlink any loose image files referenced by chapters (e.g. displacements.png)
for img in glob.glob(os.path.join(BOOK_DIR, "*.png")) + \
           glob.glob(os.path.join(BOOK_DIR, "*.jpg")) + \
           glob.glob(os.path.join(BOOK_DIR, "*.svg")):
    name = os.path.basename(img)
    link = os.path.join(SLIDES_DIR, name)
    if not os.path.exists(link):
        os.symlink(os.path.join("..", name), link)
        print(f"  symlinked: _slides/{name} -> ../{name}")

# ── Generate slide .qmd files ─────────────────────────────────────────────────
FM_RE = re.compile(r'^---\s*\n.*?\n---\s*\n', re.DOTALL)

chapters = sorted(glob.glob(os.path.join(BOOK_DIR, "ch*.qmd")))

for src_path in chapters:
    fname = os.path.basename(src_path)
    with open(src_path) as f:
        content = f.read()

    # Extract title from front matter
    title_match = re.search(r'^title:\s*"([^"]+)"', content, re.MULTILINE)
    title = title_match.group(1) if title_match else fname.replace(".qmd", "")

    # Strip front matter — _slides/_quarto.yml owns all format/execute config
    body = FM_RE.sub("", content, count=1).lstrip("\n")

    slide_content = f'---\ntitle: "{title}"\n---\n\n{body}'

    dst_path = os.path.join(SLIDES_DIR, fname)
    with open(dst_path, "w") as f:
        f.write(slide_content)
    print(f"  generated: _slides/{fname}")

print(f"Done — {len(chapters)} slide files written to _slides/")
