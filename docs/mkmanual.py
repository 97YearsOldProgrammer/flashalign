#!/usr/bin/env python3
"""Build the HTML manual from flashalign.1.

    python3 docs/mkmanual.py flashalign.1 > docs/index.html

What it does, in order:
  1. runs   mandoc -Thtml -O style=manual.css flashalign.1   (the plain conversion)
  2. reads the section headings mandoc produced (NAME, SYNOPSIS, ... and the
     "Indexing options" style sub-sections) and inserts a <nav class="toc">
     contents list in front of the manual text, so manual.css can show it as
     the sidebar. mandoc's own -O toc only works for mdoc(7) pages, and
     flashalign.1 is a man(7) page, which is why this small step exists.

Nothing else is changed: every word in the output comes from flashalign.1.
Needs Python 3 and mandoc. mandoc is preinstalled on macOS; on Debian/Ubuntu
`apt install mandoc`, on Fedora `dnf install mandoc`, Homebrew `brew install mandoc`.

You can also feed it HTML you already made:
    mandoc -Thtml -O style=manual.css flashalign.1 | python3 docs/mkmanual.py - > docs/index.html
"""
import html
import re
import shutil
import subprocess
import sys

# <h1 class="Sh" id="NAME">...</h1>   (mandoc 1.14.x)
# <h2 class="Sh" id="NAME">...</h2>   (newer mandoc)  -- the class is what matters
HEADING = re.compile(r'<h([1-6]) class="(Sh|Ss)" id="([^"]*)">(.*?)</h\1>', re.S)
CONTAINER = re.compile(r'<(?:div|main) class="manual-text">')
TAG = re.compile(r'<[^>]+>')


def heading_text(fragment):
    """'<a class="permalink" href="#x">Indexing\n  options</a>' -> 'Indexing options'"""
    return ' '.join(html.unescape(TAG.sub('', fragment)).split())


def toc_label(kind, title):
    # Section names are ALL CAPS in a man page; show them as 'Output format' in the list.
    if kind == 'Sh' and title == title.upper():
        return title[:1] + title[1:].lower()
    return title


def build_toc(page):
    sections = []  # [[id, label, [(id, label), ...]], ...]
    for m in HEADING.finditer(page):
        _, kind, hid, inner = m.groups()
        if not hid:
            continue
        label = toc_label(kind, heading_text(inner))
        if kind == 'Sh':
            sections.append([hid, label, []])
        elif sections:
            sections[-1][2].append((hid, label))
    if not sections:
        return ''

    def item(hid, label):
        return '<li><a href="#%s">%s</a>' % (html.escape(hid, quote=True), html.escape(label))

    out = ['<nav class="toc" aria-label="Contents">', '<div class="toc-title">Contents</div>', '<ul>']
    for hid, label, subs in sections:
        out.append(item(hid, label))
        if subs:
            out.append('<ul>')
            out.extend(item(h, l) + '</li>' for h, l in subs)
            out.append('</ul>')
        out.append('</li>')
    out += ['</ul>', '</nav>', '']
    return '\n'.join(out)


def add_toc(page):
    nav = build_toc(page)
    m = CONTAINER.search(page)
    if not nav or m is None:
        return page
    return page[:m.start()] + nav + page[m.start():]


def main(argv):
    if len(argv) != 2:
        sys.exit('usage: python3 docs/mkmanual.py flashalign.1 > docs/index.html')
    src = argv[1]
    if src == '-':
        page = sys.stdin.read()
    else:
        if shutil.which('mandoc') is None:
            sys.exit('mkmanual.py: mandoc not found (macOS has it; Linux: apt install mandoc)')
        run = subprocess.run(['mandoc', '-Thtml', '-O', 'style=manual.css', src],
                             capture_output=True, text=True)
        if run.returncode != 0 or not run.stdout:
            sys.exit('mkmanual.py: mandoc failed\n' + run.stderr)
        if run.stderr:
            sys.stderr.write(run.stderr)  # mandoc warnings, harmless
        page = run.stdout
    sys.stdout.write(add_toc(page))


if __name__ == '__main__':
    main(sys.argv)
