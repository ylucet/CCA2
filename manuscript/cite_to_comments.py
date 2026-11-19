#!/usr/bin/env python3
"""Convert inline [CITE: ...] annotations in a .tex file into proper LaTeX
% comments, keeping any live \\cite{...} commands active.

Usage:
    python3 cite_to_comments.py main.tex            # in place
    python3 cite_to_comments.py main.tex --dry      # show a sample, no write

Rules:
  * A CITE block starts at '[CITE:' and ends at its matching ']' (bracket
    depth is tracked, so a block may contain a full nested CITE block).
  * Any \\cite[tp]{...} inside the block is extracted and kept as live text
    at the block's position; the rest becomes comment lines.
  * If the source has a sentence period immediately after the closing ']'
    (the annotation ended the sentence), that period is moved to the end of
    the sentence text where the block stood (no doubled period if one is
    already there).
  * Everything else around the block is preserved verbatim.
"""
import re
import sys


def find_blocks(src):
    blocks = []
    i = 0
    while True:
        j = src.find('[CITE:', i)
        if j == -1:
            return blocks
        depth, k = 0, j
        while k < len(src):
            if src[k] == '[':
                depth += 1
            elif src[k] == ']':
                depth -= 1
                if depth == 0:
                    break
            k += 1
        if depth != 0:
            raise SystemExit(f'unbalanced CITE block at char {j}')
        blocks.append((j, k))
        i = k + 1


def _append_to_last_line(prefix_str, addition):
    """Append `addition` to the last non-empty line of prefix_str."""
    lines = prefix_str.split('\n')
    i = len(lines) - 1
    while i >= 0 and not lines[i].strip():
        i -= 1
    if i < 0 or addition == '':
        return prefix_str
    add = addition
    if add == '.' and lines[i].rstrip().endswith('.'):
        add = ''                     # don't double the period
    if add and not lines[i].rstrip().endswith(add[-1:]) and add != '.':
        add = ' ' + add
    lines[i] = lines[i].rstrip() + add
    return '\n'.join(lines)


def convert(src):
    blocks = find_blocks(src)
    out = []
    prev_end = 0
    for (s, e) in blocks:
        inner = src[s + 1:e]
        cits = re.findall(r'\\cite[tp]?\{[^}]*\}', inner)
        inner_clean = re.sub(r'\\cite[tp]?\{[^}]*\}', '', inner)
        tail = src[e + 1:]
        m = re.match(r'\s*\.', tail)          # sentence-final period?
        consumed = m.end() if m else 0
        # the region this block lives in ends at the end of the line holding ']'
        nl = src.find('\n', e + 1 + consumed)
        region_end = len(src) if nl == -1 else nl + 1
        suffix = src[e + 1 + consumed:nl if nl != -1 else len(src)]
        line_start = src.rfind('\n', 0, s) + 1
        indent = re.match(r'\s*', src[line_start:s]).group(0)
        lines = [ln.rstrip() for ln in inner_clean.split('\n')]
        while lines and not lines[-1].strip():
            lines.pop()
        while lines and not lines[0].strip():
            lines.pop(0)
        comment = ''.join(indent + ('% ' + ln if ln.strip() else '%') + '\n'
                          for ln in lines)
        prefix = src[line_start:s].rstrip()
        live = ' '.join(c.strip() for c in cits)
        dot = '.' if m else ''
        prefix_str = src[prev_end:line_start]
        if prefix:
            # inline: live text, citation and period stay on the prefix line
            new_prefix = prefix
            if live:
                new_prefix = prefix.rstrip() + ' ' + live.rstrip()
            if dot and not new_prefix.rstrip().endswith('.'):
                new_prefix = new_prefix.rstrip() + '.'
            out.append(prefix_str)
            out.append(new_prefix + '\n')
        else:
            # block starts its own line: append citation/period to the
            # last text line before it
            addition = live.rstrip() + dot
            out.append(_append_to_last_line(prefix_str, addition))
        out.append(comment)
        # the region covers the newline ending the ']' line; if text follows
        # ']' on that line, keep it (with its newline) -- otherwise the
        # newline is implicit at the end of the comment block
        if suffix.strip():
            out.append(suffix.rstrip() + '\n')
        prev_end = region_end
    out.append(src[prev_end:])
    return ''.join(out), len(blocks)


def main():
    args = sys.argv[1:]
    dry = '--dry' in args
    args = [a for a in args if a != '--dry']
    path = args[0] if args else 'main.tex'
    src = open(path).read()
    new, n = convert(src)
    if dry:
        import difflib
        diff = list(difflib.unified_diff(src.splitlines(), new.splitlines(),
                                         lineterm='', n=1))
        print(f'{n} blocks; diff (first 80 lines):')
        print('\n'.join(diff[:80]))
        return
    open(path, 'w').write(new)
    print(f'{n} blocks converted in {path}')


if __name__ == '__main__':
    main()
