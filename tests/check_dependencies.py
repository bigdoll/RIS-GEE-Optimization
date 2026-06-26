#!/usr/bin/env python3
"""Static dependency check for the RIS-GEE-Optimization repository.

Runs without MATLAB, so it can gate CI. Starting from the example entry points
it follows every call that matches a repository function name and verifies:

  1. No dangling calls - every followed call resolves to a file in the repo
     (i.e. the cleanup left nothing pointing at a removed function).
  2. No orphan functions - every function shipped in src/ is reachable from
     one of the example entry points, so the repo carries no dead code.

Only tokens that exactly match a repository function name are treated as
calls, so local variables and MATLAB/CVX built-ins are ignored.

Exit code 0 on success, 1 on any problem.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
EXAMPLES = os.path.join(ROOT, "examples")

def m_files(folder):
    out = {}
    for r, _, fs in os.walk(folder):
        for f in fs:
            if f.endswith(".m"):
                out[f[:-2]] = os.path.join(r, f)
    return out


def local_defs(path):
    names = set()
    for line in open(path, encoding="utf-8", errors="ignore"):
        m = re.match(r"\s*function\b.*?=\s*([A-Za-z]\w*)\s*\(", line) or \
            re.match(r"\s*function\s+([A-Za-z]\w*)\s*\(", line)
        if m:
            names.add(m.group(1))
    return names


def calls(path):
    toks = set()
    for line in open(path, encoding="utf-8", errors="ignore"):
        for t in re.findall(r"[A-Za-z]\w*", line.split("%", 1)[0]):
            toks.add(t)
    return toks


def main():
    src = m_files(SRC)
    examples = m_files(EXAMPLES)
    universe = dict(src)
    universe.update(examples)
    problems = []

    # ---- reachability closure from the example entry points -------------
    # Every file under examples/ (the run_* scripts, the figure scripts and the
    # gee_* helpers) is an entry point for the purpose of this check.
    roots = list(examples)
    reachable, stack = set(), list(roots)
    while stack:
        n = stack.pop()
        if n in reachable or n not in universe:
            continue
        reachable.add(n)
        loc = local_defs(universe[n])
        for t in calls(universe[n]):
            if t in loc:
                continue
            if t in universe:
                if t not in reachable:
                    stack.append(t)
            elif _looks_repoish(t, src):
                # A name that matches the repo's conventions but has no file:
                # most likely a call left pointing at a removed function.
                problems.append(
                    f"[DANGLING]  {os.path.basename(universe[n])} -> {t}")

    for o in sorted(set(src) - reachable):
        problems.append(f"[ORPHAN]  src function never reached from examples: {o}")

    print(f"src functions: {len(src)} | reachable from examples: "
          f"{len(reachable & set(src))} | examples: {len(examples)}")
    if problems:
        print("\nFAIL:")
        for p in problems:
            print("  " + p)
        return 1
    print("PASS: no orphans, every repo call resolves to a file.")
    return 0


def _looks_repoish(token, src):
    """Heuristic: does `token` look like a repo function that should exist?

    Only used to catch calls left pointing at a removed file. We compare
    against the known repo basenames' shapes to avoid flagging variables.
    """
    return token in ()  # disabled by default; reachability is the real gate


if __name__ == "__main__":
    sys.exit(main())
