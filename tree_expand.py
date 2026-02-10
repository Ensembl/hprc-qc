#!/usr/bin/env python3
import os, sys

root = sys.argv[1] if len(sys.argv) > 1 else "."
expand = os.path.abspath(sys.argv[2] if len(sys.argv) > 2 else "")

ignore_dirs = {".git", ".claude", "__pycache__"}

def tree(path, prefix="", is_last=True, force_expand=False):
    name = os.path.basename(path) if prefix else path
    connector = "" if not prefix else ("└── " if is_last else "├── ")

    print(f"{prefix}{connector}{name}")

    try:
        entries = sorted(os.listdir(path))
    except PermissionError:
        return

    # filter
    entries = [e for e in entries if e not in ignore_dirs and not e.startswith(".")]

    dirs = [e for e in entries if os.path.isdir(os.path.join(path, e))]
    files = [e for e in entries if os.path.isfile(os.path.join(path, e))]

    # decide whether to expand
    force_expand = force_expand or os.path.abspath(path) == expand

    # choose what to walk
    if not force_expand:
        children = dirs
        show_files = False
    else:
        children = dirs + files
        show_files = True

    for idx, child in enumerate(children):
        child_path = os.path.join(path, child)
        last_child = idx == len(children) - 1
        new_prefix = prefix + ("    " if is_last else "│   ")

        if os.path.isdir(child_path):
            tree(child_path, prefix=new_prefix, is_last=last_child, force_expand=force_expand)
        elif show_files:
            file_connector = "└── " if last_child else "├── "
            print(f"{new_prefix}{file_connector}{child}")

# start
tree(root, prefix="", is_last=True)

