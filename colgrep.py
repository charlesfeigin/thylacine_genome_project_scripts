#!/usr/bin/env python

# Copyright (c) 2016 Graham Gower <graham.gower@gmail.com>
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

import sys

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="grep a specific column of stdin, for strings in specified file")
    parser.add_argument("--file", "-f", type=str, required=True, help="list of strings to look for")
    parser.add_argument("--column", "-c", type=int, required=True, help="1-based column of input in which to search")
    parser.add_argument("--reverse", "-v", default=False, action="store_true", help="print lines that don't match")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    with open(args.file) as f:
        strings = {line.rstrip("\r\n") for line in f}

    for line in sys.stdin:
        line = line.rstrip("\r\n")
        fields = line.split()

        if len(fields) < args.column:
            continue

        if fields[args.column-1] in strings:
            if not args.reverse:
                print(line)
        else:
            if args.reverse:
                print(line)
