#! /usr/bin/env python
import sys
import argparse
import tomllib
import pprint
import json


def main():
    p = argparse.ArgumentParser()
    p.add_argument("toml_file")
    p.add_argument("-o", "--output-json", required=True)
    args = p.parse_args()

    with open(args.toml_file, "rb") as fp:
        d = tomllib.load(fp)

    # pprint.pprint(d)
    out_d = {}

    creators = []
    for author_d in d["project"]["authors"]:
        creators.append(author_d)

    out_d["creators"] = creators

    with open(args.output_json, "w") as fp:
        json.dump(out_d, fp)


if __name__ == "__main__":
    sys.exit(main())
