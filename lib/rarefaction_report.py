#!/usr/bin/env python3
"""Build the standalone HTML report from lib/rarefaction_model.py output.

Reads OUTDIR/report_data.json and writes OUTDIR/rarefaction_report.html, a
single file with the data embedded and no external dependencies other than
web fonts, so it can be opened locally or shared as is.

Usage, from the project root:
    python lib/rarefaction_report.py [--outdir results/07_rarefaction]
"""
import argparse
import json
import os

TEMPLATE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rarefaction_report.template.html")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default="results/07_rarefaction")
    args = ap.parse_args()
    with open(os.path.join(args.outdir, "report_data.json")) as fh:
        data = fh.read()
    json.loads(data)  # fail here, not in the browser, if the file is damaged
    with open(TEMPLATE) as fh:
        html = fh.read()
    # "</" inside a <script> block would end it early.
    html = html.replace("__REPORT_DATA__", data.replace("</", "<\\/"))
    out = os.path.join(args.outdir, "rarefaction_report.html")
    with open(out, "w") as fh:
        fh.write(html)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
