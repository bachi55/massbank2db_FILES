####
#
# The MIT License (MIT)
#
# Copyright 2020 Eric Bach <eric.bach@aalto.fi>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is furnished
# to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
####
import os
import sys
import argparse
import pandas as pd

if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("ifn", type=str, help="Path to the MetFrag scoring output.")
    argparser.add_argument("--use_columns", type=str, nargs="+",
                           default=["Identifier", "InChIKey1", "Score", "FragmenterScore"],
                           help="Which columns to extract from the MetFrag output")
    args = argparser.parse_args()

    if not os.path.exists(args.ifn):
        sys.exit(1)

    # Load MetFrag output keeping only the score columns and write reduced table out
    pd.read_csv(args.ifn, usecols=args.use_columns).to_csv(args.ifn, index=False)

    sys.exit(0)
