####
#
# The MIT License (MIT)
#
# Copyright 2021 Eric Bach <eric.bach@aalto.fi>
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
import sqlite3


def create_tables(conn: sqlite3.Connection):
    conn.execute("CREATE TABLE molecules("
                 "  pubchem_id          INTEGER PRIMARY KEY, "
                 "  inchi               VARCHAR NOT NULL,"
                 "  inchi2D             VARCHAR NOT NULL,"
                 "  inchikey            VARCHAR NOT NULL,"
                 "  inchikey1           VARCHAR NOT NULL,"
                 "  inchikey2           VARCHAR NOT NULL,"
                 "  inchikey3           VARCHAR NOT NULL,"
                 "  smiles_canonical    VARCHAR NOT NULL,"
                 "  smiles_isomeric     VARCHAR NOT NULL,"
                 "  monoisotopic_mass   REAL NOT NULL,"
                 "  molecular_weight    REAL NOT NULL,"
                 "  exact_mass          REAL NOT NULL)")

    conn.execute("CREATE TABLE datasets("
                 "  name        VARCHAR PRIMARY KEY "
                 ")")

