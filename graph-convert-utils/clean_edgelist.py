#!/usr/bin/python

#/*
# * @author Priyank Faldu <Priyank.Faldu@ed.ac.uk> <http://faldupriyank.com>
# *
# * Copyright 2019 The University of Edinburgh
# *
# * Licensed under the Apache License, Version 2.0 (the "License");
# * you may not use this file except in compliance with the License.
# * You may obtain a copy of the License at
# *
# *     http://www.apache.org/licenses/LICENSE-2.0
# *
# * Unless required by applicable law or agreed to in writing, software
# * distributed under the License is distributed on an "AS IS" BASIS,
# * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# * See the License for the specific language governing permissions and
# * limitations under the License.
# *
# */

import sys
import os

def main():
    #copy input file to output file as follows
    #remove any lines that is either empty or starts with comments like #,",@,%
    if len(sys.argv) != 3:
        print("Usage: python clean_edgelist.py input.el output.el")
        sys.exit()

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if not os.path.isfile(input_file):
        print("Input file {} does not exist. Exiting!".format(input_file))
        sys.exit()

    if os.path.isfile(output_file):
        print("Output file {} already exist. Please delete it before proceeding. Exiting!".format(output_file))
        sys.exit()

    with open(input_file, 'r') as inp_fp, open(output_file, 'w') as out_fp:
        line = inp_fp.readline()
        while line:
            line = line.strip()
            skip = False
            if line == "":
                skip = True
            else:
                for c in ('#', '"', '@', '%'):
                    if line[0] == c:
                        skip = True
                        break

            if not skip:
                out_fp.writelines(line+"\n")
            line = inp_fp.readline()

if __name__ == '__main__':
    main()
