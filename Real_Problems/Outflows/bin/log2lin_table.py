## This script just converts Ralph's log cooling tables 
## to linear cooling tables to be used in PLUTO

import numpy as np
import re

class TableHeader():
    def __init__(self, fname, comment_symbols=['#',]):
        self.fname = fname
        self.comment_symbol = comment_symbols
        self.num_header_lines = 0
        self.num_data_lines = 0
        

    def parse_table_header(self):

        fh = open(self.fname)

        is_header = True
        num_header_lines = 0
        while is_header:
            buf = fh.readline()
            for cs in self.comment_symbols:
                is_header = re.search('^'+cs,buf) 
                if is_header: break

            num_header_lines += 1

        self.num_header_lines = num_header_lines
        self.num_data_lines = fh.readline()

        fh.close()

        return self.num_header_lines, self.num_data_lines


class SimpleTable(TableHeader):
    def __init__(self, fname):

        self.fname = fname

        self.data = []
        self.header_parsed = False

    def parse_table_header(self, comment_symbols=['#',]):

        self.header = TableHeader(self.fname, comment_symbols)
        self.header.parse_table_header()
        self.header_parsed = True

    def parse_table_data(self, fname):

        self.data = np.loadtxt(self.fname, skiprows=self.header.num_header_lines)

        return data


fname = 'cooltable_log.dat'
t = SimpleTable(fname)
t.parse_table_header(['%', ' * [^0-9]'])
data = t.parse_table_data(fname, ['%', ' * [^0-9]'])
np.savetxt('cooltable.dat', 10**data[:,0:2], '%12.6e', '   ')



