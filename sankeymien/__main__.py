# Copyright (C) 2026 Arnaud Belcour - Univ. Grenoble Alpes, Inria, Grenoble, France Microcosme
# Univ. Grenoble Alpes, Inria, Microcosme
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import argparse
import logging
import os
import sys
import time

from sankeymien import __version__ as VERSION
from sankeymien.sankey import handle_input
from sankeymien.utils import is_valid_dir


MESSAGE = f'''
Sankey diagram for Microbial Enrichment cultures. 
'''
REQUIRES = '''
Requires: plotly, pandas and numpy.
'''

logger = logging.getLogger('sankeymien')
logger.setLevel(logging.DEBUG)


def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        'sankeymien',
        description=MESSAGE + ' For specific help on each subcommand use: esmecata {cmd} --help',
        epilog=REQUIRES
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s ' + VERSION)

    parser.add_argument(
        '-a',
        '--abundance-file',
        dest='abundance_file',
        required=True,
        help='CSV/TSV file indicating the abundance of microorganisms in samples.',
        metavar='INPUT_FILE')

    parser.add_argument(
        '-j',
        '--json-file',
        dest='json_file',
        required=True,
        help='JSON file associating samples with time steps of the enrichment cultures.',
        metavar='INPUT_FILE')

    parser.add_argument(
        '-t',
        '--taxon',
        dest='taxon',
        required=True,
        help='Name of the column in abundance file containing taxon name.',
        metavar='STRING')

    parser.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        help='Output directory path.',
        metavar='OUPUT_DIR')

    parser.add_argument(
        "-t",
        dest='abundance_threshold',
        help="Abundance threshold to show taxon name above this threshold.",
        required=False,
        type=int,
        default=0)

    parent_parser_relative = argparse.ArgumentParser(add_help=False)
    parent_parser_relative.add_argument(
        '--relative',
        dest='relative',
        help='Compute relative abundance from abundance file.',
        required=False,
        action='store_true',
        default=None)

    args = parser.parse_args()

    # If no argument print the help.
    if len(sys.argv) == 1 or len(sys.argv) == 0:
        parser.print_help()
        sys.exit(1)

    is_valid_dir(args.output)

    # add logger in file
    formatter = logging.Formatter('%(message)s')
    log_file_path = os.path.join(args.output, f'sankeymien.log')
    file_handler = logging.FileHandler(log_file_path, 'w+')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    logger.addHandler(file_handler)

    # set up the default console logger
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    logger.info("--- Generate Sankey diagram for microbial enrichment cultures ---")
    handle_input(args.abundance_file, args.json_file, args.taxon, args.output, args.abundance_threshold, args.relative)

    duration = time.time() - start_time
    logger.info("--- Total runtime %.2f seconds ---" % (duration))
    """
    logger.warning(f'--- Logs written in {log_file_path} ---')
    """

if __name__ == '__main__':
    main()
