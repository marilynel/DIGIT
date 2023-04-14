#!/bin/bash

# This script may be run as needed to remove sge directories out of the current
# working directory. Sge directories are created when batch jobs are submitted
# to the queuing system. The contents of these directories are useful in the
# event of DIGIT not running as expected. Most of the time they will not be
# needed by the user and may be safely deleted with this script.

# Usage: ./CleanUpDirectory.sh

rm -rf sge.*