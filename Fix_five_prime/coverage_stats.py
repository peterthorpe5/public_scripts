#!/usr/bin/env python
"""BAM coverage statistics using samtools idxstats and depth.

This script expects exactly five command line arguments:
 * Input BAM filename
 * Input BAI filename (via Galaxy metadata)
 * Output tabular filename
 * Max coverage depth (integer)
 * Sliding window size (integer)

This messes about with the filenames to make samtools happy, then
runs "samtools idxstats" and "samtools depth", captures and combines
the output to the desired output tabular file.

Because "samtools depth" treats the max depth a little fuzzily, this
tool tries to account for this and applies a clear max-depth cut off.
"""
import sys
import os
import subprocess
import tempfile
from collections import deque

try:
    # New in Python 3.4
    from statistics import mean
except ImportError:
    def mean(list_of_values):
        """Calculate the mean average of a list of numbers."""
        # Quick and dirty, assumes already a list not an interator
        # so don't have to worry about getting the divisor.
        # Explicit float(...) to allow for Python 2 division.
        return sum(list_of_values) / float(len(list_of_values))

if "-v" in sys.argv or "--version" in sys.argv:
    #Galaxy seems to invert the order of the two lines
    print("BAM coverage statistics v0.1.0")
    cmd = "samtools 2>&1 | grep -i ^Version"
    sys.exit(os.system(cmd))

def sys_exit(msg, error_level=1):
   """Print error message to stdout and quit with given error level."""
   sys.stderr.write("%s\n" % msg)
   sys.exit(error_level)

# TODO - Proper command line API
if len(sys.argv) == 4:
    bam_filename, bai_filename, tabular_filename = sys.argv[1:]
    max_depth = "8000"
    window_size = "0"
elif len(sys.argv) == 5:
    bam_filename, bai_filename, tabular_filename, max_depth = sys.argv[1:]
    window_size = "0"
elif len(sys.argv) == 6:
    bam_filename, bai_filename, tabular_filename, max_depth, window_size = sys.argv[1:]
else:
    sys_exit("Require 3, 4 or 5 arguments: BAM, BAI, tabular filename, [max depth, [window_size]]")

if not os.path.isfile(bam_filename):
    sys_exit("Input BAM file not found: %s" % bam_filename)
if not os.path.isfile(bai_filename):
    if bai_filename == "None":
        sys_exit("Error: Galaxy did not index your BAM file")
    sys_exit("Input BAI file not found: %s" % bai_filename)

try:
    max_depth = int(max_depth)
except ValueError:
    sys_exit("Bad argument for max depth: %r" % max_depth)
if max_depth < 0:
    sys_exit("Bad argument for max depth: %r" % max_depth)

try:
    window_size = int(window_size)
except ValueError:
    sys_exit("Bad argument for window size: %r" % window_size)
if window_size and window_size < 0:
    sys_exit("Bad argument for window size (use 0 for no window, or suggest 5 or more): %r" % window_size)

# fuzz factor to ensure can reach max_depth, e.g. region with
# many reads having a deletion present could underestimate the
# coverage by capping the number of reads considered
depth_margin = 100

#Assign sensible names with real extensions, and setup symlinks:
tmp_dir = tempfile.mkdtemp()
bam_file = os.path.join(tmp_dir, "temp.bam")
bai_file = os.path.join(tmp_dir, "temp.bam.bai")
idxstats_filename = os.path.join(tmp_dir, "idxstats.tsv")
depth_filename = os.path.join(tmp_dir, "depth.tsv")
os.symlink(os.path.abspath(bam_filename), bam_file)
os.symlink(os.path.abspath(bai_filename), bai_file)
assert os.path.isfile(bam_file), bam_file
assert os.path.isfile(bai_file), bai_file
assert os.path.isfile(bam_file + ".bai"), bam_file

class NoCoverage(Exception):
    pass

def clean_up():
    os.remove(idxstats_filename)
    os.remove(depth_filename)
    os.remove(bam_file)
    os.remove(bai_file)
    os.rmdir(tmp_dir)


def samtools_depth_opt_available():
    child = subprocess.Popen(["samtools", "depth"],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # Combined stdout/stderr in case samtools is every inconsistent
    output, tmp = child.communicate()
    assert tmp is None
    # Expect to find this line in the help text, exact wording could change:
    #    -d/-m <int>         maximum coverage depth [8000]
    return " -d/-m " in output

depth_hack = False
if not samtools_depth_opt_available():
    if max_depth + depth_margin <= 8000:
        sys.stderr.write("WARNING: The version of samtools depth installed does not "
                         "support the -d option, however, the requested max-depth "
                         "is safely under the default of 8000.\n")
        depth_hack = True
    else:
        sys_exit("The version of samtools depth installed does not support the -d option.")

# Run samtools idxstats:
cmd = 'samtools idxstats "%s" > "%s"' % (bam_file, idxstats_filename)
return_code = os.system(cmd)
if return_code:
    clean_up()
    sys_exit("Return code %i from command:\n%s" % (return_code, cmd))

# Run samtools depth:
# TODO - Parse stdout instead?
if depth_hack:
    # Using an old samtools without the -d option, but hard coded default
    # of 8000 should be fine even allowing a margin for fuzzy output
    cmd = 'samtools depth "%s" > "%s"' % (bam_file, depth_filename)
else:
    cmd = 'samtools depth -d %i "%s" > "%s"' % (max_depth + depth_margin, bam_file, depth_filename)
return_code = os.system(cmd)
if return_code:
    clean_up()
    sys_exit("Return code %i from command:\n%s" % (return_code, cmd))

def get_stats_no_window(depth_iterator, length):
    """Calculate min/max/mean coverage.

    Returns tuple (int, int, float).
    """
    # First call to iterator may return NoCoverage
    # This also makes initialising the min value easy
    try:
        ref, pos, depth = next(depth_iterator)
    except NoCoverage:
        return 0, 0, 0.0
    except StopIteration:
        raise ValueError("Internal error - was there no coverage?")
    total_cov = min_cov = max_cov = depth

    # Could check pos is strictly increasing and within 1 to length?
    for ref, pos, depth in depth_iterator:
        total_cov += depth
        min_cov = min(min_cov, depth)
        max_cov = max(max_cov, depth)

    mean_cov = total_cov / float(length)

    assert min_cov <= mean_cov <= max_cov

    return min_cov, max_cov, mean_cov


def get_stats_window(depth_iterator, length, window_size):
    """Calculate min/max/mean and min/max windowed mean.

    Assumes the depth_iterator will fill in all the implicit zero
    entries which ``samtools depth`` may omit!

    Assumes window_size < number of values in iterator!
    """
    window = deque()
    total_cov = 0
    min_cov = None
    max_cov = 0.0

    assert 1 <= window_size <= length

    prev_pos = 0
    while len(window) < window_size:
        try:
            ref, pos, depth = next(depth_iterator)
        except NoCoverage:
            return 0, 0, 0.0, 0.0, 0.0
        except StopIteration:
            raise ValueError("Not enough depth values to fill %i window" % window_size)
        prev_pos += 1
        assert pos == prev_pos, "Discontinuity in coverage values for %s position %i" % (ref, pos)
        total_cov += depth
        if min_cov is None:
            min_cov = depth
        else:
            min_cov = min(min_cov, depth)
        max_cov = max(max_cov, depth)
        window.append(depth)

    assert len(window) == window_size
    min_win = max_win = mean(window)
    for ref, pos, depth in depth_iterator:
        prev_pos += 1
        assert pos == prev_pos, "Discontinuity in coverage values for %s position %i" % (ref, pos)
        total_cov += depth
        min_cov = min(min_cov, depth)
        max_cov = max(max_cov, depth)
        window.popleft()
        window.append(depth)
        assert len(window) == window_size
        win_depth = mean(window)
        min_win = min(min_win, win_depth)
        max_win = max(max_win, win_depth)

    mean_cov = total_cov / float(length)

    assert prev_pos == length, "Missing final coverage?"
    assert len(window) == window_size
    assert min_cov <= mean_cov <= max_cov
    assert min_cov <= min_win <= max_win <= max_cov

    return min_cov, max_cov, mean_cov, min_win, max_win

def depth_iterator(handle, identifier, length):
    """Iterates over ``samtools depth`` output for coverage depths.

    Uses global variables to cache the first line of output from the
    next reference sequence. Will impute missing zero depth lines
    which ``samtools depth`` omits.
    """
    global depth_ref, depth_pos, depth_reads

    assert identifier

    if depth_ref is None:
        # Right at start of file / new contig
        line = depth_handle.readline()
        if not line:
            # Must be at the end of the file.
            # This can happen if the file contig(s) had no reads mapped
            raise NoCoverage("End of file, no reads mapped to %r" % identifier)
        depth_ref, depth_pos, depth_reads = line.rstrip("\n").split()
        depth_pos = int(depth_pos)
        depth_reads = min(max_depth, int(depth_reads))
        # Can now treat as later references where first line cached
    elif identifier != depth_ref:
        # Infer that identifier had coverage zero,
        # and so was not in the 'samtools depth'
        # output.
        raise NoCoverage("%r appears to have no coverage at all" % identifier)

    # Good, at start of expected reference
    if depth_pos > 1:
        for extra_pos in range(1, depth_pos):
            yield identifier, extra_pos, 0
    yield identifier, depth_pos, depth_reads

    last_pos = depth_pos
    depth_ref = None
    depth_pos = 0
    depth_reads = 0
    for line in depth_handle:
        ref, pos, depth = line.rstrip("\n").split()
        pos = int(pos)
        depth = min(max_depth, int(depth))
        if ref != identifier:
            # Reached the end of this identifier's coverage
            # so cache this ready for next identifier
            depth_ref, depth_pos, depth_reads = ref, pos, depth
            break
        if last_pos + 1 < pos:
            for extra_pos in range(last_pos + 1, pos):
                yield identifier, extra_pos, 0
            # print("%s has no coverage between %i and %i" % (identifier, last_pos, pos))
        yield identifier, pos, depth
        last_pos = pos

    # Reached the end of this identifier's coverage or end of file
    if last_pos < length:
        for extra_pos in range(last_pos + 1, length +1):
            yield identifier, extra_pos, 0

    # Function ready to be called again to iterator over next ref


def load_total_coverage(depth_handle, identifier, length):
    """Parse some of the 'samtools depth' output for coverages.

    Returns min_cov (int), max_cov (int), mean cov (float),
    min_window_cov (float), max_window_cov (float).

    Uses global variables to cache the first line of output from the
    next reference sequence.

    Uses global variables for the max-depth and window size.
    """
    global depth_ref, depth_pos, depth_reads

    assert identifier
    assert 0 < length

    depth_iter = depth_iterator(depth_handle, identifier, length)

    if not window_size:
        min_cov, max_cov, mean_cov = get_stats_no_window(depth_iter, length)
        # Last two values are dummy return values:
        return min_cov, max_cov, mean_cov, min_cov, max_cov

    if length < window_size:
        depths = [depth for ref, pos, depth in depth_iter]
        mean_depth = mean(depths)
        sys.stderr.write("Warning: %s length %i < window size %i (using mean coverage for window values)\n"
                         % (identifier, length, window_size))
        return min(depths), max(depths), mean_depth, mean_depth, mean_depth
    else:
        return get_stats_window(depth_iter, length, window_size)


    mean_cov = bases / float(length)

    if length < window_size:
        sys.stderr.write("Warning: %s length %i < window size %i (using mean coverage for window values)\n"
                         % (identifier, length, window_size))
        min_win_cov = mean_cov
        max_win_cov = mean_cov
    return min_cov, max_cov, mean_cov, min_win_cov, max_win_cov

# Parse and combine the output
out_handle = open(tabular_filename, "w")
if window_size:
    # Write longer header with the two extra columns
    out_handle.write("#identifer\tlength\tmapped\tplaced\tmin_cov\tmax_cov\tmean_cov\tmin_win_cov\tmax_win_cov\n")
else:
    out_handle.write("#identifer\tlength\tmapped\tplaced\tmin_cov\tmax_cov\tmean_cov\n")

idxstats_handle = open(idxstats_filename)
depth_handle = open(depth_filename)

depth_ref = None
depth_pos = 0
depth_reads = 0
err = None
for line in idxstats_handle:
    identifier, length, mapped, placed = line.rstrip("\n").split()
    length = int(length)
    mapped = int(mapped)
    placed = int(placed)
    if identifier == "*":
        min_cov = 0
        max_cov = 0
        mean_cov = 0.0
        min_win_cov = 0.0
        max_win_cov = 0.0
    else:
        min_cov, max_cov, mean_cov, min_win_cov, max_win_cov \
            = load_total_coverage(depth_handle, identifier, length)
    if window_size:
        out_handle.write("%s\t%i\t%i\t%i\t%i\t%i\t%0.2f\t%0.2f\t%0.2f\n"
                         % (identifier, length, mapped, placed,
                            min_cov, max_cov, mean_cov, min_win_cov, max_win_cov))
    else:
        # Omit the unused last two columns
        assert min_win_cov == min_cov and max_win_cov == max_cov
        out_handle.write("%s\t%i\t%i\t%i\t%i\t%i\t%0.2f\n"
                         % (identifier, length, mapped, placed,
                            min_cov, max_cov, mean_cov))
    if max_cov > max_depth:
        err = ("Using max depth %i yet saw max coverage %i for %s"
              % (max_depth, max_cov, identifier))
    if not (min_cov <= mean_cov <= max_cov):
        err = ("Min/mean/max inconsistent: "
               "expect %r <= %r < %r "
               % (min_cov, mean_cov, max_cov))
    if not (min_cov <= min_win_cov <= max_win_cov <= max_cov):
        err = ("Windowed coverage inconsistent with min/max coverage: "
               "expect %r <= %r <= %r <= %r"
                % (min_cov, min_win_cov, max_win_cov, max_cov))
    if err:
        out_handle.write("ERROR during %s: %s\n" % (identifier, err))
        idxstats_handle.close()
        depth_handle.close()
        out_handle.close()
        clean_up()
        sys_exit("ERROR during %s: %s\n" % (identifier, err))

idxstats_handle.close()
depth_handle.close()
out_handle.close()

# Remove the temp symlinks and files:
clean_up()

if depth_ref is not None:
    sys_exit("Left over output from 'samtools depth'? %r" % depth_ref)
