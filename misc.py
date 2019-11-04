
import sys
import select
import multiprocessing as mp
import joblib

import subprocess
import numpy as np
import gzip
import pickle

from queue import Empty

cachedir = '.pycache'
memory = joblib.Memory(cachedir, verbose=0)
# Usage
# @memory.cache

def save_obj(name, obj):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def save_npy(name, npy):
    np.save(name, npy)
    return

def load_npy(name):
    npy = np.load(name + ".npy")
    return npy


def shell(cmd):
    """

    Yield the output of a command

    """

    popen = subprocess.Popen(cmd,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        shell=True)

    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line

    popen.stdout.close()
    return_code = popen.wait()

    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def stdin():
    """
    Reads lines from stdin pipe and yields results one-by-one.

    Yields:
        line: The next line from stdin

    Example:
        will yield a line for each txt file found in this folder
        find . -name "*.txt" | python file.py

        you can also read a file line by line
        cat filename.txt | python file.py

    """

    while sys.stdin in select.select([sys.stdin], [], [], 0)[0]:

        line = sys.stdin.readline()

        if not line:
            yield from []
            break

        line = line.strip()
        yield line


def parallel(lines, func, args, kwargs, procs=1):
    """

    Takes iterator (list or generator) `lines` and spawns # `procs` processes, calling
    `func` with prefined arguments `args` and `kwargs`.

    Using a queue and multiprocessing to call `func` with the format

    func(line, *args, **kwargs)

    Args:
        lines: iterator to be parallelised.
        func: function to call every line.
        args: Variable length argument list for `func`.
        kwargs: Arbitrary keyword arguments for `func`.
        procs: how many processes to start.

    Returns:
        results: List of all results from the parallel call (random order).

    """

    # Start a queue with the size of processes for jobs and a result queue to
    # collect results
    q_res = mp.Queue()
    q_job = mp.Queue(maxsize=procs)

    # print lock
    iolock = mp.Lock()

    # Start the pool and await queue data
    pool = mp.Pool(procs,
        initializer=process,
        initargs=(q_job, q_res, iolock, func, args, kwargs))

    # stream the data to queue
    for line in lines:

        # halts if queue is full
        q_job.put(line)

    # stop the process and pool
    for _ in range(procs):
        q_job.put(None)

    pool.close()

    # Collect all results
    while not q_res.empty():
        yield q_res.get(block=False)

    pool.join()


def process(q, results, iolock, func, args, kwargs):
    """

    multiprocessing interface for calling

    func(x, *args, **kwargs) with `x` coming from q

    args
        q: multiprocessing queue for distributing workload.
        results: multiprocessing queue for collecting results.
        iolock: print lock.
        func: function to be called with `q` output.
        kwargs: Arbitrary keyword arguments for `func`.

    """

    kwargs["iolock"] = iolock

    while True:

        try:
            line = q.get(timeout=0.05)
        except Empty:
            continue

        if line is None: break

        result = func(line, *args, **kwargs)
        results.put(result)

    return


def get_indexes(lines, pattern):

    idxs = []

    for i, line in enumerate(lines):
        if line.find(pattern) != -1:
            idxs.append(i)

    return idxs


def get_index(lines, pattern):
    for i, line in enumerate(lines):
        if line.find(pattern) != -1:
            return i
    return None


def reverse_enum(L):
    for index in reversed(range(len(L))):
        yield index, L[index]


def get_rev_indexes(lines, patterns):

    n_patterns = len(patterns)
    i_patterns = list(range(n_patterns))

    idxs = [None]*n_patterns

    for i, line in reverse_enum(lines):

        for ip in i_patterns:

            pattern = patterns[ip]

            if line.find(pattern) != -1:
                idxs[ip] = i
                i_patterns.remove(ip)

    return idxs


def get_rev_index(lines, pattern):

    for i, line in reverse_enum(lines):
        if line.find(pattern) != -1:
            return i

    return None

