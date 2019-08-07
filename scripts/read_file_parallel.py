#!/usr/bin/env python

def readlines_chunk(fname, chunk, total_chunks):
    max_line_size = 200
    print(chunk)
    fh = open(fname)
    fsize = fh.seek(0, 2)

    beg = int(fsize / total_chunks) * chunk
    if beg:
        for pos in range(max_line_size):
            fh.seek(beg + pos)
            if fh.read(1) == '\n':
                beg += pos + 1
                break

    end = int(fsize / total_chunks) * (chunk + 1)
    for pos in range(max_line_size):
        fh.seek(end + pos)
        if fh.read(1) == '\n':
            end += pos + 1
            break
    if beg == end:
        raise Exception('ERROR: more cpus then lines... perhaps no need to'
                        ' parallelize...')
    fh.seek(beg)

    read = beg
    for line in fh:
        read += len(line)
        _ = line
        if read >= end:
            break

def main():
    print('hola')
    import sys
    from multiprocessing import Pool

    fname, ncpus = sys.argv[1:]
    ncpus = int(ncpus)

    pool = Pool(ncpus)
    procs = []
    for p in range(ncpus):
        procs.append(pool.apply_async(readlines_chunk, args=(fname, p, ncpus)))
    pool.close()
    pool.join()

if __name__ == '__main__':
    exit(main())
