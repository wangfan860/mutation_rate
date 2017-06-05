from itertools import islice
import csv



def fai_chunk(fai_path, blocksize):
    '''
    read reference index .fai file as input, and use blocksize (e.g. 1 million bases = 1000000)
    to generate interval chunk
    only use autosomes and X, Y, M
    '''
    seq_map = {}
    with open(fai_path) as handle:
        head = list(islice(handle, 25))
        for line in head:
            tmp = line.split("\t")
            seq_map[tmp[0]] = int(tmp[1])
        for seq in seq_map:
            l = seq_map[seq]
            for i in range(1, l, blocksize):
                yield (seq, i, min(i+blocksize-1, l))

def count_mutation(chr, start, end, path):
    '''
    count mutation from intervals
    '''
    reader=csv.reader(open(path,'rb'), delimiter='\t')
    counts=0
    for line in reader:
        if chr == line[0]:
            if int(start) <= int(line[1]) <= int(end):
                counts += 1
    return counts

def create_density_csv(fai_path, blocksize, input_file, output_file):
    '''
    calculate mutation density from intervals, and write to the output
    '''
    with open(output_file, 'w') as oh:
        for i, block in enumerate(fai_chunk(fai_path, blocksize)):
            mc = count_mutation(block[0], block[1], block[2], input_file)
            density = float(mc) / float(int(block[2]) + 1 - int(block[1]))
            c = '%s:%s-%s, %s, %s\n' % (block[0], block[1], block[2], mc, density)
            oh.writelines(c)

def gc_element(chr, start, end, path):
    '''
    read gc5base from ucsc data, and calculate gc contents from intervals
    '''
    pos = {}
    with open(path, 'rb') as handle:
        value = []
        for line in handle:
            if not line.startswith("variableStep"):
                tmp = line.rstrip('\n').split("\t")
                if int(start) <= int(tmp[0])& int(tmp[0]) +4 <= int(end):
                    value.append(float(tmp[1]))
        pos['%s:%s-%s' % (chr, start, end)] = float(sum(value))/20/float(end + 1 - start)
    return pos

def create_gc_csv(chr, fai_path, blocksize, input_file, output_file):
    '''
    calculate gc contents from intervals, and write to the output_file
    '''
    with open(output_file, 'w') as oh:
        for i, block in enumerate(fai_chunk(fai_path, blocksize)):
            if block[0] == chr:
                gc = gc_element(block[0], block[1], block[2], input_file)
                value = gc['%s:%s-%s' % (block[0], block[1], block[2])]
                c = '%s:%s-%s, %s\n' % (block[0], block[1], block[2], value)
                oh.writelines(c)
        oh.close()
