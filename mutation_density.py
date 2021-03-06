from itertools import islice
import csv
import io
import os

def create_bed(input):
    output = os.path.basename(input).rsplit('.', 4)[0] + '.hg19.bed'
    g = open(output, 'w')
    with open(input, 'r') as fh:
        file = fh.read().splitlines()
        for i in file:
            line = i.split("\t")
            chr = "chr" + str(line[0])
            startpoint = str(line[1])
            endpoint = str(int(line[1]) + 1)
            g.write("%s\t%s\t%s\n" % (chr, startpoint, endpoint))

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
    read gcdata from bioconductor, and calculate gc contents from intervals
    '''
    pos = {}
    with io.open(path, 'rb') as handle:
        value = []
        for line in handle:
     	    tmp = line.rstrip('\n').split(",")
	    c = True
	    if tmp[0] == chr:
		if int(tmp[1]) >= start and int(tmp[2]) <= end:
			c = True
		else:
		    c = False
	    else:
		c = False
	    if c:
		num = (int(tmp[2]) + 1 - int(tmp[1])) * float(tmp[3])
	    	value.append(num)
    pos['%s:%s-%s' % (chr, start, end)] = float(sum(value))/float(end + 1 - start)
    return pos

def create_gc_csv(chr, fai_path, blocksize, input_file, output_file):
    '''
    calculate gc contents from intervals, and write to the output_file
    '''
    with open(output_file, 'w') as oh:
	intervals = []
        for i, block in enumerate(fai_chunk(fai_path, blocksize)):
            if block[0] == chr:
		intervals.append(block)
	for block in intervals:
       	    gc = gc_element(block[0], block[1], block[2], input_file)
            value = gc['%s:%s-%s' % (block[0], block[1], block[2])]
            c = '%s:%s-%s, %s\n' % (block[0], block[1], block[2], value)
            oh.writelines(c)

def feature_element(chr, start, end, path):
    pos = {}
    with open(path, 'rb') as handle:
        value = []
        for line in handle:
            tmp = line.rstrip('\n').split("\t")
            if tmp[0] == chr:
                big_start = start
                big_end = end
                feature_start = int(tmp[3])
                feature_end = int(tmp[4])
                if big_start < feature_start and big_end < feature_start:
                    value.append(0)
                elif big_start <= feature_start and feature_start <= big_end < feature_end:
                    value.append(big_end - feature_start + 1)
                elif big_start <= feature_start and big_end >= feature_end:
                    value.append(feature_end - feature_start + 1)
                elif feature_start < big_start < feature_end and big_start <  big_end < feature_end:
                    value.append(big_end - big_start + 1)
                elif feature_start < big_start <= feature_end and big_end >=feature_end:
                    value.append(feature_end - big_start + 1)
                elif big_start > feature_end and big_end > feature_end:
                    value.append(0)
    pos['%s:%s-%s' % (chr, start, end)] = sum(value)
    return pos

def create_feature_csv(chr, fai_path, blocksize, input_file, output_file):
    '''
    calculate gene feature from intervals, and write to the output_file
    '''
    with open(output_file, 'w') as oh:
        for i, block in enumerate(fai_chunk(fai_path, blocksize)):
            if block[0] == chr:
                feature = feature_element(block[0], block[1], block[2], input_file)
                value = feature['%s:%s-%s' % (block[0], block[1], block[2])]
                c = '%s:%s-%s, %s\n' % (block[0], block[1], block[2], value)
                oh.writelines(c)
