# python librarys
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings('ignore',category=UserWarning)
    import pandas as pd

import numpy as np
import os
import re
import itertools as it
from collections import Counter
from subprocess import Popen, PIPE, STDOUT
from multiprocessing import Pool
# from __future__ import print_function
import argparse


def pp(cmd):
    return Popen(cmd,
                 shell=True,
                 stdout=PIPE,
                 stderr=STDOUT).communicate()


def select_hla_reads(sample, outdir, m, cks, function):
    """
    Select reads mapped to hla in .map file
    Example of .map format of one read:
    (http://algorithms.cnag.cat/wiki/FAQ:The_GEM_alignment_format)

    Parameters
    ----------
    sample: string with sample name
    outdir: output directory
    m: number of mismatch

    Return
    ------
    dbHLA: 

    """
    def find_gen_matches(data, pattern, m):

        if function == 'type':
            def f(x):
                num_matches = np.array(x[3].split(':'))[:(m + 1)].astype(int).sum()
                ms = ','.join(x[4].split(',')[0:num_matches])
                return pd.Series([num_matches, ms])
        else:
            def f(x):
                n_matches = np.array(x[3].split(':')).astype(int)
                num_matches = n_matches[np.nonzero(n_matches)][0]
                ms = ','.join(x[4].split(',')[0:num_matches])
                return pd.Series([num_matches, ms])

        return data[~data[3].str.contains(r'^0$')
                    & data[4].str.contains(pattern)
                    ].apply(f, axis=1)

    if function == 'type':
        pp('mkdir -p {}/hla'.format(outdir))
        pp('grep HLA {1}/map/{0}.map | grep -v HWI > {1}/hla/{0}.hmap'.format(sample, outdir))
        outfile = "{1}/hla/{0}.hmap".format(sample, outdir)
    elif function == 'express':
        pp('mkdir -p {}/xhla'.format(outdir))
        pp('grep HLA {1}/xmap/{0}.map | grep -v HWI > {1}/xhla/{0}.hmap'.format(sample, outdir))
        outfile = "{1}/xhla/{0}.hmap".format(sample, outdir)

    dmap = pd.read_table(outfile,
                         header=None,
                         chunksize=cks,
                         dtype=str,
                         na_filter=False,
                         # error_bad_lines = False,
                         delimiter='\t')

    hla_pattern = r'(HLA_[[A-Z]+1?_[0-9]{2}).+\:\:\1'

    dbHLA = pd.concat([find_gen_matches(chunk, hla_pattern, m) for chunk in dmap],
                     ignore_index=True)


    return dbHLA


def eval_match_matrix(data, sample, outdir, function):
    """
    Parameters
    ----------
    data: pd.DataFrame(index = [reads],
                        columns = [exact_mactches,locus_names])
    firt column is the number of exact matches of one read
    and second column is a string with the locus where the
    reads were mapped


    Return
    ------
    rgm: pd.DataFrame(index = [reads], columns=[2dig_HLA_names])
    the number of times on read is
        mapped into one 2 digits HLA gene
        (reads gene matches)
    """
    pattern = r'HLA_\w+_[0-9]+_[0-9]+'
    mm = data.apply(lambda x: pd.Series(re.findall(pattern, x[1]))
                     .value_counts(), axis=1)

    if "HLA_" in mm.columns[0]:
         mm.rename(columns=lambda x: x[4:], inplace=True)

    if function == 'type':
        pp('mkdir -p {}/gencounts'.format(outdir))
        mm.to_csv("{0}/gencounts/{1}.cts".format(outdir,sample))
    elif function == 'express':
        pp('mkdir -p {}/xgencounts'.format(outdir))
        mm.to_csv("{0}/xgencounts/{1}.cts".format(outdir,sample))

    return mm


def mapper(sample, outdir, m, n_threads, verbosity, function):

    if function == 'express':
        pp('mkdir -p {}/xmap'.format(outdir))
        cmd = ("gem-mapper -I {0}/index/{1}.gem ".format(outdir, sample)
               + "-p -1 {0}_1.fastq -2 {0}_2.fastq ".format(sample)
               + "-q offset-33 -d all "
               + "-o {0}/xmap/{1} -T {2}".format(outdir, sample, n_threads)
               )
    elif function == 'type': 
        pp('mkdir -p {}/map'.format(outdir))
        cmd = ("gem-mapper -I ./index/index_short.gem "
               + "-p -1 {0}_1.fastq -2 {0}_2.fastq ".format(sample)
               + "-q offset-33 -d all "
               + "-o {0}/map/{1} -m {2} -T {3}".format(outdir, sample, m, n_threads)
               )

    if verbosity:
        print cmd
    pp(cmd)
    pass


def querys_full(gtype):

    def qlist(gt, (i, j)):
        return list(map(lambda x: '&'.join(x),
                        it.product(*[['{}({} > 0)'.format(op, z)
                                      for op in ['', '~']] for z in gt[i:j]])))

    def labels(gt, (i, j)):
        return ['_'.join([k.split(' > ')[0][1:]
                          for k in q.split('&') if not k.startswith('~')])
                for q in qlist(gt.index, (i, j))]

    lC1, C1 = labels(gtype, (1, 7)), qlist(gtype, (1, 7))
    lC2, C2 = labels(gtype, (7, len(gtype))), qlist(gtype, (7, len(gtype)))

    aux = {i: j for i, j in zip(lC1 + lC2, C1 + C2)}
    aux.pop('')
    aux['nowhere'] = (C1[-1] + '&' + C2[-1])

    nlist = aux['nowhere'].split('&')
    for i in xrange(len(nlist) / 2):
        aux['no' + nlist[2 * i].split('_')[0][2:]] = (nlist[2 * i]
                                                      + '&'
                                                      + nlist[2 * i + 1])

    return aux


def querys(gtype):

    def qlist(gt, (i, j)):
        return list(map(lambda x: '&'.join(x),
                        it.product(*[['{}({} > 0)'.format(op, z)
                                      for op in ['', '~']] for z in gt[i:j]])))

    def labels(gt, (i, j)):
        return ['_'.join([k.split(' > ')[0][1:]
                          for k in q.split('&') if not k.startswith('~')])
                for q in qlist(gt.index, (i, j))]

    lC1, C1 = labels(gtype, (0, len(gtype))), qlist(gtype, (0, len(gtype)))

    aux = {i: j for i, j in zip(lC1, C1)}
    aux.pop('')
    aux['nowhere'] = (C1[-1])

    nlist = aux['nowhere'].split('&')
    for i in xrange(len(nlist) / 2):
        aux['no' + nlist[2 * i].split('_')[0][2:]] = (nlist[2 * i]
                                                      + '&'
                                                      + nlist[2 * i + 1])

    return aux


def EM_1(match_matrix, q):

    counts = {label: len(match_matrix.query(cond)) for label, cond in q.iteritems()}
    # correcting homozigozity
    hz = []
    for l in q:
        gs = l.split("_")
        if (len(gs) == 2 and gs[0][:-1] == gs[1][:-1] and
                counts[gs[0]] == 0 and
                counts[gs[1]] == 0):
            counts[gs[0]] = counts[gs[1]] = counts[l] / 2.0
            counts[l] = 0
            hz.append((l, gs[0], gs[1]))

    x = pd.Series(counts, dtype=float)
    x = x[x >= 0]
    ix = [i for i in x.index if 'no' not in i]
    x = x[ix].sort_index()
    ix = [i for i in x.index if u'_' not in i]
    tc = x[ix]
    tc0 = tc.copy() + 1.0
    tc1 = tc.copy()
    ix = [i for i in x.index if u'_' in i]
    pc = x[ix]

    pw = pd.Series(index=pc.index)
    conv = True
    while conv:
        for q in pc.index:
            pw[q] = 1.0 / tc0[q.split('_')].sum()

        spw = pw * pc
        for g in tc1.index:
            ig = pc.index.map(lambda z: g in z.split('_'))
            tc1[g] = tc[g] + (tc0[g] * spw[ig]).sum()

        if (abs(tc1 - tc0) > 1e-4).all():
            tc0 = tc1.copy()
        else:
            break

    tcounts = dict(tc1)
    # print tc2 - tc1

    if hz == []:
        return tcounts, counts
    else:
        for c, g1, g2 in hz:
            counts[c] = int(counts[g1] + counts[g2])
            counts[g1] = counts[g2] = 0
        return tcounts, counts


def get_compact_model(hit_df, to_bool=False):
    """
    turn a hit matrix dataframe (can be mapping position matrix if used with to_bool=True)
    into a smaller matrix DF that removes duplicate rows and creates the "occurence" vector
    with the number of rows the representative read represents
    """

    hit_df = hit_df.ix[hit_df.any(axis=1)]  # remove all-zero rows
    if to_bool:
        hit_df = hit_df.applymap(bool)
    occurence = {r[0]: len(r) for r in hit_df.groupby(hit_df.columns.tolist()).groups.itervalues()}
    unique_mtx = hit_df.drop_duplicates()
    return unique_mtx, occurence


def EM_2(match_matrix, q, gtype, model_id):

    counts = {label: len(match_matrix.query(cond)) for label, cond in q.iteritems()}
    # correcting homozigozity
    hz = []
    for l in q:
        gs = l.split("_")
        if (len(gs) == 2 and gs[0][:-1] == gs[1][:-1] and
                counts[gs[0]] == 0 and
                counts[gs[1]] == 0):
            counts[gs[0]] = counts[gs[1]] = counts[l] / 2.0
            counts[l] = 0
            hz.append((l, gs[0], gs[1]))

    x = pd.Series(counts, dtype=float)
    x = x[x >= 0]
    ix = [i for i in x.index if 'no' not in i]
    x = x[ix].sort_index()
    ix = [i for i in x.index if u'_' not in i]
    tc = x[ix]
    tc0 = tc.copy() + 1.0
    tc1 = tc.copy()
    ix = [i for i in x.index if u'_' in i]
    pc = x[ix]

    pw = pd.Series(index=pc.index)
    conv = True
    N = tc.sum()
    # print 'tc\n', tc, 'pc\n', pc.head()
    while conv:
        for q in pc.index:
            pw[q] = 1.0 / tc0[q.split('_')].sum()

        spw = pw * pc
        for g in tc1.index:
            ig = pc.index.map(lambda z: g in z.split('_'))
            tc1[g] = tc[g] + (tc0[g] * spw[ig]).sum()

        if (abs(tc1 - tc0) / N > 1e-3).all():
            tc0 = tc1.copy()
        else:
            break

    tc1.index = [ 'c_' + i for i in tc1.index]
    return model_id, tc1.append(gtype)


def EM_3(match_matrix, q, gtype, model_id):

    counts = {label: len(match_matrix.query(cond)) for label, cond in q.iteritems()}
    # correcting homozigozity
    hz = []
    for l in q:
        gs = l.split("_")
        if (len(gs) == 2 and gs[0][:-1] == gs[1][:-1] and
                counts[gs[0]] == 0 and
                counts[gs[1]] == 0):
            counts[gs[0]] = counts[gs[1]] = counts[l] / 2.0
            counts[l] = 0
            hz.append((l, gs[0], gs[1]))

    x = pd.Series(counts, dtype=float)
    x = x[x >= 0]
    ix = [i for i in x.index if 'no' not in i]
    x = x[ix].sort_index()
    ix = [i for i in x.index if u'_' not in i]
    tc = x[ix]

    idt = tc.index
    N = tc.sum()
    tc0 =  pd.Series(N * np.ones(len(idt))/len(idt),index=idt)
    tc1 =  pd.Series(np.ones(len(idt))/len(idt),index=idt)

    # tc0 = tc.copy() + 1.0
    # tc1 = tc.copy()
    ix = [i for i in x.index if u'_' in i]
    pc = x[ix]

    pw = pd.Series(index=pc.index)
    conv = True
    # print 'tc\n', tc, 'pc\n', pc.head()
    while conv:
        for q in pc.index:
            pw[q] = 1.0 / tc0[q.split('_')].sum()

        spw = pw * pc
        for g in tc1.index:
            ig = pc.index.map(lambda z: g in z.split('_'))
            tc1[g] = tc[g] + (tc0[g] * spw[ig]).sum()

        if (abs(tc1 - tc0) / N > 1e-3).all():
            tc0 = tc1.copy()
        else:
            break

    tc1.index = [ 'c_' + i for i in tc1.index]
    return model_id, tc1.append(gtype)


def get_top_genes(match_matrix, n):
    """
    match_matrix is a integer or boolean matrix with no not a number and no 
    duplicated gene 

    Parameters
    ----------
    match_matrix: pd.DataFrame(cols=gen_name, index=read)
    n: number of top genes to select for each HLA

    Return
    ------
    top_matches: dict((A, B, ..., DRB1): top_hlas)
    """

    top_matches = dict()
    hlas = {'abc':('A', 'B', 'C'), 'db':('DQA1', 'DQB1', 'DRB1')}
    for hla in hlas.values():
        first = True
        gs = []
        for gene in hla:
            top_matches[gene] = []

        for _ in xrange(n):
            if first:
                S = match_matrix.sum()
            else:
                idx = (match_matrix.ix[:,gs] > 0).sum(axis=1).apply(bool)
                S = match_matrix.ix[~idx,:].sum()
                #print idx.head()
            for gene in hla:
                G = [i for i in S.index if i.startswith(gene)]
                x = S[G].argmax()
                gs.append(x)
                top_matches[gene].append(x)
            first = False
        
    
    return pd.Series(top_matches).sort_index()


def remove_duplicates(match_matrix, remove_reads=False):
    """
    return a match matrix with no duplicated genes

    Parameters
    ----------
    match_matrix: any match matrix
    remove_reads=bool default False, if true return matrix with no duplicated
    genes and reads
    
    Return
    ----------
    new_match_matrix: pd.DataFrame
    """
    mm = match_matrix > 0
    mm.dropna(how='all',inplace=True)
    mg = mm.T.drop_duplicates().T
    
    if not remove_reads:
        return mg
    else:
        return mg.drop_duplicates()


def create_models(top_matches):
    """
    Parameters
    ----------
    top_matches: pd.Series(index=A,B,C.., values=genes with most matches)

    """

    hlas = {'abc':['A', 'B', 'C'], 'db':['DQA1','DQB1', 'DRB1']}
    cols = {'abc':['A1', 'A2',
                   'B1', 'B2',
                   'C1', 'C2'], 
            'db':['DQA11', 'DQA12', 'DQB11',  'DQB12', 'DRB11','DRB12']}
    qqa = {}
    out = {}
    for h in hlas:
        qa = ([list(j) for j in it.combinations_with_replacement(i, 2)][:-1]
                for i in top_matches[hlas[h]].values)

        qqa[h] = np.array([list(i) for i in it.product(*qa)])
        out[h] = pd.DataFrame(qqa[h].reshape(qqa[h].shape[0],len(cols[h])),columns=cols[h])
    return out


def run_EM(arg):
    mm, q, md, m_id = arg
    return EM_2(mm, q, md, m_id)


def get_best_model(sample, models, match_matrix, n_threads=1):

    eps = {'abc':0.02, 'db':0.04}
    md = {}
    result = {}
    for h in models.keys():
        output = {}
        aux = {}
        for i in xrange(len(models[h])):
            aux[i] = match_matrix[models[h].ix[i]].any(axis=1).sum()
        aux = pd.Series(aux).sort(ascending=False, inplace=False)
        for i in aux.index:
            j, s = EM_2(match_matrix, 
                        querys(models[h].ix[i,:]), 
                        models[h].ix[i,:], i)

            n = models[h].shape[1]
            if points(s[:n], eps[h]) != -np.inf:
                md[h] = models[h].ix[j]
                output[j] = s
                break
        result[h] = output

    out = {}
    for h in models.keys():
        out[h] = pd.DataFrame(result[h]).T
        out[h].index = [sample + '.' + str(i) for i in xrange(out[h].shape[0])]

    
    r = pd.concat(out.values(),axis = 1) 
    c = sorted([ i for i in r.columns if i.startswith('c_')])
    g = sorted([ i for i in r.columns if not i.startswith('c_')])

    return r[g + c], md


def find_duplicates(rank, match_matrix):
    genes = sorted([ i for i in rank.columns if not i.startswith('c')])
    mc = match_matrix.columns 
    x = {g: match_matrix.eq(match_matrix[rank[g]], axis=0).all() for g in genes}
    for g in genes:
        rank[g] = '/'.join(x[g][x[g]].index)

    return rank


def points(x, eps):
    y = []
    for i in xrange(len(x) / 2):
        y.append(min(x[[2* i,2 * i + 1]]) / max(x[[2* i,2 * i + 1]]))
    y = np.array(y)
    if (y > eps).all():
        return x.sum()
    else:
        return - np.inf


def read_full_hla_fasta():

    fasta = open('index/index_all_exons.fasta')
    for  hearder, group in it.groupby(fasta,lambda x: x.startswith('>')):
        if hearder:
            line = group.next()
            gene_id = line[5:-1]
            # gene_id = ('_').join(line[5:-1].split('_')[:4])
        else:
            sequence = ''.join(line.strip() for line in group)
            yield gene_id, sequence


def make_sample_index(rank, sample, outdir, n_threads, verbosity):

    g = sorted([ i for i in rank.columns if not i.startswith('c_')])
    genes = np.array([ i.split('/') for i in rank[g].values[0]]).reshape(-1)
    HLA = dict(read_full_hla_fasta())
    pp('mkdir {}/index'.format(outdir))
    fname = '{0}/index/{1}'.format(outdir,sample)
    fout = open(fname +'.fasta' ,'w')
    for gg in genes:
        fout.write('>HLA_' + gg + '\n')
        fout.write(HLA[gg]+ '\n')
    fout.close()
    cmd = 'gem-indexer -i {0}.fasta -o {0} -T {1}'.format(fname,n_threads)
    if verbosity:
        print cmd
    pp(cmd)  


def eval_expression(sample, rank, md, match_matrix_x):
    
    s = {}
    for h in md:
        j, s[h] = EM_3(match_matrix_x, querys(md[h]), md[h], 0)

    
    r = pd.DataFrame(pd.concat(s.values())).T
    r.index = rank.index 
    c = sorted([ i for i in r.columns if i.startswith('c_')])
    g = sorted([ i for i in r.columns if not i.startswith('c_')])
    ca = ['ca_' + i for i in g]
    r = r[g + c]
    r.columns = g + ca
    rr = pd.concat((rank,r.ix[:,ca]),axis=1)
    return rr

    
def genotype_pipeline(sample, outdir, m, n_threads, print_header, verbosity):

    if verbosity:
        print 'genotyping:', sample
    mapper(sample, outdir, m, n_threads, verbosity, function = 'type')

    if verbosity:
        print 'selecting hla reads'
    dbHLA = select_hla_reads(sample, outdir, m, cks = 800, function = 'type')

    if verbosity:
        print 'evaluating score matrix'
    match_matrix = (eval_match_matrix(dbHLA, sample, outdir, function = 'type') > 0)

    if verbosity:
        print 'finding best model'
    top_matches = get_top_genes(match_matrix, 3)

    top_models = create_models(top_matches)

    rank, md = get_best_model(sample, top_models, match_matrix, n_threads=n_threads)    

    rankd = find_duplicates(rank, match_matrix)

    return rank, md, match_matrix


def express_pipeline(rank, md, sample, outdir, m, n_threads, print_header, verbosity):

    if verbosity:
        print 'estimating expression:', sample

    make_sample_index(rank, sample, outdir, n_threads, verbosity)

    mapper(sample, outdir, m, n_threads, verbosity, function='express')

    if verbosity:
        print 'selecting hla reads'
    dbHLA_x = select_hla_reads(sample, outdir, m = 4, cks = 2000, function = 'express')

    match_matrix_x = (eval_match_matrix(dbHLA_x, sample, outdir, function = 'express') > 0)

    rank_x = eval_expression(sample, rank, md, match_matrix_x)

    return rank_x


def pipeline_hlatx(sample, outdir, m, n_threads, print_header, verbosity):

    rank, md, match_matrix= genotype_pipeline(sample, outdir, m, n_threads, print_header, verbosity)
    rank_x = express_pipeline(rank, md, sample, outdir, m, n_threads, print_header, verbosity)
    #rank_x = find_duplicates(rank_x, match_matrix)
    rank_x.to_csv('{0}/{0}_TX.txt'.format(outdir),
                 header=print_header,
                 mode='a')

    pass



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=' :HLAtx 4-digit typer and expression', 
                            prog='HLAtx')

    parser.add_argument('--input','-i',
                      nargs='+',
                      required=True,
                      help="Name of fastq files with fished HLA reads."
                      )

    parser.add_argument('--outdir','-o',
                      required=True,
                      help="Specifies the out directory to which all files should be written"
                      )

    parser.add_argument('--verbose','-v',
                      required=False,
                      action="store_true",
                      help="Set verbose mode on."
                      )

    parser.add_argument('--nthreads','-t',
                      type=int,
                      default=1,
                      help="Set number of parallel threads."
                      )

    parser.add_argument('--mismatch','-m',
                      type=int,
                      default=0,
                      help="Set number of mismatch."
                      )

    parser.add_argument('--printheader','-p',
                      required=False,
                      action="store_true",
                      help="Create output files with headers."
                      )

    args = parser.parse_args()

    verbosity = True if args.verbose else False
    print_header = True if args.printheader else False
    sample = args.input[0]
    outdir = args.outdir
    n_threads = args.nthreads
    m = args.mismatch

    out =  pipeline_hlatx(sample, outdir, m, n_threads, print_header, verbosity)
