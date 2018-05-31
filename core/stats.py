from astropy.table import join
import numpy as np


def rename_stats(s, post):
    for c in s.colnames:
        if c != 'age' and c != 'fe/h':
            s.rename_column(c, c+post)
    return s


def stat_attr(out, func, post, pre=None):
    meds = out.groups.aggregate(func)
    meds = rename_stats(meds, post)
    if pre is not None:
        meds = join(pre, meds, keys=['age', 'fe/h'])
    return meds


def create_stats(out):
    out = out.group_by(['age', 'fe/h'])
    stats = stat_attr(out, np.median, '_median')
    stats = stat_attr(out, np.mean, '_mean', stats)
    stats = stat_attr(out, np.std, '_std', stats)
    stats = stat_attr(out, np.min, '_min', stats)
    stats = stat_attr(out, np.max, '_max', stats)
    return stats
