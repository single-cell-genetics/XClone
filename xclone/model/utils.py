"""Utility functions for XClone.
"""
# Author: Rongting Huang

def pairwise0(iterable_arr):
    """
    s -> (s0, s1), (s2, s3), (s4, s5), ...

    Example:
    l = [1,2,3,4,5,6]
    for x, y in pairwise(l):
       print("%d + %d = %d" % (x, y, x + y))

    """
    a = iter(iterable_arr)
    return zip(a, a)

def grouped(iterable, n):
    """
    generalized
    s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...

    Example:
    for x, y in grouped(l, 2):
        print("%d + %d = %d" % (x, y, x + y))
    """
    return zip(*[iter(iterable)]*n)

def pairwise1(iterable_arr):
    """
    s -> (s0, s1), (s1, s2), (s2, s3), ...
    """
    a = iter(iterable_arr)
    b = iter(iterable_arr[1:])
    return zip(a, b)

def pairwise2(iterable_arr):
    """
    s -> (s0, s1), (s1, s2), (s2, s3), ...
    Example:
    data = [1,2,3,4,5,6]
    """
    for i,k in zip(iterable_arr[0::1], iterable_arr[1::1]):
        print(str(i), '+', str(k), '=', str(i+k))
    return 0

def class2onehot(prob_, states_class_num):
    index_ = np.argmax(prob_, axis=-1)
    one_hot_ = np.eye(states_class_num)[index_]
    return one_hot_
