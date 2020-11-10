from functools import wraps
from datetime import datetime


def decor_func(func):

    @wraps(func)    # декоратор восстанавливающий описание декорируемой функции
    def inner(*args, **kwargs):
        print(f'Begin decoration function {func.__name__}')
        print(f'Result of execution function {func.__name__} is {func(*args)}')
        print(f'End decoration function {func}')
    # inner.__name__ = func.__name__

    return inner


def decor_table(func):

    def inner(*args, **kwargs):
        print(f'<table> {func.__name__}')
        print(f'Result of execution function {func.__name__} is {func(*args)}')
        print(f'</table> {func}')

    inner.__name__ = func.__name__
    return inner

@decor_func
@decor_table
def add_num(a, b, c):
    return (a + b + c) / 2


def timer(func):

    start = datetime.now()
    @wraps(func)
    def inner(*args):

        return func(*args)

    print(datetime.now() - start)
    return inner


@timer
def generator(length):
    seq = [i ** 2 + i for i in range(length)]
    return seq


if __name__ == '__main__':

    # add_num = decor_table(decor_func(add_num))
    # print(add_num.__name__)
    # add_num(1, 3, 5)
    # print(add_num.__name__)
    #
    # add_num = decor_func(decor_table(add_num))
    # print(add_num.__name__)
    # add_num(1, 3, 5)
    # print(add_num.__name__)

    # add_num(1, 3, 5)
    # print(add_num.__name__)

    l = generator(10)
    print(l)
    print(*l)

    d = {'x': 1, 'y': 2, 'z': 3}

    def test(x, y, z):
        print(x, y, z)


    test(*d)

    pass