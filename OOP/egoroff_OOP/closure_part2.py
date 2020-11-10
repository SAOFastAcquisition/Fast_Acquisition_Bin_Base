from datetime import datetime
from time import perf_counter

def average():
    a = []

    def sequenser(number):
        nonlocal a
        a.append(number)
        print(a)
        return sum(a) / len(a)

    return sequenser


def average1():
    sum = 0
    count = 0

    def sequencer(number):
        nonlocal sum
        nonlocal count
        sum += number
        count += 1
        print(f'Summa = {sum}; count = {count}')
        return sum / count

    return sequencer


def timer():
    start = perf_counter()

    def inner_time():
        nonlocal start
        delta_time = perf_counter() - start
        start = perf_counter()
        return delta_time

    return inner_time


def add_num(a,b):
    return a + b


def counter_in(func):
    count = 0

    def inner(*args, **kwargs):
        nonlocal count
        count += 1
        print(f'Фукция {func.__name__} вызывалась {count} раз')
        return func(*args)

    return inner




if __name__ == '__main__':

    time = timer()
    t1 = time()
    t2 = time()
    t3 = time()
    t4 = time()
    t5 = time()

    print(t1, t2, t3, t4, t5)


    loc_av = average1()
    print(loc_av(3))
    print(loc_av(3))
    print(loc_av(6))

    print(datetime.now())

    f = counter_in(add_num)

    print(f(5, 10))

    print(f(5, -10))