
def main_closure(name):

    def func(var=name):
        print(f'Your name is {var}')

    return func


def adder(value):

    def inner(a):
        print(a+value)

    return inner


def counter():
    count = 0

    def inner():
        nonlocal count
        count += 1
        return count

    return inner


if __name__ == '__main__':
    a = main_closure('Ivan')
    b = main_closure('Fred')
    a()
    b()
    a('Mike')
pass

q = counter()
print(q())
print(q())
print(q())
