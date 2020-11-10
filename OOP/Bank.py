

class Bank_account:

    def __init__(self, name, balance, passport):
        self.__name = name
        self.__balance = balance
        self.__passport = passport

    # def public_print_acc(self):
    #     print(self.name, self.balance, self.passport)

    def privet_print_acc(self):
        print(self.__name, self.__balance, self.__passport)

    # def protected_print_acc(self):
    #     print(self._name, self._balance, self._passport)


account1 = Bank_account('Bob', 10, 4949494)
# account1.public_print_acc()
# print(account1._name)
account1.privet_print_acc()
print(account1.__name)

