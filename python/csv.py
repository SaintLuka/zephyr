import numpy as np


class CsvFile:

    def __init__(self, filename, variables=None):
        """
        :param filename: Имя .csv файла
        :param variables: Список переменных для чтения (массив)
        """
        file = open(filename, 'r')
        line = file.readline()
        self.variables, indices = self.read_variables(line, variables)

        values = []
        for line in file:
            nums = line.split(',')
            vals = [float(nums[i]) for i in indices]
            values.append(vals)

        self.values = np.array(values)

        file.close()

        for i, var in enumerate(self.variables):
            setattr(self, var, np.array(self.values[:, i]))

        print('Successful reading of the .csv file;')
        print('    Variables:  ', self.variables)
        print('    Data shape: ', self.values.shape)

    @staticmethod
    def read_variables(line, variables=None):
        """
        Считать список переменных из строки
        :param line: Строка с именами переменных
        :param variables: Желаемые переменные
        :return: Пара: (список переменных, их индексы в строке)
        """
        if line[0] != '#':
            message = 'Первая строка csv файла должна начинаться с решетки # ' \
                      'и содержать список имен массивов'
            raise Exception(message)

        csv_variables = line[1:-1].replace(' ', '').split(',')

        if variables is None:
            # Считать все массивы
            res_variables = csv_variables
            indices = range(len(csv_variables))
        else:
            # Выбрать массивы из списка
            res_variables = []

            indices = []

            for var in variables:
                if var in csv_variables:
                    idx = csv_variables.index(var)
                    indices.append(idx)
                    res_variables.append(var)

        indices = np.array(indices)

        return res_variables, indices
