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

    def sort1D(self):
        if not (hasattr(self, 'x')):
            raise 'Read array "x" from .csv'

        idx = np.argsort(self.x)
        for i, var in enumerate(self.variables):
            setattr(self, var, getattr(self, var)[idx])
            
        print('x.size1: ', self.x.size)
            
        arr, idx = np.unique(self.x, return_index=True)
        for i, var in enumerate(self.variables):
            setattr(self, var, getattr(self, var)[idx])
            
        print('x.size2: ', self.x.size)


    def as2D(self):
        """
        Интерпретировать данные как двумерные массивы на структурированной сетке
        Необходимо наличие массивов 'x' и 'y' в данных.
        Добавляет новые переменные xmin, xmax, ymin, ymax
        """
        if not (hasattr(self, 'x') and hasattr(self, 'y')):
            raise 'Read arrays "x" and "y" from .csv'

        print(self.x)

        ny = np.argmax(self.x > self.x[0])
        nx = self.x.size // ny

        for i, var in enumerate(self.variables):
            setattr(self, var, getattr(self, var).reshape((nx, ny)))
      

        dx = (self.x[-1, 0] - self.x[0, 0]) / (nx - 1)
        dy = (self.y[0, -1] - self.y[0, 0]) / (ny - 1)

        self.xmin = 1.0e-6 * np.round(1.0e6 * (self.x[0, 0] - 0.5 * dx))
        self.ymin = 1.0e-6 * np.round(1.0e6 * (self.y[0, 0] - 0.5 * dy))
        self.xmax = 1.0e-6 * np.round(1.0e6 * (self.x[-1, -1] + 0.5 * dx))
        self.ymax = 1.0e-6 * np.round(1.0e6 * (self.y[-1, -1] + 0.5 * dy))

