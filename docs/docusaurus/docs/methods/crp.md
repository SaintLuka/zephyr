---
sidebar_position: 4

toc_min_heading_level: 2
toc_max_heading_level: 2
---

# Составная задача Римана

## Постановка задачи

import IVP from '@site/static/img/crp/ivp.svg';

<div className="my-svg-config">
  <IVP width="100%" height="auto" />
</div>

Составная задача Римана *(Composite Riemann Problem, CRP)* — одномерная задача Коши для двухматериальной газовой динамики, в которой начальное условие задается константами в трёх областях.

$$
\mathbf{Z}(x, \, 0) = 
\begin{cases}
\mathbf{Z}^{A}_L, &\qquad\;\; x < -\delta, \\
\mathbf{Z}^{B}_L, &-\delta < x < 0, \\
\mathbf{Z}^{B}_R, &\;\;\;0 < x.
\end{cases}
$$

Параметр $\delta > 0$. Область слева ($x < -\delta$) заполнена веществом $A$, область справа — веществом $B$. 
Векторы состояния:
$$
\mathbf{Z}^{A}_L = (\rho^A_L, \, \mathbf{v}_L, \, P_L, \, 1, \, 0)^T,
$$
$$
\mathbf{Z}^{B}_L = (\rho^B_L, \, \mathbf{v}_L, \, P_L, \, 0, \, 1)^T,
$$
$$
\mathbf{Z}^{B}_R = (\rho^B_R, \, \mathbf{v}_R, \, P_R, \, 0, \, 1)^T,
$$
здесь $\rho$ — плотность, $\mathbf{v}$ — скорость, $P$ — давление. 
Последние две компоненты в векторе состояния — массовые концентрации $\mathbf{\beta}$. 
Вещества слева ($x < 0$) находятся в равновесии по скорости, давлению (и температуре).

Задача не является автомодельной. 
Точное решение состоит из каскада взаимодействий волн, которые формируются от двух исходных разрывов. 
Требуется вычислить средний поток $\mathbf{F}_{\Delta t}$ через грань $x = 0$ за время $\Delta t$.

## Определения

Для поиска приближенного решения CRP используется решатель HLLC для классической задачи Римана.
От решателя HLLC требуется полная волновая конфигурация, а не только значение потока на грани:
$$
(S_{L}, S_{C}, S_{R}, \mathbf{Q}^*_{L}, \mathbf{F}^*_{L}, \mathbf{Q}^*_{R}, \mathbf{F}^*_{R}) =
\mathrm{HLLC}(\mathbf{Q}_L, \mathbf{F}_L, \mathbf{Q}_R, \mathbf{F}_R).
$$
Важно использовать функцию от четырех аргументов:
от консервативных величин $\mathbf{Q}_{L/R}$ и потоков $\mathbf{F}_{L/R}$,
которые в общем случае не согласованы. 
Решение вдоль произвольной характеристики $S$ обозначим буквой $W$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(S; \, O_1)
\quad \text{или} \quad
\mathbf{F} = W_{\mathbf{F}}(S; \, O_1),
$$
здесь $O_1$ — точка, в которой рассматривается распад разрыва.
Функция $W$ просто выбирает величины из полной волновой конфигурации.
Для построения решения могут использоваться и другие приближенные решатели с трёхволновой конфигурацией.

Будем называть **характеристикой** луч в плоскости $(x, \, t)$, направленный вдоль возрастания времени.
Характеристика $C$ задается тремя числами $(x, \, t, \, S)$: координатами начала и наклоном.
Из точки, в которой рассматривается приближенный распад разрыва HLLC, 
выходит три характеристики $C_L$, $C_C$, $C_R$ со скоростями $S_L < S_C < S_R$ соответственно.

Пересечение двух характеристик обозначим символом $\cap$:
$$
(x^*, \, t^*) = C_1 \cap C_2, \quad C_1 = (x_1, \, t_1, \, S_1), \quad C_2 = (x_2, \, t_2, \, S_2),
$$

Интересны только пересечения характеристик «в будущем», то есть $t^* > \max(t_1, \, t_2)$. 
Если таких пересечений не существует, то полагаем $t^* = +\infty$, в этом случае значение $x^{*}$ не определено.

Также в алгоритме потребуется формула для вычисления времени пересечения характеристики $C$ с осью $x = 0$, 
которое будем обозначать как $\tau(C)$. Как и ранее, нас интересуют только пересечения в будущем:
$$
t^* = \tau(C) = \begin{cases}
t - \frac{x}{S}, & \mathrm{sgn}(x) \cdot \mathrm{sgn}(S) < 0, \\
+\infty,         & \mathrm{sgn}(x) \cdot \mathrm{sgn}(S) \ge 0.
\end{cases}
$$

---

import Oint from '@site/static/img/crp/oint.svg';

<div className="my-svg-config">
  <Oint width="100%" height="auto" />
</div>


Пусть в плоскости $(x, \, t)$ вектор консервативных величин $\mathbf{Q}$ и 
вектор потоков $\mathbf{F}$ удовлетворяют системе уравнений в слабой форме, 
то есть интеграл по любому замкнутому контуру равен нулю:
$$
\oint \mathbf{F} dt - \mathbf{Q} dx = 0.
$$

Введем обозначение для интеграла от точки $O_1$ до $O_2$:
$$
\mathbf{G}(O_1, O_2) = \int_{O_1}^{O_2} \mathbf{F} dt - \mathbf{Q} dx,
$$
Значение интеграла не зависит от пути интегрирования.

Волновые конфигурации для HLL и HLLC связаны соотношениями Ренкина—Гюгонио,
поэтому для вычисления потока $\mathbf{F}_{\Delta t}$ можно использовать интегральную формулу.
Искомый поток можно вычислить по произвольному пути в плоскости $(x, \, t)$:
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \, \mathbf{G}(O, O_{\Delta t}).
$$


## Алгоритм

import Config11 from '@site/static/img/crp/1_1.svg';
import Config12 from '@site/static/img/crp/1_2.svg';
import Config2  from '@site/static/img/crp/2.svg';
import Config3  from '@site/static/img/crp/3.svg';
import Config41 from '@site/static/img/crp/4_1.svg';
import Config42 from '@site/static/img/crp/4_2.svg';
import Config43 from '@site/static/img/crp/4_3.svg';
import Config44 from '@site/static/img/crp/4_4.svg';


Считаем разрыв в точке $O = (x_0, \, t_0) = (0, \, 0)$

$$
(S_{0L}, S_{0C}, S_{0R}, \mathbf{Q}^*_{0L}, \mathbf{F}^*_{0L}, \mathbf{Q}^*_{0R}, \mathbf{F}^*_{0R}) =
\mathrm{HLLC}(\mathbf{Q}^B_L, \mathbf{F}^B_L, \mathbf{Q}^B_R, \mathbf{F}^B_R).
$$

Также имеется характеристика $C_u = (-\delta, \, 0, \, u_L)$.

Первое взаимодействие $O_1 = (x_1, \, t_1) = C_u \cap C_{0L}$.

### 1. Нет взаимодействия, $t_1 \ge \Delta t$.

### 1.1. $u_L \le 0$ или $S_{0L} \le 0$.

<div className="my-svg-config">
  <Config11 width="100%" height="auto" />
</div>

Нет взаимодействия за время $\Delta t$. Интерфейс не пересекает грань. Поток равен потоку на грани для начального разрыва (точка $O$):
$$
\mathbf{F}_{\Delta t} = W_{\mathbf{F}}(0; \, O).
$$

<div style={{ clear: 'both' }}></div>

### 1.2. $u_L > 0$ и $S_{0L} > 0$.

<div className="my-svg-config">
  <Config12 width="100%" height="auto" />
</div>

Нет взаимодействия за время $\Delta t$. 
Интерфейс может пересекать грань за время $\Delta t$. 
Поток складывается из одной или двух компонент:
$$
\xi = \min \Big( \frac{\tau_1}{\Delta t}, \, 1 \Big),
\quad
\tau_1 = \tau(C_u).
$$
$$
\mathbf{F}_{\Delta t} = \xi \cdot \mathbf{F}^B_L + (1 - \xi) \cdot \mathbf{F}^A_L,
$$
здесь $\tau_1$ — момент, когда интерфейс пересекает грань.

<div style={{ clear: 'both' }}></div>
---

Считаем разрыв в точке $O_1 = (x_1, \, t_1)$:
$$
(S_{1L}, S_{1C}, S_{1R}, \mathbf{Q}^*_{1L}, \mathbf{F}^*_{1L}, \mathbf{Q}^*_{1R}, \mathbf{F}^*_{1R}) =
\mathrm{HLLC}(\mathbf{Q}^A_L, \mathbf{F}^A_L, \mathbf{Q}^*_{0L}, \mathbf{F}^*_{0L}).
$$

### 2. Сверхзвук, $x_1 \ge 0$.

<div className="my-svg-config">
  <Config2 width="100%" height="auto" />
</div>

Сверхзвуковое течение, точка $O_1$ справа от грани. 
Интегрируем по треугольнику $(O, O_1, O_{\Delta t})$.
$$
\mathbf{G}(O, O_1) = \mathbf{F}^{B}_{L} \cdot (t_1 - t_0) - \mathbf{Q}^{B}_{L} \cdot (x_1 - x_0).
$$
Значения на отрезке $(O_1, \, O_{\Delta t})$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(\lambda; \, O_1), \qquad \lambda = \frac{x_0 - x_1}{\Delta t - t_1}.
$$
$$
\mathbf{G}(O_1, O_{\Delta t}) = \mathbf{F} \cdot (\Delta t - t_1) - \mathbf{Q} \cdot (x_0 - x_1).
$$
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \big(
\mathbf{G}(O, O_1) + \mathbf{G}(O_1, O_{\Delta t})
\big).
$$
<div style={{ clear: 'both' }}></div>

### 3. Интерфейс движется вправо, $S_{0C} \ge 0$.

<div className="my-svg-config">
  <Config3 width="100%" height="auto" />
</div>

Интегрируем по треугольнику $(O, O_1, O_{\Delta t})$.
$$
\mathbf{G}(O, O_1) = \mathbf{F}^{*}_{0L} \cdot (t_1 - t_0) - \mathbf{Q}^{*}_{0L} \cdot (x_1 - x_0).
$$
Значения на отрезке $(O_1, \, O_{\Delta t})$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(\lambda; \, O_1), \qquad \lambda = \frac{x_0 - x_1}{\Delta t - t_1}.
$$
$$
\mathbf{G}(O_1, \Delta t) = \mathbf{F} \cdot (\Delta t - t_1) - \mathbf{Q} \cdot (x_0 - x_1).
$$
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \big(
\mathbf{G}(O, O_1) + \mathbf{G}(O_1, \Delta t)
\big).
$$

<div style={{ clear: 'both' }}></div>

### 4. Интерфейс движется влево, $S_{0C} < 0$.

Второе взаимодействие: $O_2 = (x_2, \, t_2) = C_{1R} \cap C_{0C}$.

### 4.1. Нет взаимодействия, $t_2 \ge \Delta t$.

<div className="my-svg-config">
  <Config41 width="100%" height="auto" />
</div>

Точка $O_2$ далеко слева. Поток равен потоку на грани для начального разрыва (точка $O$):
$$
\mathbf{F}_{\Delta t} = W_\mathbf{F}(0; O).
$$

<div style={{ clear: 'both' }}></div>
---

Считаем разрыв в точке $O_2 = (x_2, \, t_2)$:
$$
(S_{2L}, S_{2C}, S_{2R}, \mathbf{Q}^*_{2L}, \mathbf{F}^*_{2L}, \mathbf{Q}^*_{2R}, \mathbf{F}^*_{2R}) =
\mathrm{HLLC}(\mathbf{Q}^*_{1R}, \mathbf{F}^*_{1R}, \mathbf{Q}^*_{0R}, \mathbf{F}^*_{0R}).
$$

Третье взаимодействие $O_3 = (x_3, \, t_3) = C_{1C} \cap C_{2L}$.
<div style={{ clear: 'both' }}></div>

Найдем момент пересечения интерфейса и грани $\tau_1 = \tau(C_{1C})$.

### 4.2. $t_3 \ge \Delta t$ и $\tau_1 \le \Delta t$.

<div className="my-svg-config">
  <Config42 width="100%" height="auto" />
</div>

Экзотика. Взаимодействия $O_3$ нет или оно справа, интерфейс пересекает грань. Интегрируем по треугольнику $(O, O_1, O_{\Delta t})$.
$$
\mathbf{G}(O, O_1) = \mathbf{F}^*_{0L} \cdot (t_1 - t_0) - \mathbf{Q}^*_{0L} \cdot (x_1 - x_0).
$$
Значения на отрезке $(O_1, \, O_{\Delta t})$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(\lambda; \, O_1), \qquad \lambda = \frac{x_0 - x_1}{\Delta t - t_1}.
$$
$$
\mathbf{G}(O_1, O_{\Delta t}) = \mathbf{F} \cdot (\Delta t - t_1) - \mathbf{Q} \cdot (x_0 - x_1).
$$
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \big(
\mathbf{G}(O, O_1) + \mathbf{G}(O_1, O_{\Delta t})
\big).
$$
<div style={{ clear: 'both' }}></div>

### 4.3. $t_3 \ge \Delta t$ и $\tau_1 > \Delta t$.

<div className="my-svg-config">
  <Config43 width="100%" height="auto" />
</div>

Взаимодействия $O_3$ нет или оно далеко, интерфейс не пересекает грань.
Интегрируем по треугольнику $(O, O_2, O_{\Delta t})$.
$$
\mathbf{G}(O, O_2) = \mathbf{F}^*_{0R} \cdot (t_2 - t_0) - \mathbf{Q}^*_{0R} \cdot (x_2 - x_0).
$$
Значения на отрезке $(O_2, \, O_{\Delta t})$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(\lambda; \, O_2), \qquad \lambda = \frac{x_0 - x_2}{\Delta t - t_2}.
$$
$$
\mathbf{G}(O_2, O_{\Delta t}) = \mathbf{F} \cdot (\Delta t - t_2) - \mathbf{Q} \cdot (x_0 - x_2).
$$
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \big(
\mathbf{G}(O, O_2) + \mathbf{G}(O_2, O_{\Delta t})
\big).
$$
<div style={{ clear: 'both' }}></div>

### 4.4. Третий распад, $t_3 < \Delta t$.

Считаем разрыв в точке $O_3 = (x_3, \, t_3)$:
$$
(S_{3L}, S_{3C}, S_{3R}, \mathbf{Q}^*_{3L}, \mathbf{F}^*_{3L}, \mathbf{Q}^*_{3R}, \mathbf{F}^*_{3R}) =
\mathrm{HLLC}(\mathbf{Q}^*_{1L}, \mathbf{F}^*_{1L}, \mathbf{Q}^*_{2L}, \mathbf{F}^*_{2L}).
$$

<div className="my-svg-config">
  <Config44 width="100%" height="auto" />
</div>

Интегрируем по многоугольнику $(O, O_2, O_3, O_{\Delta t})$.
$$
\mathbf{G}(O, O_2) = \mathbf{F}^*_{0R} \cdot (t_2 - t_0) - \mathbf{Q}^*_{0R} \cdot (x_2 - x_0).class
$$
$$
\mathbf{G}(O_2, O_3) = \mathbf{F}^*_{2L} \cdot (t_3 - t_2) - \mathbf{Q}^*_{2L} \cdot (x_3 - x_2).
$$
Значения на отрезке $(O_3, O_{\Delta t})$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(\lambda; \, O_3),  \qquad \lambda = \frac{x_0 - x_3}{\Delta t - t_3}.
$$
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \big(
\mathbf{G}(O, O_2) + \mathbf{G}(O_2, O_3) + \mathbf{G}(O_3, O_{\Delta t})
\big).
$$

<div style={{ clear: 'both' }}></div>

## Модификации

Упрощения алгоритма:
- Между одинаковыми материалами можно использовать двухволновой приближенный решатель. К примеру, HLL. Это сократит число взаимодействий.
- Можно добавить ограничения на скорости $S_L \le 0$ и $S_R \ge 0$, это уменьшит число экзотических конфигураций.
- Можно рассмотреть задачу, в которой два состояния слева не находятся в равновесии.
