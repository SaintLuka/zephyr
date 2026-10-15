---
sidebar_position: 4

toc_min_heading_level: 2
toc_max_heading_level: 2
---

# Составная задача Римана

## Постановка задачи

Классическая постановка композитной задачи Римана:

Параметр $\delta > 0$. Требуется вычислить средний поток $\mathbf{F}_{\Delta t}$ через грань $x = 0$ за время $\Delta t$.
Задача не является автомодельной. Точное решение состоит из каскада взаимодействий волн, которые формируются от двух исходных разрывов.

## Некоторые определения

Для построения алгоритма приближенного решения CRP воспользуемся решателем HLLC для приближенного решения классической задачи Римана.
От решателя HLLC понадобится не только значение на грани, но и вся волновая конфигурация:
$$
(S_{L}, S_{C}, S_{R}, \mathbf{Q}^*_{L}, \mathbf{F}^*_{L}, \mathbf{Q}^*_{R}, \mathbf{F}^*_{R}) =
\mathrm{HLLC}(\mathbf{Q}_L, \mathbf{F}_L, \mathbf{Q}_R, \mathbf{F}_R).
$$
Для построения решения могут использоваться другие приближенные решатели, которые позволяют получить полную трёхволновую конфигурацию.

Характеристикой будем называть луч в плоскости $(x, t)$, направленный вдоль возрастания времени. Характеристику будем обозначать тремя числами $C = (x, \, t, \, \lambda)$. Луч из точки $(x, \, t)$ в плоскости с наклоном $\lambda = dx/dt \in \mathbb{R}$.

Пересечение двух характеристик обозначим символом $\cap$:
$$
(x^*, t^*) = C_1 \cap C_2, \quad C_1 = (x_1, t_1, S_1), \quad C_2 = (x_2, t_2, S_2),
$$

Нас интересуют только пересечения характеристик "в будущем", то есть $t^* > \max(t_1, t_2)$. Если таких пересечений не существует, то полагаем $t^* = +\infty$, в этом случае значение $x^{*}$ не определено. Такие значение считаем вполне легальными.
Формулы для вычисления пересечений:


Также в алгоритмах потребуется формула для вычисления времени пересечения характеристики с осью $x = 0$, которое будем обозначать буквой $\tau$:
$$
t_1 = \tau(C) = 
$$
Как и ранее, нас интересуют только пересечения в будущем, поэтому значение $+\infty$ актуально, если таких пересечений нет.

## Алгоритм

import Config11 from '@site/static/img/crp/1_1.svg';
import Config12 from '@site/static/img/crp/1_2.svg';
import Config2 from '@site/static/img/crp/2.svg';
import Config31 from '@site/static/img/crp/3_1.svg';
import Config32 from '@site/static/img/crp/3_2.svg';
import Config41 from '@site/static/img/crp/4_1.svg';
import Config42 from '@site/static/img/crp/4_2.svg';
import Config43 from '@site/static/img/crp/4_3.svg';
import Config44 from '@site/static/img/crp/4_4.svg';
import Config51 from '@site/static/img/crp/5_1.svg';
import Config52 from '@site/static/img/crp/5_2.svg';
import Config53 from '@site/static/img/crp/5_3.svg';
import Config54 from '@site/static/img/crp/5_4.svg';
import Config55 from '@site/static/img/crp/5_5.svg';


Считаем разрыв в точке $O = (x_0, \, t_0) = (0, \, 0)$

$$
(S_{0L}, S_{0C}, S_{0R}, \mathbf{Q}^*_{0L}, \mathbf{F}^*_{0L}, \mathbf{Q}^*_{0R}, \mathbf{F}^*_{0R}) =
\mathrm{HLLC}(\mathbf{Q}^B_L, \mathbf{F}^B_L, \mathbf{Q}^B_R, \mathbf{F}^B_R).
$$

Также имеется характеристика $C_u = (-\delta, 0, u_L)$.

Первое взаимодействие $O_1 = (x_1, \, t_1) = C_u \cap C_{0L}$.

### 1. Нет взаимодействия, $t_1 \ge \Delta t$.

### 1.1. $u_L \le 0$ или $S_L \le 0$.

<div className="my-svg-config">
  <Config11 width="100%" height="auto" />
</div>

Нет взаимодействия за время $\Delta t$. Интерфейс не пересекает грань. Поток равен потоку на грани для начального разрыва (точка $O$):
$$
\mathbf{F}_{\Delta t} = W_{\mathbf{F}}(0; \, O).
$$

<div style={{ clear: 'both' }}></div>

### 1.2. $u_L > 0$ и $S_L > 0$.

<div className="my-svg-config">
  <Config12 width="100%" height="auto" />
</div>

Нет взаимодействия за время $\Delta t$. Интерфейс может пересекать грань за время $\Delta t$. Поток складывается из одной или двух компонент:
$$
\xi = \max \Big( \frac{\tau_1}{\Delta t}, \, 1 \Big),
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

Второе взаимодействие: $O_2 = (x_2, \, t_2) = C_{1R} \cap C_{0C}$.

### 2. Сверхзвук, $t_1 < \Delta t$ и $x_1 \ge 0$.

<div className="my-svg-config">
  <Config2 width="100%" height="auto" />
</div>

Сверхзвуковое течение, точка $O_1$ справа. Интегрируем по треугольнику $(O, O_1, O_{\Delta t})$.
$$
G(O, O_1) = \mathbf{F}^{B}_{L} \cdot (x_1 - x_0) - \mathbf{Q}^{B}_{L} \cdot (t_1 - t_0).
$$
Значения на отрезке $(O_1, \, O_{\Delta t})$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(\lambda; \, O_1), \qquad \lambda = \frac{x_0 - x_1}{\Delta t - t_1}.
$$
$$
G(O_1, O_{\Delta t}) = \mathbf{F} \cdot (x_0 - x_1) - \mathbf{Q} \cdot (\Delta t - t_1).
$$
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \big(
G(O, O_1) + G(O_1, O_{\Delta t})
\big).
$$
<div style={{ clear: 'both' }}></div>


### 3. Нет второго взаимодействия, $t_2 \ge \Delta t$.

### 3.1. $t_2 = +\infty$ или ($t_2 \ge \Delta t$ и $x_2 > 0$).

<div className="my-svg-config">
  <Config31 width="100%" height="auto" />
</div>

Взаимодействия $O_2$ нет, или точка $O_3$ далеко справа.
Интегрируем по треугольнику $(O, O_1, O_{\Delta t})$.
$$
G(O, O_1) = \mathbf{F}^{*}_{0L} \cdot (x_1 - x_0) - \mathbf{Q}^{*}_{0L} \cdot (t_1 - t_0).
$$
Значения на отрезке $(O_1, \, O_{\Delta t})$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(\lambda; \, O_1), \qquad \lambda = \frac{x_0 - x_1}{\Delta t - t_1}.
$$
$$
G(O_1, \Delta t) = \mathbf{F} \cdot (x_0 - x_1) - \mathbf{Q} \cdot (\Delta t - t_1).
$$
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \big(
G(O, O_1) + G(O_1, \Delta t)
\big).
$$

<div style={{ clear: 'both' }}></div>

### 3.2. $t_2 \ge \Delta t$ и $x_2 < 0$.

<div className="my-svg-config">
  <Config32 width="100%" height="auto" />
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

Третье взаимодействие $O_3 = (x_3, \, t_3) = C_{1C} \cap S_{2L}$.
<div style={{ clear: 'both' }}></div>

### 4. Интерфейс движется вправо, $S_{0C} \ge 0$.

Точка $O_2$ расположена справа и $x_2 \ge 0$.

### 4.1. $t_3 = +\infty$ или ($t_3 \ge \Delta t$ и $x_3 \le 0$).

<div className="my-svg-config">
  <Config41 width="100%" height="auto" />
</div>

Взаимодействия $O_3$ нет, или точка далеко слева. Интегрируем по треугольнику $(O, O_2, O_{\Delta t})$.
$$
G(O, O_2) = \mathbf{F}^*_{0L} \cdot (x_2 - x_0) - \mathbf{Q}^*_{0L} \cdot (t_2 - t_0).
$$
Значения на отрезке $(O_2, \, O_{\Delta t})$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(\lambda; \, O_2), \qquad \lambda = \frac{x_0 - x_2}{\Delta t - t_2}.
$$
$$
G(O_2, O_{\Delta t}) = \mathbf{F} \cdot (x_0 - x_2) - \mathbf{Q} \cdot (\Delta t - t_2).
$$
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \big(
G(O, O_2) + G(O_2, O_{\Delta t})
\big).
$$
<div style={{ clear: 'both' }}></div>

### 4.2. $t_3 \ge \Delta t$ и $x_3 > 0$.

<div className="my-svg-config">
  <Config42 width="100%" height="auto" />
</div>

Точка $O_3$ далеко справа. Интегрируем по треугольнику $(O, O_1, O_{\Delta t})$.
$$
G(O, O_1) = \mathbf{F}^*_{0L} \cdot (x_1 - x_0) - \mathbf{Q}^*_{0L} \cdot (t_1 - t_0).
$$
Значения на отрезке $(O_1, \, O_{\Delta t})$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(\lambda; \, O_1), \qquad \lambda = \frac{x_0 - x_1}{\Delta t - t_1}.
$$
$$
G(O_1, O_{\Delta t}) = \mathbf{F} \cdot (x_0 - x_1) - \mathbf{Q} \cdot (\Delta t - t_1).
$$
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \big(
G(O, O_1) + G(O_1, O_{\Delta t})
\big).
$$
<div style={{ clear: 'both' }}></div>

### 4.3. Точка $O_3$ справа, $t_3 < \Delta t$, $x_3 \ge 0$.

<div className="my-svg-config">
  <Config43 width="100%" height="auto" />
</div>

Интегрируем по многоугольнику $(O, O_2, O_3, O_{\Delta t})$.
$$
G(O, O_2) = \mathbf{F}^*_{0L} \cdot (x_2 - x_0) - \mathbf{Q}^*_{0L} \cdot (t_2 - t_0).
$$
$$
G(O_2, O_3) = \mathbf{F}^*_{1R} \cdot (x_3 - x_2) - \mathbf{Q}^*_{1R} \cdot (t_3 - t_2).
$$
Значения на отрезке $(O_3, \, O_{\Delta t})$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(\lambda; \, O_3), \qquad \lambda = \frac{x_0 - x_3}{\Delta t - t_3}.
$$
$$
G(O_3, O_{\Delta t}) = \mathbf{F} \cdot (x_0 - x_3) - \mathbf{Q} \cdot (\Delta t - t_3).
$$
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \big(
G(O, O_2) + G(O_2, O_3) + G(O_3, O_{\Delta t})
\big).
$$
<div style={{ clear: 'both' }}></div>

### 4.4. Точка $O_3$ слева, $t_3 < \Delta t$, $x_3 < 0$.

<div className="my-svg-config">
  <Config44 width="100%" height="auto" />
</div>

Интегрируем по многоугольнику $(O, O_1, O_3, O_{\Delta t})$.
$$
G(O, O_1) = \mathbf{F}^*_{0L} \cdot (x_1 - x_0) - \mathbf{Q}^*_{0L} \cdot (t_1 - t_0).
$$
$$
G(O_1, O_3) = \mathbf{F}^*_{1R} \cdot (x_3 - x_1) - \mathbf{Q}^*_{1R} \cdot (t_3 - t_1).
$$
Значения на отрезке $(O_3, \, O_{\Delta t})$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(\lambda; \, O_3), \qquad \lambda = \frac{x_0 - x_3}{\Delta t - t_3}.
$$
$$
G(O_3, O_{\Delta t}) = \mathbf{F} \cdot (x_0 - x_3) - \mathbf{Q} \cdot (\Delta t - t_3).
$$
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \big(
G(O, O_1) + G(O_1, O_3) + G(O_3, O_{\Delta t})
\big).
$$
<div style={{ clear: 'both' }}></div>

### 5. Интерфейс движется влево, $S_{0C} < 0$.

Точка $O_2$ расположена слева и $x_2 < 0$.

Найдем момент пересечения интерфейса и грани $\tau_1 = \tau(C_{1C})$.

### 5.1. $t_3 \ge \Delta t$ и $\tau_1 \ge \Delta t$.

<div className="my-svg-config">
  <Config51 width="100%" height="auto" />
</div>

Взаимодействия $O_3$ нет или оно далеко, интерфейс не пересекает грань.
Интегрируем по треугольнику $(O, O_2, O_{\Delta t})$.
$$
G(O, O_2) = \mathbf{F}^*_{0R} \cdot (x_2 - x_0) - \mathbf{Q}^*_{0R} \cdot (t_2 - t_0).
$$
Значения на отрезке $(O_2, \, O_{\Delta t})$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(\lambda; \, O_2), \qquad \lambda = \frac{x_0 - x_2}{\Delta t - t_2}.
$$
$$
G(O_2, O_{\Delta t}) = \mathbf{F} \cdot (x_0 - x_2) - \mathbf{Q} \cdot (\Delta t - t_2).
$$
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \big(
G(O, O_2) + G(O_2, O_{\Delta t})
\big).
$$
<div style={{ clear: 'both' }}></div>

### 5.2. $\tau_1 < \Delta t$ и ($t_3 \ge \Delta t$ или $x_3 > 0$).

<div className="my-svg-config">
  <Config52 width="100%" height="auto" />
</div>

Экзотика. Взаимодействия $O_3$ нет или оно справа, интерфейс пересекает грань. Интегрируем по треугольнику $(O, O_1, O_{\Delta t})$.
$$
G(O, O_1) = \mathbf{F}^*_{0L} \cdot (x_1 - x_0) - \mathbf{Q}^*_{0L} \cdot (t_1 - t_0).
$$
Значения на отрезке $(O_1, \, O_{\Delta t})$:
$$
(\mathbf{Q}, \, \mathbf{F}) = W(\lambda; \, O_1), \qquad \lambda = \frac{x_0 - x_1}{\Delta t - t_1}.
$$
$$
G(O_1, O_{\Delta t}) = \mathbf{F} \cdot (x_0 - x_1) - \mathbf{Q} \cdot (\Delta t - t_1).
$$
$$
\mathbf{F}_{\Delta t} = \frac{1}{\Delta t} \big(
G(O, O_1) + G(O_1, O_{\Delta t})
\big).
$$
<div style={{ clear: 'both' }}></div>

### 5.3.

<div className="my-svg-config">
  <Config53 width="100%" height="auto" />
</div>

<div style={{ clear: 'both' }}></div>

### 5.4.

<div className="my-svg-config">
  <Config54 width="100%" height="auto" />
</div>

<div style={{ clear: 'both' }}></div>

### 5.5.

<div className="my-svg-config">
  <Config55 width="100%" height="auto" />
</div>

<div style={{ clear: 'both' }}></div>


## Модификации

Упрощения алгоритма:
- Между одинаковыми матеаралами можно использовать двухволновой приближенный решатель. К примеру, HLL. Это значительно сокращает число конфигураций.
- Можно добавить ограничения на скорости $S_L \le 0$ и $S_R \ge 0$, это также сократит число взаимодействий.
