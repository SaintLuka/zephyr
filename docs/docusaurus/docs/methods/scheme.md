---
sidebar_position: 2
description: "Схема предиктор-корректор и MUSCL"
---

# Схемы второго порядка

## Уравнения Эйлера

Рассматривается система уравнений Эйлера в консервативной форме:

$$
\frac{\partial \mathbf{q}}{\partial t} + \frac{\partial \mathbf{f}}{\partial x} + \frac{\partial \mathbf{g}}{\partial y} + \frac{\partial \mathbf{h}}{\partial z} = 0,
$$

где

$$
\mathbf{q} = \begin{pmatrix}\rho \\ \rho u \\ \rho v \\ \rho w \\ \rho E\end{pmatrix},\quad
\mathbf{f} = \begin{pmatrix}\rho u \\ \rho u^2+P \\ \rho uv \\ \rho uw \\ u(\rho E+P)\end{pmatrix},\quad
\mathbf{g} = \begin{pmatrix}\rho v \\ \rho vu \\ \rho v^2+P \\ \rho vw \\ v(\rho E+P)\end{pmatrix},\quad
\mathbf{h} = \begin{pmatrix}\rho w \\ \rho wu \\ \rho wv \\ \rho w^2+P \\ w(\rho E+P)\end{pmatrix}.
$$

Здесь $\rho$ — плотность, 
$\mathbf{v} = (u, v, w)^T$ — скорость, 
$P$ — давление, 
$e$ — удельная внутренняя энергия, 
$E = e + \frac{1}{2}\mathbf{v}^2$ — удельная полная энергия. 
Введём вектор примитивных переменных $\mathbf{z} = (\rho, u, v, w, P)^T$.

## Расчётная схема

Используется явная конечно-объёмная схема второго порядка по времени и пространству. 
Для повышения точности по пространству применяется реконструкция градиента примитивных 
переменных и интерполяция на грани; по времени — двухстадийная схема предиктор–корректор[^Rodionov]. 

Расчётная сетка состоит из ячеек произвольного вида (многоугольников/многогранников); 
адаптивные ячейки рассматриваются как ячейки общего вида.

---

Этапы одного шага по времени:

1. **Ограничение по Куранту.** Вычисление $\Delta t$, удовлетворяющего условию устойчивости Куранта.
2. **Реконструкция градиента.** В каждой ячейке восстанавливается $\nabla \mathbf{z}_i$ 
   (обычный или [лимитированный градиент](slope/#22-ячейко-ориентированная-интерполяция)).
3. **Предиктор.** Вычисляются значения на полушаге:
   $$
   \mathbf{q}^{n+1/2}_i = \mathbf{q}^n_i - \frac{\Delta t}{2V_i}\sum_{j\in\nu(i)} \mathcal{R}^T \mathbf{f}(\mathcal{R}\mathbf{z}^n_{ij}) S_{ij}, \qquad \mathbf{q}^{n+1/2}_i \to \mathbf{z}^{n+1/2}_i,
   $$
   где $\mathcal{R}$ — матрица поворота в локальную систему координат грани $\sigma_{ij}$, 
   $\mathbf{z}^n_{ij}$ — интерполированные примитивные переменные на грани, $S_{ij}$ — площадь грани. 
   Заметим, что на данном этапе поток вычисляется по простой формуле от значений внутри ячейки.
4. **Корректор.** Вычисляются значения на новом временном слое:
   $$
   \mathbf{q}^{n+1}_i = \mathbf{q}^n_i - \frac{\Delta t}{V_i}\sum_{j\in\nu(i)} \mathcal{R}^T \mathbf{F}(\mathcal{R}\mathbf{z}^{n+1/2}_{ij}, \mathcal{R}\mathbf{z}^{n+1/2}_{ji}) S_{ij}, \qquad \mathbf{q}^{n+1}_i \to \mathbf{z}^{n+1}_i.
   $$
   Здесь $\mathbf{F}(\mathbf{z}_L, \mathbf{z}_R)$ — численный поток (метод Годунова, Русанова, Роу, HLL, HLLC и другие).
5. **Сеточная адаптация.** В соответствии с [некоторым критерием](criterion#пороговый-критерий-адаптации) 
   выставляются флаги адаптации (`SPLIT`, `COARSE`, `NONE`), затем выполняется адаптация сетки.


[^Rodionov]: *A.V. Rodionov.* Methods of increasing the accuracy in Godunov’s scheme. 
            USSR Comput. Math. Math. Phys. 27 (6), 164–169 (1987).