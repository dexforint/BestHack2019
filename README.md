# BestHack2019
Solutions of BestHack hackathon
Команда: satana
Капитан: Ларичев Дмитрий Юрьевич
Email: dexforint@mail.ru

# Требования:
Python 3
numpy
pandas
scipy
argparse (стандартная библиотека python 3)

# Алгоритм проверки решения:
1) Поместите в папку с solution.py файлы F.csv и Wind.csv
2) Запустите в консоли программу с аргументами:
  python ./solution.py -V0 [V0] -H0 [H0] -m [m] [-a [a1,a2,a3,...] -R [R] -dt [dt]]
  
  # Обязательные параметры
    -V0 - начальная скорость
    -H0 - начальная высота
    -m - масса
  
  # Необязательные параметры
    -a - углы через запятую без пробелов в радианах (например, -a 0,1.07,3.14 ) [0]
    -R - радиус шара [1]
    -dt - временной шаг (чем меньше, тем больше точность и время выполнения) [0.01]
  
# В итоге у вас для каждого заданного угла выведутся искомые координаты, а в папке появятся файлы с зависимостями координат и скорости от времени
