<p align="center">
  <img src="docs/logo.svg" width="400" alt="non-parametric-dsp Logo">
</p>

[ [English](docs/README.eng.md) ]

non-parametric-dsp - C++ библиотека для цифровой обработки нестационарных сигналов

# Оглавление
- [Функционал](#Функционал)
- [Начало работы](#начало-работы)
- [Credits](#credits)
- [Лицензия](#лицензия)
- [Контакты](#контакты)

# Функционал
non-parametric-dsp - C++ библиотека для цифровой обработки нестационарных сигналов, предоставляющая следующий функционал:
- аппроксимация сигналов различными методами
- численное дифференцирование и интегрирование (в том числе и дробных степеней)
- широкий спектр алгоритмов фильтрации
- вычисление мгновенных частот и амплитуд по Гильберту и по Тихонову
- фазовращение сигналов
- гибко настраиваемая эмпирическая модовая декомпозиция
- частотно-временной анализ
- токенизация сигналов для обучения нейросетей
non-parametric-dsp - исследовательский проект, предоставляющий большое количество экспериментальных и уникальных подходов;

Для более ознакомления с функционалом non-parametric-dsp пожалуйста, ознакомьтесь с примерами из директории examples

> [!IMPORTANT]
> Для использования в своих проектах, пожалуйста, используйте векту “for_linkage”; Она предоставляет headeronly библиотеку без встроенных возможностей визуализации

## Начало работы
Для начала работы вам потребуется установить следующий набор зависимостей:

#### <img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/git/git-original.svg" width="24" height="24" alt="Git Logo"/> Git

#### <img src="https://xmake.io/assets/img/logo.svg" width="24" height="24" alt="XMake Logo"/> XMake

#### <img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/cplusplus/cplusplus-original.svg" width="24" height="24" alt="C++ Logo"/> C++ Компиляторы:
- Clang
- GCC
- MSVC

#### <img src="https://images.sftcdn.net/images/t_app-icon-m/p/4f6f9692-96da-11e6-9846-00163ed833e7/2948317542/gnuplot-gnuplot-logo.png" width="24" height="24" alt="Gnuplot Logo"/> Gnuplot

Далее выполните следующий набор команд:
```bash
git clone https://github.com/Dimandragon/non-parametric-dsp.git
cd non-parametric-dsp
git submodule update --init --recursive
./build_matplot.sh
xmake
```
И, наконец, попробуйте выполнить примеры
```bash
xmake r [имя_примера]
```
Названия примеров вы можете посмотреть в файле xmake.lua; все примеры являются “binary” таргетами системы сборки xmake


# Credits
non-parametric-dsp использует следующий проекты:
- [pocketfft](https://github.com/mreineck/pocketfft) - быстрая легковесная C++ реализации быстрого преобразования Фурье для работы с произвольными размерами входных данных
- [icecream](https://github.com/renatoGarcia/icecream-cpp) - C++ библиотека для простого и удобного формативрования вывода
- [matplotplusplus](https://github.com/alandefreitas/matplotplusplus) - C++ библиотека для визуализации данных
- [gnuplot](http://www.gnuplot.info/) - консольный инструмент для визуализации данных
- [boost](https://www.boost.org/) - большая коллекция разнообразных C++ библиотек
- [alglib](https://www.alglib.net/) - обширная библиотека численных методов
- [xmake](https://github.com/xmake-io/xmake) - удобная декларативная кроссплатформенная система сборки, написанная на lua

# Лицензия
[MIT license](LICENSE)

---

# Контакты
[![Email](https://img.shields.io/badge/Email-D14836?style=flat&logo=gmail&logoColor=white)](mailto:dkuznetsov071105@gmail.com)
```
dkuznetsov071105@gmail.com
```

[![Telegram](https://img.shields.io/badge/Telegram-2CA5E0?style=flat&logo=telegram&logoColor=white)](https://t.me/diman_botan)
```
@diman_botan
```