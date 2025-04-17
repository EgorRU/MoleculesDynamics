# MoleculesDynamics

**MoleculesDynamics** — это приложение для моделирования молекулярной динамики с 3D-визуализацией. Поддерживает шесть численных методов интегрирования и предоставляет сравнение точности по среднеквадратичному отклонению (MSE).

## Возможности

- **6 численных методов интегрирования:**
  - Метод Рунге-Кутта 4 порядка
  - Метод Верле
  - Скоростной метод Верле
  - Leapfrog (Метод с перескоками)
  - Метод Бимана-Шофилда
  - Метод предиктор-корректор

- **Визуализация**
  - 6 интерактивных 3D-графиков (по одному на каждый метод)
  - Отображение положения молекул в реальном времени
  - Подписи и графики обновляются на каждом шаге моделирования

- **Настраиваемые параметры:**
  - Количество молекул: 2 – 1000
  - Плотность ρ: 0.1 – 5.0
  - Масса: 0.1 – 1000
  - Потенциал Леннарда-Джонса:
    - ε (глубина ямы): 0.1 – 1000
    - σ (радиус взаимодействия): 0.1 – 1000
  - Начальные скорости Δ: 0.01 – 2
  - Шаг по времени dt: 0.000001 – 0.01
  - Количество шагов моделирования: 100 – 100000

- **Работа с MSE:**
  - Выбор эталонного метода для сравнения
  - Расчёт MSE на каждом шаге и усреднение по ходу симуляции
  - Отображение MSE для каждого метода

- **Дополнительные функции:**
  - Отражающие граничные условия
  - Пауза и продолжение симуляции
  - Сброс положения камер графиков
  - Генерация новых начальных условий или перезапуск с сохранёнными

##  Интерфейс

### Правая панель управления:

- Поля ввода параметров с пояснением действия при изменении:
  - **Генерация+перезапуск**: генерирует новые начальные условия
  - **Перезапуск**: перезапускает симуляцию с текущими условиями
  - **Без действия**: изменение не влияет на симуляцию

- **Выбор метода-эталона для MSE**
- **Отображение MSE для всех методов**
- **Кнопки управления:**
  - "Сгенерировать новые начальные условия"
  - "Начать моделирование заново на старых начальных условиях"
  - "Пауза/Продолжить"
  - "Сброс камер для графиков"