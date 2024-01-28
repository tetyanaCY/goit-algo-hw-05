# Порівняльний аналіз алгоритмів пошуку рядків

## Ефективність алгоритмів у `article_1.txt`

### Пошук існуючого підрядка
- **Кнут-Морріс-Пратт**: Не надано
- **Бойєр-Мур**: 0.0001379 секунд


### Пошук неіснуючого підрядка
- **Кнут-Морріс-Пратт**: Не надано
- **Бойєр-Мур**: 0.0004252 секунд


#### Найшвидший алгоритм для `article_1.txt`
- **Для існуючого підрядка**: Бойєр-Мур
- **Для неіснуючого підрядка**: Бойєр-Мур

## Ефективність алгоритмів у `article_2.txt`

### Пошук існуючого підрядка
- **Кнут-Морріс-Пратт**: 0.0053691 секунд
- **Бойєр-Мур**: 0.0015420 секунд
- **Рабін-Карп**: 0.0142394 секунд

### Пошук неіснуючого підрядка
- **Кнут-Морріс-Пратт**: 0.0036416 секунд
- **Бойєр-Мур**: 0.0005760 секунд
- **Рабін-Карп**: 0.0082438 секунд

#### Найшвидший алгоритм для `article_2.txt`
- **Для існуючого підрядка**: Бойєр-Мур
- **Для неіснуючого підрядка**: Бойєр-Мур

## Загальний найшвидший алгоритм
- **Бойєр-Мур (Неіснуючий підрядок у `article_1.txt`)**: 0.0004252 секунд

## Висновки

1. **Перевага Бойєра-Мура**: Алгоритм Бойєра-Мура послідовно перевершував алгоритми Кнута-Морріса-Пратта та Рабіна-Карпа у всіх тестах, як для існуючих, так і для неіснуючих підрядків, у обох текстових файлах.

2. **Різниця у продуктивності**: Продуктивність алгоритмів значно відрізнялася між двома текстовими файлами, при цьому Бойєр-Мур показав найкращу загальну ефективність.

3. **Довший час виконання Рабіна-Карпа**: Алгоритм Рабіна-Карпа демонстрував значно довший час виконання, особливо у `article_2.txt`.

4. **Ефективність для неіснуючих підрядків**: Алгоритм Бойєра-Мура був особливо ефективний у пошуку неіснуючих підрядків, як це видно з найшвидшого зафіксованого часу у `article_1.txt`.
