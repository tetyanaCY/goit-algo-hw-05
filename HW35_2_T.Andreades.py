def binary_search(arr, target):
    left, right = 0, len(arr) - 1
    iterations = 0
    upper_bound = None

    while left <= right:
        iterations += 1
        mid = left + (right - left) // 2

        if arr[mid] == target:
            # Якщо значення знайдено, верхня межа це знайдене значення
            upper_bound = arr[mid]
            break
        elif arr[mid] < target:
            left = mid + 1
        else:
            right = mid - 1
            upper_bound = arr[mid]  # Оновлюємо верхню межу

    if upper_bound is None and left < len(arr):
        upper_bound = arr[left]  # Встановлюємо верхню межу на найменший елемент праворуч

    return iterations, upper_bound

# Тестування функції:
arr = [1.5, 2.3, 4.7, 6.1, 7.8, 9.4]
print(binary_search(arr, 5))  # Наприклад, пошук значення 5 у відсортованому масиві
